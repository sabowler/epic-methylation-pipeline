#!/usr/bin/env Rscript
# Unsupervised Analysis of Illumina EPIC v2 Methylation Array Data
# Scott A. Bowler — Ndhlovu Lab, Weill Cornell Medicine
# 2026-02-13

# Installation (run once if packages not installed):
# Uncomment the block below on first use, then comment it out again.
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("minfi",
#                        "IlluminaHumanMethylationEPICv2manifest",
#                        "IlluminaHumanMethylationEPICv2anno.20a1.hg38"))
# install.packages(c("ggplot2", "pheatmap", "RColorBrewer",
#                    "gridExtra", "ggrepel", "umap"))

suppressPackageStartupMessages({
  library(minfi)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(ggrepel)
})

cat("EPIC v2 Methylation Array — Unsupervised Analysis\n\n")

# Set working directory
setwd("/path/to/your/data")
base_dir   <- getwd()
output_dir <- file.path(base_dir, "methylation_analysis_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output directory: %s\n\n", output_dir))


# Step 1. Sample sheet and IDAT loading
cat("Step 1: Creating sample sheet and reading .idat files...\n")

sample_dirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
sample_dirs <- sort(sample_dirs[grepl("^ACCN\\d+$", sample_dirs)])
cat(sprintf("Found %d ACCN directories\n", length(sample_dirs)))

sample_sheet <- data.frame()
for (dir in sample_dirs) {
  dir_path  <- file.path(base_dir, dir)
  grn_files <- list.files(dir_path, pattern = "_Grn\\.idat$", full.names = FALSE)
  if (length(grn_files) > 0) {
    base_idat <- sub("_Grn\\.idat$", "", grn_files[1])
    red_file  <- paste0(base_idat, "_Red.idat")
    if (file.exists(file.path(dir_path, red_file))) {
      sample_sheet <- rbind(sample_sheet, data.frame(
        Sample_Name = dir,
        Basename    = base_idat,
        Directory   = dir,
        stringsAsFactors = FALSE
      ))
    }
  }
}
cat(sprintf("Found %d complete sample pairs (Grn + Red)\n", nrow(sample_sheet)))
print(sample_sheet)

sample_sheet$Basename_full <- file.path(base_dir, sample_sheet$Directory,
                                        sample_sheet$Basename)

cat("\nReading IDAT files (this may take a few minutes)...\n")
RGset <- read.metharray(basenames = sample_sheet$Basename_full,
                        verbose = TRUE, force = TRUE)
sampleNames(RGset)          <- sample_sheet$Sample_Name
annotation(RGset)["array"]  <- "IlluminaHumanMethylationEPICv2"
annotation(RGset)["annotation"] <- "20a1.hg38"

cat(sprintf("\nLoaded %d samples | Array: %s | Annotation: %s | Probes: %s\n\n",
            ncol(RGset), annotation(RGset)["array"],
            annotation(RGset)["annotation"], nrow(RGset)))


# Step 2. Quality control
cat("Step 2: Quality Control...\n")
detP <- detectionP(RGset)
qc_metrics <- data.frame(
  Sample         = sampleNames(RGset),
  Mean_detP      = colMeans(detP),
  Failed_probes  = colSums(detP > 0.01),
  Percent_failed = colSums(detP > 0.01) / nrow(detP) * 100
)

pdf(file.path(output_dir, "qc_report.pdf"), width = 12, height = 8)
barplot(colMeans(detP),
        las = 2, cex.names = 0.7,
        main = "Mean Detection P-values per Sample",
        ylab = "Mean Detection P-value",
        col = ifelse(colMeans(detP) > 0.01, "red", "darkgreen"))
abline(h = 0.01, col = "red", lty = 2, lwd = 2)
legend("topleft", legend = c("Good (p < 0.01)", "Poor (p > 0.01)"),
       fill = c("darkgreen", "red"))

barplot(qc_metrics$Percent_failed,
        names.arg = qc_metrics$Sample,
        las = 2, cex.names = 0.7,
        main = "Percentage of Failed Probes per Sample",
        ylab = "% Failed Probes (p > 0.01)",
        col = ifelse(qc_metrics$Percent_failed > 5, "red", "darkgreen"))
abline(h = 5, col = "red", lty = 2, lwd = 2)

qcReport(RGset, sampNames = sampleNames(RGset), pdf = NULL)
dev.off()

write.csv(qc_metrics, file.path(output_dir, "qc_metrics.csv"), row.names = FALSE)

cat("QC Summary:\n")
print(qc_metrics)
cat(sprintf("\nSamples with >5%% failed probes: %d\n\n",
            sum(qc_metrics$Percent_failed > 5)))


# Step 3. Preprocessing and normalization
cat("Step 3: Preprocessing and normalization...\n")

cat("  Running NOOB background correction...\n")
MSet <- preprocessNoob(RGset)

cat("  Calculating Beta and M values...\n")
beta <- getBeta(MSet)
M    <- getM(MSet)

failed_probes <- rowSums(detP > 0.01) > ncol(detP) * 0.5
cat(sprintf("  Removing %d failed probes\n", sum(failed_probes)))
beta_clean <- beta[!failed_probes, ]
M_clean    <- M[!failed_probes, ]

cat("  Filtering SNP and cross-reactive probes...\n")
tryCatch({
  beta_clean <- rmSNPandCH(beta_clean, what = "both", rmcrosshyb = TRUE)
  M_clean    <- M_clean[rownames(beta_clean), ]
  cat(sprintf("  Remaining probes after filtering: %d\n", nrow(beta_clean)))
}, error = function(e) {
  cat("  Warning: Could not remove SNP/cross-reactive probes — requires EPIC v2 annotation package\n")
})

cat(sprintf("\nFinal dataset: %d probes x %d samples\n\n", nrow(beta_clean), ncol(beta_clean)))

save(beta_clean, M_clean, sample_sheet,
     file = file.path(output_dir, "processed_data.RData"))


# Step 4. Principal component analysis
cat("Step 4: Principal Component Analysis...\n")

var_probes <- order(apply(M_clean, 1, var), decreasing = TRUE)[1:min(10000, nrow(M_clean))]
M_pca      <- M_clean[var_probes, ]

pca_result  <- prcomp(t(M_pca), scale. = TRUE, center = TRUE)
var_explained <- summary(pca_result)$importance[2, ] * 100

pca_data <- data.frame(
  Sample = colnames(M_pca),
  PC1    = pca_result$x[, 1],
  PC2    = pca_result$x[, 2],
  PC3    = pca_result$x[, 3],
  PC4    = pca_result$x[, 4]
)
write.csv(pca_data, file.path(output_dir, "pca_coordinates.csv"), row.names = FALSE)

pdf(file.path(output_dir, "pca_plots.pdf"), width = 12, height = 10)

p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, alpha = 0.7, color = "steelblue") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "PCA: PC1 vs PC2",
       x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_explained[2])) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(pca_data, aes(x = PC3, y = PC4, label = Sample)) +
  geom_point(size = 4, alpha = 0.7, color = "coral") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "PCA: PC3 vs PC4",
       x = sprintf("PC3 (%.1f%% variance)", var_explained[3]),
       y = sprintf("PC4 (%.1f%% variance)", var_explained[4])) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

scree_data     <- data.frame(PC = paste0("PC", 1:min(10, length(var_explained))),
                              Variance = var_explained[1:min(10, length(var_explained))])
scree_data$PC  <- factor(scree_data$PC, levels = scree_data$PC)

p3 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  geom_line(aes(group = 1), color = "red", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1); print(p2); print(p3)
grid.arrange(p1, p2, p3, ncol = 1)
dev.off()

cat("PCA plots saved\n\n")


# Step 5. Hierarchical clustering
cat("Step 5: Hierarchical Clustering...\n")

dist_samples <- dist(t(M_pca))
hc           <- hclust(dist_samples, method = "complete")

pdf(file.path(output_dir, "hierarchical_clustering.pdf"), width = 14, height = 8)
plot(hc,
     main = "Hierarchical Clustering Dendrogram (Most Variable Probes)",
     xlab = "Samples", ylab = "Distance", cex = 0.8, hang = -1)

for (k in 2:5) {
  plot(hc,
       main = sprintf("Hierarchical Clustering — %d Clusters", k),
       xlab = "Samples", ylab = "Distance", cex = 0.8, hang = -1)
  rect.hclust(hc, k = k, border = 2:6)
}
dev.off()

cluster_assignments <- data.frame(Sample = sample_sheet$Sample_Name)
for (k in 2:5) {
  clusters <- cutree(hc, k = k)
  cluster_assignments[paste0("k", k)] <- clusters[sample_sheet$Sample_Name]
}
write.csv(cluster_assignments,
          file.path(output_dir, "cluster_assignments.csv"), row.names = FALSE)

cat("Hierarchical clustering complete\n\n")


# Step 6. Heatmap of most variable probes
cat("Step 6: Heatmap of most variable probes...\n")

n_top_probes   <- min(1000, nrow(beta_clean))
top_var_probes <- order(apply(beta_clean, 1, var), decreasing = TRUE)[1:n_top_probes]
beta_heatmap   <- beta_clean[top_var_probes, ]

annotation_col <- data.frame(
  Cluster_k2 = factor(cutree(hc, k = 2)),
  Cluster_k3 = factor(cutree(hc, k = 3)),
  Cluster_k4 = factor(cutree(hc, k = 4))
)
rownames(annotation_col) <- colnames(beta_heatmap)

pdf(file.path(output_dir, "heatmap_top_variable_probes.pdf"), width = 12, height = 14)
pheatmap(beta_heatmap,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE, show_colnames = TRUE,
         annotation_col = annotation_col,
         main = sprintf("Top %d Most Variable Probes", n_top_probes),
         fontsize_col = 8, border_color = NA)
dev.off()

cat("Heatmap created\n\n")


# Step 7. Sample correlation matrix
cat("Step 7: Sample correlation matrix...\n")

cor_matrix <- cor(beta_clean, method = "pearson")

pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"), width = 12, height = 11)
pheatmap(cor_matrix,
         color = colorRampPalette(c("yellow", "orange", "red"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_distance_rows = as.dist(1 - cor_matrix),
         clustering_distance_cols = as.dist(1 - cor_matrix),
         clustering_method = "complete",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "Sample-to-Sample Correlation Heatmap",
         fontsize = 8, border_color = "grey90",
         annotation_col = annotation_col,
         annotation_row = annotation_col)
dev.off()

cat("Correlation heatmap created\n\n")


# Step 8. Multidimensional scaling
cat("Step 8: Multidimensional Scaling (MDS)...\n")

mds_result <- cmdscale(dist_samples, k = 4, eig = TRUE)
mds_var    <- mds_result$eig / sum(mds_result$eig) * 100

mds_data <- data.frame(
  Sample = colnames(M_pca),
  MDS1   = mds_result$points[, 1],
  MDS2   = mds_result$points[, 2],
  MDS3   = mds_result$points[, 3],
  MDS4   = mds_result$points[, 4]
)

pdf(file.path(output_dir, "mds_plots.pdf"), width = 12, height = 5)

p1 <- ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Sample)) +
  geom_point(size = 4, alpha = 0.7, color = "purple") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "MDS: Dimension 1 vs 2",
       x = sprintf("MDS1 (%.1f%% variance)", mds_var[1]),
       y = sprintf("MDS2 (%.1f%% variance)", mds_var[2])) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(mds_data, aes(x = MDS3, y = MDS4, label = Sample)) +
  geom_point(size = 4, alpha = 0.7, color = "darkgreen") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "MDS: Dimension 3 vs 4",
       x = sprintf("MDS3 (%.1f%% variance)", mds_var[3]),
       y = sprintf("MDS4 (%.1f%% variance)", mds_var[4])) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("MDS plots created\n\n")


# Step 9. UMAP
cat("Step 9: UMAP projection...\n")

tryCatch({
  library(umap)
  umap_result <- umap(t(M_pca), n_neighbors = min(15, ncol(M_pca) - 1))
  umap_data   <- data.frame(
    Sample = colnames(M_pca),
    UMAP1  = umap_result$layout[, 1],
    UMAP2  = umap_result$layout[, 2]
  )

  pdf(file.path(output_dir, "umap_plot.pdf"), width = 10, height = 8)
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, label = Sample)) +
    geom_point(size = 4, alpha = 0.7, color = "darkred") +
    geom_text_repel(size = 3, max.overlaps = 20) +
    labs(title = "UMAP Projection", x = "UMAP1", y = "UMAP2") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(p)
  dev.off()

  write.csv(umap_data, file.path(output_dir, "umap_coordinates.csv"), row.names = FALSE)
  cat("UMAP complete\n\n")

}, error = function(e) {
  cat("umap package not available — skipping. Install with: install.packages('umap')\n\n")
})


# Step 10. Summary report
cat("Step 10: Writing summary report...\n")

summary_text <- sprintf(
"METHYLATION ARRAY ANALYSIS SUMMARY
===================================
Date:                  %s
Array type:            EPIC v2
Samples:               %d
Total probes (raw):    %s
Probes after QC:       %d

QC
--
Samples with >5%% failed probes:  %d
Mean detection p-value range:     %.4f - %.4f

PCA
---
PC1: %.2f%%
PC2: %.2f%%
PC3: %.2f%%
PC4: %.2f%%
",
Sys.Date(),
ncol(beta_clean),
nrow(RGset),
nrow(beta_clean),
sum(qc_metrics$Percent_failed > 5),
min(qc_metrics$Mean_detP),
max(qc_metrics$Mean_detP),
var_explained[1],
var_explained[2],
var_explained[3],
var_explained[4])

cat(summary_text)
writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))

cat(sprintf("\nAnalysis complete. Results saved to: %s\n", output_dir))
