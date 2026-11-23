# nolint start
setwd(this.path::here())

library(tidyverse)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(heplots)
library(vegan)
library(MVN)

# ==========================================
# 1. Load Data & Align Sample IDs
# ==========================================
# Script resides in Project/Script; popgen_out resides in Project/White-Space/popgen_out
data_dir <- file.path("..", "White-Space", "popgen_out")

# Read covariance matrix, SNP weights, admixture Q (K=2), and BAM sample list
C <- as.matrix(read.table(file.path(data_dir, "pcangsd_output.cov")))
W <- as.matrix(read.table(file.path(data_dir, "pcangsd_output.weights")))
Q2 <- as.matrix(read.table(file.path(data_dir, "pcangsd_output.admix.2.Q")))
sites_kept <- scan(file.path(data_dir, "pcangsd_output.sites"), what = integer(), quiet = TRUE)

bams <- read.table(file.path(data_dir, "bam_list.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(bams) <- "id"
bams$id <- gsub(pattern = "^.*/", replacement = "", bams$id)
bams$id <- gsub(pattern = "\\.sorted\\.bam$|\\.bam$", replacement = "", bams$id)

rownames(C) <- bams$id
colnames(C) <- bams$id
rownames(Q2) <- bams$id
colnames(Q2) <- c("Cluster1", "Cluster2")

cat("=== Loaded PCAngsd Summary ===\n")
cat("Samples (n):", nrow(C), "\n")
cat("Total sites in mask:", length(sites_kept), "| Sites retained:", sum(sites_kept == 1), "\n")
cat("SNP weights dimensions:", nrow(W), "x", ncol(W), "\n\n")

# ==========================================
# 2. Eigen-Decomposition & Inferred Clusters
# ==========================================
# Derive K=2 cluster assignments from PCAngsd NMF admixture proportions
# (Used as inferred genetic clusters for downstream MVN, Box's M, and PERMANOVA)
inferred_cluster <- factor(ifelse(Q2[, "Cluster1"] >= Q2[, "Cluster2"], "Cluster_1", "Cluster_2"))

e <- eigen(C)
eigvals <- e$values
totvar <- sum(eigvals[eigvals > 0])
var_exp <- round((eigvals / totvar) * 100, 2)

cat("=== Eigenvalues & Percent Variance Explained ===\n")
eig_df <- data.frame(
  PC = paste0("PC", seq_along(eigvals)),
  Eigenvalue = round(eigvals, 6),
  Pct_Variance = var_exp,
  Cum_Pct_Variance = round(cumsum(var_exp), 2)
)
print(eig_df)
cat("\n")

# Extract top PCs for multivariate modeling (PC1-PC3)
pcadat <- as.data.frame(e$vectors[, 1:3])
colnames(pcadat) <- c("PC1", "PC2", "PC3")
pcadat$id <- bams$id
pcadat$Cluster <- inferred_cluster
pcadat$Cluster1_Prob <- round(Q2[, "Cluster1"], 4)
pcadat$Cluster2_Prob <- round(Q2[, "Cluster2"], 4)

cat("=== Sample Coordinates & Inferred Admixture Clusters ===\n")
print(pcadat)
cat("\n")

# ==========================================
# 3. SNP Loadings (Weights) & PCA Scatterplots
# ==========================================
# Inspect top SNP weights from PCAngsd (PC1 weights)
snp_weights_df <- data.frame(
  Site_Index = seq_len(nrow(W)),
  PC1_Weight = W[, 1],
  Abs_PC1_Weight = abs(W[, 1])
)
top_snps <- snp_weights_df %>%
  arrange(desc(Abs_PC1_Weight)) %>%
  slice_head(n = 20)

cat("=== Top 20 SNP Loadings (Weights) on PC1 ===\n")
print(top_snps)
cat("\n")

# Plot A: PCA Scatterplot (PC1 vs PC2) colored by inferred K=2 cluster
p_pca_12 <- ggplot(pcadat, aes(x = PC1, y = PC2, color = Cluster, label = id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text(vjust = -1.1, size = 4, fontface = "bold", show.legend = FALSE) +
  labs(
    title = "PCA of Oryza sativa Japonica Accessions (PCAngsd)",
    subtitle = "Colored by PCAngsd Admixture Cluster (K = 2)",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Inferred Cluster"
  ) +
  scale_color_manual(values = c("Cluster_1" = "#1b9e77", "Cluster_2" = "#d95f02")) +
  theme_bw() +
  theme(
    text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "top"
  )
print(p_pca_12)

# Plot B: PCA Scatterplot (PC1 vs PC3)
p_pca_13 <- ggplot(pcadat, aes(x = PC1, y = PC3, color = Cluster, label = id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text(vjust = -1.1, size = 4, fontface = "bold", show.legend = FALSE) +
  labs(
    title = "PCA of Oryza sativa Japonica Accessions (PC1 vs PC3)",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC3 (", var_exp[3], "%)"),
    color = "Inferred Cluster"
  ) +
  scale_color_manual(values = c("Cluster_1" = "#1b9e77", "Cluster_2" = "#d95f02")) +
  theme_bw() +
  theme(
    text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "top"
  )
print(p_pca_13)

# Plot C: Top 20 SNP Weights (Loadings) Barplot
p_loadings <- ggplot(top_snps, aes(x = reorder(factor(Site_Index), Abs_PC1_Weight), y = PC1_Weight, fill = PC1_Weight > 0)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2b8cbe", "FALSE" = "#e34a33"), labels = c("Negative", "Positive"), name = "Sign") +
  labs(
    title = "Top 20 SNP Weights (Loadings) Contributing to PC1",
    x = "Filtered SNP Site Index",
    y = "PCAngsd SNP Weight (PC1)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold"),
    legend.position = "top"
  )
print(p_loadings)

# Plot D: Admixture Barplot (K=2)
admix_long <- pcadat %>%
  select(id, Cluster, Cluster1_Prob, Cluster2_Prob) %>%
  pivot_longer(cols = c(Cluster1_Prob, Cluster2_Prob), names_to = "Ancestry", values_to = "Probability")

p_admix <- ggplot(admix_long, aes(x = Probability, y = factor(id), fill = Ancestry)) +
  geom_col(color = "white", linewidth = 0.3) +
  facet_grid(Cluster ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Cluster1_Prob" = "#1b9e77", "Cluster2_Prob" = "#d95f02"), labels = c("Cluster 1", "Cluster 2")) +
  labs(
    title = "PCAngsd Admixture Proportions (K = 2)",
    x = "Ancestry Proportion",
    y = "Accession (SRR)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold"),
    legend.position = "top"
  )
print(p_admix)

# ==========================================
# 4. Multivariate Normality (MVN) & Box's M
# ==========================================
# Evaluate MVN on PC1-PC2 (guarded with tryCatch for small n = 11 sample size)
cat("=== Multivariate Normality (MVN) on PC1 & PC2 ===\n")
mvn_res <- tryCatch(
  MVN::mvn(data = pcadat[, c("PC1", "PC2")], mvnTest = "mardia"),
  error = function(e) {
    cat("MVN Mardia test warning/error:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(mvn_res)) {
  print(mvn_res$multivariateNormality)
  print(mvn_res$univariateNormality)
}
cat("\n")

# Evaluate Box's M test for homogeneity of covariance matrices across inferred clusters
# Requires >= 2 clusters and sufficient observations per cluster (n_i > p)
cat("=== Box's M Test for Homogeneity of Covariance Matrices ===\n")
cluster_counts <- table(pcadat$Cluster)
print(cluster_counts)

if (length(cluster_counts) >= 2 && all(cluster_counts > 2)) {
  boxm_res <- tryCatch(
    heplots::boxM(pcadat[, c("PC1", "PC2")], group = pcadat$Cluster),
    error = function(e) {
      cat("Box's M test could not be computed (singular covariance or small group size):", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(boxm_res)) {
    print(boxm_res)
  }
} else {
  cat("Skipping Box's M: at least one inferred cluster has <= 2 samples (requires n_i > p = 2).\n")
}
cat("\n")

# ==========================================
# 5. PERMANOVA & Multivariate Dispersion (betadisper)
# ==========================================
# Convert genetic covariance matrix C to Euclidean distance matrix (or use PC coordinates)
# Using Euclidean distance across all eigenvectors / standardized covariance distance
dist_mat <- dist(scale(e$vectors[, 1:3]), method = "euclidean")

cat("=== PERMANOVA (adonis2) Across Inferred K=2 Clusters ===\n")
if (length(cluster_counts) >= 2) {
  set.seed(42)
  perm_res <- vegan::adonis2(dist_mat ~ Cluster, data = pcadat, permutations = 999)
  print(perm_res)
  cat("\n")

  cat("=== Homogeneity of Multivariate Dispersion (betadisper) ===\n")
  disp_res <- vegan::betadisper(dist_mat, group = pcadat$Cluster)
  print(disp_res)
  print(anova(disp_res))
  print(vegan::permutest(disp_res, permutations = 999))
} else {
  cat("Skipping PERMANOVA and betadisper: all samples assigned to a single cluster.\n")
}
# nolint end