
library(ggplot2)
library(dplyr)

setwd("Your own pathway")

# --- Parameters -------------------------------------------------------------
group_levels <- c("Control", "ADHD")
group_colors <- c("Control" = "#4DBBD5", "ADHD" = "#E64B35")

# --- Load data --------------------------------------------------------------
expr_data <- read.csv("expression_matrix_for_limma.csv", row.names = 1, check.names = FALSE)
sample_info <- read.csv("sample_info.csv", row.names = 1, check.names = FALSE)

expr_data <- as.matrix(expr_data)
storage.mode(expr_data) <- "numeric"

common_samples <- intersect(colnames(expr_data), rownames(sample_info))
if (length(common_samples) < 2) stop("Fewer than 2 matched samples.")

expr_data <- expr_data[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]
sample_info$Group <- factor(sample_info$Group, levels = group_levels)

# --- PCA: top 5000 variable genes -------------------------------------------
gene_var <- apply(expr_data, 1, var, na.rm = TRUE)
expr_pca <- expr_data[gene_var > 0, , drop = FALSE]
top_idx <- order(gene_var[gene_var > 0], decreasing = TRUE)[1:min(5000, sum(gene_var > 0))]
expr_pca <- expr_pca[top_idx, , drop = FALSE]

pca_result <- prcomp(t(expr_pca), center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- sample_info[pca_df$Sample, "Group"]

var_exp <- summary(pca_result)$importance[2, 1:2] * 100

# --- Plot -------------------------------------------------------------------
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    stat_ellipse(level = 0.95, alpha = 0.12, geom = "polygon", show.legend = FALSE) +
    geom_point(size = 3.2, alpha = 0.9, stroke = 0.3) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    labs(
        title = "Principal Component Analysis",
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.35, linetype = "dashed"),
        panel.grid.minor = element_blank()
    )

ggsave("PCA_Plot.pdf", p_pca, width = 7, height = 6.5)
ggsave("PCA_Plot.tiff", p_pca, width = 7, height = 6.5, dpi = 600, compression = "lzw")
cat("✓ PCA plot saved.\n")
