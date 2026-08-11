
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("ggplot2", "ggrepel", "ggpubr", "dplyr", "scales")
bioc_pkgs <- c("limma")

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

setwd("Your own")

# --- Parameters -------------------------------------------------------------
p_cutoff <- 0.05
logfc_cutoff <- 0.4
group_levels <- c("Control", "ADHD")
group_colors <- c("Control" = "#4DBBD5", "ADHD" = "#E64B35")

# --- Load & filter ----------------------------------------------------------
expr_data <- read.csv("expression_matrix_for_limma.csv", row.names = 1, check.names = FALSE)
sample_info <- read.csv("sample_info.csv", row.names = 1, check.names = FALSE)

expr_data <- as.matrix(expr_data)
storage.mode(expr_data) <- "numeric"

common_samples <- intersect(colnames(expr_data), rownames(sample_info))
if (length(common_samples) < 2) stop("Fewer than 2 matched samples.")

expr_data <- expr_data[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]

if (!"Group" %in% colnames(sample_info)) stop("sample_info must contain 'Group'.")
if (!all(group_levels %in% unique(sample_info$Group))) stop("Group must contain Control and Micro.")

keep <- rowSums(expr_data > 1, na.rm = TRUE) >= ceiling(ncol(expr_data) / 2)
expr_filtered <- expr_data[keep, , drop = FALSE]
cat("Genes:", nrow(expr_data), "-> retained for limma:", nrow(expr_filtered), "\n")

# --- Limma ------------------------------------------------------------------
group <- factor(sample_info$Group, levels = group_levels)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

weights <- arrayWeights(expr_filtered, design)

fit <- lmFit(expr_filtered, design, weights = weights)
contrast_matrix <- makeContrasts(Micro - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE)

results <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
results$Gene <- rownames(results)

deg <- results[results$P.Value < p_cutoff & abs(results$logFC) > logfc_cutoff, , drop = FALSE]

write.csv(results, "limma_all_results.csv", row.names = FALSE)
write.csv(deg, "limma_significant_DEGs.csv", row.names = FALSE)
cat("Significant DEGs:", nrow(deg), "\n")

# --- Volcano plot -----------------------------------------------------------
results$Label <- "Ns"
results$Label[results$logFC > logfc_cutoff & results$P.Value < p_cutoff] <- "Up"
results$Label[results$logFC < -logfc_cutoff & results$P.Value < p_cutoff] <- "Down"

top_labels <- bind_rows(
    results %>% filter(Label == "Up") %>% arrange(P.Value) %>% slice_head(n = 10),
    results %>% filter(Label == "Down") %>% arrange(P.Value) %>% slice_head(n = 10)
)

p_vol <- ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Label)) +
    geom_point(alpha = 0.65, size = 1.6) +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", Ns = "grey82")) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_text_repel(data = top_labels, aes(label = Gene), max.overlaps = 20, size = 3, min.segment.length = 0) +
    labs(title = "Volcano Plot: Micro vs Control", x = expression(log[2]~Fold~Change), y = expression(-log[10](P)), color = NULL) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey93", linewidth = 0.35, linetype = "dashed"),
        panel.grid.minor = element_blank()
    )

ggsave("Volcano_Plot.pdf", p_vol, width = 8, height = 7)
ggsave("Volcano_Plot.tiff", p_vol, width = 8, height = 7, dpi = 600, compression = "lzw")
cat("✓ Volcano plot saved.\n")

# --- Target-gene boxplots ---------------------------------------------------
target_genes <- c("NFE2L2", "IL6", "IL6ST", "IL6R", "PIK3CA","PIK3CB", "AKT1", "NQO1",
                  "GCLC", "TNS1", "HMOX1", "KEAP1", "SOD1", "SOD2","DTHD1", "PIK3IP1")

for (gene in target_genes) {
    if (!gene %in% rownames(expr_data)) { cat("⚠ Gene", gene, "not found.\n"); next }
    
    plot_data <- data.frame(
        Expression = as.numeric(expr_data[gene, ]),
        Group = factor(sample_info$Group, levels = group_levels)
    )
    
    y_range <- diff(range(plot_data$Expression, na.rm = TRUE))
    if (!is.finite(y_range) || y_range == 0) y_range <- 1
    max_y <- max(plot_data$Expression, na.rm = TRUE)
    
    p_box <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
        stat_boxplot(geom = "errorbar", width = 0.2, linewidth = 0.45) +
        geom_boxplot(width = 0.5, alpha = 0.9, outlier.shape = NA, linewidth = 0.6) +
        geom_jitter(width = 0.15, size = 2, shape = 21, color = "black", fill = "white", alpha = 0.8) +
        stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5,
                           label.y = max_y + y_range * 0.10, size = 6) +
        stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5,
                           label.y = max_y + y_range * 0.22, size = 3.2) +
        scale_fill_manual(values = group_colors) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
        labs(title = gene, y = "Log2 TPM", x = NULL) +
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
            axis.text = element_text(color = "black", size = 12),
            axis.title.y = element_text(face = "bold"),
            axis.line = element_line(linewidth = 0.8),
            legend.position = "none"
        )
    
    ggsave(paste0("Boxplot_Single_", gene, ".pdf"), p_box, width = 4, height = 5)
    cat("✓ Boxplot saved:", gene, "\n")
}

cat("\n=== Differential analysis completed ===\n")
