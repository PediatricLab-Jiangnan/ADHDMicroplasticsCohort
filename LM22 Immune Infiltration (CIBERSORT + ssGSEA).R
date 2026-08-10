if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer")
bioc_pkgs <- c("GSVA", "BiocParallel", "limma")

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

setwd("Your own pathway")

# --- Parameters -------------------------------------------------------------
expression_is_log2_tpm <- TRUE
group_levels <- c("Control", "Micro")
group_colors <- c("Control" = "#4DBBD5", "Micro" = "#E64B35")
cell_order_mode <- "LM22"   # "LM22" | "P_value" | "P_adjust"
show_ns <- TRUE
filter_by_cibersort_p <- FALSE
ssgsea_top_genes_per_cell <- 80
ssgsea_min_gene_set_size <- 10

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

# ============================================================================
# 4.1 CIBERSORT
# ============================================================================
cat(">>> Preparing CIBERSORT input...\n")

lm22_expr <- expr_data
lm22_expr <- lm22_expr[!is.na(rownames(lm22_expr)) & trimws(rownames(lm22_expr)) != "", , drop = FALSE]

if (expression_is_log2_tpm) lm22_expr <- 2^lm22_expr - 1
lm22_expr[!is.finite(lm22_expr)] <- 0
lm22_expr[lm22_expr < 0] <- 0

lm22_df <- data.frame(GeneSymbol = toupper(trimws(rownames(lm22_expr))), lm22_expr, check.names = FALSE) %>%
    group_by(GeneSymbol) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

lm22_values <- as.matrix(lm22_df[, -1, drop = FALSE])
lm22_df <- lm22_df[rowSums(lm22_values, na.rm = TRUE) > 0, , drop = FALSE]

write.table(lm22_df, "CIBERSORT_LM22_mixture.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

if (!file.exists("CIBERSORT.R")) stop("CIBERSORT.R not found.")
if (!file.exists("LM22.txt")) stop("LM22.txt not found.")

source("CIBERSORT.R")
cibersort_result <- CIBERSORT(sig_matrix = "LM22.txt", mixture_file = "CIBERSORT_LM22_mixture.txt", perm = 1000, QN = FALSE)
write.csv(cibersort_result, "CIBERSORT_LM22_Result.csv", row.names = TRUE)

# --- Parse CIBERSORT results ------------------------------------------------
cibersort_df <- as.data.frame(cibersort_result, check.names = FALSE, stringsAsFactors = FALSE)
cibersort_df$Sample <- rownames(cibersort_df)

lm22_matrix <- as.matrix(read.delim("LM22.txt", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
storage.mode(lm22_matrix) <- "numeric"
immune_cell_names <- colnames(lm22_matrix)
immune_cell_columns <- intersect(immune_cell_names, colnames(cibersort_df))

if (length(immune_cell_columns) == 0) stop("No LM22 cell columns found in CIBERSORT output.")

# Optional: filter by CIBERSORT sample P-value
if (filter_by_cibersort_p) {
    p_candidates <- c("P-value", "P.value", "P_Value", "Pvalue", "p.value")
    p_col <- intersect(p_candidates, colnames(cibersort_df))
    if (length(p_col) > 0) {
        p_col <- p_col[1]
        before <- nrow(cibersort_df)
        cibersort_df <- cibersort_df[!is.na(cibersort_df[[p_col]]) & cibersort_df[[p_col]] < 0.05, , drop = FALSE]
        cat("CIBERSORT P-filter:", before, "->", nrow(cibersort_df), "samples.\n")
    }
}

group_df <- data.frame(Sample = rownames(sample_info), Group = sample_info$Group, stringsAsFactors = FALSE)
cibersort_df$Sample <- as.character(cibersort_df$Sample)
group_df$Sample <- as.character(group_df$Sample)

cibersort_plot_data <- cibersort_df %>%
    select(Sample, all_of(immune_cell_columns)) %>%
    left_join(group_df, by = "Sample")

if (any(is.na(cibersort_plot_data$Group))) stop("Sample name mismatch between CIBERSORT and sample_info.")

cibersort_plot_data$Group <- factor(cibersort_plot_data$Group, levels = group_levels)
cibersort_plot_data[immune_cell_columns] <- lapply(cibersort_plot_data[immune_cell_columns], function(x) as.numeric(as.character(x)))

write.csv(cibersort_plot_data, "CIBERSORT_LM22_Fraction_With_Group.csv", row.names = FALSE)

cibersort_long <- cibersort_plot_data %>%
    pivot_longer(cols = all_of(immune_cell_columns), names_to = "ImmuneCell", values_to = "Fraction") %>%
    filter(!is.na(Fraction), !is.na(Group))

write.csv(cibersort_long, "CIBERSORT_LM22_Long_Format.csv", row.names = FALSE)

# --- CIBERSORT statistics (raw + BH) ----------------------------------------
cibersort_statistics <- cibersort_long %>%
    group_by(ImmuneCell) %>%
    summarise(
        Control_n = sum(Group == group_levels[1]), Micro_n = sum(Group == group_levels[2]),
        Control_mean = mean(Fraction[Group == group_levels[1]], na.rm = TRUE),
        Micro_mean = mean(Fraction[Group == group_levels[2]], na.rm = TRUE),
        Difference = Micro_mean - Control_mean,
        P_value = tryCatch(wilcox.test(Fraction ~ Group, exact = FALSE)$p.value, error = function(e) NA_real_),
        fraction_min = min(Fraction, na.rm = TRUE), fraction_max = max(Fraction, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        P_adjust = p.adjust(P_value, method = "BH"),
        Raw_significance = case_when(is.na(P_value) ~ "", P_value < 0.001 ~ "***", P_value < 0.01 ~ "**", P_value < 0.05 ~ "*", TRUE ~ "ns"),
        BH_significance = case_when(is.na(P_adjust) ~ "", P_adjust < 0.001 ~ "***", P_adjust < 0.01 ~ "**", P_adjust < 0.05 ~ "*", TRUE ~ "ns"),
        fraction_range = ifelse(is.finite(fraction_max - fraction_min) & (fraction_max - fraction_min) > 0, fraction_max - fraction_min, 0.05),
        label_y = fraction_max + 0.18 * fraction_range
    )

write.csv(cibersort_statistics, "CIBERSORT_LM22_Group_Statistics_Raw_and_BH.csv", row.names = FALSE)

# --- Cell ordering ----------------------------------------------------------
if (cell_order_mode == "P_value") {
    immune_cell_order <- cibersort_statistics %>% arrange(P_value) %>% pull(ImmuneCell)
} else if (cell_order_mode == "P_adjust") {
    immune_cell_order <- cibersort_statistics %>% arrange(P_adjust) %>% pull(ImmuneCell)
} else {
    immune_cell_order <- immune_cell_names[immune_cell_names %in% immune_cell_columns]
}

cibersort_long$ImmuneCell <- factor(cibersort_long$ImmuneCell, levels = immune_cell_order)
cibersort_statistics$ImmuneCell <- factor(cibersort_statistics$ImmuneCell, levels = immune_cell_order)
cibersort_statistics <- cibersort_statistics %>% arrange(ImmuneCell)

# --- CIBERSORT facet plot function ------------------------------------------
make_cibersort_facet_plot <- function(label_column, use_violin = TRUE, subtitle_text, output_prefix) {
    label_data <- cibersort_statistics
    label_data$PlotLabel <- label_data[[label_column]]
    if (!show_ns) label_data$PlotLabel[label_data$PlotLabel == "ns"] <- ""
    
    p <- ggplot(cibersort_long, aes(x = Group, y = Fraction, fill = Group))
    if (use_violin) p <- p + geom_violin(trim = FALSE, width = 0.90, alpha = 0.38, color = NA, na.rm = TRUE)
    
    p <- p +
        geom_boxplot(width = ifelse(use_violin, 0.25, 0.58), alpha = 0.88, outlier.shape = NA, color = "black", linewidth = 0.48, na.rm = TRUE) +
        geom_jitter(aes(color = Group), width = 0.13, height = 0, size = 1.35, alpha = 0.70, show.legend = FALSE, na.rm = TRUE) +
        geom_text(data = label_data, aes(x = 1.5, y = label_y, label = PlotLabel), inherit.aes = FALSE, size = 4.0, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ ImmuneCell, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = group_colors, drop = FALSE) +
        scale_color_manual(values = group_colors, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
        labs(title = "CIBERSORT LM22 Immune-Cell Infiltration", subtitle = subtitle_text, x = NULL, y = "Estimated immune-cell fraction", fill = NULL) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey35", margin = margin(b = 10)),
              strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.4),
              strip.text = element_text(face = "bold", size = 9),
              axis.text.x = element_text(color = "black", size = 9),
              axis.text.y = element_text(color = "black", size = 8),
              axis.title.y = element_text(face = "bold", size = 11),
              axis.line = element_line(color = "black", linewidth = 0.5),
              legend.position = "top",
              panel.spacing = unit(0.75, "lines"))
    
    ggsave(paste0(output_prefix, ".pdf"), p, width = 13, height = 15, limitsize = FALSE)
    ggsave(paste0(output_prefix, ".tiff"), p, width = 13, height = 15, dpi = 600, compression = "lzw", limitsize = FALSE)
    return(p)
}

make_cibersort_facet_plot("Raw_significance", TRUE, "Wilcoxon rank-sum test; stars = unadjusted P values (exploratory)", "CIBERSORT_LM22_Facet_Violin_Boxplot_RawP")
make_cibersort_facet_plot("BH_significance", TRUE, "Wilcoxon rank-sum test; stars = BH-adjusted P values", "CIBERSORT_LM22_Facet_Violin_Boxplot_BH")
make_cibersort_facet_plot("Raw_significance", FALSE, "Wilcoxon rank-sum test; stars = unadjusted P values (exploratory)", "CIBERSORT_LM22_Facet_Boxplot_RawP")
make_cibersort_facet_plot("BH_significance", FALSE, "Wilcoxon rank-sum test; stars = BH-adjusted P values", "CIBERSORT_LM22_Facet_Boxplot_BH")

# --- Article-style horizontal boxplot ---------------------------------------
make_cibersort_article_plot <- function(label_column, subtitle_text, output_prefix) {
    label_data <- cibersort_statistics
    label_data$PlotLabel <- label_data[[label_column]]
    if (!show_ns) label_data$PlotLabel[label_data$PlotLabel == "ns"] <- ""
    
    p <- ggplot(cibersort_long, aes(x = ImmuneCell, y = Fraction, fill = Group)) +
        geom_boxplot(position = position_dodge(width = 0.78), width = 0.66, alpha = 0.82, outlier.shape = NA, linewidth = 0.45) +
        geom_point(aes(color = Group), position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.78), size = 1.30, alpha = 0.58, show.legend = FALSE) +
        geom_text(data = label_data, aes(x = ImmuneCell, y = label_y, label = PlotLabel), inherit.aes = FALSE, size = 3.7, fontface = "bold") +
        scale_fill_manual(values = group_colors, drop = FALSE) +
        scale_color_manual(values = group_colors, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
        labs(title = "CIBERSORT LM22 Immune-Cell Infiltration", subtitle = subtitle_text, x = NULL, y = "Estimated immune-cell fraction", fill = NULL) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
              plot.subtitle = element_text(hjust = 0.5, size = 10.5, color = "grey35"),
              axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.title.y = element_text(face = "bold", size = 12),
              legend.position = "top",
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              plot.margin = margin(10, 15, 10, 10))
    
    ggsave(paste0(output_prefix, ".pdf"), p, width = 20, height = 9, limitsize = FALSE)
    ggsave(paste0(output_prefix, ".tiff"), p, width = 20, height = 9, dpi = 600, compression = "lzw", limitsize = FALSE)
    return(p)
}

make_cibersort_article_plot("Raw_significance", "Wilcoxon rank-sum test; stars = unadjusted P values (exploratory)", "CIBERSORT_LM22_Article_Boxplot_RawP")
make_cibersort_article_plot("BH_significance", "Wilcoxon rank-sum test; stars = BH-adjusted P values", "CIBERSORT_LM22_Article_Boxplot_BH")

cat("CIBERSORT nominal P < 0.05:", nrow(cibersort_statistics %>% filter(!is.na(P_value), P_value < 0.05)), "\n")
cat("CIBERSORT BH P < 0.05:", nrow(cibersort_statistics %>% filter(!is.na(P_adjust), P_adjust < 0.05)), "\n")

# ============================================================================
# 4.2 ssGSEA
# ============================================================================
cat("\n>>> Running LM22-ssGSEA...\n")

expr_ssgsea <- expr_data
if (!expression_is_log2_tpm) expr_ssgsea <- log2(pmax(expr_ssgsea, 0) + 1)
expr_ssgsea[!is.finite(expr_ssgsea)] <- NA_real_
expr_ssgsea <- expr_ssgsea[rowSums(is.na(expr_ssgsea)) == 0, , drop = FALSE]

rownames(expr_ssgsea) <- toupper(trimws(rownames(expr_ssgsea)))
rownames(lm22_matrix) <- toupper(trimws(rownames(lm22_matrix)))

expr_ssgsea <- limma::avereps(expr_ssgsea, ID = rownames(expr_ssgsea))
ssgsea_keep <- apply(expr_ssgsea, 1, function(x) sd(x, na.rm = TRUE) > 0)
ssgsea_keep[is.na(ssgsea_keep)] <- FALSE
expr_ssgsea <- expr_ssgsea[ssgsea_keep, , drop = FALSE]

# Build gene sets from LM22 weights
signature_lm22 <- list()
signature_lm22_info <- list()

for (i in seq_len(ncol(lm22_matrix))) {
    cell_type <- colnames(lm22_matrix)[i]
    weights_i <- lm22_matrix[, i]
    names(weights_i) <- rownames(lm22_matrix)
    positive_weights <- weights_i[is.finite(weights_i) & weights_i > 0]
    positive_weights <- sort(positive_weights, decreasing = TRUE)
    candidate_genes <- names(positive_weights)
    matched_genes <- candidate_genes[candidate_genes %in% rownames(expr_ssgsea)]
    final_genes <- head(matched_genes, ssgsea_top_genes_per_cell)
    
    if (length(final_genes) >= ssgsea_min_gene_set_size) {
        signature_lm22[[cell_type]] <- final_genes
        signature_lm22_info[[cell_type]] <- data.frame(
            Cell_Type = cell_type,
            Positive_Genes_in_LM22 = length(candidate_genes),
            Matched_Genes = length(matched_genes),
            Used_Genes = length(final_genes),
            Top5_Genes = paste(head(final_genes, 5), collapse = ", "),
            stringsAsFactors = FALSE
        )
    }
}

if (length(signature_lm22) == 0) stop("No LM22 cell gene sets met the minimum gene requirement.")

signature_lm22_info_df <- bind_rows(signature_lm22_info)
write.csv(signature_lm22_info_df, "LM22_ssGSEA_Gene_Set_Info.csv", row.names = FALSE)
cat("ssGSEA retained", length(signature_lm22), "LM22 immune-cell gene sets.\n")
print(signature_lm22_info_df)

# Run ssGSEA
if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    ssgsea_parameter <- GSVA::ssgseaParam(exprData = expr_ssgsea, geneSets = signature_lm22,
                                          minSize = ssgsea_min_gene_set_size, maxSize = 500, alpha = 0.25, normalize = TRUE)
    ssgsea_score <- GSVA::gsva(ssgsea_parameter, verbose = TRUE, BPPARAM = BiocParallel::SerialParam())
} else {
    ssgsea_score <- GSVA::gsva(expr = expr_ssgsea, gset.idx.list = signature_lm22, method = "ssgsea",
                               min.sz = ssgsea_min_gene_set_size, max.sz = 500, tau = 0.25, ssgsea.norm = TRUE,
                               verbose = TRUE, parallel.sz = 1)
}

ssgsea_score <- as.matrix(ssgsea_score)[, rownames(sample_info), drop = FALSE]

write.csv(ssgsea_score, "LM22_ssGSEA_Score_CellBySample.csv", row.names = TRUE)
write.csv(t(ssgsea_score), "LM22_ssGSEA_Score_SampleByCell.csv", row.names = TRUE)
save(ssgsea_score, signature_lm22, signature_lm22_info_df, sample_info, file = "LM22_ssGSEA_Results.RData")

# Long format
ssgsea_long <- as.data.frame(t(ssgsea_score), check.names = FALSE)
ssgsea_long$Sample <- rownames(ssgsea_long)
ssgsea_group_df <- data.frame(Sample = rownames(sample_info), Group = sample_info$Group, stringsAsFactors = FALSE)

ssgsea_long <- ssgsea_long %>%
    left_join(ssgsea_group_df, by = "Sample") %>%
    pivot_longer(cols = -c(Sample, Group), names_to = "Cell_Type", values_to = "Score")

ssgsea_long$Group <- factor(ssgsea_long$Group, levels = group_levels)
ssgsea_cell_order <- immune_cell_names[immune_cell_names %in% rownames(ssgsea_score)]
ssgsea_long$Cell_Type <- factor(ssgsea_long$Cell_Type, levels = ssgsea_cell_order)
write.csv(ssgsea_long, "LM22_ssGSEA_Long_Format.csv", row.names = FALSE)

# ssGSEA statistics
ssgsea_statistics <- ssgsea_long %>%
    group_by(Cell_Type) %>%
    summarise(
        Control_n = sum(Group == group_levels[1]), Micro_n = sum(Group == group_levels[2]),
        Control_mean = mean(Score[Group == group_levels[1]], na.rm = TRUE),
        Micro_mean = mean(Score[Group == group_levels[2]], na.rm = TRUE),
        Difference = Micro_mean - Control_mean,
        P_value = tryCatch(wilcox.test(Score ~ Group, exact = FALSE)$p.value, error = function(e) NA_real_),
        score_min = min(Score, na.rm = TRUE), score_max = max(Score, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        P_adjust = p.adjust(P_value, method = "BH"),
        Raw_significance = case_when(is.na(P_value) ~ "", P_value < 0.001 ~ "***", P_value < 0.01 ~ "**", P_value < 0.05 ~ "*", TRUE ~ "ns"),
        BH_significance = case_when(is.na(P_adjust) ~ "", P_adjust < 0.001 ~ "***", P_adjust < 0.01 ~ "**", P_adjust < 0.05 ~ "*", TRUE ~ "ns"),
        score_range = ifelse(is.finite(score_max - score_min) & (score_max - score_min) > 0, score_max - score_min, 0.05),
        label_y = score_max + 0.18 * score_range
    )

write.csv(ssgsea_statistics, "LM22_ssGSEA_Group_Statistics_Raw_and_BH.csv", row.names = FALSE)

# --- ssGSEA heatmap ---------------------------------------------------------
ssgsea_scaled <- t(scale(t(ssgsea_score)))
ssgsea_scaled[ssgsea_scaled > 3] <- 3
ssgsea_scaled[ssgsea_scaled < -3] <- -3

annotation_col <- data.frame(Group = sample_info$Group, row.names = rownames(sample_info))
annotation_colors <- list(Group = group_colors)

pdf("LM22_ssGSEA_Heatmap.pdf", width = 14, height = 9)
pheatmap(ssgsea_scaled, show_colnames = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "LM22 Immune-Cell Enrichment by ssGSEA", fontsize_row = 9, border_color = NA)
dev.off()

# --- ssGSEA cell-cell correlation -------------------------------------------
ssgsea_cell_cor <- cor(t(ssgsea_score), method = "spearman", use = "pairwise.complete.obs")
pdf("LM22_ssGSEA_Cell_Correlation.pdf", width = 12, height = 11)
pheatmap(ssgsea_cell_cor, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-1, 1, length.out = 101), main = "LM22 Cell-Cell Correlation (Spearman)", fontsize = 8, border_color = "grey90")
dev.off()

# --- ssGSEA distribution plot -----------------------------------------------
p_dist <- ggplot(ssgsea_long, aes(x = Cell_Type, y = Score, fill = Cell_Type)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(ssgsea_cell_order))) +
    labs(title = "LM22 Immune-Cell Score Distribution", x = NULL, y = "ssGSEA score") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold", size = 9),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

ggsave("LM22_ssGSEA_All_Cell_Distribution.pdf", p_dist, width = 18, height = 8, limitsize = FALSE)
ggsave("LM22_ssGSEA_All_Cell_Distribution.tiff", p_dist, width = 18, height = 8, dpi = 600, compression = "lzw", limitsize = FALSE)

# --- ssGSEA article-style boxplot -------------------------------------------
make_ssgsea_article_plot <- function(label_column, subtitle_text, output_prefix) {
    label_data <- ssgsea_statistics
    label_data$PlotLabel <- label_data[[label_column]]
    if (!show_ns) label_data$PlotLabel[label_data$PlotLabel == "ns"] <- ""
    
    p <- ggplot(ssgsea_long, aes(x = Cell_Type, y = Score, fill = Group)) +
        geom_boxplot(position = position_dodge(width = 0.80), width = 0.68, alpha = 0.82, outlier.shape = NA, linewidth = 0.45) +
        geom_point(aes(color = Group), position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.80), size = 1.30, alpha = 0.55, show.legend = FALSE) +
        geom_text(data = label_data, aes(x = Cell_Type, y = label_y, label = PlotLabel), inherit.aes = FALSE, size = 3.7, fontface = "bold") +
        scale_fill_manual(values = group_colors, drop = FALSE) +
        scale_color_manual(values = group_colors, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
        labs(title = "LM22 Immune-Cell Enrichment: Control vs Micro", subtitle = subtitle_text, x = NULL, y = "ssGSEA score", fill = NULL) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.title.y = element_text(face = "bold", size = 12),
              legend.position = "top",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
              plot.subtitle = element_text(hjust = 0.5, size = 10.5, color = "grey35"),
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              plot.margin = margin(10, 15, 10, 10))
    
    ggsave(paste0(output_prefix, ".pdf"), p, width = 20, height = 9, limitsize = FALSE)
    ggsave(paste0(output_prefix, ".tiff"), p, width = 20, height = 9, dpi = 600, compression = "lzw", limitsize = FALSE)
    return(p)
}

make_ssgsea_article_plot("Raw_significance", "Wilcoxon rank-sum test; stars = unadjusted P values (exploratory)", "LM22_ssGSEA_Article_Boxplot_RawP")
make_ssgsea_article_plot("BH_significance", "Wilcoxon rank-sum test; stars = BH-adjusted P values", "LM22_ssGSEA_Article_Boxplot_BH")

cat("\nssGSEA nominal P < 0.05:", nrow(ssgsea_statistics %>% filter(!is.na(P_value), P_value < 0.05)), "\n")
cat("ssGSEA BH P < 0.05:", nrow(ssgsea_statistics %>% filter(!is.na(P_adjust), P_adjust < 0.05)), "\n")
cat("\n=== Immune infiltration analysis completed ===\n")
