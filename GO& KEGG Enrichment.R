
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("ggplot2", "dplyr", "stringr", "scales")
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

setwd("Your Own")

bubble_top_n <- 20

# --- Helper functions -------------------------------------------------------
prepare_bubble_data <- function(enrich_df, top_n = 20, wrap_width = 46) {
    if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
    df <- enrich_df %>%
        filter(!is.na(Description), !is.na(GeneRatio), !is.na(Count), !is.na(pvalue)) %>%
        arrange(pvalue, desc(Count)) %>% slice_head(n = top_n)
    if (nrow(df) == 0) return(NULL)
    df$GeneRatioNum <- vapply(strsplit(df$GeneRatio, "/", fixed = TRUE),
                                function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    df$NegLog10P <- -log10(pmax(df$pvalue, .Machine$double.xmin))
    df$DescriptionWrap <- stringr::str_wrap(df$Description, width = wrap_width)
    lev <- df$DescriptionWrap[order(df$GeneRatioNum)]
    df$DescriptionWrap <- factor(df$DescriptionWrap, levels = unique(lev))
    return(df)
}

make_bubble_plot <- function(df, plot_title) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    ggplot(df, aes(x = GeneRatioNum, y = DescriptionWrap)) +
        geom_segment(aes(x = 0, xend = GeneRatioNum, y = DescriptionWrap, yend = DescriptionWrap),
                     color = "grey84", linewidth = 0.65, lineend = "round") +
        geom_point(aes(size = Count, fill = NegLog10P), shape = 21, color = "black", stroke = 0.35, alpha = 0.95) +
        scale_fill_gradientn(colors = c("#3B4CC0", "#6F92F3", "#B9D0F9", "#F7B89C", "#D73027"),
                             name = expression(-log[10](P))) +
        scale_size_continuous(range = c(4.5, 11), breaks = scales::pretty_breaks(n = 4), name = "Gene count") +
        scale_x_continuous(labels = scales::label_number(accuracy = 0.01),
                           expand = expansion(mult = c(0.01, 0.12))) +
        labs(title = plot_title, x = "Gene ratio", y = NULL) +
        guides(size = guide_legend(order = 1, override.aes = list(fill = "grey80", alpha = 1)),
               fill = guide_colorbar(order = 2, barheight = grid::unit(40, "mm"), barwidth = grid::unit(5, "mm"))) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 11) +
        theme(
            axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.x = element_line(color = "black", linewidth = 0.7),
            axis.text.x = element_text(color = "black", size = 10),
            axis.text.y = element_text(color = "black", size = 10, lineheight = 0.95, hjust = 1),
            axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 8)),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 14)),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9),
            legend.box = "vertical",
            plot.margin = margin(10, 20, 10, 10)
        )
}

# --- Load DEGs & convert IDs ------------------------------------------------
deg <- read.csv("limma_significant_DEGs.csv", check.names = FALSE)
if (nrow(deg) == 0) stop("No significant DEGs found.")

gene_ids <- bitr(unique(deg$Gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_ids <- gene_ids[!duplicated(gene_ids$ENTREZID), , drop = FALSE]

# --- GO BP ------------------------------------------------------------------
cat(">>> Running GO enrichment...\n")
ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.5, qvalueCutoff = 0.5, readable = TRUE)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    go_all <- as.data.frame(ego)
    write.csv(go_all, "GO_BP_Relaxed_Results.csv", row.names = FALSE)
    go_df <- prepare_bubble_data(go_all, bubble_top_n, wrap_width = 46)
    p_go <- make_bubble_plot(go_df, "GO Biological Process Enrichment")
    if (!is.null(p_go)) {
        go_height <- max(7.2, 0.36 * nrow(go_df) + 2)
        ggsave("GO_BP_Beautiful_Bubble.pdf", p_go, width = 10.5, height = go_height, limitsize = FALSE)
        ggsave("GO_BP_Beautiful_Bubble.tiff", p_go, width = 10.5, height = go_height, dpi = 600, compression = "lzw", limitsize = FALSE)
        cat("✓ GO bubble plot saved.\n")
    }
} else {
    cat("⚠ No GO terms found.\n")
}

# --- KEGG -------------------------------------------------------------------
cat(">>> Running KEGG enrichment...\n")
kk <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.5, qvalueCutoff = 0.5)

if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kk_all <- as.data.frame(kk)
    write.csv(kk_all, "KEGG_Relaxed_Results.csv", row.names = FALSE)
    kk_df <- prepare_bubble_data(kk_all, bubble_top_n, wrap_width = 46)
    p_kegg <- make_bubble_plot(kk_df, "KEGG Pathway Enrichment")
    if (!is.null(p_kegg)) {
        kegg_height <- max(7.2, 0.36 * nrow(kk_df) + 2)
        ggsave("KEGG_Beautiful_Bubble.pdf", p_kegg, width = 10.5, height = kegg_height, limitsize = FALSE)
        ggsave("KEGG_Beautiful_Bubble.tiff", p_kegg, width = 10.5, height = kegg_height, dpi = 600, compression = "lzw", limitsize = FALSE)
        cat("✓ KEGG bubble plot saved.\n")
    }
} else {
    cat("⚠ No KEGG pathways found.\n")
}

cat("\n=== GO/KEGG enrichment completed ===\n")
