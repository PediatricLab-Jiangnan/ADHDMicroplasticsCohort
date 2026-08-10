if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr")
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

# --- Load full Limma results ------------------------------------------------
results <- read.csv("limma_all_results.csv", check.names = FALSE)

gene_list <- results$logFC
names(gene_list) <- results$Gene
gene_list <- gene_list[is.finite(gene_list) & !is.na(names(gene_list)) & names(gene_list) != ""]
gene_list <- sort(gene_list, decreasing = TRUE)

# --- Convert to Entrez ID ---------------------------------------------------
ids <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

gsea_map <- data.frame(SYMBOL = names(gene_list), logFC = as.numeric(gene_list)) %>%
    inner_join(ids, by = "SYMBOL") %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
    ungroup()

gsea_list <- gsea_map$logFC
names(gsea_list) <- gsea_map$ENTREZID
gsea_list <- sort(gsea_list, decreasing = TRUE)

# --- Run GSEA ---------------------------------------------------------------
cat(">>> Running GSEA...\n")
gsea_res <- gseKEGG(geneList = gsea_list, organism = "hsa", minGSSize = 10, maxGSSize = 500,
                    pvalueCutoff = 1, verbose = FALSE)

if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_res <- setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    gsea_all <- as.data.frame(gsea_res)
    write.csv(gsea_all, "GSEA_KEGG_All_Results.csv", row.names = FALSE)
    
    gsea_plot_data <- gsea_all %>% filter(pvalue < 0.2) %>% arrange(pvalue)
    num_plots <- nrow(gsea_plot_data)
    cat("Generating plots for", num_plots, "GSEA pathways...\n")
    
    if (num_plots > 0) {
        pdf("GSEA_All_Pathways_Plots.pdf", width = 10, height = 7)
        for (i in seq_len(num_plots)) {
            pathway_id <- gsea_plot_data$ID[i]
            pathway_desc <- gsea_plot_data$Description[i]
            pathway_nes <- round(gsea_plot_data$NES[i], 2)
            pathway_p <- format(gsea_plot_data$pvalue[i], digits = 3, scientific = TRUE)
            
            p_gsea <- gseaplot2(gsea_res, geneSetID = pathway_id,
                                title = paste0(pathway_desc, "\nNES = ", pathway_nes, "; P = ", pathway_p),
                                color = "firebrick", subplots = 1:3)
            print(p_gsea)
        }
        dev.off()
        cat("✓ GSEA plots saved.\n")
    }
} else {
    cat("⚠ GSEA analysis found no pathways.\n")
}

cat("\n=== GSEA completed ===\n")
