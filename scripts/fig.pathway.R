# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(gridExtra)
library(reshape2)

# 创建输出目录
output_dir <- "../results/output_figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 设置工作目录 - 建议从项目根目录运行
# setwd("C:\\Users\\hao55\\Desktop\\MR\\Herpes_zoster")

# Create lists of significant genes from each omics layer
methylation_genes <- c("MRPL24", "SLC25A13", "MCL1", "MOCS1", "MCCD1")
expression_genes <- c("CPLX1", "CTS5", "HORMAD1", "HLA-B", "HLA-H", "HLA-W")
protein_genes <- c("DTYMK", "SNAP29", "NFU1", "ACP6", "MECR", "UNG")

# Combine all unique genes
all_genes <- unique(c(methylation_genes, expression_genes, protein_genes))

# 保存图表为PNG和PDF格式的辅助函数
save_plot_both_formats <- function(plot_obj, base_name, width = 10, height = 8) {
  # 设置文件路径
  png_file <- file.path(output_dir, paste0(base_name, ".png"))
  pdf_file <- file.path(output_dir, paste0(base_name, ".pdf"))
  
  # 保存PNG版本（高分辨率）
  ggsave(png_file, plot_obj, width = width, height = height, dpi = 300)
  cat("创建PNG文件:", png_file, "\n")
  
  # 保存PDF版本
  ggsave(pdf_file, plot_obj, width = width, height = height)
  cat("创建PDF文件:", pdf_file, "\n")
}

# Combine into a complete list with source information
all_sig_genes <- data.frame(
  gene = c(methylation_genes, expression_genes, protein_genes),
  source = c(rep("methylation", length(methylation_genes)),
             rep("expression", length(expression_genes)),
             rep("protein", length(protein_genes)))
)

# Add effect direction (from your results)
all_sig_genes$direction <- ifelse(all_sig_genes$gene %in% 
                                    c("MRPL24", "SLC25A13", "DTYMK", "SNAP29", "NFU1"), "risk", "protective")

# Convert gene symbols to Entrez IDs
gene_list <- bitr(all_genes, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Run GO enrichment analysis
go_enrich <- enrichGO(gene = gene_list$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.2,
                      qvalueCutoff = 0.5)

# Process and plot mitochondrial GO terms
if(nrow(go_enrich@result) > 0) {
  # Filter for mitochondrial terms
  mito_results <- go_enrich@result[grepl("mitochond|respir|electron|ATP|oxida", 
                                         go_enrich@result$Description, 
                                         ignore.case = TRUE), ]
  
  if(nrow(mito_results) > 0) {
    # Sort by p-value
    mito_results <- mito_results[order(mito_results$p.adjust), ]
    
    # Take top terms if there are many
    top_terms <- head(mito_results, 15)
    
    # Create a manual dotplot
    p1 <- ggplot(top_terms, aes(x = Count, y = reorder(Description, -p.adjust), 
                                size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      theme_bw() +
      labs(title = "Mitochondrial Pathways Associated with Herpes Zoster Risk",
           x = "Gene Count", y = "GO Term") +
      theme(axis.text.y = element_text(size = 10))
    
    # 保存第一个图表为PNG和PDF
    save_plot_both_formats(p1, "plot1_mitochondrial_pathways", width = 8, height = 10)
  } else {
    cat("No mitochondrial GO terms were found.\n")
  }
}

# Create gene lists for each omics layer
gene_lists <- list(
  Methylation = gene_list$ENTREZID[gene_list$SYMBOL %in% methylation_genes],
  Expression = gene_list$ENTREZID[gene_list$SYMBOL %in% expression_genes],
  Protein = gene_list$ENTREZID[gene_list$SYMBOL %in% protein_genes]
)

# Run compareCluster analysis
comp_go <- compareCluster(gene_lists, 
                          fun = "enrichGO",
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.2)

# Check if we have a valid compareClusterResult object
if(class(comp_go)[1] == "compareClusterResult" && nrow(comp_go@compareClusterResult) > 0) {
  # Create dotplot for omics comparison
  p2 <- dotplot(comp_go, showCategory = 5, 
                title = "Pathway Enrichment Across Omics Layers")
  
  # 保存第二个图表为PNG和PDF
  save_plot_both_formats(p2, "plot2_omics_comparison", width = 9, height = 12)
  
  # 过滤线粒体相关通路
  mito_compare_results <- comp_go@compareClusterResult[
    grepl("mitochond|respir|electron|ATP|oxida", 
          comp_go@compareClusterResult$Description, 
          ignore.case = TRUE), ]
  
  if(nrow(mito_compare_results) > 0) {
    # 准备热图数据
    plot_data <- dcast(mito_compare_results, Description ~ Cluster, 
                       value.var = "p.adjust", fill = 1)
    
    # 转换为长格式用于ggplot
    plot_data_long <- melt(plot_data, id.vars = "Description", 
                           variable.name = "Omics_Layer", value.name = "p.adjust")
    
    # 转换p值用于可视化（-log10）
    plot_data_long$neg_log10_p <- -log10(plot_data_long$p.adjust)
    plot_data_long$neg_log10_p[plot_data_long$neg_log10_p > 5] <- 5  # 设置上限
    
    # 创建热图
    p3 <- ggplot(plot_data_long, aes(x = Omics_Layer, y = Description, fill = neg_log10_p)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "navy", 
                          name = "-log10(p-value)") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Mitochondrial Pathways Across Omics Layers",
           x = "Omics Layer", y = "Pathway")
    
    # 保存第三个图表为PNG和PDF
    save_plot_both_formats(p3, "plot3_mitochondrial_heatmap", width = 10, height = 8)
  }
}

# 组合所有可用图表为一个大图
plots_to_combine <- list()
if(exists("p1")) plots_to_combine[[length(plots_to_combine)+1]] <- p1
if(exists("p2")) plots_to_combine[[length(plots_to_combine)+1]] <- p2
if(exists("p3")) plots_to_combine[[length(plots_to_combine)+1]] <- p3

if(length(plots_to_combine) > 0) {
  combined_plot <- do.call(grid.arrange, c(plots_to_combine, ncol=1))
  
  # 保存组合图表为PNG和PDF
  save_plot_both_formats(combined_plot, "all_mitochondrial_plots", 
                         width = 10, height = 7*length(plots_to_combine))
} else {
  cat("没有成功生成任何图表\n")
}

cat("\n所有图表已成功保存为PNG和PDF格式，存储在", output_dir, "文件夹中。\n") 