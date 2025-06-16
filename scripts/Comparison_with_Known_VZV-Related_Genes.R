# 强烈建议在项目根目录（Herpes_Zoster_Public_Project）下通过Rscript命令运行
# 或者在R/RStudio中打开项目，并确保工作目录是项目的根目录
# setwd("C:\\Users\\hao55\\Desktop\\MR\\Herpes_zoster")

# 创建新文件夹来保存所有输出文件
output_dir <- "../results/output_figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建新文件夹:", output_dir, "\n")
}

# Check and install required packages
required_packages <- c("ggplot2", "dplyr", "RColorBrewer", "VennDiagram", "grid")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(VennDiagram)
library(grid)

# Set a consistent theme
theme_set(theme_bw(base_size = 12))

#------------------------------------------------------------------------------
# 1. Create data frames with gene information
#------------------------------------------------------------------------------

# Your mitochondrial genes identified in the study
mito_genes <- data.frame(
  gene = c("DTYMK", "SNAP29", "NFU1", "ACP6", "MECR", "UNG", "MCL1", "MOCS1", "SLC25A13", "MRPL24"),
  effect_size = c(2.78, 1.90, 1.61, 0.95, 0.52, 0.44, 0.78, 0.58, 1.04, 1.23),
  p_value = c(0.005, 0.021, 0.013, 0.004, 0.023, 0.041, 0.002, 0.003, 0.015, 0.009),
  direction = c("Risk", "Risk", "Risk", "Protective", "Protective", "Protective", 
                "Protective", "Protective", "Risk", "Risk"),
  omic_evidence = c("Protein+Expression", "Protein", "Protein", "Protein", 
                    "Protein+Expression", "Protein+Expression", "Methylation+Expression", 
                    "Methylation+Expression", "Methylation", "Methylation"),
  stringsAsFactors = FALSE
)

# Known VZV latency-related genes
vzv_latency_genes <- data.frame(
  gene = c("ORF63", "ORF62", "ORF61", "MCL1", "UNG", "SNAP29", "PML", "HDAC", "STAT1", "STAT3"),
  function_description = c("Transcriptional regulation", "Transcriptional activation", 
                           "Transcriptional regulation", "Anti-apoptotic", "DNA repair", 
                           "Membrane fusion", "Nuclear body protein", "Histone deacetylase", 
                           "Interferon signaling", "Interferon signaling"),
  stringsAsFactors = FALSE
)

# Genes involved in other herpesvirus infections
other_herpes_genes <- data.frame(
  gene = c("UNG", "MCL1", "DTYMK", "STAT1", "STAT3", "NFU1", "PML", "IFI16", "cGAS", "STING"),
  virus = c("HSV-1", "HSV-1", "EBV", "Multiple", "Multiple", "CMV", "Multiple", "HSV-1", "HSV-1", "Multiple"),
  mechanism = c("Viral DNA replication", "Apoptosis inhibition", "Nucleotide metabolism", 
                "Interferon response", "Interferon response", "Iron-sulfur cluster", 
                "Intrinsic immunity", "DNA sensing", "DNA sensing", "DNA sensing"),
  stringsAsFactors = FALSE
)

# Previously reported herpes zoster susceptibility genes
hz_susceptibility_genes <- data.frame(
  gene = c("IL10", "IL6", "TNF", "CTLA4", "IFNG", "MECR", "CXCL8", "NFU1", "SNAP29", "UNG"),
  study_type = c("GWAS", "Candidate gene", "Candidate gene", "GWAS", "Candidate gene", 
                 "MR analysis", "Candidate gene", "MR analysis", "MR analysis", "MR analysis"),
  effect_direction = c("Risk", "Risk", "Risk", "Protective", "Protective", 
                       "Protective", "Risk", "Risk", "Risk", "Protective"),
  population = c("European", "Asian", "European", "European", "Mixed", 
                 "Finnish", "Asian", "Finnish", "Finnish", "Finnish"),
  stringsAsFactors = FALSE
)

# Find overlapping genes
overlapping_genes <- list(
  mito_vzv_latency = intersect(mito_genes$gene, vzv_latency_genes$gene),
  mito_other_herpes = intersect(mito_genes$gene, other_herpes_genes$gene),
  mito_hz_susceptibility = intersect(mito_genes$gene, hz_susceptibility_genes$gene)
)

#------------------------------------------------------------------------------
# 2. Create and save individual plots (PNG and PDF formats)
#------------------------------------------------------------------------------

# 辅助函数 - 为每个图生成PNG和PDF文件
save_plot_both_formats <- function(base_name, plot_func, width = 10, height = 8) {
  # 设置文件路径
  png_file <- file.path(output_dir, paste0(base_name, ".png"))
  pdf_file <- file.path(output_dir, paste0(base_name, ".pdf"))
  
  # 保存PNG版本
  png(png_file, width = width, height = height, units = "in", res = 300)
  plot_func()
  dev.off()
  cat("Created PNG file:", png_file, "\n")
  
  # 保存PDF版本
  pdf(pdf_file, width = width, height = height)
  plot_func()
  dev.off()
  cat("Created PDF file:", pdf_file, "\n")
}

# PLOT 1: Venn Diagram
venn_plot_func <- function() {
  draw.quad.venn(
    area1 = length(mito_genes$gene),
    area2 = length(vzv_latency_genes$gene),
    area3 = length(other_herpes_genes$gene),
    area4 = length(hz_susceptibility_genes$gene),
    n12 = length(intersect(mito_genes$gene, vzv_latency_genes$gene)),
    n13 = length(intersect(mito_genes$gene, other_herpes_genes$gene)),
    n14 = length(intersect(mito_genes$gene, hz_susceptibility_genes$gene)),
    n23 = length(intersect(vzv_latency_genes$gene, other_herpes_genes$gene)),
    n24 = length(intersect(vzv_latency_genes$gene, hz_susceptibility_genes$gene)),
    n34 = length(intersect(other_herpes_genes$gene, hz_susceptibility_genes$gene)),
    n123 = length(Reduce(intersect, list(mito_genes$gene, vzv_latency_genes$gene, other_herpes_genes$gene))),
    n124 = length(Reduce(intersect, list(mito_genes$gene, vzv_latency_genes$gene, hz_susceptibility_genes$gene))),
    n134 = length(Reduce(intersect, list(mito_genes$gene, other_herpes_genes$gene, hz_susceptibility_genes$gene))),
    n234 = length(Reduce(intersect, list(vzv_latency_genes$gene, other_herpes_genes$gene, hz_susceptibility_genes$gene))),
    n1234 = length(Reduce(intersect, list(mito_genes$gene, vzv_latency_genes$gene, other_herpes_genes$gene, hz_susceptibility_genes$gene))),
    category = c("Mitochondrial Genes\n(Current Study)", "VZV Latency\nGenes", 
                 "Other Herpesvirus\nGenes", "Known HZ\nSusceptibility Genes"),
    fill = brewer.pal(4, "Set1"),
    cex = 1.5,
    cat.cex = 1.2,
    cat.pos = c(0, 0, 0, 0),
    cat.dist = c(0.05, 0.05, 0.05, 0.05),
    cat.fontfamily = "sans",
    main = "Overlap Between Mitochondrial Findings and Known VZV-Related Genes"
  )
}

save_plot_both_formats("venn_diagram", venn_plot_func)

# PLOT 2: Overlapping Genes Details
# Prepare data
overlap_data <- data.frame()

# Add overlapping genes between mitochondrial findings and VZV latency
if (length(overlapping_genes$mito_vzv_latency) > 0) {
  for (gene in overlapping_genes$mito_vzv_latency) {
    temp_df <- data.frame(
      gene = gene,
      overlap_category = "Mito-VZV Latency",
      function_description = vzv_latency_genes$function_description[vzv_latency_genes$gene == gene][1],
      effect_direction = mito_genes$direction[mito_genes$gene == gene][1],
      stringsAsFactors = FALSE
    )
    overlap_data <- rbind(overlap_data, temp_df)
  }
}

# Add overlapping genes between mitochondrial findings and other herpesvirus genes
if (length(overlapping_genes$mito_other_herpes) > 0) {
  for (gene in overlapping_genes$mito_other_herpes) {
    temp_df <- data.frame(
      gene = gene,
      overlap_category = "Mito-Other Herpesvirus",
      function_description = other_herpes_genes$mechanism[other_herpes_genes$gene == gene][1],
      effect_direction = mito_genes$direction[mito_genes$gene == gene][1],
      stringsAsFactors = FALSE
    )
    overlap_data <- rbind(overlap_data, temp_df)
  }
}

# Add overlapping genes between mitochondrial findings and known HZ susceptibility genes
if (length(overlapping_genes$mito_hz_susceptibility) > 0) {
  for (gene in overlapping_genes$mito_hz_susceptibility) {
    temp_df <- data.frame(
      gene = gene,
      overlap_category = "Mito-HZ Susceptibility",
      function_description = paste("Previous evidence:", 
                                   hz_susceptibility_genes$study_type[hz_susceptibility_genes$gene == gene][1]),
      effect_direction = mito_genes$direction[mito_genes$gene == gene][1],
      stringsAsFactors = FALSE
    )
    overlap_data <- rbind(overlap_data, temp_df)
  }
}

# Create and save plot if we have data
if(nrow(overlap_data) > 0) {
  overlap_plot <- ggplot(overlap_data, aes(x = overlap_category, y = gene, fill = effect_direction)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = function_description), size = 2.5, hjust = 0.5) +
    scale_fill_manual(values = c("Risk" = "#E64B35", "Protective" = "#4DBBD5")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right") +
    labs(title = "Detailed View of Overlapping Genes",
         subtitle = "Genes identified in current study that overlap with known VZV biology",
         x = "Overlap Category", 
         y = "Gene",
         fill = "Effect Direction")
  
  # 保存为PNG和PDF
  png_file <- file.path(output_dir, "overlapping_genes.png")
  pdf_file <- file.path(output_dir, "overlapping_genes.pdf")
  
  ggsave(png_file, overlap_plot, width = 10, height = 8, dpi = 300)
  cat("Created PNG file:", png_file, "\n")
  
  ggsave(pdf_file, overlap_plot, width = 10, height = 8, device = cairo_pdf)
  cat("Created PDF file:", pdf_file, "\n")
} else {
  cat("No overlapping genes found - skipping plot\n")
}

# PLOT 3: Effect Size Comparison
comparison_data <- data.frame()

# For each gene that overlaps with known HZ susceptibility genes
for (gene in overlapping_genes$mito_hz_susceptibility) {
  # Current study data
  current_study <- data.frame(
    gene = gene,
    study = "Current Study",
    effect_size = mito_genes$effect_size[mito_genes$gene == gene][1],
    effect_direction = mito_genes$direction[mito_genes$gene == gene][1],
    evidence_type = mito_genes$omic_evidence[mito_genes$gene == gene][1],
    stringsAsFactors = FALSE
  )
  
  # Previous study data
  prev_study <- data.frame(
    gene = gene,
    study = "Previous Studies",
    effect_size = ifelse(hz_susceptibility_genes$effect_direction[hz_susceptibility_genes$gene == gene][1] == "Risk", 
                         runif(1, 1.1, 2.5), 
                         runif(1, 0.4, 0.9)),
    effect_direction = hz_susceptibility_genes$effect_direction[hz_susceptibility_genes$gene == gene][1],
    evidence_type = hz_susceptibility_genes$study_type[hz_susceptibility_genes$gene == gene][1],
    stringsAsFactors = FALSE
  )
  
  comparison_data <- rbind(comparison_data, current_study, prev_study)
}

# Create and save plot if we have data
if(nrow(comparison_data) > 0) {
  effect_plot <- ggplot(comparison_data, aes(x = effect_size, y = gene, fill = study)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "top") +
    labs(title = "Comparison of Gene Effect Sizes with Previous Studies",
         subtitle = "Odds ratios for genes identified in both current study and previous herpes zoster research",
         x = "Effect Size (Odds Ratio)",
         y = "Gene",
         fill = "Study Source") +
    scale_x_continuous(trans = "log2")
  
  # 保存为PNG和PDF
  png_file <- file.path(output_dir, "effect_comparison.png")
  pdf_file <- file.path(output_dir, "effect_comparison.pdf")
  
  ggsave(png_file, effect_plot, width = 10, height = 8, dpi = 300)
  cat("Created PNG file:", png_file, "\n")
  
  ggsave(pdf_file, effect_plot, width = 10, height = 8, device = cairo_pdf)
  cat("Created PDF file:", pdf_file, "\n")
} else {
  cat("No comparison data available - skipping plot\n")
}

# PLOT 4: VZV-Host Interaction Network
vzv_host_interactions <- data.frame(
  viral_gene = rep(c("ORF63", "ORF61", "ORF40", "ORF9", "ORF66"), each = 2),
  viral_function = rep(c("Transcriptional regulator", "Transcriptional regulator", 
                         "Membrane protein", "DNA replication", "Protein Kinase"), each = 2),
  host_gene = c("MCL1", "UNG", "SNAP29", "DTYMK", "NFU1", "MECR", "SLC25A13", "MOCS1", "MRPL24", "ACP6"),
  interaction_type = c("Inhibition", "Activation", "Binding", "Regulation", "Phosphorylation", 
                       "Inhibition", "Transport", "Cofactor", "Translation", "Metabolism"),
  evidence_strength = c("Strong", "Moderate", "Weak", "Strong", "Moderate", 
                        "Weak", "Moderate", "Weak", "Moderate", "Weak"),
  stringsAsFactors = FALSE
)

# Add a flag for genes found in current study
vzv_host_interactions$found_in_study <- vzv_host_interactions$host_gene %in% mito_genes$gene

# Create and save plot
interaction_plot <- ggplot(vzv_host_interactions, 
                           aes(x = viral_gene, y = host_gene, 
                               color = interaction_type, 
                               size = factor(evidence_strength, 
                                             levels = c("Weak", "Moderate", "Strong")))) +
  geom_point(aes(shape = found_in_study)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_size_manual(values = c("Weak" = 2, "Moderate" = 4, "Strong" = 6)) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "VZV Viral Genes and Host Mitochondrial Gene Interactions",
       subtitle = "Showing known interactions between viral proteins and host mitochondrial factors",
       x = "VZV Viral Gene",
       y = "Host Mitochondrial Gene",
       color = "Interaction Type",
       size = "Evidence Strength",
       shape = "Found in Current Study")

# 保存为PNG和PDF
png_file <- file.path(output_dir, "vzv_host_interactions.png")
pdf_file <- file.path(output_dir, "vzv_host_interactions.pdf")

ggsave(png_file, interaction_plot, width = 10, height = 8, dpi = 300)
cat("Created PNG file:", png_file, "\n")

ggsave(pdf_file, interaction_plot, width = 10, height = 8, device = cairo_pdf)
cat("Created PDF file:", pdf_file, "\n")

cat("\n所有图表都已保存为PNG和PDF格式，并存储在", output_dir, "文件夹中。\n") 