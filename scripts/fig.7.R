# Load necessary libraries
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(dplyr)

# 创建输出目录（如果不存在）
output_dir <- "../results/output_figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# Set seed for reproducibility
set.seed(123)

# Create a sample dataset based on your findings
# Define the genes we want to include
genes <- c(
  "MRPL24", "SLC25A13", "MCL1", "MOCS1", "DTYMK", 
  "NFU1", "MECR", "UNG", "CPLX1", "HLA-B", 
  "HORMAD1", "SETDB1", "MCCD1", "RAB2A", "SNAP29"
)

# Define gene functions - using a more organized categorization
gene_functions <- c(
  "Ribosome", "Transport", "Apoptosis", "Cofactor", "Nucleotide",
  "Fe-S", "Lipid", "DNA repair", "Exocytosis", "Immunity",
  "Cell cycle", "Epigenetic", "Metabolism", "Trafficking", "Membrane fusion"
)

# Create a data frame with gene names and functions
gene_info <- data.frame(
  Gene = genes,
  Function = gene_functions,
  stringsAsFactors = FALSE
)

# Create data structures for each omics level
# For methylation data (values represent -log10(p-value) * sign(effect))
methylation_data <- matrix(NA, nrow = length(genes), ncol = 3)
colnames(methylation_data) <- c("CpG1", "CpG2", "CpG3")
rownames(methylation_data) <- genes

# For expression data
expression_data <- matrix(NA, nrow = length(genes), ncol = 3)
colnames(expression_data) <- c("Blood", "Skin_NE", "Skin_SE")
rownames(expression_data) <- genes

# For protein data
protein_data <- matrix(NA, nrow = length(genes), ncol = 1)
colnames(protein_data) <- c("Abundance")
rownames(protein_data) <- genes

# For colocalization evidence
pph4_data <- matrix(NA, nrow = length(genes), ncol = 1)
colnames(pph4_data) <- c("PPH4")
rownames(pph4_data) <- genes

# Fill in data based on your results
# Methylation data
# Positive values = risk increasing, Negative values = protective
methylation_data["MRPL24", "CpG1"] <- 2.1   # Significant, risk increasing
methylation_data["MRPL24", "CpG2"] <- 1.4   # Non-significant
methylation_data["SLC25A13", "CpG1"] <- 1.8  # Significant, risk increasing
methylation_data["MCL1", "CpG1"] <- -2.3    # Significant, protective
methylation_data["MCL1", "CpG2"] <- -1.9    # Significant, protective
methylation_data["MOCS1", "CpG1"] <- -2.1   # Significant, protective
methylation_data["MCCD1", "CpG1"] <- -3.2   # Significant, protective

# Expression data
expression_data["MRPL24", "Skin_NE"] <- -0.9  # Non-significant, protective
expression_data["MRPL24", "Skin_SE"] <- -1.2  # Non-significant, protective
expression_data["SLC25A13", "Skin_NE"] <- -0.8  # Non-significant, protective
expression_data["CPLX1", "Skin_NE"] <- -2.1   # Significant, protective
expression_data["CPLX1", "Skin_SE"] <- -1.8   # Significant, protective
expression_data["HLA-B", "Skin_NE"] <- 1.3    # Significant, risk increasing
expression_data["HLA-B", "Skin_SE"] <- 1.2    # Significant, risk increasing
expression_data["HORMAD1", "Skin_NE"] <- -1.5  # Significant, protective
expression_data["HORMAD1", "Skin_SE"] <- -1.4  # Significant, protective
expression_data["SETDB1", "Skin_NE"] <- 1.4    # Significant, risk increasing
expression_data["RAB2A", "Skin_SE"] <- 1.6     # Significant, risk increasing

# Protein data
protein_data["DTYMK", "Abundance"] <- 1.0    # Significant, risk increasing
protein_data["NFU1", "Abundance"] <- 0.5     # Significant, risk increasing
protein_data["MECR", "Abundance"] <- -0.7    # Significant, protective
protein_data["UNG", "Abundance"] <- -0.8     # Significant, protective
protein_data["SNAP29", "Abundance"] <- 0.6   # Significant, risk increasing

# Colocalization evidence (PPH4 values)
pph4_data["MRPL24", "PPH4"] <- 0.62
pph4_data["SLC25A13", "PPH4"] <- 0.59
pph4_data["MCL1", "PPH4"] <- 0.23
pph4_data["MOCS1", "PPH4"] <- 0.58
pph4_data["CPLX1", "PPH4"] <- 0.71
pph4_data["HLA-B", "PPH4"] <- 0.65
pph4_data["HORMAD1", "PPH4"] <- 0.54
pph4_data["MCCD1", "PPH4"] <- 0.48

# Create significance markers for the visualization
# In a real implementation, you would base this on your actual significance thresholds
methylation_sig <- matrix(FALSE, nrow = length(genes), ncol = 3)
colnames(methylation_sig) <- colnames(methylation_data)
rownames(methylation_sig) <- genes
methylation_sig["MRPL24", "CpG1"] <- TRUE
methylation_sig["SLC25A13", "CpG1"] <- TRUE
methylation_sig["MCL1", "CpG1"] <- TRUE
methylation_sig["MCL1", "CpG2"] <- TRUE
methylation_sig["MOCS1", "CpG1"] <- TRUE
methylation_sig["MCCD1", "CpG1"] <- TRUE

expression_sig <- matrix(FALSE, nrow = length(genes), ncol = 3)
colnames(expression_sig) <- colnames(expression_data)
rownames(expression_sig) <- genes
expression_sig["CPLX1", "Skin_NE"] <- TRUE
expression_sig["CPLX1", "Skin_SE"] <- TRUE
expression_sig["HLA-B", "Skin_NE"] <- TRUE
expression_sig["HLA-B", "Skin_SE"] <- TRUE
expression_sig["HORMAD1", "Skin_NE"] <- TRUE
expression_sig["HORMAD1", "Skin_SE"] <- TRUE
expression_sig["SETDB1", "Skin_NE"] <- TRUE
expression_sig["RAB2A", "Skin_SE"] <- TRUE

protein_sig <- matrix(FALSE, nrow = length(genes), ncol = 1)
colnames(protein_sig) <- colnames(protein_data)
rownames(protein_sig) <- genes
protein_sig["DTYMK", "Abundance"] <- TRUE
protein_sig["NFU1", "Abundance"] <- TRUE
protein_sig["MECR", "Abundance"] <- TRUE
protein_sig["UNG", "Abundance"] <- TRUE
protein_sig["SNAP29", "Abundance"] <- TRUE

# Combine data for visualization
combined_data <- cbind(expression_data, methylation_data, protein_data)
combined_sig <- cbind(expression_sig, methylation_sig, protein_sig)

# Calculate evidence count for each gene (for potential sorting or annotation)
evidence_count <- rowSums(!is.na(combined_data) & abs(combined_data) > 0, na.rm = TRUE)

# Create a multi-omics status indicator 
# This shows at how many omics levels each gene has evidence
multi_omics_status <- rep(0, length(genes))
names(multi_omics_status) <- genes

# Check expression evidence
has_expression_evidence <- rowSums(!is.na(expression_data) & abs(expression_data) > 0, na.rm = TRUE) > 0
# Check methylation evidence
has_methylation_evidence <- rowSums(!is.na(methylation_data) & abs(methylation_data) > 0, na.rm = TRUE) > 0
# Check protein evidence
has_protein_evidence <- rowSums(!is.na(protein_data) & abs(protein_data) > 0, na.rm = TRUE) > 0

# Sum up the number of omics levels with evidence
for (i in 1:length(genes)) {
  multi_omics_status[i] <- sum(has_expression_evidence[i], has_methylation_evidence[i], has_protein_evidence[i])
}

# Define improved color mapping function for the heatmap
# Use a more perceptually uniform and colorblind-friendly palette
col_fun <- colorRamp2(
  c(-4, -2, 0, 2, 4), 
  c("#053061", "#4393c3", "#f7f7f7", "#d6604d", "#67001f")  # ColorBrewer RdBu palette
)

# Better column names for improved readability
better_column_names <- c(
  "Blood", "Non-Sun\nExposed\nSkin", "Sun\nExposed\nSkin",  # Expression
  "CpG1", "CpG2", "CpG3",  # Methylation
  "Protein\nAbundance"  # Protein
)
names(better_column_names) <- colnames(combined_data)

# Create annotation for evidence type with improved spacing and clarity
evidence_type <- c(
  rep("Expression", ncol(expression_data)), 
  rep("Methylation", ncol(methylation_data)), 
  rep("Protein", ncol(protein_data))
)

# Define a more harmonious color palette for evidence types
evidence_colors <- c(
  "Expression" = "#4575b4",  # Blue
  "Methylation" = "#d73027", # Red
  "Protein" = "#4daf4a"      # Green
)

# Create improved column annotation with better spacing
column_anno <- HeatmapAnnotation(
  Evidence = evidence_type,
  col = list(Evidence = evidence_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  gap = unit(2, "mm"),
  height = unit(5, "mm")
)

# Create a more visually appealing functional annotation
# Using a color-blind friendly palette for gene functions
n_functions <- length(unique(gene_info$Function))
function_colors <- c(
  "Ribosome" = "#E69F00", 
  "Transport" = "#56B4E9", 
  "Apoptosis" = "#009E73", 
  "Cofactor" = "#F0E442", 
  "Nucleotide" = "#0072B2",
  "Fe-S" = "#D55E00", 
  "Lipid" = "#CC79A7", 
  "DNA repair" = "#999999", 
  "Exocytosis" = "#E69F00", 
  "Immunity" = "#56B4E9",
  "Cell cycle" = "#009E73", 
  "Epigenetic" = "#F0E442", 
  "Metabolism" = "#0072B2", 
  "Trafficking" = "#D55E00", 
  "Membrane fusion" = "#CC79A7"
)

function_anno <- rowAnnotation(
  Function = gene_info$Function,
  col = list(Function = function_colors),
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    Function = list(
      title = "Functional Category",
      at = unique(gene_info$Function),
      labels = unique(gene_info$Function),
      nrow = ceiling(n_functions/3)
    )
  ),
  width = unit(5, "mm"),
  annotation_name_side = "bottom",
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold")
)

# Add an annotation for multi-omics evidence
multi_omics_anno <- rowAnnotation(
  "Multi-omics\nEvidence" = anno_barplot(
    multi_omics_status, 
    bar_width = 0.8,
    gp = gpar(fill = c("1" = "#fee8c8", "2" = "#fdbb84", "3" = "#e34a33")[as.character(multi_omics_status)]),
    axis = FALSE,
    width = unit(10, "mm")
  ),
  annotation_name_side = "bottom",
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 10)
)

# Create the colocalization evidence heatmap with improved visualization
pph4_heatmap <- Heatmap(
  pph4_data,
  name = "PPH4",
  col = colorRamp2(c(0, 0.5, 0.7, 1), c("#FFFFFF", "#FFFFAA", "#FFAA00", "#FF0000")),
  show_row_names = FALSE,
  width = unit(15, "mm"),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  
  # Add cell values and horizontal guides to show threshold
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(pph4_data[i, j])) {
      # Add value text
      grid.text(sprintf("%.2f", pph4_data[i, j]), x, y, gp = gpar(fontsize = 8))
      
      # Add a subtle horizontal line at PPH4 = 0.5 threshold
      if (pph4_data[i, j] >= 0.5 && pph4_data[i, j] < 0.51) {
        grid.lines(
          c(x - width/2, x + width/2), 
          c(y, y), 
          gp = gpar(col = "#00000044", lty = 2)
        )
      }
    }
  },
  
  # Add a border to separate from main heatmap
  border = TRUE,
  
  # Add better legends
  heatmap_legend_param = list(
    title = "Colocalization\nEvidence (PPH4)",
    at = c(0, 0.5, 0.7, 1),
    labels = c("0", "0.5", "0.7", "1"),
    legend_height = unit(50, "mm"),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

# Group genes by evidence pattern for better visualization
# This helps identify genes with consistent effects across omics levels
# Prioritize genes with multi-omics evidence

# Define a custom order based on evidence patterns and strength
# This creates a more storytelling arrangement of genes
gene_order <- order(
  -multi_omics_status,  # First by number of omics levels
  -rowSums(combined_sig, na.rm = TRUE),  # Then by total significance
  decreasing = TRUE  # Descending order
)

# Create the main heatmap with improved visualization
main_heatmap <- Heatmap(
  combined_data,
  name = "Effect",
  col = col_fun,
  na_col = "#EEEEEE",  # Light gray for missing values (distinct from white = no effect)
  row_names_side = "left",
  column_names_side = "top",
  
  # Use improved column names
  column_labels = better_column_names,
  
  # Add the enhanced annotations
  top_annotation = column_anno,
  left_annotation = function_anno,
  right_annotation = multi_omics_anno,
  
  # Improve text elements
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  
  # Add asterisks for significant findings with improved visibility
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (combined_sig[i, j]) {
      grid.text("*", x, y, gp = gpar(fontsize = 14, fontface = "bold", col = "black"))
    }
  },
  
  # Improved grouping and organization
  column_split = evidence_type,
  cluster_columns = FALSE,
  
  # Use custom row order instead of clustering
  cluster_rows = FALSE,
  row_order = gene_order,
  
  # Add borders between evidence types for clarity
  column_gap = unit(5, "mm"),
  border = TRUE,
  
  # Add better legends
  heatmap_legend_param = list(
    title = "Effect",
    at = c(-4, -2, 0, 2, 4),
    labels = c("-4", "-2", "0", "2", "4"),
    legend_height = unit(50, "mm"),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

# Combine the heatmaps with improved spacing
combined_heatmap <- main_heatmap + pph4_heatmap

# Set the plot title with improved positioning and style
title <- "Multi-Omic View of Mitochondrial Gene Associations with Herpes Zoster Risk"
subtitle <- "Mendelian Randomization Analysis of FinnGen R12 Data"

# 定义输出文件路径
png_file <- file.path(output_dir, "Multi_Omic_Heatmap.png")
pdf_file <- file.path(output_dir, "Multi_Omic_Heatmap.pdf")

# 首先保存为PNG格式
png(png_file, width = 10, height = 8, units = "in", res = 300)

# 绘制热图
draw(combined_heatmap, 
     column_title = title,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     column_title_side = "top",
     padding = unit(c(1, 1, 2, 1), "cm"),  # top, right, bottom, left padding
     heatmap_legend_side = "right",
     legend_grouping = "original")

# 添加副标题
grid.text(subtitle, 
          x = unit(0.5, "npc"), 
          y = unit(0.97, "npc"),
          gp = gpar(fontsize = 12, fontface = "italic"))

# 添加解释说明
grid.text("Values represent -log10(p-value) × sign(effect): Red = increased risk, Blue = protective", 
          x = unit(0.5, "npc"), 
          y = unit(0.02, "npc"),
          gp = gpar(fontsize = 9))

# 添加显著性指示符说明
grid.text("* FDR-significant finding", 
          x = unit(0.85, "npc"), 
          y = unit(0.04, "npc"),
          gp = gpar(fontsize = 10, fontface = "bold"))

# 关闭PNG设备
dev.off()
cat("创建PNG文件:", png_file, "\n")

# 然后保存为PDF格式
pdf(pdf_file, width = 10, height = 8)

# 再次绘制相同的热图
draw(combined_heatmap, 
     column_title = title,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     column_title_side = "top",
     padding = unit(c(1, 1, 2, 1), "cm"),  # top, right, bottom, left padding
     heatmap_legend_side = "right",
     legend_grouping = "original")

# 添加副标题
grid.text(subtitle, 
          x = unit(0.5, "npc"), 
          y = unit(0.97, "npc"),
          gp = gpar(fontsize = 12, fontface = "italic"))

# 添加解释说明
grid.text("Values represent -log10(p-value) × sign(effect): Red = increased risk, Blue = protective", 
          x = unit(0.5, "npc"), 
          y = unit(0.02, "npc"),
          gp = gpar(fontsize = 9))

# 添加显著性指示符说明
grid.text("* FDR-significant finding", 
          x = unit(0.85, "npc"), 
          y = unit(0.04, "npc"),
          gp = gpar(fontsize = 10, fontface = "bold"))

# 关闭PDF设备
dev.off()
cat("创建PDF文件:", pdf_file, "\n")

# 通知用户处理完成
cat("\n热图已成功保存为PNG和PDF格式，存储在", output_dir, "文件夹中。\n") 