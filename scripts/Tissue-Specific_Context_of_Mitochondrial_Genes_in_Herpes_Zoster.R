# 强烈建议在项目根目录（Herpes_Zoster_Public_Project）下通过Rscript命令运行
# 或者在R/RStudio中打开项目，并确保工作目录是项目的根目录
# setwd("C:\\Users\\hao55\\Desktop\\MR\\Herpes_zoster")

# 创建输出目录
output_dir <- "../results/output_figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 加载必要的库
library(ggplot2)
library(patchwork)  # 用于组合图表
library(dplyr)

# 第一部分: 创建解剖学背景图
#-----------------------------------------
tissue_data <- data.frame(
  x = c(1, 2, 3),
  y = c(1, 1, 1),
  tissue = c("Dorsal Root\nGanglia", "Peripheral\nNerve", "Skin"),
  role = c("Viral Latency", "Viral Transport", "Rash Formation"),
  genes = c("MCL1, UNG", "MECR, MOCS1", "DTYMK")
)

anatomical_plot <- ggplot(tissue_data, aes(x = x, y = y)) +
  # 创建组织节点
  geom_point(size = 20, color = "skyblue", alpha = 0.4) +
  # 添加组织标签
  geom_text(aes(label = tissue), fontface = "bold", size = 4.5) +
  # 添加功能描述
  geom_text(aes(label = role, y = y - 0.2), size = 3.5) +
  # 添加基因注释
  geom_text(aes(label = genes, y = y - 0.4), size = 3.2, fontface = "italic", color = "darkblue") +
  # 添加连接箭头
  geom_segment(aes(x = 1.2, xend = 1.8, y = 1, yend = 1), 
               arrow = arrow(length = unit(0.2, "cm")), size = 0.7) +
  geom_segment(aes(x = 2.2, xend = 2.8, y = 1, yend = 1), 
               arrow = arrow(length = unit(0.2, "cm")), size = 0.7) +
  # 移除背景元素
  theme_void() +
  # 添加解释性标题
  labs(title = "Herpes Zoster Progression Through Tissues",
       subtitle = "Key mitochondrial genes associated with each stage") +
  # 设置图表边界
  xlim(0.5, 3.5) +
  ylim(0.3, 1.5)

# 第二部分: 使用ggplot2创建热图
#-----------------------------------------
# 创建热图的示例数据 - 替换为你的实际表达数据
expression_data <- data.frame(
  gene = rep(c("DTYMK", "MCL1", "MECR", "UNG", "MOCS1"), each = 2),
  tissue = rep(c("Sun-Exposed Skin", "Non-Sun-Exposed Skin"), 5),
  expression = c(
    1.8, 1.5,  # DTYMK
    1.2, 0.9,  # MCL1
    -0.5, -0.8,  # MECR (表示下调的负值)
    -0.7, -1.0,  # UNG
    -0.3, -0.6   # MOCS1
  ),
  significant = c(
    TRUE, TRUE,    # DTYMK
    TRUE, FALSE,   # MCL1
    TRUE, TRUE,    # MECR
    FALSE, TRUE,   # UNG
    TRUE, FALSE    # MOCS1
  )
)

# 创建自定义ggplot热图
expression_heatmap <- ggplot(expression_data, 
                             aes(x = tissue, y = gene, fill = expression)) +
  # 创建瓦片
  geom_tile(color = "white", size = 0.5) +
  # 添加显著性标记（星号）
  geom_text(aes(label = ifelse(significant, "*", "")), 
            color = "black", size = 6) +
  # 使用发散色标（红色表示上调，蓝色表示下调）
  scale_fill_gradient2(
    low = "blue3", 
    mid = "white", 
    high = "red3", 
    midpoint = 0,
    name = "Log2 Fold\nChange"
  ) +
  # 改善外观
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Mitochondrial Gene Expression in Skin Tissues",
    x = "",
    y = ""
  )

# 第三部分: 创建MR结果的森林图
#-----------------------------------------
# 你的MR结果 - 替换为你的实际数据
mr_results <- data.frame(
  gene = c("DTYMK", "SNAP29", "NFU1", "ACP6", "MECR", "UNG"),
  odds_ratio = c(2.78, 1.90, 1.61, 0.95, 0.52, 0.44),
  lower_ci = c(1.35, 1.10, 1.11, 0.91, 0.30, 0.20),
  upper_ci = c(5.70, 3.28, 2.34, 0.98, 0.91, 0.97),
  p_value = c(0.005, 0.021, 0.013, 0.004, 0.023, 0.041)
)

mr_plot <- ggplot(mr_results, 
                  aes(x = odds_ratio, y = reorder(gene, -odds_ratio), 
                      xmin = lower_ci, xmax = upper_ci,
                      color = odds_ratio < 1)) +
  # 在OR = 1处添加参考线
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
  # 添加误差线和点
  geom_errorbarh(height = 0.2, size = 0.7) +
  geom_point(size = 3.5) +
  # 添加p值注释
  geom_text(aes(label = paste0("p = ", p_value), x = ifelse(odds_ratio > 1, upper_ci * 1.1, lower_ci * 0.9)), 
            hjust = ifelse(mr_results$odds_ratio > 1, 0, 1),
            size = 3.5) +
  # 设置风险增加与保护效应的颜色
  scale_color_manual(values = c("firebrick", "steelblue"), 
                     labels = c("Risk Increasing", "Protective"),
                     name = "Effect Direction") +
  # 使用对数比例尺表示比值比
  scale_x_continuous(trans = "log10") +
  # 改善外观
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "Mitochondrial Gene Associations with Herpes Zoster Risk",
    subtitle = "Mendelian Randomization Results from FinnGen R12",
    x = "Odds Ratio (95% CI, log scale)", 
    y = ""
  )

# 第四部分: 创建通路富集摘要
#-----------------------------------------
# 示例通路数据 - 替换为你的实际通路分析结果
pathway_data <- data.frame(
  pathway = c(
    "Oxidative stress response",
    "Neuronal apoptotic signaling",
    "Mitochondrial fusion",
    "Cytochrome c release",
    "ATP biosynthesis"
  ),
  negative_log_p = c(3.5, 3.2, 2.9, 2.7, 2.4),
  gene_count = c(12, 8, 6, 7, 9)
)

pathway_plot <- ggplot(pathway_data, 
                       aes(x = reorder(pathway, negative_log_p), 
                           y = negative_log_p,
                           size = gene_count)) +
  # 添加点
  geom_point(color = "purple4", alpha = 0.8) +
  # 添加显著性阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  # 改善外观
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "Mitochondrial Pathways Associated with Herpes Zoster",
    x = "",
    y = "-log10(p-value)",
    size = "Gene Count"
  )

# 第五部分: 使用patchwork组合所有图表
#-----------------------------------------
combined_figure <- (anatomical_plot) / 
  (expression_heatmap | pathway_plot) / 
  (mr_plot) +
  plot_layout(heights = c(1, 1.5, 1.5)) +
  plot_annotation(
    title = "Tissue-Specific Context of Mitochondrial Genes in Herpes Zoster",
    subtitle = "Integrating findings from Mendelian Randomization analysis of FinnGen R12",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# 显示组合图
print(combined_figure)

# 定义PNG和PDF输出文件路径
png_file <- file.path(output_dir, "herpes_zoster_mitochondrial_genes.png")
pdf_file <- file.path(output_dir, "herpes_zoster_mitochondrial_genes.pdf")

# 保存PNG版本
ggsave(png_file, combined_figure, width = 12, height = 14, dpi = 300)
cat("创建PNG文件:", png_file, "\n")

# 保存PDF版本
ggsave(pdf_file, combined_figure, width = 12, height = 14)
cat("创建PDF文件:", pdf_file, "\n")

# 保存单独的组件图表（可选）
individual_plots <- list(
  anatomical = anatomical_plot,
  heatmap = expression_heatmap,
  pathway = pathway_plot,
  mr = mr_plot
)

# 如果需要单独的组件图表，取消下面的注释
for (name in names(individual_plots)) {
  png_component <- file.path(output_dir, paste0("component_", name, ".png"))
  pdf_component <- file.path(output_dir, paste0("component_", name, ".pdf"))
  
  ggsave(png_component, individual_plots[[name]], width = 8, height = 6, dpi = 300)
  ggsave(pdf_component, individual_plots[[name]], width = 8, height = 6)
  
  cat("创建组件图表:", name, "(PNG和PDF)\n")
}

cat("\n所有图表已成功保存为PNG和PDF格式，存储在", output_dir, "文件夹中。\n") 