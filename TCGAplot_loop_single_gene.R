# Load necessary package
library(TCGAplot)
library(beepr)

# ====== USER SETTINGS ======
gene <- "ATAD2"  # Change this to your gene of interest
output_dir <- "TCGA_plots_ATAD2"  # Output folder

# ====== SETUP ======
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Helper function to save plot using png() with inches and 300 DPI
save_plot <- function(filename, width_in, height_in, plot_expr) {
  filepath <- file.path(output_dir, filename)
  png(filename = filepath, width = width_in, height = height_in, units = "in", res = 300)
  print(plot_expr)
  dev.off()
}

# ====== PLOT FUNCTIONS ======

# 1. Pan-cancer boxplot
save_plot("pan_boxplot.png", 10, 5.0, {
  pan_boxplot(gene = gene, palette = "jco", legend = "right")
})

# 2. Paired tumor-normal boxplot
save_plot("paired_boxplot.png", 10.5, 5.0, {
  pan_paired_boxplot(gene, palette = "jco", legend = "right", method = "wilcox.test")
})

# 3. Tumor-only boxplot
save_plot("tumor_boxplot.png", 10, 5.0, {
  pan_tumor_boxplot(gene)
})

# 4. TMB radar
save_plot("TMB_radar.png", 8, 8, {
  gene_TMB_radar(gene, method = "pearson")
})

# 5. MSI radar
save_plot("MSI_radar.png", 8, 8, {
  gene_MSI_radar(gene)
})

# 6. Immune checkpoint heatmap
save_plot("checkpoint_heatmap.png", 8, 3.5, {
  gene_checkpoint_heatmap(gene)
})

# 7. Chemokine heatmap
save_plot("chemokine_heatmap.png", 8, 8, {
  gene_chemokine_heatmap(gene, method = "pearson")
})

# 8. Receptor heatmap
save_plot("receptor_heatmap.png", 8.2, 4.5, {
  gene_receptor_heatmap(gene)
})

# 9. Immune stimulator heatmap
save_plot("immustimulator_heatmap.png", 8, 9, {
  gene_immustimulator_heatmap(gene)
})

# 10. Immune inhibitor heatmap
save_plot("immuinhibitor_heatmap.png", 8, 5.8, {
  gene_immuinhibitor_heatmap(gene)
})

# 11. Immune cell heatmap
save_plot("immucell_heatmap.png", 9, 5.2, {
  gene_immucell_heatmap(gene, method = "pearson")
})

# 12. Immune score heatmap
save_plot("immunescore_heatmap.png", 8, 2.5, {
  gene_immunescore_heatmap(gene, method = "pearson")
})

# 13. Immune score triangle
save_plot("immunescore_triangle.png", 9, 4, {
  gene_immunescore_triangle(gene, method = "pearson")
})

# 14. Forest plot
save_plot("forest_plot.png", 8.2, 8, {
  pan_forest(gene)
})

message("✅ All plots saved to: ", output_dir)
