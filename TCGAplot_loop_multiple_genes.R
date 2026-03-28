# Load necessary packages
library(TCGAplot)
library(beepr)
library(parallel)
library(doParallel)
library(foreach)
library(progress)

# ====== USER SETTINGS ======
genes <- c("gene1", "gene2", "gene3") #Repleace with official gene symbols

# ====== PARALLEL SETUP ======
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent | Gene :current/:total | ETA: :eta",
  total = length(genes), clear = FALSE, width = 60
)

# Log file
log_file <- "failed_genes.log"
write("", file = log_file)

# ====== PARALLEL LOOP ======
results <- foreach(gene = genes,
                   .packages = c("TCGAplot")) %dopar% {
                     
                     tryCatch({
                       
                       pb$tick()
                       
                       # Output directory
                       output_dir <- paste0("TCGA_plots_", gene)
                       if (!dir.exists(output_dir))
                         dir.create(output_dir, recursive = TRUE)
                       
                       # Helper function: only run plot if the file does not yet exist
                       save_plot <- function(filename, width_in, height_in, plot_expr) {
                         filepath <- file.path(output_dir, filename)
                         
                         # Skip if plot already exists
                         if (file.exists(filepath)) {
                           return("SKIPPED")
                         }
                         
                         png(filename = filepath, width = width_in, height = height_in,
                             units = "in", res = 300)
                         print(plot_expr)
                         dev.off()
                         return("CREATED")
                       }
                       
                       # ==== ALL PLOT FUNCTIONS WITH AUTO-SKIP ====
                       
                       save_plot("pan_boxplot.png", 10, 5.0, {
                         pan_boxplot(gene = gene, palette = "jco", legend = "right")
                       })
                       
                       save_plot("paired_boxplot.png", 10.5, 5.0, {
                         pan_paired_boxplot(gene, palette = "jco", legend = "right",
                                            method = "wilcox.test")
                       })
                       
                       save_plot("tumor_boxplot.png", 10, 5.0, {
                         pan_tumor_boxplot(gene)
                       })
                       
                       save_plot("TMB_radar.png", 8, 8, {
                         gene_TMB_radar(gene, method = "pearson")
                       })
                       
                       save_plot("MSI_radar.png", 8, 8, {
                         gene_MSI_radar(gene)
                       })
                       
                       save_plot("checkpoint_heatmap.png", 8, 3.5, {
                         gene_checkpoint_heatmap(gene)
                       })
                       
                       save_plot("chemokine_heatmap.png", 8, 8, {
                         gene_chemokine_heatmap(gene, method = "pearson")
                       })
                       
                       save_plot("receptor_heatmap.png", 8.2, 4.5, {
                         gene_receptor_heatmap(gene)
                       })
                       
                       save_plot("immustimulator_heatmap.png", 8, 9, {
                         gene_immustimulator_heatmap(gene)
                       })
                       
                       save_plot("immuinhibitor_heatmap.png", 8, 5.8, {
                         gene_immuinhibitor_heatmap(gene)
                       })
                       
                       save_plot("immucell_heatmap.png", 9, 5.2, {
                         gene_immucell_heatmap(gene, method = "pearson")
                       })
                       
                       save_plot("immunescore_heatmap.png", 8, 2.5, {
                         gene_immunescore_heatmap(gene, method = "pearson")
                       })
                       
                       save_plot("immunescore_triangle.png", 9, 4, {
                         gene_immunescore_triangle(gene, method = "pearson")
                       })
                       
                       save_plot("forest_plot.png", 8.2, 8, {
                         pan_forest(gene)
                       })
                       
                       return(list(gene = gene, status = "OK"))
                       
                     }, error = function(e) {
                       
                       msg <- paste0("FAILED: ", gene, " | ", e$message, "\n")
                       cat(msg, file = log_file, append = TRUE)
                       message(msg)
                       
                       return(list(gene = gene, status = "ERROR"))
                     })
                   }

# Stop workers
stopCluster(cl)

beep(3)
message("🎉 ALL GENES PROCESSED WITH PARALLEL EXECUTION!")