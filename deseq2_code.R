library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(data.table)

# -----------------------------------
# Your cancers
# -----------------------------------

cancers <- "all"   # OR c("COAD","READ","STAD","ESCA")

data_dir <- "Z:/DATA/Guru/datasets/TCGA_raw_counts/"

if (length(cancers) == 1 && cancers == "all") {
  
  files <- list.files(
    path = data_dir,
    pattern = "^raw_counts_.*\\.csv$",
    full.names = FALSE
  )
  
  # Remove merged file if present
  files <- files[!grepl("merged_rawcounts", files)]
  
  # Extract cancer codes
  cancers <- gsub("raw_counts_|\\.csv", "", files)
}

# -----------------------------------
# LOOP THROUGH EACH CANCER
# -----------------------------------

for (cancer in cancers)
{
  cat("Processing:", cancer, "(", which(cancers==cancer), "/", length(cancers), ")\n")
  
  # -----------------------------------
  # 1. LOAD CSV
  # -----------------------------------
  
  file_name <- paste0(data_dir, "raw_counts_", cancer, ".csv")
  df <- as.data.frame(fread(file_name))
  
  # -----------------------------------
  # 2. REMOVE gene_type COLUMN
  # -----------------------------------
  
  df <- df[df$gene_type == "protein_coding",]
  
  df <- df[, -2]   # removes 2nd column (gene_type)
  
  rownames(df) <- df$Gene
  
  # Remove Gene column before DESeq
  df$Gene <- NULL
  
  # -----------------------------------
  # 4. CREATE METADATA FROM TCGA BARCODES
  # -----------------------------------
  
  sample_names <- colnames(df)
  
  sample_type_code <- substr(sample_names, 14, 15)
  
  condition <- ifelse(sample_type_code %in% c("01","02","06"),
                      "Tumor",
                      ifelse(sample_type_code == "11", "Normal", NA))
  
  metadata <- data.frame(
    sample = sample_names,
    condition = condition
  )
  
  rownames(metadata) <- metadata$sample
  
  # Remove NA samples if any
  keep <- !is.na(metadata$condition)
  metadata <- metadata[keep, ]
  df <- df[, rownames(metadata)]
  
  # -----------------------------------
  # CHECK SAMPLE COUNTS
  # -----------------------------------
  
  sample_counts <- table(metadata$condition)
  
  cat("Tumor samples:", sample_counts["Tumor"], "\n")
  cat("Normal samples:", sample_counts["Normal"], "\n")
  
  # Ensure at least a few normal samples exist
  if (!("Normal" %in% names(sample_counts)) || sample_counts["Normal"] < 3) {
    cat("Skipping", cancer, "- insufficient Normal samples\n\n")
    next
  }
  
  # Ensure at least a few tumor samples exist
  if (!("Tumor" %in% names(sample_counts)) || sample_counts["Tumor"] < 3) {
    cat("Skipping", cancer, "- insufficient Tumor samples\n\n")
    next
  }
  
  # -----------------------------------
  # 5. SET NORMAL AS REFERENCE
  # -----------------------------------
  
  metadata$condition <- factor(metadata$condition)
  metadata$condition <- relevel(metadata$condition, ref = "Normal")
  
  # -----------------------------------
  # 6. CREATE DESEQ OBJECT
  # -----------------------------------
  
  cat("Counts matrix dimensions:", dim(df), "\n")
  cat("Metadata dimensions:", dim(metadata), "\n")
  
  dds <- DESeqDataSetFromMatrix(
    countData = df,
    colData = metadata,
    design = ~ condition
  )
  
  # -----------------------------------
  # 7. RUN DESEQ
  # -----------------------------------
  
  dds <- DESeq(dds)
  
  # -----------------------------------
  # 8. GET RESULTS
  # -----------------------------------
  
  #res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
  res <- results(dds)
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  res_df <- res_df %>%
    arrange(padj)
  
  # -----------------------------------
  # 9. SAVE RESULTS CSV
  # -----------------------------------
  
  fwrite(
    res_df,
    paste0("DESeq2_results_", cancer, ".csv"),
    row.names = FALSE
  )
  
  # -----------------------------------
  # 10. VST TRANSFORMATION
  # -----------------------------------
  
  vsd <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vsd)
  
  # Save VST matrix
  fwrite(as.data.frame(vst_matrix),
         paste0("VST_matrix_", cancer, ".csv"))
}

cat("Analysis complete for all cancers.\n")