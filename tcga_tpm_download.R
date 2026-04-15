# ==== Downloading the files ======
library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)

# TCGA cancer projects
tcga_cancers <- c(
  "TCGA-COAD", "TCGA-ESCA", "TCGA-READ", "TCGA-STAD"
)

# Output directory
outdir <- "TCGA_STAR_TPM_values"
dir.create(outdir, showWarnings = FALSE)

for (proj in tcga_cancers) {
  
  message("Processing ", proj, " ...")
  
  query <- GDCquery(
    project = proj,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query)
  
  se <- GDCprepare(query, summarizedExperiment = FALSE)
  
  #code to extract only unstranded and keep gene names as 1st column and remove the rest, and also the last few rows which aren't genes
  # ================== Extract unstranded tpm ==================
  
  # 1. Keep only gene rows (remove N_* summary rows)
  #    These rows have empty gene_name or gene_id starting with "N_"
  se_genes <- se[gene_name != "" & !grepl("^N_", gene_id)]
  
  # 2. Identify unstranded raw count columns ONLY
  #    (exclude fpkm and other assays)
  unstranded_cols <- grep("^tpm_unstranded_", colnames(se_genes), value = TRUE)
  
  # 3. Subset to gene_name + unstranded tpm
  tpm_dt <- se_genes[
    ,
    c("gene_name", "gene_type", unstranded_cols),
    with = FALSE
  ]
  
  # Set gene_name as first column name explicitly
  setnames(tpm_dt, "gene_name", "Gene")
  
  # 4. Optional: remove duplicated gene names (keep first)
  #raw_counts_dt <- raw_counts_dt[!duplicated(Gene)]
  
  # I am keeping the row with max mean expression of a duplicate, because the duplicated row has mostly 0 and the actual has the values.
  # i.e, which row has a higher mean.
  # Identify numeric columns
  num_cols <- setdiff(colnames(tpm_dt), c("Gene", "gene_type"))
  
  # Compute mean expression per row
  tpm_dt[, mean_expr := rowMeans(.SD), .SDcols = num_cols]
  
  # Keep highest per gene
  tpm_dt <- tpm_dt[
    order(Gene, -mean_expr)
  ][
    , .SD[1], by = Gene
  ][
    , mean_expr := NULL
  ]
  
  old_names <- colnames(tpm_dt)
  
  # Keep "Gene" as-is, strip "unstranded_" from sample columns
  new_names <- old_names
  new_names[-1] <- sub("^tpm_unstranded_", "", old_names[-1])
  
  setnames(tpm_dt, old = old_names, new = new_names)
  
  # Short cancer code (ACC, BRCA, etc.)
  cancer_code <- sub("TCGA-", "", proj)
  
  # Write CSV using fwrite (much faster than write.csv)
  fwrite(
    tpm_dt,
    file = file.path(outdir, paste0("tpm_", cancer_code, ".csv"))
  )

  # Cleanup
  rm(query, se, tpm_dt, unstranded_cols, new_names, old_names, se_genes)
  gc()
}

message("All TCGA TPM values CSV files generated successfully.")
