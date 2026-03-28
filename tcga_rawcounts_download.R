library(TCGAbiolinks)
library(SummarizedExperiment)

# 1) Create a query
query <- GDCquery(
  project      = "TCGA-ACC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 2) Download raw count files
GDCdownload(query)

# 3) Prepare as a RangedSummarizedExperiment
se <- GDCprepare(query)

# 4) Extract raw count matrix (genes x samples)
counts <- assay(se, "unstranded")

# ===== loop =====
# ==== Downloading the files ======
library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)

# TCGA cancer projects
tcga_cancers <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

# Output directory
#outdir <- "TCGA_STAR_raw_counts"
#dir.create(outdir, showWarnings = FALSE)

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
  # ================== Extract unstranded raw counts ==================
  
  # 1. Keep only gene rows (remove N_* summary rows)
  #    These rows have empty gene_name or gene_id starting with "N_"
  se_genes <- se[gene_name != "" & !grepl("^N_", gene_id)]
  
  # 2. Identify unstranded raw count columns ONLY
  #    (exclude fpkm/tpm and other assays)
  unstranded_cols <- grep("^unstranded_", colnames(se_genes), value = TRUE)
  
  # 3. Subset to gene_name + unstranded counts
  raw_counts_dt <- se_genes[
    ,
    c("gene_name", unstranded_cols),
    with = FALSE
  ]
  
  # Set gene_name as first column name explicitly
  setnames(raw_counts_dt, "gene_name", "Gene")
  
  # 4. Optional: remove duplicated gene names (keep first)
  #raw_counts_dt <- raw_counts_dt[!duplicated(Gene)]
  
  # I am summing up the counts of duplicate gene rows because the duplicated row has mostly 0 and the actual has the values.
  # Aggregate duplicate genes by summing all counts
  raw_counts_dt <- raw_counts_dt[, lapply(.SD, sum), by = Gene]
  
  old_names <- colnames(raw_counts_dt)
  
  # Keep "Gene" as-is, strip "unstranded_" from sample columns
  new_names <- old_names
  new_names[-1] <- sub("^unstranded_", "", old_names[-1])
  
  setnames(raw_counts_dt, old = old_names, new = new_names)
  
  # Short cancer code (ACC, BRCA, etc.)
  cancer_code <- sub("TCGA-", "", proj)
  
  # Write CSV using fwrite (much faster than write.csv)
  # fwrite(
  #   raw_counts_dt,
  #   file = file.path(outdir, paste0("raw_counts_", cancer_code, ".csv"))
  # )
  fwrite(
    raw_counts_dt,
    file = paste0("raw_counts_", cancer_code, ".csv")
  )
  
  # Cleanup
  rm(query, se, raw_counts_dt, unstranded_cols, new_names, old_names, se_genes)
  gc()
}

message("All TCGA raw count CSV files generated successfully.")