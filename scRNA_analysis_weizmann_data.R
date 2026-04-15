library(Seurat)
library(Matrix)
library(tidyverse)

# Set path to your files
matrix_file <- "raw_files/Data_Lung/Data_Song2019_Lung/Exp_data_UMIcounts.mtx"
genes_file <- "raw_files/Data_Lung/Data_Song2019_Lung/Genes.txt"
barcodes_file <- "raw_files/Data_Lung/Data_Song2019_Lung/Cells.csv"

# 1. Load count matrix
counts <- readMM(matrix_file)

# 2. Load gene names
genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- genes$V1  # or genes$V2 if gene symbols are in 2nd column

# 3. Load barcodes + metadata
barcodes <- read.csv(barcodes_file, header = TRUE, stringsAsFactors = FALSE)
colnames(counts) <- barcodes$cell_name  # must match the names in barcodes

# 4. Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = "Song et al. 2019", min.cells = 3, min.features = 200)

# 5. Add metadata
# Match by cell name
metadata <- barcodes %>% column_to_rownames("cell_name")
seurat_obj <- AddMetaData(seurat_obj, metadata)

# 6. Add percent mitochondrial content
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 7. QC visualization
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "complexity"), ncol = 4)

# 8. Filtering (adjust thresholds as needed)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize, find variable genes, scale, PCA, clustering, UMAP, etc.
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)


seurat_obj$cell_type[seurat_obj$cell_type == ""] <- "Unknown"
seurat_obj$cell_type <- factor(seurat_obj$cell_type)  # Clean up factor levels
DimPlot(seurat_obj, group.by = "cell_type", label = TRUE, repel = TRUE)

# Annotation using SingleR
library(SingleR)
library(celldex)

# Load Seurat clusters
seurat_clusters <- Idents(seurat_obj)

# Human Primary Cell Atlas
ref_hpa <- celldex::HumanPrimaryCellAtlasData()
pred_hpa <- SingleR(test = GetAssayData(seurat_obj, slot = "data"),
                    ref = ref_hpa, labels = ref_hpa$label.main)
seurat_obj$SingleR_HPA <- pred_hpa$labels
seurat_obj$ClusterAnnotation_HPA <- tapply(pred_hpa$labels, seurat_clusters, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})[as.character(seurat_clusters)]

# Blueprint + ENCODE
ref_blueprint <- celldex::BlueprintEncodeData()
pred_blueprint <- SingleR(test = GetAssayData(seurat_obj, slot = "data"),
                          ref = ref_blueprint, labels = ref_blueprint$label.main)
seurat_obj$SingleR_Blueprint <- pred_blueprint$labels
seurat_obj$ClusterAnnotation_Blueprint <- tapply(pred_blueprint$labels, seurat_clusters, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})[as.character(seurat_clusters)]

# Monaco Immune
ref_monaco <- celldex::MonacoImmuneData()
pred_monaco <- SingleR(test = GetAssayData(seurat_obj, slot = "data"),
                       ref = ref_monaco, labels = ref_monaco$label.main)
seurat_obj$SingleR_Monaco <- pred_monaco$labels
seurat_obj$ClusterAnnotation_Monaco <- tapply(pred_monaco$labels, seurat_clusters, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})[as.character(seurat_clusters)]

# Database of Immune Cell Expression
ref_dice <- celldex::DatabaseImmuneCellExpressionData()
pred_dice <- SingleR(test = GetAssayData(seurat_obj, slot = "data"),
                     ref = ref_dice, labels = ref_dice$label.main)
seurat_obj$SingleR_DICE <- pred_dice$labels
seurat_obj$ClusterAnnotation_DICE <- tapply(pred_dice$labels, seurat_clusters, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})[as.character(seurat_clusters)]

# Plot single-cell annotations from different references
DimPlot(seurat_obj, group.by = "SingleR_HPA", label = TRUE, repel = TRUE)
DimPlot(seurat_obj, group.by = "SingleR_Blueprint", label = TRUE, repel = TRUE)
DimPlot(seurat_obj, group.by = "SingleR_Monaco", label = TRUE, repel = TRUE)
DimPlot(seurat_obj, group.by = "SingleR_DICE", label = TRUE, repel = TRUE)

# Plot cluster-level annotations
DimPlot(seurat_obj, group.by = "ClusterAnnotation_HPA", label = TRUE, repel = TRUE)
DimPlot(seurat_obj, group.by = "ClusterAnnotation_Blueprint", label = TRUE, repel = TRUE)
DimPlot(seurat_obj, group.by = "ClusterAnnotation_Monaco", label = TRUE, repel = TRUE)
DimPlot(seurat_obj, group.by = "ClusterAnnotation_DICE", label = TRUE, repel = TRUE)

saveRDS(seurat_obj, "annotated_rds/song2019_lung.rds")

# Get the list of genes
genes <- rownames(seurat_obj)

# Print the first few genes
head(genes)

# If you want to save the list to a file
write.csv(genes, "genes_song2019_lung.csv", row.names = FALSE)