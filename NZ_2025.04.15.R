install.packages("Seurat")
install.packages("Matrix")
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("JASPAR2020", "TFBSTools"))


install.packages(c("Seurat", "SeuratData"))

# Data load
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(Azimuth)
library(matrixStats)

library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
library(hdf5r)


library(Seurat)
library(dplyr)
library(Matrix)
setwd('D:/R_Genomix/S2/unser/data1/')


seurat_object1 <- readRDS("10963_Fcol_5GEX.rds")
seurat_object1 <- CreateSeuratObject(counts = seurat_object1)
print(seurat_object1)
#DimPlot(seurat_object1)


seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^MT-")

VlnPlot(seurat_object1, features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object1), 10)

plot1 <- FeatureScatter(seurat_object1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

seurat_object1 <- FindVariableFeatures(seurat_object1)

all.genes <- rownames(seurat_object1)
seurat_object1 <- NormalizeData(seurat_object1, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object1 <- ScaleData(seurat_object1, features = all.genes)
seurat_object1 <- RunPCA(seurat_object1, features = VariableFeatures(object = seurat_object1))

VizDimLoadings(seurat_object1, dims = 1:2, reduction = "pca")

# Testing integration methods
ElbowPlot(seurat_object1)
seurat_object1 <- FindNeighbors(seurat_object1, dims = 1:11)
seurat_object1 <- FindClusters(seurat_object1, resolution = 0.5)
seurat_object1 <- RunUMAP(seurat_object1, dims = 1:11)

DimPlot(seurat_object1, reduction = "umap")

# Відсіяти клітини, де частка мітохондріальних генів <= 10%
seurat_object1mt10 <- subset(seurat_object1, subset = percent.mt <= 10)

top10mt10 <- head(VariableFeatures(seurat_object1mt10), 10)

plot1mt10 <- FeatureScatter(seurat_object1mt10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2mt10 <- FeatureScatter(seurat_object1mt10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1mt10 + plot2mt10

seurat_object1mt10 <- FindVariableFeatures(seurat_object1mt10)

all.genes <- rownames(seurat_object1mt10)
seurat_object1mt10 <- NormalizeData(seurat_object1mt10, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object1mt10 <- ScaleData(seurat_object1mt10, features = all.genes)
seurat_object1mt10 <- RunPCA(seurat_object1mt10, features = VariableFeatures(object = seurat_object1mt10))

VizDimLoadings(seurat_object1mt10, dims = 1:2, reduction = "pca")

ElbowPlot(seurat_object1mt10)
seurat_object1mt10 <- FindNeighbors(seurat_object1mt10, dims = 1:5)
seurat_object1mt10 <- FindClusters(seurat_object1mt10, resolution = 0.5)
seurat_object1mt10 <- RunUMAP(seurat_object1mt10, dims = 1:5)

DimPlot(seurat_object1mt10, reduction = "umap")
ElbowPlot(seurat_object1mt10)
seurat_object1mt10_11 <- FindNeighbors(seurat_object1mt10, dims = 1:11)
seurat_object1mt10_11 <- FindClusters(seurat_object1mt10, resolution = 0.5)
seurat_object1mt10_11 <- RunUMAP(seurat_object1mt10, dims = 1:11)

DimPlot(seurat_object1mt10_11, reduction = "umap")

seurat_object1mt10_11 <- FindVariableFeatures(seurat_object1mt10_11)
head(VariableFeatures(seurat_object1mt10_11))
markers1mt10_11 <- FindAllMarkers(seurat_object1mt10_11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers1mt10_11)
DotPlot(seurat_object1mt10_11, features = markers1mt10_11$gene) + RotatedAxis()

duplicated_genes <- markers1mt10_11$gene[duplicated(markers1mt10_11$gene)]
print(duplicated_genes)
unique_genes <- unique(markers1mt10_11$gene)
DotPlot(seurat_object1mt10_11, features = unique_genes) + RotatedAxis()
Idents(seurat_object1mt10_11) <- factor(Idents(seurat_object1mt10_11), levels = unique(Idents(seurat_object1mt10_11)))

unique_markers1 <- markers1mt10_11[!duplicated(markers1mt10_11$gene), ]
genes_by_cluster1 <- split(unique_markers1$gene, unique_markers1$cluster)
print(genes_by_cluster1[["1"]])

str(genes_by_cluster1)
print(genes_by_cluster1[["1"]])
DotPlot(seurat_object1mt10_11, features = genes_by_cluster1[["1"]]) + RotatedAxis()

duplicated_levels <- markers1mt10_11$gene[duplicated(markers1mt10_11$gene)]
print(duplicated_levels)
markers1mt10_11$gene <- make.unique(as.character(markers1mt10_11$gene))
DotPlot(seurat_object1mt10_11, features = markers1mt10_11$gene) + RotatedAxis()

#знайти де дані для обрахунку
slotNames(seurat_object1mt10_11@assays$RNA)
names(seurat_object1mt10_11@assays$RNA@layers)

#Calculate Fold Change:
expression_data <- seurat_object1mt10_11@assays$RNA@layers$data

markers1mt10_11$fold_change <- apply(expression_data[markers1mt10_11$gene, ], 1, function(x) max(x) / min(x + 1e-10))

slotNames(seurat_object1mt10_11@assays$RNA)
names(seurat_object1mt10_11@assays$RNA@layers)
#Check if All Genes Exist in expression_data: Verify that the genes in markers1mt10_11$gene are present in the row names of expression_data:
missing_genes <- setdiff(markers1mt10_11$gene, rownames(expression_data))
print(missing_genes)
#Filter markers1mt10_11 for Existing Genes: Remove genes from markers1mt10_11 that are not present in expression_data
markers1mt10_11 <- markers1mt10_11[markers1mt10_11$gene %in% rownames(expression_data), ]
#Recalculate Fold Change: Once the markers1mt10_11 object only contains genes present in expression_data, retry your fold change calculation
markers1mt10_11$fold_change <- apply(expression_data[markers1mt10_11$gene, ], 1, function(x) max(x) / min(x + 1e-10))
#Filter Markers with Fold Change > 3
filtered_markers <- markers1mt10_11[markers1mt10_11$fold_change > 1, ]
#Create DotPlot with Filtered Markers
DotPlot(seurat_object1mt10_11, features = filtered_markers$gene) + RotatedAxis()

missing_genes <- setdiff(filtered_markers$gene, rownames(seurat_object1mt10_11@assays$RNA@layers$data))
print(missing_genes)
DefaultAssay(seurat_object1mt10_11) <- "RNA"
DotPlot(seurat_object1mt10_11, features = filtered_markers$gene) + RotatedAxis()
rlang::last_trace()

missing_genes <- setdiff(filtered_markers$gene, rownames(seurat_object1mt10_11@assays$RNA@layers$data))
print(missing_genes)
filtered_markers <- filtered_markers[filtered_markers$gene %in% rownames(seurat_object1mt10_11@assays$RNA@layers$data), ]
DotPlot(seurat_object1mt10_11, features = filtered_markers$gene) + RotatedAxis()


