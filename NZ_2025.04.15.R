library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(Azimuth)
library(matrixStats)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
library(hdf5r)
library(Matrix)
setwd('D:/R_Genomix/S2/unser/data1/')

seurat_object1 <- readRDS("10963_Fcol_5GEX.rds")
seurat_object1 <- CreateSeuratObject(counts = seurat_object1)
print(seurat_object1)

# Переклад шару "counts" до "data"
DefaultAssay(seurat_object1) <- "RNA"
seurat_object1 <- SetAssayData(seurat_object1, slot = "data", new.data = GetAssayData(seurat_object1, slot = "counts"))
print(seurat_object1)
#DimPlot(seurat_object1)

seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^MT-")

VlnPlot(seurat_object1, features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)

# Відсіяти клітини, де частка мітохондріальних генів <= 5%
seurat_object1mt5 <- subset(seurat_object1, subset = percent.mt <= 5)
print(seurat_object1mt5)

# Identify the 10 most highly variable genes
# спершу знайти гени, що варіюють, і далі визначити 10, що найбільше
seurat_object1mt5 <- FindVariableFeatures(seurat_object1mt5)
top10mt5 <- head(VariableFeatures(seurat_object1mt5), 10)
print(top10mt5)


plot1mt5 <- FeatureScatter(seurat_object1mt5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2mt5 <- FeatureScatter(seurat_object1mt5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1mt5 + plot2mt5


all.genes <- rownames(seurat_object1mt5)
seurat_object1mt5 <- NormalizeData(seurat_object1mt5, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object1mt5 <- ScaleData(seurat_object1mt5, features = all.genes)
seurat_object1mt5 <- RunPCA(seurat_object1mt5, features = VariableFeatures(object = seurat_object1mt5))

VizDimLoadings(seurat_object1mt5, dims = 1:2, reduction = "pca")

ElbowPlot(seurat_object1mt5)
#згідно графіку, лише перші 4 компоненти найбільш вагомі, тому далі лише 4
seurat_object1mt5 <- FindNeighbors(seurat_object1mt5, dims = 1:4)
seurat_object1mt5 <- FindClusters(seurat_object1mt5, resolution = 0.5)
seurat_object1mt5 <- RunUMAP(seurat_object1mt5, dims = 1:4)

DimPlot(seurat_object1mt5, reduction = "umap")

saveRDS(seurat_object1mt5, file = "seurat_object1mt5_cluster.rds")

#знаходження маркерів для кожного кластеру
cluster0.markers <- FindMarkers(seurat_object1mt5, ident.1 = 0)
head(cluster0.markers, n = 5)

cluster1.markers <- FindMarkers(seurat_object1mt5, ident.1 = 1)
head(cluster1.markers, n = 5)

cluster2.markers <- FindMarkers(seurat_object1mt5, ident.1 = 2)
head(cluster2.markers, n = 5)

cluster3.markers <- FindMarkers(seurat_object1mt5, ident.1 = 3)
head(cluster3.markers, n = 5)

cluster4.markers <- FindMarkers(seurat_object1mt5, ident.1 = 4)
head(cluster4.markers, n = 5)

cluster5.markers <- FindMarkers(seurat_object1mt5, ident.1 = 5)
head(cluster5.markers, n = 5)

cluster6.markers <- FindMarkers(seurat_object1mt5, ident.1 = 6)
head(cluster6.markers, n = 5)

cluster7.markers <- FindMarkers(seurat_object1mt5, ident.1 = 7)
head(cluster7.markers, n = 5)

cluster8.markers <- FindMarkers(seurat_object1mt5, ident.1 = 8)
head(cluster8.markers, n = 5)

cluster9.markers <- FindMarkers(seurat_object1mt5, ident.1 = 9)
head(cluster9.markers, n = 5)

# find all markers for all clusters
all_markers <- FindAllMarkers(seurat_object1mt5)
head(all_markers) 

#DoHeatmap(seurat_object1mt5, features = top10mt5$gene) 

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_object1mt5, features = top10$gene) 

saveRDS(seurat_object1mt5, file = "seurat_object1mt5_cluster9.rds")


#перейменування кластерів відподно до клітин в них

# Перевизначити ідентифікатори, паралельно обєднавши кластери 0 та 1
seurat_object1mt5 <- RenameIdents(seurat_object1mt5,
                              '0' = 'T cells', 
                              '1' = 'T cells',
                              '2' = 'B cells',
                              '3' = 'Nerve and smooth muscle cells',
                              '4' = 'Glandular epithelium',
                              '5' = 'Endothelium',
                              '6' = 'Prostate cancer stromal cells',
                              '7' = 'Tumor-associated macrophages',
                              '8' = 'Tumor-associated monocytes',
                              '9' = 'Cancer-associated fibroblasts (CAFs)'
                              )


table(Idents(seurat_object1mt5))
DimPlot(seurat_object1mt5, reduction = "umap")
saveRDS(seurat_object1mt5, file = "seurat_object1mt5_cluster_NameClusters.rds")

#подивитися ген CD31 або ще названий PECAM1 (ендотеліальний маркер, але не у статті) у всіх кластерах
FeaturePlot(seurat_object1mt5, features = "PECAM1")
VlnPlot(seurat_object1mt5, features = "PECAM1", pt.size = 0.1)

#замінити ген CD31 // PECAM1 на ген CD200 як маркер ендотелію
FeaturePlot(seurat_object1mt5, features = "CD200")
VlnPlot(seurat_object1mt5, features = "CD200", pt.size = 0.1)

#маркер на фібробласти PDPN Podoplanin
FeaturePlot(seurat_object1mt5, features = "PDPN")
VlnPlot(seurat_object1mt5, features = "PDPN", pt.size = 0.1)
