# Датасет https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206426
# Title "A study on the radiosensitivity of radiation-induced lung injury at the acute phase based on single-cell transcriptomics"
# GEO accession: Series GSE206426.
# Organism 	Mus musculus

setwd("D:/Courses/2025_Genomics_R/1_Project/GSE206426")

library(Seurat)

# Якщо файли, такі як barcodes.tsv, відсутні, можна виконати попередню обробку даних scRNA-seq, використовуючи альтернативні підходи:
# Завантаження даних у форматі матриці каунтів, можна завантажити її безпосередньо в Seurat:

counts <- as.matrix(read.table("RILI_scRNAseq_UMIcount.txt", header = TRUE, row.names = 1))

seurat_obj <- CreateSeuratObject(counts = counts, project = "scRNAseq")

## Функції для нормалізації та фільтрації даних, як у стандартному аналізі:
seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000) # Ідентифікація змінних генів

# Візуалізація змінних генів
VariableFeaturePlot(seurat_obj)

# Масштабування даних
seurat_obj <- ScaleData(seurat_obj)

# Виконання PCA (головних компонент)
seurat_obj <- RunPCA(seurat_obj)

# Візуалізація PCA
DimPlot(seurat_obj, reduction = "pca")

# Визначення кластерів
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Виконання UMAP для візуалізації кластерів
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Визначення маркерів для кожного кластеру
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Перегляд маркерів для конкретного кластеру
head(cluster_markers[cluster_markers$cluster == 1, ])

# Візуалізація маркерів
FeaturePlot(seurat_obj, features = c("Fbln1", "Clec3b"))

######################################
#Встановлення пакета - scMCA (Single-cell Mouse Cell Atlas) для визначення типів клітин у даних scRNA-seq
install.packages("devtools")
install.packages("usethis")

library(usethis)
library(devtools)
install.packages("shinythemes")
install_github("ggjlab/scMCA")

#Завантаження бібліотеки
library(scMCA)

###визначення типів клітин у даних scRNA-seq мишей за допомогою пакета scMCA після кластеризації з використанням Seurat:

# Отримання матриці експресії генів
expression_matrix <- as.matrix(seurat_obj@assays$RNA@data)
#######error!!!!


# Побудова теплової карти для маркерів
library(dplyr)
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Збереження результатів top10 як .csv файл
write.csv(top10, file = "top10_markers.csv", row.names = FALSE)

DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()

### Масштабування всіх генів: За замовчуванням функція ScaleData масштабує лише змінні гени (Variable Features). Щоб масштабувати всі гени, виконайте наступне:
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

### Перевірка наявності генів: Переконайтеся, що всі гени з top10$gene присутні у даних:
intersect_genes <- intersect(top10$gene, rownames(seurat_obj[["RNA"]])) 

### Оновлення списку генів: Якщо деякі гени відсутні, можна оновити список для теплової карти:
top10_filtered <- top10[top10$gene %in% rownames(seurat_obj[["RNA"]]), ]
DoHeatmap(seurat_obj, features = top10_filtered$gene) + NoLegend()

### Функціональний аналіз генів. Пакети clusterProfiler для аналізу шляхів і функцій:
# Встановлення та завантаження пакетів:
BiocManager::install("org.Mm.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(org.Mm.eg.db)
library(clusterProfiler)

gene_list <- cluster_markers$gene

# Аналіз для біологічних процесів (BP - збагачення генів у біологічних процесах)
enrich_BP <- enrichGO(
  gene = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05
)
dotplot(enrich)

# Аналіз для молекулярних функцій MF - молекулярні функції
enrich_MF <- enrichGO(
  gene = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  pvalueCutoff = 0.05
)
dotplot(enrich_MF)


# Аналіз для клітинних компонентів CC - клітинні компоненти
enrich_CC <- enrichGO(
  gene = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pvalueCutoff = 0.05
)
dotplot(enrich_CC)
