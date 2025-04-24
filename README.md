This project was made during the DAAD course of sc-RNAseq analysis (Team 7: Anna Yatsyshyna, Dmytro Aleksandrovych, Mariia Grandova, Zarina Nidoieva)
With the aim to study radisensitivity and radioresistance of cells, we have analysed several datasets:
1.GSE206426 - mouse lung samples in control group and after irradiation,
2.10963_Fcol_5GEX and 10997_PCaFcol_5GEX: irradiated human prostate tumors
Key stages of the analysis:
1.Data preprocessing and quality:
  Cell filtering: Removal of low-quality cells (low read counts, high mitochondrial gene abundance).
  Normalization (to eliminate variations caused by technical factors).
  Identification of variable genes (i.e. those that contribute most to differences between cells).
2. Data dimensionality reduction:
  Using PCA, UMAP, to visualize multidimensional data and reduce complexity to 2-3 dimensions.
3.Clustering:
  Grouping cells based on similarities in their expression profile to identify cell populations.
4.Cluster annotation:
  Identification of cell types in each cluster by matching expressed genes with known biomarkers.
5.Differential expression analysis:
  Search for genes whose expression differs between clusters or conditions (irradiated, control).
6.Functional analysis:
  Interpretation of biological functions using tools such as Gene Ontology (GO) or Pathway Enrichment.
7.Visualization of results:
  Construction of heat maps, gene expression plots, clustering and cell distribution.




