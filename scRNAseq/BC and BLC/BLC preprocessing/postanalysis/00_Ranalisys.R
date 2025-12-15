# conda activate monocle3

library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(hdf5r)

file <- "../analysis/cellranger/count/sample-HRI040891/outs/filtered_feature_bc_matrix.h5"

# load the count matrices
gene_names <- rownames(Seurat::Read10X_h5(file))
counts <- Seurat::Read10X_h5(file, use.names = F)

adj.matrix <- Read10X_h5(file)
srat <- CreateSeuratObject(adj.matrix, project = "HRI040891") 
adj.matrix <- NULL
srat

meta <- srat@meta.data
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

# Add mitocondrial and ribosomal fraction of reads
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")
srat[[]]
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) 
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")

# Normalization and dimensionality reduction
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat), 10)
top10 
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
VizDimLoadings(srat, dims = 1:9, reduction = "pca")
DimHeatmap(srat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(srat, reduction = "pca")
ElbowPlot(srat)

srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.5)
srat <- RunUMAP(srat, dims = 1:10, verbose = F)
table(srat@meta.data$seurat_clusters)
DimPlot(srat,label.size = 4,repel = T,label = T)

# Visualise gene expression
FeaturePlot(srat, features = "percent.mt") 
FeaturePlot(srat, features = "nFeature_RNA")
FeaturePlot(srat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

# Cellcycle
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)
table(srat[[]]$Phase)

# SCTransform normalization and clustering
srat <- SCTransform(srat, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = F)
table(srat[[]]$seurat_clusters)
DimPlot(srat, label = T)

all.markers <- FindAllMarkers(srat, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
table(all.markers$cluster)

top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers
write.table(as.data.frame(AverageExpression(object = srat, return.seurat = F)$RNA), sep = "\t", col.names = T, row.names = T, file = "AverageExpression.tsv")

# Cell type annotation using SingleR
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(DietSeurat(srat))
sce
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
table(hpca.main$pruned.labels)
table(hpca.fine$pruned.labels)
srat@meta.data$hpca.main   <- hpca.main$pruned.labels
srat@meta.data$hpca.fine   <- hpca.fine$pruned.labels

srat <- SetIdent(srat, value = "hpca.main")
DimPlot(srat, label = T , repel = T, label.size = 3) + NoLegend()

srat <- SetIdent(srat, value = "hpca.fine")
DimPlot(srat, label = T , repel = T, label.size = 3) + NoLegend()

# Work only with fibroblasts
srat_fibroblasts <- subset(x = srat, subset = hpca.main == "Fibroblasts")

srat_fibroblasts <- NormalizeData(srat_fibroblasts)
srat_fibroblasts <- FindVariableFeatures(srat_fibroblasts, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat_fibroblasts), 10)
top10 
plot1 <- VariableFeaturePlot(srat_fibroblasts)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(srat_fibroblasts)
srat_fibroblasts <- ScaleData(srat_fibroblasts, features = all.genes)
srat_fibroblasts <- RunPCA(srat_fibroblasts, features = VariableFeatures(object = srat_fibroblasts))
VizDimLoadings(srat_fibroblasts, dims = 1:9, reduction = "pca")
DimHeatmap(srat_fibroblasts, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(srat_fibroblasts, reduction = "pca")
ElbowPlot(srat_fibroblasts)

srat_fibroblasts <- FindNeighbors(srat_fibroblasts, dims = 1:10)
srat_fibroblasts <- FindClusters(srat_fibroblasts, resolution = 0.5)
srat_fibroblasts <- RunUMAP(srat_fibroblasts, dims = 1:10, verbose = F)
table(srat_fibroblasts@meta.data$seurat_clusters)
DimPlot(srat_fibroblasts,label.size = 4,repel = T,label = T)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat_fibroblasts <- CellCycleScoring(srat_fibroblasts, s.features = s.genes, g2m.features = g2m.genes)
table(srat_fibroblasts[[]]$Phase)

srat_fibroblasts <- SCTransform(srat_fibroblasts, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
srat_fibroblasts <- RunPCA(srat_fibroblasts, verbose = F)
srat_fibroblasts <- RunUMAP(srat_fibroblasts, dims = 1:30, verbose = F)
srat_fibroblasts <- FindNeighbors(srat_fibroblasts, dims = 1:30, verbose = F)
srat_fibroblasts <- FindClusters(srat_fibroblasts, verbose = F)
table(srat_fibroblasts[[]]$seurat_clusters)
DimPlot(srat_fibroblasts, label = T)

all.markers_fibroblasts <- FindAllMarkers(srat_fibroblasts, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers_fibroblasts)
table(all.markers_fibroblasts$cluster)

top3_markers_fibroblasts <- as.data.frame(all.markers_fibroblasts %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers_fibroblasts
write.table(as.data.frame(AverageExpression(object = srat_fibroblasts, return.seurat = F)$RNA), sep = "\t", col.names = T, row.names = T, file = "AverageExpression_fibroblasts.tsv")

# Work only with clusters 11, 19, 20, 24 and 25
srat <- SetIdent(srat, value = "seurat_clusters")
srat_select <- subset(x = srat, ident = c(11, 19, 20, 24, 25))

srat_select <- NormalizeData(srat_select)
srat_select <- FindVariableFeatures(srat_select, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat_select), 10)
top10 
plot1 <- VariableFeaturePlot(srat_select)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(srat_select)
srat_select <- ScaleData(srat_select, features = all.genes)
srat_select <- RunPCA(srat_select, features = VariableFeatures(object = srat_select))
VizDimLoadings(srat_select, dims = 1:9, reduction = "pca")
DimHeatmap(srat_select, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(srat_select, reduction = "pca")
ElbowPlot(srat_select)

srat_select <- FindNeighbors(srat_select, dims = 1:10)
srat_select <- FindClusters(srat_select, resolution = 0.5)
srat_select <- RunUMAP(srat_select, dims = 1:10, verbose = F)
table(srat_select@meta.data$seurat_clusters)
DimPlot(srat_select,label.size = 4,repel = T,label = T)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat_select <- CellCycleScoring(srat_select, s.features = s.genes, g2m.features = g2m.genes)
table(srat_select[[]]$Phase)

srat_select <- SCTransform(srat_select, method = "glmGamPoi", ncells = 8824, 
                                vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
srat_select <- RunPCA(srat_select, verbose = F)
srat_select <- RunUMAP(srat_select, dims = 1:30, verbose = F)
srat_select <- FindNeighbors(srat_select, dims = 1:30, verbose = F)
srat_select <- FindClusters(srat_select, verbose = F)
table(srat_select[[]]$seurat_clusters)
DimPlot(srat_select, label = T)

all.markers_select <- FindAllMarkers(srat_select, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers_select)
table(all.markers_select$cluster)

top3_markers_select <- as.data.frame(all.markers_select %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers_select
write.table(as.data.frame(AverageExpression(object = srat_select, return.seurat = F)$RNA), sep = "\t", col.names = T, row.names = T, file = "AverageExpression_clusters-11-19-20-24-25.tsv")



