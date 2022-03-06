
#####################################SingleCellData_Analysis##########################################################

--Step by step process of analyzing single cell data in "SEURAT"

=======================================================================================================================

library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
setwd("C:/Users/xxx/OneDrive/Documents/ISB/Data.....")
discovery.data <- Read10X(data.dir ="Discovery/")
discovery <- CreateSeuratObject(counts = discovery.data, project = "T2D", min.cells = 3, min.features = 200)
discovery
discovery[["percent.mt"]] <- PercentageFeatureSet(discovery, pattern = "^MT-")
discovery[["percent.mt"]]
VlnPlot(discovery, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
discovery <- subset(discovery, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
discovery
View(discovery)
View(discovery.data)
discovery <- CreateSeuratObject(counts = discovery.data, project = "T2D", min.cells = 3, min.features = 200)
discovery
VlnPlot(discovery, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
discovery[["percent.mt"]] <- PercentageFeatureSet(discovery, pattern = "^MT-")
VlnPlot(discovery, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
discovery <- subset(discovery, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
discovery <- NormalizeData(discovery, normalization.method = "LogNormalize", scale.factor = 10000)
discovery <- FindVariableFeatures(discovery, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(discovery), 10)
plot1 <- VariableFeaturePlot(discovery)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(discovey)
discovery <- ScaleData(discovery, features = all.genes)
all.genes <- rownames(discovery)
discovery <- ScaleData(discovery, features = all.genes)
discovery <- RunPCA(discovery, features = VariableFeatures(object = discovery))
discovery<- RunTSNE(discovery, dims = 1:15)
DimHeatmap(discovery, dims = 1, cells = 500, balanced = TRUE)
discovery <- JackStraw(discovery, num.replicate = 100)
discovery<- RunPCA(discovery, features = VariableFeatures(object = discovery))
DimPlot(discovery, reduction = "pca")
DimHeatmap(discovery, dims = 1, cells = 500, balanced = TRUE)
discovery <- JackStraw(discovery, num.replicate = 100)
discovery <- ScoreJackStraw(discovery, dims = 1:10)
ElbowPlot(discovery)
discovery <- FindNeighbors(discovery, dims = 1:15)
discovery <- FindClusters(discovery, resolution = 0.5)
head(Idents(discovery, 5))
discovery.markers <- FindAllMarkers(discovery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
discovery.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
VlnPlot(discovery, features = c("MMP12"), slot = "counts", log =TRUE)
FeaturePlot(discovery, features = c("MMP12", "CD9"))
top2000_discoveryMarkers <- discovery.markers %>% group_by(cluster) %>% top_n(n = 2000, wt = avg_logFC)
write.csv(file = "top200_discoveryMarkers.csv", top2000_discoveryMarkers)
saveRDS(discovery, file="discovery_4_14_2020.rds")
DimPlot(discovery)
FeaturePlot(discovery, features = c("GNG11", "SEPW1"))

discovery.metadata <- read.table(file = "GSE129363_Discovery_Cohort_CellAnnotation.txt", row.names = 1)
discovery_meta <- AddMetaData(object = discovery, metadata = discovery.metadata)
V3 <- FetchData(object = discovery_meta, vars = 'V3')
VlnPlot(discovery_meta, features = c("ITGAM", "ITGAX", "PTPRC","CD74", "CD14", "CD68", "FLT3"))
V2 <- FetchData(object = discovery_meta, vars = 'V2')
V4 <- FetchData(object = discovery_meta, vars = 'V4')
saveRDS(discovery, file="discovery_4_14_2020.rds")

DimPlot(discovery_meta, group.by = "V2", label = TRUE)
DimPlot(discovery_meta, group.by = "V3", label = TRUE)
t2d_subset <- subset(discovery_meta, subset = (V3 == "Diabetic"))
DimPlot(t2d_subset)
Nont2d_subset <- subset(discovery_meta, subset = (V3 == "NonDiabetic"))
####
Discovery_MFDCs_subset <- subset(discovery_meta, idents = c('0', '1', '2', '3', '4', '7', '8'), invert = TRUE)
Discovery_MFDCs <- NormalizeData(Discovery_MFDCs_subset)
Discovery_MFDCs <- FindVariableFeatures(Discovery_MFDCs, selection.method = "vst", nfeatures = 2000)
Discovery_MFDCs_top100 <- head(VariableFeatures(Discovery_MFDCs), 100)
write.csv(file = "Discovery_MFDCs_top100_Features.csv", Discovery_MFDCs_top100)
####
FeaturePlot(discovery_meta, features = c("ITGAM", "ITGAX", "PTPRC","CD74", "CD14", "CD68", "FLT3"))
FeaturePlot(discovery_meta, features = c("CD20", "CD19", "IGHA1" ))
plot1 <- VariableFeaturePlot(Discovery_MFDCs)
plot2 <- LabelPoints(plot = plot1, points = Discovery_MFDCs_top100, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(Discovery_MFDCs)
Discovery_MFDCs <- ScaleData(Discovery_MFDCs, features = all.genes)
VlnPlot(Discovery_MFDCs, features = c("nFeature_RNA", "nCount_RNA"))
Discovery_MFDCs<- RunPCA(Discovery_MFDCs, features = VariableFeatures(object = Discovery_MFDCs))
Discovery_MFDCs <- FindNeighbors(Discovery_MFDCs, dims = 1:15)
Discovery_MFDCs <- FindClusters(Discovery_MFDCs, resolution = 0.5)
Discovery_MFDCsmarkers <- FindAllMarkers(Discovery_MFDCs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Discovery_MFDCsmarkers <- Discovery_MFDCsmarkers %>% group_by(cluster) %>% top_n(n = 2000, wt = avg_logFC)
write.csv(file = "Discovery_MFDCsmarkers.csv", Discovery_MFDCsmarkers)
DimPlot(Discovery_MFDCs)
DimPlot(Discovery_MFDCs, label = TRUE)

##################Patientwise data###
Idents(object = DiscoveryMeta) <- 'V2'
levels(x =DiscoveryMeta)
Ate11VAT<- subset(DiscoveryMeta, subset = (V2 == "Ate-11-VAT"))
Ate11SAT<- subset(DiscoveryMeta, subset = (V2 == "Ate-11-SAT"))
Ate11<-merge(x = Ate11VAT, y = Ate11SAT, add.cell.ids = NULL, merge.data = TRUE, project = "Ate11")
Ate11 <-NormalizeData(Ate11)
Ate11 <-ScaleData(Ate11)
Ate11 <-FindVariableFeatures(Ate11, selection.method = "vst", nfeatures = 2000)
Ate11 <-RunPCA(Ate11, features = VariableFeatures(object = Ate11))
Ate11 <-RunTSNE(Ate11)
Ate11 <- FindNeighbors(Ate11, dims = 1:15)
Ate11 <- FindClusters(Ate11, resolution = 0.5)
VlnPlot(Ate11, features = c("ITGAM", "ITGAX", "PTPRC","CD74", "CD14", "CD68", "FLT3", "MMP12"))
FeaturePlot(Ate11, features = c("ITGAM", "ITGAX", "PTPRC","CD74", "CD14", "CD68", "FLT3", "MMP12"))
DimPlot(Ate11 , reduction = "tsne")
saveRDS(file="Ate11.rds", Ate11)
##Ate11<-readRDS(file ="Ate11.rds")

##################Cluster_RNA_Avg###

Ate11_clusterAverages <- AverageExpression(Ate11)
write.csv(file = "Ate11_clusterAverages.csv", Ate11_clusterAverages)
#########ConditionsCombined###########
t2d_mfdc<- readRDS(file="t2d_mfdc_subset_4_16_2020.rds")
nont2d_mfdc<- readRDS(file="Nont2d_mfdc_subset_4_16_2020.rds")
t2d_mfdc$condition <- "T2D"
nont2d_mfdc$condition <- "NonT2D"
immune.anchors <- FindIntegrationAnchors(object.list = list(t2d_mfdc, nont2d_mfdc), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DefaultAssay(immune.combined) <- "RNA"
clustcells_7.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "condition", verbose = FALSE)
head(clustcells_7.markers.markers)
DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("red", "blue"), dot.scale = 8,
split.by = "condition") + RotatedAxis()

###########Discovery_combined#########
> setwd("C:/Users/pavel/OneDrive/Documents/ISB/Data.....")
> discovery_meta<- readRDS(file="discovery_meta_4_14_2020.rds")
> t2d_subset <- subset(discovery_meta, subset = (V3 == "Diabetic"))
> Nont2d_subset <- subset(discovery_meta, subset = (V3 == "NonDiabetic"))
> t2d_subset$condition <- "T2D"
> Nont2d_subset$condition <- "NonT2D"
> discovety.anchors <- FindIntegrationAnchors(object.list = list(t2d_subset, Nont2d_subset), dims = 1:20)
discovery.combined <- IntegrateData(anchorset = discovety.anchors, dims = 1:20)
DefaultAssay(discovery.combined) <- "integrated"
discovery.combined <- ScaleData(discovery.combined, verbose = FALSE)
discovery.combined <- RunPCA(discovery.combined, npcs = 30, verbose = FALSE)
discovery.combined <- RunTSNE(discovery.combined, reduction = "pca", dims = 1:20)
discovery.combined <- FindNeighbors(discovery.combined, reduction = "pca", dims = 1:20)
discovery.combined <- FindClusters(discovery.combined, resolution = 0.5)
DefaultAssay(discovery.combined) <- "RNA"
DotPlot(discovery.combined, features = rev(markers.to.plot), cols = c("red", "blue"), dot.scale = 8,
split.by = "condition") + RotatedAxis()
 DotPlot(discovery.combined, features = rev(markers.to.plot), cols = c("red", "blue"), dot.scale = 8,
+ split.by = "condition", group.by = "condition") + RotatedAxis()
clustcells_1.markers <- FindConservedMarkers(discovery.combined, ident.1 = 1, grouping.var = "condition", verbose = FALSE)
head(clustcells_1.markers.markers)
############Cluster Names##########
dc_clusterRenamed <-RenameIdents(dc, `0` = "P2", `1` = "P1",  `2` ="IR-ATMs", `3` = "P3", `4` = "P6", `5` = "CD4 T cells", `6` = "APCs", `7` =  "NK cells", `8` =  "P5",  `9` = "E1", `10` =  "P4", `11` = "CD14 high cells", `12` = "E2", `13` ="B cells", `14` = "CD8 T cells")
> DimPlot(dc_clusterRenamed)
> DimPlot(dc_clusterRenamed, label = TRUE)
########## FC- condition##########
foldChanges<-FindMarkers(discovery_comb, ident.1 = 'T2D', ident.2='NonT2D', min.pct = 0.05)

#########Identify differential expressed genes across conditions###cluster vs cluster plot
> cluster1t.cells <- subset(immune.combined, idents = '1')
> Idents(cluster1t.cells) <- "condition"
> avg.cluster1t.cells <- log1p(AverageExpression(cluster1t.cells, verbose = FALSE)$RNA)
> avg.cluster1t.cells$gene <- rownames(avg.cluster1t.cells)
> p1 <- ggplot(avg.cluster1t.cells, aes(T2D, NonT2D)) + geom_point() + ggtitle("Cluster 1 Cells")
> p1 <- LabelPoints(plot = p1, points = markers.to.plot, repel = FALSE)
> plot(p1)

###################################################################################################################

How to install "SingCellaR" for Single Cell SingleCell Data Analysis

====================================================================================================================

```{r}
{install.packages('devtools')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
library(devtools)
install_github('supatt-lab/SingCellaR',ref='master', repos = BiocManager::repositories())}
```


library(SingCellaR)
library(ggplot2)
library(survminer)
library(circlize)
library(gridExtra)
library(ggpubr)

#########################################################################################################################

How to "read RDS" file and "convert" to Ann data objec (h5ad)

=========================================================================================================================

########### AnnData/H5AD files to SeuratDisk- intall in R ########

if (!requireNamespace("remotes", quietly = TRUE)) {
install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

######### Need this library and load it ################
library(SeuratDisk)

###########################This is what you need ######################

#########################Converting from Seurat to AnnData via h5Seurat##########
ab<-readRDS(file="200AB_projected.rds")
ab
SaveH5Seurat(ab, filename = "200AB.h5Seurat")
Convert("200AB.h5Seurat", dest = "h5ad")

####saved file in AnnData## 200AB.h5ad ######################

#####################Additional info #############################

# This creates a copy of this .h5ad object reformatted into .h5seurat
Convert("filename.h5ad", "filename.h5seurat")

# This .d5seurat object can then be read in manually
data <- LoadH5Seurat("filename.h5seurat")

####or to load RNA ####
data <- LoadH5Seurat("filename.h5seurat", assays = "RNA")

##SaveH5Seurat(data_2, overwrite = TRUE)
> SaveH5Seurat(data_2, filename = "filename_data.h5Seurat")
##data_2 <- Connect("filename_data.h5Seurat")