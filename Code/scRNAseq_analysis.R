library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(SingleR)
library(enrichR)
library(CellChat)
library(SingleCellExperiment)
library(SeuratWrappers)
library(tidyverse)
library(monocle3)
library(celldex)
library(ggpubr)


#read and save all samples-------------------------------------------------------------------------------------------------
BMMC_D1T1 <- readRDS('path..../GSM4138872_scRNA_BMMC_D1T1.rds')
BMMC_D1T2 <- readRDS('path..../GSM4138873_scRNA_BMMC_D1T2.rds')
CD34_D2T1 <- readRDS('path..../GSM4138874_scRNA_CD34_D2T1.rds')
CD34_D3T1 <- readRDS('path..../GSM4138875_scRNA_CD34_D3T1.rds')

#Create the Seurat Object for each sample-----------------------------------------------------------------------------------
BMMC_D1T1_so <- CreateSeuratObject(counts = BMMC_D1T1, project = "BMMC_D1T1", min.cells = 3, min.features = 200)
BMMC_D1T2_so <- CreateSeuratObject(counts = BMMC_D1T2, project = "BMMC_D1T2", min.cells = 3, min.features = 200)
CD34_D2T1_so <- CreateSeuratObject(counts = CD34_D2T1, project = "CD34_D2T1", min.cells = 3, min.features = 200)
CD34_D3T1_so <- CreateSeuratObject(counts = CD34_D3T1, project = "CD34_D3T1", min.cells = 3, min.features = 200)

#Create the metadata for all sample depends on assignment fiel-----------------------------------------------------------
BMMC_D1T1_so$Donor <- "D1"
BMMC_D1T1_so$Replictae <- "T1"
BMMC_D1T1_so$Sex <- "F"


BMMC_D1T2_so$Donor <- "D1"
BMMC_D1T2_so$Replictae <- "T2"
BMMC_D1T2_so$Sex <- "F"


CD34_D2T1_so$Donor <- "D2"
CD34_D2T1_so$Replictae <- "T1"
CD34_D2T1_so$Sex <- "M"


CD34_D3T1_so$Donor <- "D3"
CD34_D3T1_so$Replictae <- "T1"
CD34_D3T1_so$Sex <- "F"


#tring to calculate the percentage where is the MT or not and we cant find any percent---------------------------------------
BMMC_D1T1_so[["percent.mt"]] <- PercentageFeatureSet(object = BMMC_D1T1_so, pattern = "^MT-")
BMMC_D1T2_so[["percent.mt"]] <- PercentageFeatureSet(object = BMMC_D1T2_so, pattern = "^MT-")
CD34_D2T1_so[["percent.mt"]] <- PercentageFeatureSet(object = CD34_D2T1_so, pattern = "^MT-")
CD34_D3T1_so[["percent.mt"]] <- PercentageFeatureSet(object = CD34_D3T1_so, pattern = "^MT-")

VlnPlot(BMMC_D1T1_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(BMMC_D1T2_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CD34_D2T1_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CD34_D3T1_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Plot the Scattering for each sample-------------------------------------------------------------------------------------
FeatureScatter(BMMC_D1T1_so,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+ geom_smooth(method='lm')
FeatureScatter(BMMC_D1T2_so,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+ geom_smooth(method='lm')
FeatureScatter(CD34_D2T1_so,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+ geom_smooth(method='lm')
FeatureScatter(CD34_D3T1_so,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+ geom_smooth(method='lm')

#Filteration-------------------------------------------------------------------------------------------------------------
BMMC_D1T1_so <- subset(BMMC_D1T1_so,subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
BMMC_D1T2_so <- subset(BMMC_D1T2_so,subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
CD34_D2T1_so <- subset(CD34_D2T1_so,subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
CD34_D3T1_so <- subset(CD34_D3T1_so,subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

#Normalization-----------------------------------------------------------------------------------------------------------
BMMC_D1T1_so <- NormalizeData(BMMC_D1T1_so, normalization.method = "LogNormalize", scale.factor = 10000) 
BMMC_D1T2_so <- NormalizeData(BMMC_D1T2_so, normalization.method = "LogNormalize", scale.factor = 10000)
CD34_D2T1_so <- NormalizeData(CD34_D2T1_so, normalization.method = "LogNormalize", scale.factor = 10000)
CD34_D3T1_so <- NormalizeData(CD34_D3T1_so, normalization.method = "LogNormalize", scale.factor = 10000)

#preprocessing -------------------------------------------------------------------------------------------------------
BMMC_D1T1_so <- FindVariableFeatures(object = BMMC_D1T1_so)
BMMC_D1T1_so <- ScaleData(object =BMMC_D1T1_so)
BMMC_D1T1_so <- RunPCA(object = BMMC_D1T1_so)
ElbowPlot(BMMC_D1T1_so)
BMMC_D1T1_so <- FindNeighbors(object =  BMMC_D1T1_so, dims = 1:20)
BMMC_D1T1_so <- FindClusters(object = BMMC_D1T1_so)
BMMC_D1T1_so <- RunUMAP(object = BMMC_D1T1_so, dims = 1:20)

BMMC_D1T2_so <- FindVariableFeatures(object = BMMC_D1T2_so )
BMMC_D1T2_so  <- ScaleData(object =BMMC_D1T2_so )
BMMC_D1T2_so  <- RunPCA(object = BMMC_D1T2_so )
ElbowPlot(BMMC_D1T2_so )
BMMC_D1T2_so  <- FindNeighbors(object =  BMMC_D1T2_so , dims = 1:20)
BMMC_D1T2_so  <- FindClusters(object = BMMC_D1T2_so )
BMMC_D1T2_so  <- RunUMAP(object = BMMC_D1T2_so , dims = 1:20)

CD34_D2T1_so <- FindVariableFeatures(object = CD34_D2T1_so )
CD34_D2T1_so  <- ScaleData(object =CD34_D2T1_so )
CD34_D2T1_so  <- RunPCA(object = CD34_D2T1_so )
ElbowPlot(CD34_D2T1_so )
CD34_D2T1_so  <- FindNeighbors(object =  CD34_D2T1_so , dims = 1:20)
CD34_D2T1_so  <- FindClusters(object = CD34_D2T1_so )
CD34_D2T1_so  <- RunUMAP(object = CD34_D2T1_so , dims = 1:20)

CD34_D3T1_so <- FindVariableFeatures(object = CD34_D3T1_so )
CD34_D3T1_so  <- ScaleData(object =CD34_D3T1_so )
CD34_D3T1_so  <- RunPCA(object = CD34_D3T1_so )
ElbowPlot(CD34_D3T1_so )
CD34_D3T1_so  <- FindNeighbors(object =  CD34_D3T1_so , dims = 1:20)
CD34_D3T1_so  <- FindClusters(object = CD34_D3T1_so )
CD34_D3T1_so  <- RunUMAP(object = CD34_D3T1_so , dims = 1:20)


#no ground truth to find the optimal pk value -----------------------------------------------------------------------------
BMM_D1T1_so_Sweep_res <- paramSweep_v3(BMMC_D1T1_so, PCs = 1:20, sct = FALSE)
BMM_D1T1_so_Sweep_summ <- summarizeSweep(BMM_D1T1_so_Sweep_res, GT = FALSE)
PK_BMM_D1T1_so <- find.pK(BMM_D1T1_so_Sweep_summ)

ggplot(PK_BMM_D1T1_so, aes(pK, BCmetric, group = 1)) +
  geom_point()+
  geom_line()


BMMC_D1T2_so_Sweep_res <- paramSweep_v3(BMMC_D1T2_so, PCs = 1:20, sct = FALSE)
BMMC_D1T2_so_Sweep_summ <- summarizeSweep(BMMC_D1T2_so_Sweep_res, GT = FALSE)
PK_BMMC_D1T2_so <- find.pK(BMMC_D1T2_so_Sweep_summ)

ggplot(PK_BMMC_D1T2_so, aes(pK, BCmetric, group = 1)) +
  geom_point()+
  geom_line()


CD34_D2T1_so_Sweep_res <- paramSweep_v3(CD34_D2T1_so, PCs = 1:20, sct = FALSE)
CD34_D2T1_so_Sweep_summ <- summarizeSweep(CD34_D2T1_so_Sweep_res, GT = FALSE)
PK_CD34_D2T1_so<- find.pK(CD34_D2T1_so_Sweep_summ)

ggplot(PK_CD34_D2T1_so, aes(pK, BCmetric, group = 1)) +
  geom_point()+
  geom_line()


CD34_D3T1_so_Sweep_res <- paramSweep_v3(CD34_D3T1_so, PCs = 1:20, sct = FALSE)
CD34_D3T1_so_Sweep_summ <- summarizeSweep(CD34_D3T1_so_Sweep_res, GT = FALSE)
PK_CD34_D3T1_so<- find.pK(CD34_D3T1_so_Sweep_summ)

ggplot(PK_CD34_D3T1_so, aes(pK, BCmetric, group = 1)) +
  geom_point()+
  geom_line()


#homotyic finding --------------------------------------------------------------------------------------------------------
optimal_pk_BMMC_D1T1 <- PK_BMM_D1T1_so %>%
  filter(BCmetric == max(BCmetric))  %>%
  select(pK)
optimal_pk_BMMC_D1T1 <- as.numeric(as.character((optimal_pk_BMMC_D1T1[[1]])))
an_BMMC_D1T1_so <- BMMC_D1T1_so@meta.data$seurat_clusters
homo_BMMC_D1T1_so <- modelHomotypic(an_BMMC_D1T1_so)
exp_BMMC_D1T1 <- round(0.75*nrow(BMMC_D1T1_so@meta.data))
exp_doublet_BMMC_D1T1 <- round(exp_BMMC_D1T1*(1-homo_BMMC_D1T1_so))

optimal_pk_BMMC_D1T2 <- PK_BMMC_D1T2_so %>%
  filter(BCmetric == max(BCmetric))  %>%
  select(pK)
optimal_pk_BMMC_D1T2 <- as.numeric(as.character((optimal_pk_BMMC_D1T2[[1]])))
an_BMMC_D1T2_so <- BMMC_D1T2_so@meta.data$seurat_clusters
homo_BMMC_D1T2_so <- modelHomotypic(an_BMMC_D1T2_so)
exp_BMMC_D1T2 <- round(0.75*nrow(BMMC_D1T2_so@meta.data))
exp_doublet_BMMC_D1T2 <- round(exp_BMMC_D1T2*(1-homo_BMMC_D1T2_so))

optimal_pk_CD34_D2T1 <- PK_CD34_D2T1_so %>%
  filter(BCmetric == max(BCmetric))  %>%
  select(pK)
optimal_pk_CD34_D2T1 <- as.numeric(as.character((optimal_pk_CD34_D2T1[[1]])))
an_CD34_D2T1_so <- CD34_D2T1_so@meta.data$seurat_clusters
homo_CD34_D2T1_so <- modelHomotypic(an_CD34_D2T1_so)
exp_CD34_D2T1 <- round(0.75*nrow(CD34_D2T1_so@meta.data))
exp_doublet_CD34_D2T1 <- round(exp_CD34_D2T1*(1-homo_CD34_D2T1_so))

optimal_pk_CD34_D3T1 <- PK_CD34_D3T1_so %>%
  filter(BCmetric == max(BCmetric))  %>%
  select(pK)
optimal_pk_CD34_D3T1 <- as.numeric(as.character((optimal_pk_CD34_D3T1[[1]])))
an_CD34_D3T1_so <- CD34_D3T1_so@meta.data$seurat_clusters
homo_CD34_D3T1_so <- modelHomotypic(an_CD34_D3T1_so)
exp_CD34_D3T1 <- round(0.75*nrow(CD34_D3T1_so@meta.data))
exp_doublet_CD34_D3T1 <- round(exp_CD34_D3T1*(1-homo_CD34_D3T1_so))

# Doubletfinder runing --------------------------------------------------------------------------------------------------
BMMC_D1T1_so <- doubletFinder_v3(BMMC_D1T1_so,
                                 PCs = 1:20,
                                 pN =0.25,
                                 pK = optimal_pk_BMMC_D1T1,
                                 nExp = exp_BMMC_D1T1,
                                 reuse.pANN = FALSE, sct = FALSE)

BMMC_D1T2_so <- doubletFinder_v3(BMMC_D1T2_so,
                                 PCs = 1:20,
                                 pN =0.25,
                                 pK = optimal_pk_BMMC_D1T2,
                                 nExp = exp_BMMC_D1T2,
                                 reuse.pANN = FALSE, sct = FALSE)

CD34_D2T1_so <- doubletFinder_v3(CD34_D2T1_so,
                                 PCs = 1:20,
                                 pN =0.25,
                                 pK = optimal_pk_CD34_D2T1,
                                 nExp = exp_CD34_D2T1,
                                 reuse.pANN = FALSE, sct = FALSE)

CD34_D3T1_so <- doubletFinder_v3(CD34_D3T1_so,
                                 PCs = 1:20,
                                 pN =0.25,
                                 pK = optimal_pk_CD34_D3T1,
                                 nExp = exp_CD34_D3T1,
                                 reuse.pANN = FALSE, sct = FALSE)

#visualization ---------------------------------------------------------------------------------------------------------
DimPlot(BMMC_D1T1_so, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_4276")
DimPlot(BMMC_D1T2_so, reduction = 'umap', group.by = "DF.classifications_0.25_0.28_4330")
DimPlot(CD34_D2T1_so, reduction = 'umap', group.by = "DF.classifications_0.25_0.18_1159")
DimPlot(CD34_D3T1_so, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_4026")


#merging data without Batch correction ---------------------------------------------------------------------------------
merged_lis_so <- merge(CD34_D2T1_so, y = c(CD34_D3T1_so, BMMC_D1T2_so,BMMC_D1T1_so))


merged_lis_so <- NormalizeData(object = merged_lis_so)
merged_lis_so <- FindVariableFeatures(object = merged_lis_so )
merged_lis_so  <- ScaleData(object = merged_lis_so )
merged_lis_so  <- RunPCA(object = merged_lis_so )
ElbowPlot(merged_lis_so )
merged_lis_so  <- FindNeighbors(object =  merged_lis_so , dims = 1:20)
merged_lis_so  <- FindClusters(object = merged_lis_so, resolution = 0.26 )
merged_lis_so <- RunUMAP(object = merged_lis_so , dims =1:20)


#Batch effect using seurat merging --------------------------------------------------------------------------------------
Integrated_so <-SplitObject(merged_lis_so, split.by = "orig.ident")
Integrated_so <- lapply(X = Integrated_so, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Integrated_so)

anchors <- FindIntegrationAnchors(object.list = Integrated_so, anchor.features = features)

seurat.Inegrated <- IntegrateData(anchorset = anchors)


seurat.Inegrated <- ScaleData(seurat.Inegrated, verbose = FALSE)
seurat.Inegrated <- RunPCA(seurat.Inegrated, npcs = 30, verbose = FALSE)
seurat.Inegrated <- FindNeighbors(seurat.Inegrated, reduction = "pca", dims = 1:20)
seurat.Inegrated <- FindClusters(seurat.Inegrated, resolution = 0.26)
seurat.Inegrated <- RunUMAP(seurat.Inegrated, reduction = "pca", dims = 1:20)


DimPlot(merged_lis_so, reduction = 'umap', group.by ="orig.ident")+
  DimPlot(seurat.Inegrated, reduction = 'umap', group.by ="orig.ident")


DimPlot(merged_lis_so, reduction = 'umap', group.by ="Sex")+
  DimPlot(seurat.Inegrated, reduction = 'umap', group.by ="Sex")


DimPlot(merged_lis_so, reduction = 'umap', group.by ="Replictae")+
  DimPlot(seurat.Inegrated, reduction = 'umap', group.by ="Replictae")


DimPlot(merged_lis_so, reduction = 'umap', group.by ="Donor")+
  DimPlot(seurat.Inegrated, reduction = 'umap', group.by ="Donor")

#Dim reduction -----------------------------------------------------------------------------------------------------------
DR_seurat.Inegrated <- FindVariableFeatures(seurat.Inegrated, selection.method = "vst", nfeatures = 2000)
DR_top10 <- head(VariableFeatures(DR_seurat.Inegrated), 10)
plot1 <- VariableFeaturePlot(DR_seurat.Inegrated)
plot2 <- LabelPoints(plot = plot1, points = DR_top10, repel = TRUE)
plot1 + plot2

DR_seurat.Inegrated <- RunPCA(object = seurat.Inegrated,features = VariableFeatures(object = seurat.Inegrated))
ElbowPlot(DR_seurat.Inegrated)

VizDimLoadings(DR_seurat.Inegrated, dims = 1:2, reduction = "pca")
DimPlot(DR_seurat.Inegrated, reduction = "pca")
DimHeatmap(DR_seurat.Inegrated, dims = 1, cells = 500, balanced = TRUE)
DR_seurat.Inegrated <- JackStraw(DR_seurat.Inegrated, num.replicate = 100)
DR_seurat.Inegrated <- ScoreJackStraw(DR_seurat.Inegrated, dims = 1:20)
JackStrawPlot(DR_seurat.Inegrated, dims = 1:15)
DR_seurat.Inegrated <- RunUMAP(DR_seurat.Inegrated, dims = 1:20)
DimPlot(DR_seurat.Inegrated, reduction = 'umap')

#clusstering -------------------------------------------------------------------------------------------------------------
DR_seurat.Inegrated <- FindNeighbors(DR_seurat.Inegrated, reduction = "pca", dims = 1:30)
DR_seurat.Inegrated <- FindClusters(DR_seurat.Inegrated, resolution = 0.24)
DR_seurat.Inegrated <- RunUMAP(DR_seurat.Inegrated, dims = 1:10)

clusterd_int_seurat <- DR_seurat.Inegrated

# visualize data -----------------------------------------------------------------------------------------------------------
DimPlot(DR_seurat.Inegrated, reduction ='umap',group.by = 'orig.ident')+
  DimPlot(DR_seurat.Inegrated, reduction ='umap',group.by = 'seurat_clusters',label = TRUE)


#Automatic annotation -----------------------------------------------------------------------------------------------------
Human_mapping <- HumanPrimaryCellAtlasData()

sc_Data <- as.SingleCellExperiment(DietSeurat(clusterd_int_seurat))
Annotation_Auto <- SingleR(test = sc_Data, ref = Human_mapping, assay.type.test = 1,labels = Human_mapping$label.main)
table(Annotation_Auto$pruned.labels)

clusterd_int_seurat@meta.data$annotation_auto <- Annotation_Auto$pruned.labels
clusterd_int_seurat <- SetIdent(clusterd_int_seurat, value = "Annotation_Auto")
DimPlot(clusterd_int_seurat, reduction = "umap", group.by = "annotation_auto"  , label = T, repel = T, label.size = 3)+
  DimPlot(DR_seurat.Inegrated, reduction ='umap',group.by = 'seurat_clusters',label = TRUE)

#manual annotation --------------------------------------------------------------------------------------------------------
Markers <- FindAllMarkers(DR_seurat.Inegrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


HSC <- FeaturePlot(object = DR_seurat.Inegrated,
                   features = c("CD34", "CD38", "Sca1", "Kit"),
                   sort.cell = TRUE,
                   min.cutoff = 'q10',
                   label = TRUE,
                   repel = TRUE)


LMPP <- FeaturePlot(object = DR_seurat.Inegrated,
                    features = c("CD38", "CD52", "CSF3R", "ca1", "Kit", "CD34", "Flk2"),
                    sort.cell = TRUE,
                    min.cutoff = 'q10',
                    label = TRUE,
                    repel = TRUE)

CLP <- FeaturePlot(object = DR_seurat.Inegrated,
                   features = c("IL7R"),
                   sort.cell = TRUE,
                   min.cutoff = 'q10',
                   label = TRUE,
                   repel = TRUE)



GMP <- FeaturePlot(object = DR_seurat.Inegrated,
                   features = c("ELANE"),
                   sort.cell = TRUE,
                   min.cutoff = 'q10',
                   label = TRUE,
                   repel = TRUE)


Common_Progenitor <- FeaturePlot(object = DR_seurat.Inegrated,
                                 features = c("IL3","GM-CSF","M-CSF"),
                                 sort.cell = TRUE,
                                 min.cutoff = 'q10',
                                 label = TRUE,
                                 repel = TRUE)



B_Cell <- FeaturePlot(object = DR_seurat.Inegrated,
                      features = c("CD19","CD20","CD38"),
                      sort.cell = TRUE,
                      min.cutoff = 'q10',
                      label = TRUE,
                      repel = TRUE)


Pre_B <- FeaturePlot(object = DR_seurat.Inegrated,
                     features = c("CD19","CD34"),
                     sort.cell = TRUE,
                     min.cutoff = 'q10',
                     label = TRUE,
                     repel = TRUE)



Plasma <- FeaturePlot(object = DR_seurat.Inegrated,
                      features = c("SDC1","IGHA1","IGLC1","MZB1","JCHAIN"),
                      sort.cell = TRUE,
                      min.cutoff = 'q10',
                      label = TRUE,
                      repel = TRUE)





CD8 <- FeaturePlot(object = DR_seurat.Inegrated,
                   features = c("CD3D","CD3E","CD8A","CD8B"),
                   sort.cell = TRUE,
                   min.cutoff = 'q10',
                   label = TRUE,
                   repel = TRUE)


CD4 <- FeaturePlot(object = DR_seurat.Inegrated,
                   features = c("CD3D","CD3E","CD4"),
                   sort.cell = TRUE,
                   min.cutoff = 'q10',
                   label = TRUE,
                   repel = TRUE)



NK <- FeaturePlot(object = DR_seurat.Inegrated,
                  features = c("FCGR3A","NCAM1","NKG7","KLRB1"),
                  sort.cell = TRUE,
                  min.cutoff = 'q10',
                  label = TRUE,
                  repel = TRUE)


Erythrocytes <- FeaturePlot(object = DR_seurat.Inegrated,
                            features = c("GATA1","HBB","HBA1","HBA2"),
                            sort.cell = TRUE,
                            min.cutoff = 'q10',
                            label = TRUE,
                            repel = TRUE)


pDC<- FeaturePlot(object = DR_seurat.Inegrated,
                  features = c("IRF8","IRF4","IRF7"),
                  sort.cell = TRUE,
                  min.cutoff = 'q10',
                  label = TRUE,
                  repel = TRUE)



cDC <- FeaturePlot(object = DR_seurat.Inegrated,
                   features = c("CD1C","CD207","ITGAM","NOTCH2","SIRPA"),
                   sort.cell = TRUE,
                   min.cutoff = 'q10',
                   label = TRUE,
                   repel = TRUE)


CD14 <- FeaturePlot(object = DR_seurat.Inegrated,
                    features = c("CD14","CCL3","CCL4","IL1B"),
                    sort.cell = TRUE,
                    min.cutoff = 'q10',
                    label = TRUE,
                    repel = TRUE)



CD16 <- FeaturePlot(object = DR_seurat.Inegrated,
                    features = c("FCGR3A","CD68","S100A12"),
                    sort.cell = TRUE,
                    min.cutoff = 'q10',
                    label = TRUE,
                    repel = TRUE)



Basophils <- FeaturePlot(object = DR_seurat.Inegrated,
                         features = c("GATA2"),
                         sort.cell = TRUE,
                         min.cutoff = 'q10',
                         label = TRUE,
                         repel = TRUE)

plot(HSC)
plot(LMPP)
plot(CLP)
plot(GMP)
plot(Common_Progenitor)
plot(B_Cell)
plot(Pre_B)
plot(Plasma)
plot(CD8)
plot(CD4)
plot(NK)
plot(Erythrocytes)
plot(pDC)
plot(cDC)
plot(CD14)
plot(CD16)
plot(Basophils)

#Assign the markers -------------------------------------------------------------------------------------------------------
New.cluster.ids <- c("T-Cells","CD14","GMP","Plasma","Basophils","cDC","Pre-Cell","Erythrocytes","CLP","CD8","B-Cell","CD4","LMPP","NK","CD16")
names(New.cluster.ids) <- levels(DR_seurat.Inegrated)
DR_seurat.Inegrated <- RenameIdents(DR_seurat.Inegrated,New.cluster.ids) 

DimPlot(object = DR_seurat.Inegrated, reduction = "umap",  label = TRUE )+
  DimPlot(clusterd_int_seurat, reduction = "umap", group.by = "annotation_auto"  , label = T, repel = T, label.size = 3)


# gene-expression -------------------------------------------------------------------------------------------------------
VlnPlot(DR_seurat.Inegrated, features = c ("CD19"), group.by = "seurat_clusters")
VlnPlot(DR_seurat.Inegrated, features = c ("NKG7"), group.by = "seurat_clusters")
VlnPlot(DR_seurat.Inegrated, features = c ("CD3D"), group.by = "seurat_clusters")
FeaturePlot(DR_seurat.Inegrated, features = c("CD19","NKG7","CD3D"), min.cutoff = "q9")




# Volcano plot --------------------------------------------------------------------------------------------------------------

Markers_B_T <-  FindMarkers(object =  DR_seurat.Inegrated, ident.1 = c("CD8","CD4","T-Cells"), ident.2 = "B-Cell",  min.pct = 0.25)
head(Markers_B_T, n = 5) 

Markers_B_T$expressed <- "NO"
dim(Markers_B_T) 

Markers_B_T$expressed[Markers_B_T$avg_log2FC > 0.25 & Markers_B_T$p_val_adj < 0.01] <- "Up"
Markers_B_T$expressed[Markers_B_T$avg_log2FC < 0.25 & Markers_B_T$p_val_adj < 0.01] <- "DOWN"

Markers_B_T$labels <- NA
ggplot(data = Markers_B_T,aes(x=avg_log2FC , y= -log10(p_val), col = expressed, label = labels))+
  geom_point()+
  theme_minimal()+
  geom_text()+
  scale_color_manual(values = c('blue','black','red'))+
  geom_vline(xintercept = c(-0.2,0.2),col='red')+
  geom_hline(yintercept = -log10(0.00001),col='red')+
  theme(text = element_text(size = 20))

Markers_B_T$rn <- row.names(Markers_B_T)
arrange(Markers_B_T, p_val)   
head(arrange(Markers_B_T,p_val),10)$p_val

thresh = head(arrange(Markers_B_T,p_val),10)$p_val[10]
thresh
Markers_B_T$labels[Markers_B_T$p_val<= thresh] <- (Markers_B_T$rn[Markers_B_T$p_val <= thresh])
head(arrange(Markers_B_T,p_val),10)


ggplot(data = Markers_B_T,aes(x=avg_log2FC , y= -log10(p_val), col = expressed, label = labels))+
  geom_point()+
  theme_minimal()+
  geom_text()+
  scale_color_manual(values = c('blue','black','red'))+
  geom_vline(xintercept = c(-0.2,0.2),col='red')+
  geom_hline(yintercept = -log10(0.00001),col='red')+
  theme(text = element_text(size = 20))


# volcano plot CD4 and CD14-----------------------------------------------------------------------------------------------

Markers_4_14 <-  FindMarkers(object =  DR_seurat.Inegrated, ident.1 = c("CD8","CD4"), ident.2 = c("CD14","CD16"),  min.pct = 0.25)
head(Markers_4_14, n = 5) 

Markers_4_14$expressed <- "NO"
dim(Markers_4_14) 

Markers_4_14$expressed[Markers_4_14$avg_log2FC > 0.25 & Markers_4_14$p_val_adj < 0.01] <- "Up"
Markers_4_14$expressed[Markers_4_14$avg_log2FC < 0.25 & Markers_4_14$p_val_adj < 0.01] <- "DOWN"

Markers_4_14$labels <- NA
ggplot(data = Markers_4_14,aes(x=avg_log2FC , y= -log10(p_val), col = expressed, label = labels))+
  geom_point()+
  theme_minimal()+
  geom_text()+
  scale_color_manual(values = c('blue','black','red'))+
  geom_vline(xintercept = c(-0.2,0.2),col='red')+
  geom_hline(yintercept = -log10(0.00001),col='red')+
  theme(text = element_text(size = 20))


Markers_4_14$rn <- row.names(Markers_4_14)
arrange(Markers_4_14, p_val)   
head(arrange(Markers_4_14,p_val),10)$p_val

thresh = head(arrange(Markers_4_14,p_val),10)$p_val[10]
thresh
Markers_4_14$labels[Markers_4_14$p_val<= thresh] <- (Markers_4_14$rn[Markers_4_14$p_val <= thresh])
head(arrange(Markers_4_14,p_val),10)


ggplot(data = Markers_4_14,aes(x=avg_log2FC , y= -log10(p_val), col = expressed, label = labels))+
  geom_point()+
  theme_minimal()+
  geom_text()+
  scale_color_manual(values = c('blue','black','red'))+
  geom_vline(xintercept = c(-0.2,0.2),col='red')+
  geom_hline(yintercept = -log10(0.00001),col='red')+
  theme(text = element_text(size = 20))




