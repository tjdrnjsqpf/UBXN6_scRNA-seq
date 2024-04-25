library(Seurat)
library(dplyr)

data_UN <- Read10X(data.dir = "~$PATH/PBMC_UN")
data_LPS <- Read10X(data.dir = "~$PATH/PBMC_LPS")

# Create SeuratObject
data_UN.obj <- CreateSeuratObject(counts = data_UN,project = "UN")
data_LPS.obj <- CreateSeuratObject(counts = data_LPS,project = "LPS")

# Estimating ratio of mitochondrial RNAs
data_UN.obj[["percent.MTA"]] <- PercentageFeatureSet(data_UN.obj,pattern="MT.AT")
data_UN.obj[["percent.MTC"]] <- PercentageFeatureSet(data_UN.obj,pattern="MT.C")
data_UN.obj[["percent.MTND"]] <- PercentageFeatureSet(data_UN.obj,pattern="MT.ND")
data_UN.obj[["percent.mt"]] <- data_UN.obj$percent.MTA + data_UN.obj$percent.MTC +data_UN.obj$percent.MTND

data_LPS.obj[["percent.MTA"]] <- PercentageFeatureSet(data_LPS.obj,pattern="MT.AT")
data_LPS.obj[["percent.MTC"]] <- PercentageFeatureSet(data_LPS.obj,pattern="MT.C")
data_LPS.obj[["percent.MTND"]] <- PercentageFeatureSet(data_LPS.obj,pattern="MT.ND")
data_LPS.obj[["percent.mt"]] <- data_LPS.obj$percent.MTA + data_LPS.obj$percent.MTC +data_LPS.obj$percent.MTND

# Filtering 
data_UN.obj <- subset(data_UN.obj, subset= nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 10)
data_LPS.obj <- subset(data_LPS.obj, subset= nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 10)

Merged.obj <- merge(data_UN.obj, y=c(data_LPS.obj), add.cell.ids=c("UN","LPS"),project="Merged")


# SCTransform
Merged.obj <- SCTransform(Merged.obj,vars.to.regress = "percent.mt")


#FindVariableFeatures
Merged.obj <- FindVariableFeatures(Merged.obj)


#Scaling the data
all.genes <- row.names(Merged.obj)
Merged.obj <- ScaleData(Merged.obj, features = all.genes)


# PCA
Merged.obj <- RunPCA(Merged.obj, features = VariableFeatures(object = Merged.obj))


# Cluster the cells
Merged.obj <- FindNeighbors(Merged.obj, dims = 1:10)
Merged.obj <- FindClusters(Merged.obj, resolution = 0.05)

# Run non-linear dimensional reduction (UMAP/tSNE)

Merged.obj <- RunUMAP(Merged.obj, dims = 1:10)
DimPlot(Merged.obj, reduction = "umap")


# tSNE
Merged.obj <- RunTSNE(Merged.obj, dims = 1:3)
DimPlot(Merged.obj, reduction = "tsne",group.by = "orig.ident")
