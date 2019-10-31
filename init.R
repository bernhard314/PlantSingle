require(Seurat)
# rename files to matrix.mtx, genes.tsv and barcodes.tsv
EGEOD121619normalised<-Read10X("./Data/raw/")
Annotation<-read.table("./Data/ExpDesign-E-GEOD-121619.tsv",header = TRUE,sep="\t",row.names = 1)
EGEOD121619.seurat = CreateSeuratObject(counts = EGEOD121619normalised)
EGEOD121619.seurat<-AddMetaData(EGEOD121619.seurat,Annotation)
EGEOD121619.seurat<-NormalizeData(EGEOD121619.seurat,normalization.method = "LogNormalize")
EGEOD121619.seurat<-FindVariableFeatures(EGEOD121619.seurat)
EGEOD121619.seurat<-ScaleData(EGEOD121619.seurat,features = rownames(EGEOD121619.seurat)) 
EGEOD121619.seurat<-RunPCA(object = EGEOD121619.seurat)
DimPlot(EGEOD121619.seurat, reduction = "pca",group.by = "Factor.Value.temperature.")
EGEOD121619.seurat<-RunUMAP(EGEOD121619.seurat,dims = 1:10)
EGEOD121619.seurat<-AddMetaData(EGEOD121619.seurat,Annotation)
DimPlot(EGEOD121619.seurat, reduction = "umap",group.by = "Factor.Value.temperature.")
EGEOD121619.seurat_list<-SplitObject(EGEOD121619.seurat,split.by ="Factor.Value.temperature." )
for (i in names(EGEOD121619.seurat_list)) {
  EGEOD121619.seurat_list[[i]] <- SCTransform(EGEOD121619.seurat_list[[i]], verbose = TRUE)
}

for (i in names(EGEOD121619.seurat_list)) {
  EGEOD121619.seurat_list[[i]] <- RunPCA(EGEOD121619.seurat_list[[i]], verbose = TRUE)
}


gc()
Integration.features <- SelectIntegrationFeatures(object.list = EGEOD121619.seurat_list, nfeatures = 3000)
gc()

EGEOD121619.seurat_list <- PrepSCTIntegration(object.list = EGEOD121619.seurat_list, anchor.features = Integration.features)

Integration.anchors <- FindIntegrationAnchors(object.list = EGEOD121619.seurat_list, normalization.method = "SCT",  anchor.features = Integration.features,reduction = "cca")
# change rpca to cca
EGEOD121619.integrated <- IntegrateData(anchorset = Integration.anchors, normalization.method = "SCT")

EGEOD121619.integrated <- RunPCA(object = EGEOD121619.integrated, verbose = FALSE)
EGEOD121619.integrated <- RunUMAP(object = EGEOD121619.integrated, dims = 1:30)
EGEOD121619.integrated<-FindNeighbors(EGEOD121619.integrated,verbose = TRUE,dims=1:30)
EGEOD121619.integrated<-FindClusters(EGEOD121619.integrated,verbose = TRUE,dims=1:30,resolution = .2) #change resolution for more clusters



DimPlot(EGEOD121619.integrated,group.by =  "Factor.Value.temperature.")
DimPlot(EGEOD121619.integrated)
EGEOD121619.markers <- FindAllMarkers(EGEOD121619.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For 
globalDifExp<-FindMarkers(EGEOD121619.seurat,group.by = "Factor.Value.temperature.",ident.1 = "22 degree celsius",ident.2 = "38 degree celsius",test.use = "t")
write.csv2(globalDifExp,file="globalDE.csv")
