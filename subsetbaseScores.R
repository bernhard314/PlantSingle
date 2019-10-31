
SubsetCells<- colnames(EGEOD121619.seurat)[EGEOD121619.seurat$ModScoreCluster%in% c(11,12)]
EGEOD121619.seurat_sub<- subset(EGEOD121619.seurat,cells=SubsetCells)
EndodermDifExp<-FindMarkers(EGEOD121619.seurat_sub,group.by = "Factor.Value.temperature.",ident.1 = "22 degree celsius",ident.2 = "38 degree celsius",test.use = "t")
write.csv2(EndodermDifExp,file="EndodermDifExp.csv")

SubsetCells<- colnames(EGEOD121619.seurat)[EGEOD121619.seurat$ModScoreCluster%in% c(13,14,15,16)]
EGEOD121619.seurat_sub<- subset(EGEOD121619.seurat,cells=SubsetCells)
SteleDifExp<-FindMarkers(EGEOD121619.seurat_sub,group.by = "Factor.Value.temperature.",ident.1 = "22 degree celsius",ident.2 = "38 degree celsius",test.use = "t")
write.csv2(SteleDifExp,file="SteleDifExp.csv")


SubsetCells<- colnames(EGEOD121619.seurat)[EGEOD121619.seurat$ModScoreCluster%in% c(8)]
EGEOD121619.seurat_sub<- subset(EGEOD121619.seurat,cells=SubsetCells)
Cortex2DifExp<-FindMarkers(EGEOD121619.seurat_sub,group.by = "Factor.Value.temperature.",ident.1 = "22 degree celsius",ident.2 = "38 degree celsius",test.use = "t")
write.csv2(Cortex2DifExp,file="Cortex2DifExp.csv")



SubsetCells<- colnames(EGEOD121619.seurat)[EGEOD121619.seurat$ModScoreCluster%in% c(7,10)]
EGEOD121619.seurat_sub<- subset(EGEOD121619.seurat,cells=SubsetCells)
Cortex7_10DifExp<-FindMarkers(EGEOD121619.seurat_sub,group.by = "Factor.Value.temperature.",ident.1 = "22 degree celsius",ident.2 = "38 degree celsius",test.use = "t")
write.csv2(Cortex7_10DifExp,file="Cortex7_10DifExp.csv")



