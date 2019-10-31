library(magrittr)
library(dplyr)

clusterFiles<-list.files("./Data/cluster/",full.names = TRUE)
IDs<-c(1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,8,9)
clusterFiles<-clusterFiles[order(IDs)]
temp<-read.csv("./Data/cluster/Cluster1.tsv",sep="\t",skip=1)
temp<-cbind(temp,cluster=1)
CellTypeAnnotation<-temp
for(i in 2:length(clusterFiles)){
temp<-read.csv(clusterFiles[i],sep="\t",skip=1)
temp<-cbind(temp,cluster=i)
colnames(temp)<-colnames(CellTypeAnnotation)
CellTypeAnnotation<-rbind(CellTypeAnnotation,temp)  
}

Celltypes<-c("Cluster 1: Non-hair cells I","Cluster 2: Non-hair cells II","Cluster 3: Non-hair cells III","Cluster 4: Non-hair cells IV","Cluster 5: Columella I","Cluster 6: Columella II",
             "Cluster 7: Cortex I","Cluster 8: Cortex II","Cluster 9: Cortex_Hair cells I","Cluster 10: Cortex_Hair cells II",
             "Cluster 11: Endodermis I","Cluster 12: Endodermis II","Cluster 13: Stele I","Cluster 14: Stele II","Cluster 15: Stele III","Cluster 16: Stele IV","Cluster 17: Hair cells")

Top3Makerlist<-CellTypeAnnotation  %>% group_by(cluster) %>% top_n(n = 4, wt = Average.Difference)
Top3Makerlist<-as.data.frame(Top3Makerlist)
for(i in 1:17){ 
  pdf(file=paste(Celltypes[[i]],".pdf",sep=""))
   print(FeaturePlot(EGEOD121619.seurat, features = as.character(Top3Makerlist[Top3Makerlist$cluster==i,1])))
   dev.off()
}

for (i in 1:17) EGEOD121619.integrated<-AddModuleScore(EGEOD121619.integrated,features = list(drop(as.character(CellTypeAnnotation[CellTypeAnnotation$cluster==i,1]))),name = paste("Cluster",i,"fullList",sep=""))
  
ModScores<-grep("fullList",colnames(EGEOD121619.integrated[[]]),value=TRUE)

for(i in 1:17){ 
  jpeg(filename = paste(Celltypes[[i]],"fullList.jpg",sep=""))
  print(FeaturePlot(EGEOD121619.integrated, features = ModScores[i]))
  dev.off()
}

ModScoreM<-EGEOD121619.integrated[[]][,ModScores]

 ModScoreM<-EGEOD121619.integrated[[]][,ModScores]
 ModScoreM<-t(t(ModScoreM)/apply(ModScoreM,2,function(v) quantile(v,.99)))
 boxplot(ModScoreM)
ModScoreCluster<- apply(ModScoreM,1,which.max)
EGEOD121619.seurat<-AddMetaData(EGEOD121619.seurat,ModScoreCluster,col.name = "ModScoreCluster")
EGEOD121619.integrated<-AddMetaData(EGEOD121619.integrated,ModScoreCluster,col.name = "ModScoreCluster")
 




