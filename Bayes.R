library(pheatmap)
library(gplots)

#set working directory
setwd("~/Desktop/ABC_GCB_classifier/")
n<-50 #top n differentially expressed genes to be used in the predictor


################################################################################################

scale.normalization<-function(matrix)
{
  scale.norm<-log(t(apply(matrix,1,scale))+10)
  return(scale.norm)
}

p.subgroup<-function(pred,ABC.pred,GCB.pred)
{
  #  pred<-ABC.pred
  mean.ABC<-mean(ABC.pred)
  sd.ABC<-sd(ABC.pred)
  mean.GCB<-mean(GCB.pred)
  sd.GCB<-sd(GCB.pred)
  p.ABC<-dnorm(pred,mean.ABC,sd.ABC)/(dnorm(pred,mean.ABC,sd.ABC)+dnorm(pred,mean.GCB,sd.GCB))
  p.GCB<-1-p.ABC  
  p<-rbind(p.ABC,p.GCB)
  rownames(p)<-c("p.ABC","p.GCB")
  return(p)
}

###############################################################################################

## loading labeled microarray 
staudt<-read.table("staudt_dlbcl_pnas_2008.txt",row.names=1,sep="\t",header=T)
label<-read.table("staudt_dlbcl_pnas_2008.GCB_ABC_labels.txt",row.names=1,sep="\t",header=F)
colnames(staudt)<-substr(colnames(staudt),nchar(colnames(staudt))-11,nchar(colnames(staudt)))
rownames(label)<-substr(rownames(label),nchar(rownames(label))-11,nchar(rownames(label)))

## scale normalization of microarray
staudt.scale<-scale.normalization(staudt)
colnames(staudt.scale)<-colnames(staudt)

## divide train and test set
ABC.scale<-staudt.scale[,label=="ABC "]
ABC.train<-ABC.scale[,1:37]
ABC.test<-ABC.scale[,38:74]
GCB.scale<-staudt.scale[,label=="GCB "]
GCB.train<-GCB.scale[,1:36]
GCB.test<-GCB.scale[,37:71]

## select predictors
p.value.scale<-vector()
t.statistic.scale<-vector()
for(i in 1:dim(ABC.train)[1])
{
  #print(i)
  t.test.result<-t.test(ABC.train[i,],GCB.train[i,])
  p.value.scale[i]<-t.test.result$p.value 
  t.statistic.scale[i]<-t.test.result$statistic 
}
names(p.value.scale)<-rownames(ABC.train)
names(t.statistic.scale)<-rownames(ABC.train)
FDR<-p.adjust(p.value.scale,"BH")
t.statistic<-t.statistic.scale[rank(FDR)<=n]
ABC.pred<-t.statistic %*% ABC.train[rank(FDR)<=n,]
GCB.pred<-t.statistic %*% GCB.train[rank(FDR)<=n,]

## validate on train and test set
ABC.test.LPS<-t.statistic %*% ABC.test[rank(FDR)<=n,]
GCB.test.LPS<-t.statistic %*% GCB.test[rank(FDR)<=n,]
p.ABC.train<-p.subgroup(ABC.pred,ABC.pred,GCB.pred)
p.GCB.train<-p.subgroup(GCB.pred,ABC.pred,GCB.pred)
p.ABC.test<-p.subgroup(ABC.test.LPS,ABC.pred,GCB.pred)
p.GCB.test<-p.subgroup(GCB.test.LPS,ABC.pred,GCB.pred)

## load RNA-seq data (Please modify this part to accomodate with your dataset: DLBCL is the to-be-classified gene X sample matrix)
file<-read.csv('FPKM_DLBCL.csv')
file<-file[!duplicated(file[,2]),]
rownames(file)<-file[,2]
DLBCL<-file

## prepare predictors
pred<-names(FDR[rank(FDR)<=n])[names(FDR[rank(FDR)<=n]) %in% rownames(DLBCL)]
DLBCL.pred<-scale.normalization(DLBCL[pred,])
colnames(DLBCL.pred)<-colnames(DLBCL)
DLBCL.LPS<-t.statistic[pred] %*% DLBCL.pred
ABC.pred.test<-t.statistic[pred] %*% ABC.train[pred,]
GCB.pred.test<-t.statistic[pred] %*% GCB.train[pred,]

test.subgroup<-p.subgroup(DLBCL.LPS,ABC.pred.test,GCB.pred.test)
sort.DLBCL<-sort(test.subgroup[1,],decreasing=T)

## prepare heatmap annotation
color<-c('pink','blue','grey')
annotation_colors=list(ClassifiedLabel=c(ABC=color[1],GCB=color[2],unclassified=color[3]))
list<-c(rep("unclassified",dim(DLBCL)[2]))
list[test.subgroup[1,]>0.9]<-"ABC"
list[test.subgroup[2,]>0.9]<-"GCB"
names(list)<-colnames(test.subgroup)
annotation_col<-data.frame(ClassifiedLabel=list)

## plot heatmap 
pdf(paste0(getwd(),"/RNASeq-heatmap.pdf"),onefile=F,width=8,height=6)
pheatmap(as.matrix(DLBCL.pred[,names(sort.DLBCL)]),cluster_cols=F,col=greenred(75),border_color = NA,scale="row",fontsize_row=6,fontsize_col=5,annotation_col=annotation_col,annotation_colors = annotation_colors)
dev.off()

## plot probability
pdf(paste0(getwd(),"/probability-plot.pdf"),height=3,width=5)
plot(sort.DLBCL,axes=F,type="n",xlab="",ylab="Probability")
axis(1,labels=F)
axis(2)
lines(sort.DLBCL,col="dark red",lwd=3)
lines(1-sort.DLBCL,col="dark blue",lwd=3)
abline(h=0.9,lty=3,col="grey",lwd=3)
legend("right",legend=c("ABC Probability","GCB Probability"),col=c("dark red","dark blue"),lwd=c(3,3),cex=0.8)
dev.off()

## output classification label
write.table(list,file=paste0(getwd(),"/classification-label.txt"),sep="\t",col.names=F,quote=F)
