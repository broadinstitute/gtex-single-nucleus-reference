
library(hdf5r)
library(Seurat)
library(ggplot2)
library(dplyr)
library(randomForest)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(ggbeeswarm)


annotations_sub=fread("~/2021-04-05.gtex_metadata_obs.csv")
annotations=annotations_sub%>%subset(!channel%in%c("skin_CST_GTEX-1CAMR","skin_TST_GTEX-15EOM","skin_TST_GTEX-1CAMR","breast_EZ_GTEX-1R9PN","skin_NST_GTEX-1CAMR"))


#### load prostate nuclei ####

prostate_nuc=annotations_sub%>%subset(tissue=="prostate")%>%subset(prep%in%c("TST"))
names(prostate_nuc)[1]="cellnames"
prostate_nuc$cellnames=gsub("-prostate","-1",prostate_nuc$cellnames)

prostate_nuc=prostate_nuc%>%select(cellnames,prep,`Broad cell type`,`Granular cell type`,`Tissue composition`,individual)

colnames(prostate_nuc)[3]="Broad_cell_type"
colnames(prostate_nuc)[4]="granular_cell_type"
colnames(prostate_nuc)[5]="compartment"

load("prostate_TST_sub_fpr1.RObj")
prostate_TST_sub_sub=subset(prostate_TST_sub,idents=levels(prostate_TST_sub@active.ident)[c(1:10,12)])
prostate=readRDS("prostate_cells.RDS")
prostate_sub=subset(prostate,ident=c("Epithelial cell (basal)","Epithelial cell (luminal)","Fibroblast","Vascular Endothelial","SM","Epithelial cell (Hillock)","Pericyte","Epithelial cell (club)","Macrophage","LEC","Neuroendocrine"))
prostate_sub<- FindVariableFeatures(prostate_sub, selection.method = "vst", nfeatures = 3000)

#######################################################################################################################################################
##################################################################2D PROSTATE##########################################################################
#######################################################################################################################################################

##### train classifier on TST nuc ####

training.set = c()
test.set=c()
training.label = c()
test.label=c()

#genes.use=intersect(prostate_TST_sub_sub@assays$RNA@var.features,prostate@assays$RNA@var.features)

genes.use=intersect(prostate_TST_sub_sub@assays$RNA@var.features,prostate_sub@assays$RNA@var.features) #1372

I=grep("^MT-",genes.use)
II=grep("^RPL",genes.use)
III=grep("^RPS",genes.use)
genes.use[c(I,II,III)]

genes.use=genes.use[-c(I,II,III)]

for (i in (levels(prostate_TST_sub_sub@active.ident))){
  cells.in.clust = WhichCells(prostate_TST_sub_sub,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(prostate_TST_sub_sub@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_prostate_nuc_01cb = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 


test.predict = predict(rf_prostate_nuc_01cb,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)

#### predict on cells ####


prostate_cells=prostate_sub

scprostate.rf = prostate_cells@assays$RNA@data[genes.use, ] #1383 genes
scprostate.rf = t(scale(t(scprostate.rf), center=TRUE, scale=TRUE))
scprostate.rf[is.na(scprostate.rf)] = 0
scprostate.ident = factor(prostate_cells@active.ident) 

scprostate.predict <- predict(rf_prostate_nuc_01cb,t(scprostate.rf))

Conf_test = table(scprostate.ident, scprostate.predict[names(scprostate.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))



Conf_test_M_sum=data.frame(scprostate.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="scprostate.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="snuccelltype"


Conf_test_M$snuccelltype=as.character(Conf_test_M$snuccelltype)



Conf_test_M$scprostate.ident <- factor(Conf_test_M$scprostate.ident, levels = rev(c("Epithelial cell (basal)","Epithelial cell (luminal)","Epithelial cell (club)","Epithelial cell (Hillock)","Fibroblast","Pericyte","SM","Vascular Endothelial","LEC","Macrophage","Neuroendocrine")))
Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = c("Epithelial cell (basal)","Epithelial cell (luminal)","Epithelial cell (club)","Epithelial cell (Hillock)","Fibroblast","Pericyte/SMC","Myocyte (smooth muscle)","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Neuroendocrine"))



pdf("2D",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=snuccelltype,y=scprostate.ident))+xlab("Nuclei/Predicted Labels")+ylab("Cells/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()




#########################################################################################################################################################
##################################################################S11H PROSTATE##########################################################################
#########################################################################################################################################################

##### train on cells ######

prostate_cells=prostate_sub
training.set = c()
test.set=c()
training.label = c()
test.label=c()


for (i in (levels(prostate_cells@active.ident))){
  cells.in.clust = WhichCells(prostate_cells,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(prostate_cells@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_prostate_cells_01cb = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_prostate_cells_01cb,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)


#### predict on nuclei ####

scprostate.rf = prostate_TST_sub_sub@assays$RNA@data[genes.use, ] #1246 genes
scprostate.rf = t(scale(t(scprostate.rf), center=TRUE, scale=TRUE))
scprostate.rf[is.na(scprostate.rf)] = 0
snuccelltype= factor(prostate_TST_sub_sub@active.ident) 

scprostate.predict <- predict(rf_prostate_cells_01cb,t(scprostate.rf))

Conf_test = table(snuccelltype, scprostate.predict[names(snuccelltype)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))



Conf_test_M_sum=data.frame(snuccelltype=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="snuccelltype")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="scprostate.ident"


Conf_test_M$scprostate.ident=as.character(Conf_test_M$scprostate.ident)



Conf_test_M$scprostate.ident <- factor(Conf_test_M$scprostate.ident, levels = c("Epithelial cell (basal)","Epithelial cell (luminal)","Epithelial cell (club)","Epithelial cell (Hillock)","Fibroblast","Pericyte","SM","Vascular Endothelial","LEC","Macrophage","Neuroendocrine"))


Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = rev(c("Epithelial cell (basal)","Epithelial cell (luminal)","Epithelial cell (club)","Epithelial cell (Hillock)","Fibroblast","Pericyte/SMC","Myocyte (smooth muscle)","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Neuroendocrine")))

pdf("S11H",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=scprostate.ident,y=snuccelltype))+ylab("Nuclei/Test set")+xlab("Cells/Predicted Labels")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()






