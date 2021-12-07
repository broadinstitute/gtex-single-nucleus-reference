

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


###### load cells ######
load("~/skin_c_cells.RObj")
skin_c_sub=subset(skin_c, idents=levels(skin_c@active.ident)[-7])
###### load nuclei######
load(paste0("~/skin_CSTTST_sub_fpr1.RObj"))

#### annotations ####
skin_nuc=annotations%>%subset(tissue=="skin")%>%subset(prep%in%c("CST","TST"))
names(skin_nuc)[1]="cellnames"
skin_nuc$cellnames=gsub("-skin","-1",skin_nuc$cellnames)
skin_nuc=skin_nuc%>%select(cellnames,prep,`Broad cell type`,`Granular cell type`,`Tissue composition`,individual)
colnames(skin_nuc)[3]="Broad_cell_type"
colnames(skin_nuc)[4]="granular_cell_type"
colnames(skin_nuc)[5]="compartment"

###################################################################################################################################################
##################################################################2C SKIN##########################################################################
###################################################################################################################################################

##### train classifier on CST nuc ####

skin_CST_sub_sub=subset(skin_CST_sub,ident=levels(skin_CST_sub@active.ident)[-c(2,4,6,7,13,15)])
skin_c_sub=subset(skin_c,ident=levels(skin_c@active.ident)[-c(1,7,11,12,16)])
skin_CST_sub_sub <- FindVariableFeatures(skin_CST_sub_sub, selection.method = "vst", nfeatures = 3000)
skin_c_sub <- FindVariableFeatures(skin_c_sub, selection.method = "vst", nfeatures = 3000)
skin_cells=skin_c_sub
training.set = c()
test.set=c()
training.label = c()
test.label=c()

genes.use=intersect(skin_c_sub@assays$RNA@var.features,skin_CST_sub_sub@assays$RNA@var.features)

I=grep("^MT-",genes.use)
II= grep("^RPS",genes.use)
III= grep("^RPL",genes.use)
genes.use=genes.use[-c(I,II,III)] #1075 

skin_common=genes.use
saveRDS(skin_common,paste0(work.d,"2021-01-18.rf_skin_commgenes.RDS"))

for (i in (levels(skin_CST_sub_sub@active.ident))){
  cells.in.clust = WhichCells(skin_CST_sub_sub,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(skin_CST_sub_sub@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_skin_nuc_01cb = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_skin_nuc_01cb,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)


#### predict on cells ####
skin_cells=skin_c_sub
scskin.rf = skin_cells@assays$RNA@data[genes.use, ]
scskin.rf = t(scale(t(scskin.rf), center=TRUE, scale=TRUE))
scskin.rf[is.na(scskin.rf)] = 0
scskin.ident = factor(skin_cells@active.ident) 

scskin.predict <- predict(rf_skin_nuc_01cb,t(scskin.rf))

Conf_test = table(scskin.ident, scskin.predict[names(scskin.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))



Conf_test_M_sum=data.frame(scskin.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="scskin.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="snuccelltype"


Conf_test_M$snuccelltype=as.character(Conf_test_M$snuccelltype)



Conf_test_M$scskin.ident <- factor(Conf_test_M$scskin.ident, levels = rev(c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Melanocyte","sweat gland","Fibroblast","Pericyte","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (Langerhans)","Immune (DC)","Immune (T cell)")))
Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Melanocyte","Sweat gland cell","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (Langerhans)","Immune (T cell)"))


pdf("2C.pdf",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=snuccelltype,y=scskin.ident))+xlab("Nuclei/Predicted Labels")+ylab("Cells/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()

###################################################################################################################################################
##################################################################S11F SKIN##########################################################################
###################################################################################################################################################

######### tRAIN on cells ########


skin_cells=skin_c_sub
training.set = c()
test.set=c()
training.label = c()
test.label=c()


for (i in (levels(skin_cells@active.ident))){
  cells.in.clust = WhichCells(skin_cells,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(skin_cells@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_skin_cells_01cb = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_skin_cells_01cb,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)


#### predict on nuclei ####

scskin.rf = skin_CST_sub_sub@assays$RNA@data[genes.use, ] #1246 genes
scskin.rf = t(scale(t(scskin.rf), center=TRUE, scale=TRUE))
scskin.rf[is.na(scskin.rf)] = 0
snuccelltype= factor(skin_CST_sub_sub@active.ident) 

scskin.predict <- predict(rf_skin_cells_01cb,t(scskin.rf))

Conf_test = table(snuccelltype, scskin.predict[names(snuccelltype)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))



Conf_test_M_sum=data.frame(snuccelltype=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="snuccelltype")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="scskin.ident"


Conf_test_M$scskin.ident=as.character(Conf_test_M$scskin.ident)


Conf_test_M$scskin.ident <- factor(Conf_test_M$scskin.ident, levels = c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Melanocyte","sweat gland","Fibroblast","Pericyte","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (Langerhans)","Immune (DC)","Immune (T cell)"))


Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = rev(c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Melanocyte","Sweat gland cell","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (Langerhans)","Immune (T cell)")))


pdf("S11F.pdf",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=scskin.ident,y=snuccelltype))+ylab("Nuclei/Test set")+xlab("Cells/Predicted Labels")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()




