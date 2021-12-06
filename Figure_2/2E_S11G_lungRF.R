


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


####################################################################################################################################################
##################################################################2E LUNG##########################################################################
####################################################################################################################################################


#### load lung nuclei ####

lung_nuc=annotations_sub%>%subset(tissue=="lung")%>%subset(prep%in%c("TST"))
names(lung_nuc)[1]="cellnames"
lung_nuc$cellnames=gsub("-lung","-1",lung_nuc$cellnames)

lung_nuc=lung_nuc%>%select(cellnames,prep,`Broad cell type`,`Granular cell type`,`Tissue composition`,individual)

colnames(lung_nuc)[3]="Broad_cell_type"
colnames(lung_nuc)[4]="granular_cell_type"
colnames(lung_nuc)[5]="compartment"

load("~/lung_TST_sub_fpr1.RObj")

lung_c=readRDS("~/lung_cells.rds")


###### training on nuclei broad cell types#######

lung_TST_sub <- FindVariableFeatures(lung_TST_sub, selection.method = "vst", nfeatures = 3000)

lung_sub <- FindVariableFeatures(lung_sub, selection.method = "vst", nfeatures = 3000)
lung_subells=lung_sub
training.set = c()
test.set=c()
training.label = c()
test.label=c()

genes.use=intersect(lung_sub@assays$originalexp@var.features,lung_TST_sub@assays$RNA@var.features)#1227

I=grep("^MT-",genes.use)
II= grep("^RPS",genes.use)
III= grep("^RPL",genes.use)
#genes.use=genes.use[-c(I,II,III)] #1374

for (i in (levels(lung_TST_sub@active.ident))){
  cells.in.clust = WhichCells(lung_TST_sub,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(lung_TST_sub@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_lung_nuc_broad = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_lung_nuc_broad,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)


######## predict on cells ###########

genes.use=rownames(rf_lung_nuc_broad$importance)
lung_sub_subells=lung_sub
sclung.rf = lung_sub_subells@assays$originalexp@data[genes.use, ]
sclung.rf=as.matrix(sclung.rf)
sclung.rf = t(scale(t(sclung.rf), center=TRUE, scale=TRUE))
sclung.rf[is.na(sclung.rf)] = 0
Idents(lung_sub_subells) <- "Broad.cell.type"
sclung.ident = factor(lung_sub_subells@active.ident) 

sclung.predict <- predict(rf_lung_nuc_broad,t(sclung.rf))

Conf_test = table(sclung.ident, sclung.predict[names(sclung.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))

names(Conf_test_M)[1]="sclung.ident"

Conf_test_M_sum=data.frame(sclung.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="sclung.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="snuccelltype"


Conf_test_M$snuccelltype=as.character(Conf_test_M$snuccelltype)

Conf_test_M$sclung.ident=as.character(Conf_test_M$sclung.ident)

Conf_test_M$sclung.ident <- factor(Conf_test_M$sclung.ident, levels = rev(c("Epithelial cell (alveolar type I)","Epithelial cell (alveolar type II)","Epithelial cell (basal)","Epithelial cell (club)","Epithelial cell (ciliated)","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (alveolar macrophage)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)")))


Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = c("Epithelial cell (alveolar type I)", "Epithelial cell (alveolar type II)","Epithelial cell (basal)","Epithelial cell (club)","Epithelial cell (ciliated)","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (alveolar macrophage)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)"))




pdf("2E.pdf",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=snuccelltype,y=sclung.ident))+xlab("Nuclei/Predicted Labels")+ylab("Cells/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()

####################################################################################################################################################
##################################################################S11G LUNG#########################################################################
####################################################################################################################################################


######### build RF on cells ##########
training.set = c()
test.set=c()
training.label = c()
test.label=c()
for (i in (levels(lung_sub_subells@active.ident))){
  cells.in.clust = WhichCells(lung_sub_subells,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(lung_sub_subells@assays$originalexp@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_lung_cell_broad = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_lung_nuc_broad,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)



######## predict on nuclei ###########

lung_sub_subells=lung_TST_sub
snlung.rf = lung_sub_subells@assays$RNA@data[genes.use, ]
snlung.rf=as.matrix(snlung.rf)
snlung.rf = t(scale(t(snlung.rf), center=TRUE, scale=TRUE))
snlung.rf[is.na(snlung.rf)] = 0
snlung.ident = factor(lung_sub_subells@active.ident) 

snlung.predict <- predict(rf_lung_cell_broad,t(snlung.rf))

Conf_test = table(snlung.ident, snlung.predict[names(snlung.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))

names(Conf_test_M)[2]="sclung.ident"

Conf_test_M_sum=data.frame(snlung.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="snlung.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

Conf_test_M$sclung.ident=as.character(Conf_test_M$sclung.ident)

Conf_test_M$snlung.ident=as.character(Conf_test_M$snlung.ident)

Conf_test_M$sclung.ident <- factor(Conf_test_M$sclung.ident, levels = c("Epithelial cell (alveolar type I)","Epithelial cell (alveolar type II)","Epithelial cell (basal)","Epithelial cell (club)","Epithelial cell (ciliated)","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (alveolar macrophage)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)"))

Conf_test_M$snlung.ident <- factor(Conf_test_M$snlung.ident, levels = rev(c("Epithelial cell (alveolar type I)", "Epithelial cell (alveolar type II)","Epithelial cell (basal)","Epithelial cell (club)","Epithelial cell (ciliated)","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (alveolar macrophage)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)")))




pdf("S11G.pdf",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=sclung.ident,y=snlung.ident))+xlab("Cells training/Predicted Labels")+ylab("Nuclei/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()



