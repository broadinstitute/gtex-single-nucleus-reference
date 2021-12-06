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

################################################################################################################################################
##################################################################S13O SKIN##########################################################################
################################################################################################################################################

###### load cells ######
load("~/skin_c_cells.RObj")
skin_c_sub=subset(skin_c, idents=levels(skin_c@active.ident)[-7])
###### load nuclei######
load(paste0("~/skin_CSTTST_sub_fpr1.RObj"))

skin_nuc=annotations_sub%>%subset(tissue=="skin")%>%subset(prep%in%c("CST","TST"))
names(skin_nuc)[1]="cellnames"
skin_nuc$cellnames=gsub("-skin","-1",skin_nuc$cellnames)
skin_nuc=skin_nuc%>%select(cellnames,prep,`Broad cell type`,`Granular cell type`,`Tissue composition`,individual)
colnames(skin_nuc)[3]="Broad_cell_type"
colnames(skin_nuc)[4]="granular_cell_type"
colnames(skin_nuc)[5]="compartment"
df=data.frame(cellnames=rownames(skin_CST_sub@meta.data))
df=df%>%left_join(skin_nuc,by="cellnames")
df_sub=df%>%select(cellnames,granular_cell_type)
colnames(df_sub)[2]="granular"
rownames(df_sub)=df_sub$cellnames
skin_CST_sub=AddMetaData(skin_CST_sub,metadata=df_sub)
Idents(skin_CST_sub) <- "granular"

###### training on nuclei #######

skin_CST_sub <- FindVariableFeatures(skin_CST_sub, selection.method = "vst", nfeatures = 3000)
skin_c <- FindVariableFeatures(skin_c, selection.method = "vst", nfeatures = 3000)
skin_cells=skin_c
training.set = c()
test.set=c()
training.label = c()
test.label=c()

genes.use=intersect(skin_c@assays$RNA@var.features,skin_CST_sub@assays$RNA@var.features)

I=grep("^MT-",genes.use)
II= grep("^RPS",genes.use)
III= grep("^RPL",genes.use)
genes.use=genes.use[-c(I,II,III)] #1101

for (i in (levels(skin_CST_sub@active.ident))){
  cells.in.clust = WhichCells(skin_CST_sub,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(skin_CST_sub@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_skin_nuc_g = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_skin_nuc_g,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)

######## predict on cells ###########

skin_c_subells=skin_c_sub

genes.use=rownames(rf_skin_nuc_g$importance)


scskin.rf = skin_c_subells@assays$RNA@data[genes.use, ]
scskin.rf = t(scale(t(scskin.rf), center=TRUE, scale=TRUE))
scskin.rf[is.na(scskin.rf)] = 0
scskin.ident = factor(skin_c_subells@active.ident) 

scskin.predict <- predict(rf_skin_nuc_g,t(scskin.rf))

Conf_test = table(scskin.ident, scskin.predict[names(scskin.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))



Conf_test_M_sum=data.frame(scskin.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="scskin.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="snuccelltype"


Conf_test_M$snuccelltype=as.character(Conf_test_M$snuccelltype)



Conf_test_M$scskin.ident <- factor(Conf_test_M$scskin.ident, levels = rev(c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Melanocyte","sweat gland","Fibroblast","Pericyte","Mesenchymal","SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (Langerhans)","Immune (DC)","Immune (T cell)","Preadipocyte/Fibroblast","Remove","Neuro?")))
Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Epithelial cell (cornified keratinocyte)","Epithelial cell (mature keratinocyte)","Melanocyte","Sweat gland cell","Sebaceous gland cell (CYP4F3+)","Sebaceous gland cell (KRT79 hi)","Sebaceous gland cell (SEC14L6+)","Fibroblast","Pericyte/SMC I","Pericyte/SMC II","Endothelial cell (vascular) I","Endothelial cell (vascular) II","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (Langerhans)","Immune (T cell)","Immune (mast cell)","Adipocyte","Unknown (ATP1B3-GJB6)","Unknown (SFN-GJA1)"))


pdf(paste0("S13O.pdf"),10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=snuccelltype,y=scskin.ident))+xlab("Nuclei/Predicted Labels")+ylab("Cells/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()



####################################################################################################################################################
##################################################################S13G PROSTATE##########################################################################
####################################################################################################################################################
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
prostate_sub=subset(prostate,idents=levels(prostate@active.ident)[-19])


df=data.frame(cellnames=rownames(prostate_TST_sub@meta.data))
df=df%>%left_join(prostate_nuc,by="cellnames")
df_sub=df%>%select(cellnames,granular_cell_type)
colnames(df_sub)[2]="granular"

rownames(df_sub)=df_sub$cellnames
prostate_TST_sub=AddMetaData(prostate_TST_sub,metadata=df_sub)

Idents(prostate_TST_sub) <- "granular"

###### training on nuclei #######

prostate_TST_sub <- FindVariableFeatures(prostate_TST_sub, selection.method = "vst", nfeatures = 3000)

prostate_sub <- FindVariableFeatures(prostate_sub, selection.method = "vst", nfeatures = 3000)
prostate_subells=prostate_sub
training.set = c()
test.set=c()
training.label = c()
test.label=c()

genes.use=intersect(prostate_sub@assays$RNA@var.features,prostate_TST_sub@assays$RNA@var.features)

I=grep("^MT-",genes.use)
II= grep("^RPS",genes.use)
III= grep("^RPL",genes.use)
genes.use=genes.use[-c(I,II,III)] #1374

for (i in (levels(prostate_TST_sub@active.ident))){
  cells.in.clust = WhichCells(prostate_TST_sub,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(prostate_TST_sub@assays$RNA@data[genes.use,])
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))

rf_prostate_nuc_g = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_prostate_nuc_g,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)


######## predict on cells ###########


genes.use=rownames(rf_prostate_nuc_g$importance)

prostate_sub_subells=prostate_sub
scprostate.rf = prostate_sub_subells@assays$RNA@data[genes.use, ]
scprostate.rf = t(scale(t(scprostate.rf), center=TRUE, scale=TRUE))
scprostate.rf[is.na(scprostate.rf)] = 0
scprostate.ident = factor(prostate_sub_subells@active.ident) 

scprostate.predict <- predict(rf_prostate_nuc_g,t(scprostate.rf))

Conf_test = table(scprostate.ident, scprostate.predict[names(scprostate.ident)])

Conf_test_M=melt(Conf_test,id=rownames(Conf_test))


names(Conf_test_M)[1]="scprostate.ident"

Conf_test_M$scprostate.ident=as.character(Conf_test_M$scprostate.ident)

I=which(Conf_test_M$scprostate.ident%in%c("FOXI1+","Stress_BE","Club_BE","SEMG1+","Mixed_Epithelial","Stress","Cycling"))
Conf_test_M$scprostate.ident[I]="Epithelial_Unknown"

Conf_test_M_agg=Conf_test_M%>%group_by(scprostate.ident,Var2)%>%summarize(value_agg=sum(value))

Conf_test_M_sum=data.frame(scprostate.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
I=which(Conf_test_M_sum$scprostate.ident%in%c("FOXI1+","Stress_BE","Club_BE","SEMG1+","Mixed_Epithelial","Stress","Cycling"))
Conf_test_M_sum$scprostate.ident[I]="Epithelial_Unknown"

Conf_test_M_sum_agg=Conf_test_M_sum%>%group_by(scprostate.ident)%>%summarize(sums_agg=sum(sums))
names(Conf_test_M_agg)[2]="snuccelltype"
Conf_test_M_agg=Conf_test_M_agg%>%left_join(Conf_test_M_sum_agg,by="scprostate.ident")
Conf_test_M_agg=Conf_test_M_agg%>%mutate(prop=(value_agg/sums_agg))


Conf_test_M_agg$snuccelltype=as.character(Conf_test_M_agg$snuccelltype)

Conf_test_M_agg$scprostate.ident <- factor(Conf_test_M_agg$scprostate.ident, levels = rev(c("Epithelial cell (basal)","Epithelial cell (luminal)","Stress_luminal","Epithelial cell (club)","Epithelial cell (Hillock)","Epithelial_Unknown","Fibroblast","Pericyte","SM","Vascular Endothelial","LEC","Macrophage","Other","Neuroendocrine")))


Conf_test_M_agg$snuccelltype <- factor(Conf_test_M_agg$snuccelltype, levels = c("Epithelial cell (basal I)", "Epithelial cell (basal II)","Epithelial cell (luminal)","Epithelial cell (club)","Epithelial cell (Hillock)","Fibroblast","Pericyte/SMC","Myocyte (smooth muscle)","Endothelial cell (vascular) I","Endothelial cell (vascular) II","Endothelial cell (vascular) III","Endothelial cell (lymphatic)","Immune (macrophage I)","Immune (macrophage II)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)","Immune (neutrophil)","Neuroendocrine","Schwann cell"))


pdf(paste0("S13G.pdf"),10,10,useDingbats=FALSE)
Conf_test_M_agg%>%ggplot(aes(x=snuccelltype,y=scprostate.ident))+xlab("Nuclei/Predicted Labels")+ylab("Cells/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()

####################################################################################################################################################
##################################################################S13W LUNG##########################################################################
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

###### training on nuclei granular cell types#######
df=data.frame(cellnames=rownames(lung_TST_sub@meta.data))
df=df%>%left_join(lung_nuc,by="cellnames")
df_sub=df%>%select(cellnames,granular_cell_type)
colnames(df_sub)[2]="granular"

rownames(df_sub)=df_sub$cellnames
lung_TST_sub=AddMetaData(lung_TST_sub,metadata=df_sub)

Idents(lung_TST_sub) <- "granular"

lung_TST_sub <- FindVariableFeatures(lung_TST_sub, selection.method = "vst", nfeatures = 3000)

lung_c <- FindVariableFeatures(lung_c, selection.method = "vst", nfeatures = 3000)

training.set = c()
test.set=c()
training.label = c()
test.label=c()

genes.use=intersect(lung@assays$originalexp@var.features,lung_TST_sub@assays$RNA@var.features)#1183

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

rf_lung_nuc_g = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
test.predict = predict(rf_lung_nuc_g,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)

###### granular prediction on cells ######
lung_sub_subells=lung
sclung.rf = lung_sub_subells@assays$originalexp@data[genes.use, ]
sclung.rf=as.matrix(sclung.rf)
sclung.rf = t(scale(t(sclung.rf), center=TRUE, scale=TRUE))

sclung.rf[is.na(sclung.rf)] = 0
Idents(lung_sub_subells) <- "Broad.cell.type"
sclung.ident = factor(lung_sub_subells@active.ident) 

sclung.predict <- predict(rf_lung_nuc_g,t(sclung.rf))

Conf_test = table(sclung.ident, sclung.predict[names(sclung.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))

names(Conf_test_M)[1]="sclung.ident"

Conf_test_M_sum=data.frame(sclung.ident=(rownames(Conf_test)),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="sclung.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="snuccelltype"



Conf_test_M$snuccelltype=as.character(Conf_test_M$snuccelltype)

Conf_test_M$sclung.ident=as.character(Conf_test_M$sclung.ident)

Conf_test_M$snuccelltype <- factor(Conf_test_M$snuccelltype, levels = c("Epithelial cell (alveolar type I)","Epithelial cell (alveolar type II)","Epithelial cell (basal)","Epithelial cell (club)","Epithelial cell (ciliated)","Fibroblast I","Fibroblast II","Pericyte/SMC","Endothelial cell (vascular) I","Endothelial cell (vascular) II","Endothelial cell (vascular) III","Endothelial cell (lymphatic)","Immune (macrophage)","Immune (macrophage activated)","Immune (alveolar macrophage)","Immune (alveolar macrophage activated)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)"))

Conf_test_M$sclung.ident <- factor(Conf_test_M$sclung.ident, levels = rev(c("Epithelial cell (alveolar type I)","Epithelial cell (alveolar type II)","Epithelial cell (ionocyte)","Serous","Epithelial cell (basal)","Epithelial cell (goblet)","Mucous","Epithelial cell (club)","Epithelial cell (ciliated)","Mesothelial","Fibroblast","Lipofibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (alveolar macrophage)","Immune (B cell)","Immune (NK cell)","Immune (T cell)","Immune (mast cell)","Plasmacytoid Dendritic","Platelet/Megakaryocyte","Neuroendocrine")))




pdf("S13W.pdf",10,10,useDingbats=FALSE)
Conf_test_M%>%ggplot(aes(x=snuccelltype,y=sclung.ident))+xlab("Nuclei Training set/Predicted Granular Labels")+ylab("Cells/Test set")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()



