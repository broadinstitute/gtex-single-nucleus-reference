
#######################################################################################################################################
############################################# DISSOCIATION SCORES #############################################
#######################################################################################################################################
diss_score=NULL

library(GSA)
x=GSA.read.gmt("~/vanOudenaarden_dissociation_human.gmt")

#### lung ####

load("~/lung_TST_sub_fpr1.RObj")
lung_c=readRDS("~/lung_cells.rds")


lung_cell=AddModuleScore(lung_c,features=list(x$genesets[[1]]),name="DISS")
lung_nuc=AddModuleScore(lung_TST_sub,features=list(x$genesets[[1]]),name="DISS")

dfi=data.frame(group="frozen",celltype=lung_nuc@active.ident,score=lung_nuc@meta.data$DISS1)
dfi=rbind(dfi,data.frame(group="fresh",celltype=lung_cell@active.ident,score=lung_cell@meta.data$DISS1))

df$celltype=as.character(df$celltype)

I=which(df$celltype=="AT1")
df$celltype[I]="Epithelial cell (alveolar type I)"

I=which(df$celltype=="AT2")
df$celltype[I]="Epithelial cell (alveolar type II)"

I=which(df$celltype=="Club")
df$celltype[I]="Epithelial cell (club)"

I=which(df$celltype=="Mast")
df$celltype[I]="Immune (mast cell)"

I=which(df$celltype=="Ciliated")
df$celltype[I]="Epithelial cell (ciliated)"

I=which(df$celltype=="Immune (NK)")
df$celltype[I]="Immune (NK cell)"

I=which(df$celltype=="Immune (macrophage)")
df$celltype[I]="Immune (DC/macrophage)"

I=which(df$celltype=="Basal")
df$celltype[I]="Epithelial cell (basal)"

I=which(df$celltype=="B.cells")
df$celltype[I]="Immune (B cell)"

I=which(df$celltype=="Pericytes")
df$celltype[I]="Pericyte/SMC"

I=which(df$celltype=="LEC")
df$celltype[I]="Endothelial cell (lymphatic)"

I=which(df$celltype=="Endothelial")
df$celltype[I]="Endothelial cell (vascular)"

dfi_sub=dfi%>%subset(celltype%in%unique(lung_nuc@active.ident))


pdf("S14B.pdf",w=15,h=10,useDingbats=FALSE)
dfi_sub%>%ggplot(aes(x=celltype,y=score,fill=group))+geom_boxplot()+scale_fill_manual(values=c("#0072ce","#dc9cbf"))+theme(axis.text.x = element_text(angle=90,hjust=1))
dev.off()

df$celltype=factor(df$celltype,levels=c("Epithelial cell (alveolar type I)","Epithelial cell (alveolar type II)","Epithelial cell (ciliated)","Epithelial cell (basal)","Epithelial cell (club)","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (alveolar macrophage)","Immune (mast cell)","Immune (T cell)","Immune (NK cell)","Immune (B cell)"))

dfi_sub$celltype=factor(dfi_sub$celltype,levels=unique(dfi_sub$celltype))

lung_pv=NULL
for (celltype in levels(dfi_sub$celltype)){
  print(celltype)
  dat=dfi_sub%>%subset(celltype=celltype)
  pv=wilcox.test(dat$score ~ dat$group)$p.value
  lung_pv=rbind(lung_pv,data.frame(tissue="lung",celltype=celltype,pval=pv))
}


#### prostate ####

prostate=readRDS("~/prostate_cells.RDS")
table(prostate@active.ident)
prostate_sub=subset(prostate,ident=c("Epithelial cell (basal)","Epithelial cell (luminal)","Fibroblast","Vascular Endothelial","SM","Epithelial cell (Hillock)","Pericyte","Epithelial cell (club)","Macrophage","LEC","Neuroendocrine"))


load("~/prostate_TST_sub_fpr1.RObj")
prostate_TST_sub_sub=subset(prostate_TST_sub,idents=levels(prostate_TST_sub@active.ident)[c(1:10,12)])


prostate_cell=AddModuleScore(prostate_sub,features=list(x$genesets[[1]]),name="DISS")
prostate_nuc=AddModuleScore(prostate_TST_sub_sub,features=list(x$genesets[[1]]),name="DISS")

df=data.frame(group="frozen",celltype=prostate_nuc@active.ident,score=prostate_nuc@meta.data$DISS1)
df=rbind(df,data.frame(group="fresh",celltype=prostate_cell@active.ident,score=prostate_cell@meta.data$DISS1))

df$celltype=as.character(df$celltype)

I=which(df$celltype=="Luminal Epithelial")
df$celltype[I]="Epithelial cell (luminal)"

I=which(df$celltype=="Club Epithelial")
df$celltype[I]="Epithelial cell (club)"

I=which(df$celltype=="Hillock Epithelial")
df$celltype[I]="Epithelial cell (Hillock)"

I=which(df$celltype=="LEC")
df$celltype[I]="Endothelial cell (lymphatic)"

I=which(df$celltype=="Basal Epithelia")
df$celltype[I]="Epithelial cell (basal)"

I=which(df$celltype=="SM")
df$celltype[I]="Myocyte (smooth muscle)"

I=which(df$celltype=="Vascular Endothelial")
df$celltype[I]="Endothelial cell (vascular)"

I=which(df$celltype=="Macrophage")
df$celltype[I]="Immune (DC/macrophage)"

I=which(df$celltype=="Macrophage")
df$celltype[I]="Immune (DC/macrophage)"

I=which(df$celltype=="Pericyte")
df$celltype[I]="Pericyte/SMC"

df$celltype=factor(df$celltype,levels=c("Epithelial cell (basal)","Epithelial cell (luminal)","Epithelial cell (club)","Epithelial cell (Hillock)","Fibroblast","Pericyte/SMC","Myocyte (smooth muscle)","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Neuroendocrine"))

pdf("S14A.pdf",w=15,h=10,useDingbats=FALSE)
df%>%ggplot(aes(x=celltype,y=score,fill=group))+geom_boxplot()+scale_fill_manual(values=c("#0072ce","#dc9cbf"))+theme(axis.text.x = element_text(angle=90,hjust=1))
dev.off()




prostate_pv=NULL
for (celltype in levels(df$celltype)){
  print(celltype)
  dat=df%>%subset(celltype=celltype)
  pv=wilcox.test(dat$score ~ dat$group)$p.value
  prostate_pv=rbind(prostate_pv,data.frame(tissue="prostate",celltype=celltype,pval=pv))
}



#### skin ####

load("~/skin_c_cells.RObj")
load(paste0("~/skin_CSTTST_sub_fpr1.RObj"))

skin_cell=AddModuleScore(skin_c_sub,features=list(x$genesets[[1]]),name="DISS")
skin_nuc=AddModuleScore(skin_CST_sub_sub,features=list(x$genesets[[1]]),name="DISS")

df=data.frame(group="frozen",celltype=skin_nuc@active.ident,score=skin_nuc@meta.data$DISS1)
df=rbind(df,data.frame(group="fresh",celltype=skin_cell@active.ident,score=skin_cell@meta.data$DISS1))

df$celltype=as.character(df$celltype)

I=which(df$celltype=="sweat gland")
df$celltype[I]="Sweat gland cell"

I=which(df$celltype=="Pericyte")
df$celltype[I]="Pericyte/SMC"

I=which(df$celltype=="Immune (DC)")
df$celltype[I]="Immune (DC/macrophage)"

df$celltype=factor(df$celltype,levels=c("Epithelial cell (basal keratinocyte)","Epithelial cell (suprabasal keratinocyte)","Melanocyte","Sweat gland cell","Fibroblast","Pericyte/SMC","Endothelial cell (vascular)","Endothelial cell (lymphatic)","Immune (DC/macrophage)","Immune (Langerhans)","Immune (T cell)"))

pdf("2F.pdf",w=15,h=10,useDingbats=FALSE)
df%>%ggplot(aes(x=celltype,y=score,fill=group))+geom_boxplot()+scale_fill_manual(values=c("#0072ce","#dc9cbf"))+theme(axis.text.x = element_text(angle=90,hjust=1))
dev.off()


skin_pv=NULL
for (celltype in levels(df$celltype)){
  print(celltype)
  dat=df%>%subset(celltype=celltype)
  pv=wilcox.test(dat$score ~ dat$group)$p.value
  skin_pv=rbind(skin_pv,data.frame(tissue="skin",celltype=celltype,pval=pv))
}

diss_score=NULL
diss_score=rbind(diss_score,data.frame(tissue="lung",df))
diss_score=rbind(diss_score,data.frame(tissue="skin",df))
diss_score=rbind(diss_score,data.frame(tissue="prostate",df))

tissue_pv=NULL
tissue_pv=rbind(tissue_pv,lung_pv)
tissue_pv=rbind(tissue_pv,skin_pv)
tissue_pv=rbind(tissue_pv,prostate_pv)


