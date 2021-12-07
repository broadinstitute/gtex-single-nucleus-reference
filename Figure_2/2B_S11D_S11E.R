library(hdf5r)
library(Seurat)
library(ggplot2)
library(dplyr)
library(randomForest)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(ggbeeswarm)
install.packages("ggh4x")
library(ggh4x)

annotations_sub=fread("~/2021-04-05.gtex_metadata_obs.csv")
annotations=annotations_sub%>%subset(!channel%in%c("skin_CST_GTEX-1CAMR","skin_TST_GTEX-15EOM","skin_TST_GTEX-1CAMR","breast_EZ_GTEX-1R9PN","skin_NST_GTEX-1CAMR"))


tcols=c("#603813","#A86D27","#ffc372","#a01b68","#416ab1","#ebd4d4","#8ac28e")

tissdiss=annotations%>%select(V1,tissue,channel,`Tissue composition`,prep)
names(tissdiss)[4]="tisscomp"
names(tissdiss)[1]="barcodes"
tissdiss$tisscomp=factor(tissdiss$tisscomp,levels=c("Adipose","Endothelial","Epithelial","Fibroblast","Immune","Muscle","Other"))


pdf("2B.pdf",w=15,30,useDingbats=FALSE)
tissdiss%>%ggplot(aes(x=channel,y=1,fill=tisscomp))+geom_col(position="fill")+scale_fill_manual(values = tcols)+facet_wrap(~tissue,scales="free_x") +theme(axis.text.x = element_text(angle = 45,hjust=1,size=20)) +ylab("")+xlab("")
dev.off()

tissdiss_num=tissdiss%>%group_by(channel,prep,tissue,tisscomp)%>%summarize(num_cells=n())
tissdiss_total=tissdiss%>%group_by(channel,prep,tissue)%>%summarize(num_total=n())

tissdiss_total_sub=tissdiss_total%>%subset(num_total>=30)

pdf("S11D.pdf",w=42,h=20,useDingbats=FALSE)
tissdiss%>%arrange(tissue,prep,channel)%>%ggplot(aes(x=channel,y=1,fill=tisscomp))+geom_bar(position="fill")+scale_fill_manual(values = tcols)+facet_nested(~tissue+prep,space="free",scales="free_x",nest_line = TRUE) +theme(axis.text.y = element_text(size=20),axis.text.x = element_text(angle = 45,hjust=1,size=20)) +ylab("")+xlab("")+ theme(legend.text=element_text(size=20),legend.title=element_blank(),strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 20))
dev.off()

pdf("S11E",w=30,h=10,useDingbats=FALSE)
tissdiss%>%arrange(tissue,prep)%>%ggplot(aes(x=prep,y=1,fill=tisscomp))+geom_bar(stat="identity",position="fill")+scale_fill_manual(values = tcols)+facet_nested(~tissue+prep,space="free",scales="free_x",nest_line = TRUE) +theme(axis.text.y = element_text(size=20),axis.text.x = element_text(angle = 45,hjust=1,size=20)) +ylab("")+xlab("")+ theme(legend.text=element_text(size=20),legend.title=element_blank(),strip.background =element_rect(fill="white"),strip.text.x = element_text(size = 20))
dev.off()
