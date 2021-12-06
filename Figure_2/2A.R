library(hdf5r)
library(dplyr)
library(ggplot2)
theme_set(theme_cowplot())
library(data.table)
annotations_sub=fread("~/2021-04-05.gtex_metadata_obs.csv")
theme_set(theme_cowplot())
library(ggbeeswarm)

annotations=annotations_sub%>%subset(!channel%in%c("skin_CST_GTEX-1CAMR","skin_TST_GTEX-15EOM","skin_TST_GTEX-1CAMR","breast_EZ_GTEX-1R9PN","skin_NST_GTEX-1CAMR"))

all_buffer=annotations%>%group_by(channel,`Broad cell type`)%>%summarise(num_cells=n())
all_tot=annotations%>%group_by(channel)%>%summarize(num_tot=n())

names(all_buffer)[2]="celltype"
all_buffer=all_buffer%>%left_join(all_tot,by="channel")
all_buffer=all_buffer%>%mutate(p=(num_cells/num_tot))
all_buffer_ent=all_buffer%>%group_by(channel)%>%summarise(entropy=(-sum(p*log2(p))))

df=annotations%>%select(channel,prep,tissue)
df=unique(df)
names(df)[1]="channel"
all_buffer_ent=all_buffer_ent%>%left_join(df,by="channel")

pdf(paste0("panelA.pdf"),w=14,10,useDingbats=FALSE)
all_buffer_ent%>%ggplot(aes(x=tissue,y=entropy,fill=prep))+geom_quasirandom(bandwidth=0.05,dodge.width=0.6,shape = 21,size=4,color="black")  +scale_fill_manual(values = c("#79a339","#a7cee2","#e1c7a7","#24467c")) +theme(axis.text.x = element_text(angle = 45,hjust=1,size=20)) +ylab("Shannon's entropy for celltype diversity")+xlab("")
dev.off()



all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%summarize(var=var(entropy))
# # A tibble: 8 x 2
# tissue                  var
# <chr>                 <dbl>
# 1 breast              0.585  
# 2 esophagusmucosa     0.0520 
# 3 esophagusmuscularis 0.0183 
# 4 heart               0.00785
# 5 lung                0.0268 
# 6 prostate            0.382  
# 7 skeletalmuscle      0.0963 
# 8 skin                0.243   
# > mean(0.0520,0.0183,0.00785,0.0268,0.382,0.0963)
# [1] 0.052

varlist = list(lung=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="lung"))$entropy,prostate=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="prostate"))$entropy,skin=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="skin"))$entropy,skeletalmuscle=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="skeletalmuscle"))$entropy,heart=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="heart"))$entropy,breast=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="breast"))$entropy,esophagusmucosa=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="esophagusmucosa"))$entropy,esophagusmuscularis=(all_buffer_ent%>%subset(prep%in%c("CST","TST","NST"))%>%group_by(tissue)%>%subset(tissue=="esophagusmuscularis"))$entropy)

sbp_list=list(skin=varlist$skin, breast=varlist$breast, prostate=varlist$prostate)
fligner.test(sbp_list)

rest_list=list(lung=varlist$lung, skeletalmuscle=varlist$skeletalmuscle, esophagusmuscularis=varlist$esophagusmuscularis,esophagusmucosa=varlist$esophagusmucosa,heart=varlist$heart)

two_list=list(sbp=sbp_list,rest=rest_list)

I= which(all_buffer_ent$tissue%in%c("breast","prostate","skin"))
group[I]="sbp"
group[-I]="rest"
all_buffer_ent_data=data.frame(all_buffer_ent,group=group)
fligner.test(entropy~group,data=all_buffer_ent_data)

# > fligner.test(entropy~group,data=all_buffer_ent_data)
# 
# Fligner-Killeen test of homogeneity of variances
# 
# data:  entropy by group
# Fligner-Killeen:med chi-squared = 3.4722, df = 1, p-value = 0.06241


all_buffer_ent$prep=factor(all_buffer_ent$prep,levels=c("CST","NST","TST","EZ"))

library(lme4)
library(lmerTest)

formula="entropy ~ prep + (1|tissue)"
m=lmer(formula,data=all_buffer_ent)
summary(m)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method [
# lmerModLmerTest]
# Formula: formula
#    Data: all_buffer_ent
# 
# REML criterion at convergence: 145
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.82942 -0.43705  0.05346  0.60323  2.13106 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  tissue   (Intercept) 0.2518   0.5018  
#  Residual             0.2213   0.4704  
# Number of obs: 90, groups:  tissue, 8
# 
# Fixed effects:
#             Estimate Std. Error       df t value        Pr(>|t|)    
# (Intercept)  2.39569    0.20195  9.82347  11.863 0.0000003840889 ***
# prepNST     -0.11433    0.13736 78.89176  -0.832           0.408    
# prepTST      0.03165    0.13747 78.94046   0.230           0.818    
# prepEZ      -1.08666    0.14281 78.95928  -7.609 0.0000000000502 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#         (Intr) prpNST prpTST
# prepNST -0.333              
# prepTST -0.330  0.488       
# prepEZ  -0.319  0.471  0.472


