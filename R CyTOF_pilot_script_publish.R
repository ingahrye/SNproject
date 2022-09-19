
library(ggplot2)
library(reshape2)
library(gsubfn)
library(randomForest)
library(preprocessCore)
library(e1071)
library(plotrix)
library(spatstat)
library("splines")
library("survival")
library("survminer")
library("KMsurv")
library("plyr")

setwd("C:/Users/ingah/Dropbox/CyTOF_pilot/Rscripts_files/FilesToGitHub_publication")

###load files
fracSub=read.table("../FilesToGitHub_publication/Cellpopulations_abundance.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")

Treg =read.table("../FilesToGitHub_publication/ExpressionCell_median_Treg_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD4=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD4_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD8=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD8_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD8RA=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD8RA_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD8RO=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD8RO_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD4RO=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD4RO_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD4RA=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD4RA_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
NK=read.table("../FilesToGitHub_publication/ExpressionCell_median_NKcells_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
TFH=read.table("../FilesToGitHub_publication/ExpressionCell_median_TFH_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
DNT=read.table("../FilesToGitHub_publication/ExpressionCell_median_DNTcell_Ab.txt",header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
APC=read.table("../FilesToGitHub_publication/ExpressionCell_median_APC_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
NaiveB=read.table("../FilesToGitHub_publication/ExpressionCell_median_naiveBcell_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
MemB=read.table("../FilesToGitHub_publication/ExpressionCell_median_memBcell_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
Plasmablast=read.table("../FilesToGitHub_publication/ExpressionCell_median_plasma_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
GCB=read.table("../FilesToGitHub_publication/ExpressionCell_median_GCB_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
NKdim=read.table("../FilesToGitHub_publication/ExpressionCell_NKdim_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
NKbright=read.table("../FilesToGitHub_publication/ExpressionCell_median_NKbright_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
TCRgd=read.table("../FilesToGitHub_publication/ExpressionCell_median_TCRgd_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD8RAROdp=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD8RAROdp_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CD4CD8dp=read.table("../FilesToGitHub_publication/ExpressionCell_median_CD4CD8dp_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
Tcells=read.table("../FilesToGitHub_publication/ExpressionCell_median_Tcells_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
Bcells=read.table("../FilesToGitHub_publication/ExpressionCell_median_Bcells_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
NonBT=read.table("../FilesToGitHub_publication/ExpressionCell_median_nonBTcells_Ab.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")


####tumor cells 
TClivemean=read.table("../FilesToGitHub_publication/Pilot_mean_TClive_n52.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
TCdeadmean=read.table("../FilesToGitHub_publication/Pilot_mean_TCdead_n52.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
TClivemed=read.table("../FilesToGitHub_publication/Pilot_median_TClive_n52.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
TCdeadmed=read.table("../FilesToGitHub_publication/Pilot_median_TCdead_n52.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")


clin=read.table("../FilesToGitHub_publication/ClinicalFile.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
clin$ID.int=as.integer(clin$ID)   
clin$ID.chr=clin$ID   

clin$TC=rep(NA, nrow(clin))
clin$TC[clin$Tumor.cells=="Less than 50"]="tc"
clin$TC[clin$Tumor.cells=="With tumor"]="tc"
clin$TC[clin$Tumor.cells=="only dead tumor and few"]="tc"
clin$TC[clin$Tumor.cells=="No tumor"]="no tc"
clin$TC<-factor(clin$TC, levels=c("tc", "no tc"))

clin$TC_group=rep(NA, nrow(clin))
clin$TC_group[clin$Tumor.cells=="Less than 50"]="deadLess50"
clin$TC_group[clin$Tumor.cells=="With tumor"]="More50"
clin$TC_group[clin$Tumor.cells=="only dead tumor and few"]="deadLess50"
clin$TC_group[clin$Tumor.cells=="No tumor"]="noTC"
clin$TC_group<-factor(clin$TC_group, levels=c("noTC","deadLess50", "More50"))

clin$TC_CyTOF=rep(NA, nrow(clin))
clin$TC_CyTOF[clin$TC_CyTOF_percent<0.02]="TCneg"
clin$TC_CyTOF[clin$TC_CyTOF_percent>=0.02]="TCpos"
clin$TC_CyTOF<-factor(clin$TC_CyTOF, levels=c("TCneg", "TCpos"))

##split the samples into groups according to TC/LC ratio
clin$TC_CyTOF_split=rep(NA, nrow(clin))
clin$TC_CyTOF_split[clin$TC_CyTOF_percent<0.02]="TCneg"
clin$TC_CyTOF_split[clin$TC_CyTOF_percent>=0.02&clin$TC_CyTOF_percent<=0.13]="TClow"
clin$TC_CyTOF_split[clin$TC_CyTOF_percent>0.13]="TChigh"
clin$TC_CyTOF_split<-factor(clin$TC_CyTOF_split, levels=c("TCneg", "TClow", "TChigh"))

clin$ERper_pat_split=rep(NA, nrow(clin))
clin$ERper_pat_split[clin$ER_per_pat<1]="ERneg"
clin$ERper_pat_split[clin$ER_per_pat>=1&clin$ER_per_pat<=10]="ERlow"
clin$ERper_pat_split[clin$ER_per_pat>=10&clin$ER_per_pat<50]="ERint"
clin$ERper_pat_split[clin$ER_per_pat>=50]="ERhigh"
clin$ERper_pat_split<-factor(clin$ERper_pat_split, levels=c("ERneg", "ERlow","ERint", "ERhigh"))

clin$KI67per_pat_split=rep(NA, nrow(clin))
clin$KI67per_pat_split[clin$Ki67_per_pat<=30]="KI67low"
clin$KI67per_pat_split[clin$Ki67_per_pat>30&clin$Ki67_per_pat<=70]="KI67med"
clin$KI67per_pat_split[clin$Ki67_per_pat>70]="KI67high"
clin$KI67per_pat_split<-factor(clin$KI67per_pat_split, levels=c("KI67low", "KI67med","KI67high"))

clin$pNed2=rep(NA, nrow(clin))
clin$pNed2[clin$pNed==0]="neg"
clin$pNed2[clin$pNed=="mi"]="mi"
clin$pNed2[clin$pNed>=1 & clin$pNed<=3]="macro"
clin$pNed2<-factor(clin$pNed2, levels=c("neg", "mi","macro"))


clin$metStatus=rep(NA, nrow(clin))
clin$metStatus[clin$metastaticStatus==0]="no"
clin$metStatus[clin$metastaticStatus==1]="yes"
clin$metStatus<-factor(clin$metStatus, levels=c("no", "yes"))



clin$PT_size2=rep(NA, nrow(clin))
clin$PT_size2[clin$PT_size_pat<18]="<18"
clin$PT_size2[clin$PT_size_pat>17&clin$PT_size_pat<30]="18-30"
clin$PT_size2[clin$PT_size_pat>=30 & clin$PT_size_pat<=60]="30-60"
clin$PT_size2<-factor(clin$PT_size2, levels=c("<18", "18-30","30-60"))




fracSub$CD4RORA=(fracSub$CD4RO)/(fracSub$CD4RA)
fracSub$CD8RORA=(fracSub$CD8RO)/(fracSub$CD8RA)

data=rbind(cbind(fracSub, subfrac="fracSub"))
           
sNames=colnames(data[c(2:28)]) #
data_all=melt(data, measure.vars=sNames)

data_all$ID_int=as.integer(data_all$ID)


dataClin=merge(data_all, clin, by.x="ID", by.y="ID.chr", all.x=T)

dataClin$type<-  factor(dataClin$type , levels = c("Sn", "SNmet", "LKmet", "ctr"))
dataClin$variable<-  factor(dataClin$variable , levels = c("Bcell", "Tcell","Non_BT","APC","NK","GCBcell","Plasmacells",
                                                           "naiveBcell","memBcell","CD4","CD8","CD4CD8","CD4CD8dp","dnTcell","gdTcell",
                                                           "Treg","TFH","CD4RO","CD4RA","CD4RORA","CD4RAROdp","CD8RA","CD8RO",
                                                           "CD8RORA","CD8RAROdp", "Nkdim","Nkbright"))

dataClin$run<-  factor(dataClin$run , levels = c( "run1","run2","run3","run4","run5","run6","run7","run8","run9","run10", "run11"))
dataClin$pNed<-  factor(dataClin$pNed , levels = c( "0","mi","1","2","3"))
dataClin$value=as.integer(dataClin$value)

#Figure: abundance cell populations all samples
ggplot(data=subset(dataClin, type%in%c("LKmet", "Sn","SNmet")&variable!="NA"), aes(ID,value)) +
  geom_jitter(position=position_jitter(width=0.1, height = .2),
              size=0.5) + ggtitle("abundance cell populations all samples")+
  facet_wrap(~variable, scales="fixed") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5))+theme_bw()
  


##figure: abundance cell populations split on axillary status
ggplot(data=subset(dataClin, type%in%c("LKmet", "Sn","SNmet")&variable!="NA"), aes(type,value)) +
  geom_boxplot(outlier.shape=NA,aes(fill=type),position = position_dodge(width=0.9)) + scale_fill_brewer(palette="Greens")+
  geom_jitter(aes(fill=type),position=position_jitter(width=0.1, height = .2),
              size=1) + ggtitle("figure: abundance cell populations split on axillary status")+
facet_wrap(~variable, scales="fixed") +
theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5)) + theme_bw()


#figure 1_abundance Split on TC cytof
ggplot(data=subset(dataClin, type%in%c("LKmet", "Sn","SNmet")&subfrac=="fracSub"), aes(TC_CyTOF_split,value)) +
  geom_boxplot(outlier.shape=NA,aes(fill=TC_CyTOF_split),position = position_dodge(width=0.9)) + scale_fill_brewer(palette="Reds")+
  geom_jitter(aes(fill=TC_CyTOF_split),position=position_jitter(width=0.1, height = .2),
              size=1) + ggtitle("all CD45 abundance( fraction of subtype")+
  facet_wrap(~variable, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5)) + theme_bw()


#### wilcox fracsub=subfrac three groups

N.wilcox = length(unique(dataClin$variable))
data.wilcox = subset(dataClin,type!="ctr" &variable!="NA")[, c("variable", "type", "value")]
res.wilcox = data.frame(
  variable    = unique(data.wilcox$variable), 
  Sn.SNmet    = rep(NA, N.wilcox),
  Sn.LKmet    = rep(NA, N.wilcox),
  SNmet.LKmet = rep(NA, N.wilcox),
  stringsAsFactors = F)

for(i in 1:N.wilcox) {
  ct = res.wilcox$variable[i]
  d.Sn    = subset(data.wilcox, variable==ct & type=="Sn")$value
  d.SNmet = subset(data.wilcox, variable==ct & type=="SNmet")$value
  d.LKmet = subset(data.wilcox, variable==ct & type=="LKmet")$value
  res.wilcox$Sn.SNmet[i] = wilcox.test(d.Sn, d.SNmet)$p.value
  res.wilcox$Sn.LKmet[i] = wilcox.test(d.Sn, d.LKmet)$p.value
  res.wilcox$SNmet.LKmet[i] = wilcox.test(d.SNmet, d.LKmet)$p.value
  
}


gplots::heatmap.2(-log10(res.wilcox.mat), trace="none", margins=c(15,15) )
write.table(res.wilcox, "wilcox.pairwise.txt", sep="\t")

##by bonferronix40, memory and Naive B cells are significantly different between grade 2 and 3.
#memory Bcells are sig lower in grade 3 tumors compared to grade2


##### split by ER status (pathology) neg and high(>50%)

#kruskal-Wallis test (tester tre grupper, nb; bonferroni på 3 grupper)

> kruskal.test(value ~ variable, data=data.wilcox)$p.value




##alternative regression plot + coefficients
regressions = data.frame()
for (v in levels(dataClin$variable)) {
  cytof = subset(dataClin, 
                 type %in% c("LKmet", "Sn","SNmet") &
                   variable == v)
  
  lm.r=lm(value~log2(TC_CyTOF_percent_corrigated), data=cytof)
  cor.coeff = cor(cytof$value, log2(cytof$TC_CyTOF_percent_corrigated),method="spearman") #legg til method="spearman"
  
  regressions = rbind(regressions, data.frame(
    variable    = v,
    intercept   = summary(lm.r)$coefficient[1,1],
    slope       = summary(lm.r)$coefficient[2,1],
    correlation = cor.coeff,
    slope.p     = summary(lm.r)$coefficient[2,4]  
  ))
}

regressions = regressions[order(regressions$slope.p), ]
regressions$description=sprintf("%s (cor=%0.2f, p=%0.7f)", 
                                regressions$variable, regressions$correlation, regressions$slope.p)
annotation = regressions[, c("variable", "description")]
dataClin2 = merge(dataClin, annotation, by.x="variable", by.y="variable")
dataClin2$description = ordered(dataClin$description, regressions$description)

##figure regression (min 1149 SNmet)
ggplot(data=subset(dataClin2, type%in%c("LKmet", "Sn","SNmet")&subfrac=="fracSub"), aes(log2(TC_CyTOF_percent_corrigated),value)) +
  geom_point(aes(color=TC_CyTOF, size=1),position = position_dodge(width=0.9),size=1) +
  ggtitle("Regression_TC% vs subfrac")+geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~description, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=12, vjust=0.5))+ theme_bw()



####################################
##
ggplot(data=subset(dataClin, type%in%c("LKmet", "Sn","SNmet")&subfrac=="fracSub"), aes(type,value)) +
  geom_boxplot(outlier.shape=NA,aes(fill=type),position = position_dodge(width=0.9)) + scale_fill_brewer(palette="Dark2")+
  geom_point(aes(fill=type),position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.1),
             size=2) + ggtitle("overview all markers subfraction")+
  facet_wrap(~variable, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5))+theme_bw()

N.wilcox = length(unique(dataClin$variable))
data.wilcox = subset(dataClin, subfrac=="fracSub" &type!="ctr")[, c("variable", "type", "value")]
res.wilcox = data.frame(
  variable    = unique(data.wilcox$variable), 
  Sn.SNmet    = rep(NA, N.wilcox),
  Sn.LKmet    = rep(NA, N.wilcox),
  SNmet.LKmet = rep(NA, N.wilcox),
  stringsAsFactors = F)

for(i in 1:N.wilcox) {
  ct = res.wilcox$variable[i]
  d.Sn    = subset(data.wilcox, variable==ct & type=="Sn")$value
  d.SNmet = subset(data.wilcox, variable==ct & type=="SNmet")$value
  d.LKmet = subset(data.wilcox, variable==ct & type=="LKmet")$value
  res.wilcox$Sn.SNmet[i] = wilcox.test(d.Sn, d.SNmet)$p.value
  res.wilcox$Sn.LKmet[i] = wilcox.test(d.Sn, d.LKmet)$p.value
  res.wilcox$SNmet.LKmet[i] = wilcox.test(d.SNmet, d.LKmet)$p.value
  
}

gplots::heatmap.2(-log10(res.wilcox.mat), trace="none", margins=c(15,15) )
write.table(res.wilcox, "wilcox.pairwise.txt", sep="\t")


                                               
                                                        
######################################################################################################


####expression Subpop_merge all subPop

data_expSubPop=rbind(cbind(APC, celltype="APC"),
                     cbind(Bcells, celltype="Bcells"),
                     cbind(CD4, celltype="CD4"),
                     cbind(CD4CD8dp, celltype="CD4CD8dp"),
                     cbind(CD4RA, celltype="CD4RA"),
                     cbind(CD4RO, celltype="CD4RO"),
                     cbind(CD8, celltype="CD8"),
                     cbind(CD8RA, celltype="CD8RA"),
                     cbind(CD8RAROdp, celltype="CD8RAROdp"),
                     cbind(CD8RO, celltype="CD8RO"),
                     cbind(DNT, celltype="DNT"),
                     cbind(TCRgd, celltype="TCRgd"),          #check this file
                     cbind(MemB, celltype="memB"),
                     cbind(NaiveB, celltype="naiveB"),
                     cbind(NK, celltype="NK"),
                     cbind(NKbright, celltype="NKbright"),
                     cbind(NKdim, celltype="NKdim"),
                     cbind(NonBT, celltype="NonBT"),
                     cbind(Plasmablast, celltype="plasma"),
                     cbind(Tcells, celltype="Tcell"),
                     cbind(TFH, celltype="TFH"),
                     cbind(Treg, celltype="Treg"))
        


###arcsine transform
data_expSubPopAh<-data_expSubPop[,c(-1,-39)]/5 #remove ID and celltype
data_expSubPopAh<-asinh(data_expSubPopAh)
data_expSubPopAh$ID<-data_expSubPop$ID
data_expSubPopAh$celltype<-data_expSubPop$celltype
##

sNames=colnames(data_expSubPopAh[c(1:37)]) #
data_expSubPop_all=melt(data_expSubPopAh, measure.vars=sNames)

ExprClin=merge(data_expSubPop_all, clin, by.x="ID", by.y="ID.chr", all.x=T)

ExprClin$type<-  factor(ExprClin$type , levels = c("Sn", "SNmet", "LKmet", "ctr"))


################################regression for tigit on the different subpopulations
regressions = data.frame()
for (v in levels(ExprClin$celltype)) {
  cytof = subset(ExprClin, 
                 type %in% c("LKmet", "Sn","SNmet") &
                   variable=="TIGIT"&
                   celltype == v)
  
  lm.r=lm(value~log2(TC_CyTOF_percent_corrigated), data=cytof)
  cor.coeff = cor(cytof$value, log2(cytof$TC_CyTOF_percent_corrigated),method="spearman") #legg til method="spearman"
  
  regressions = rbind(regressions, data.frame(
    celltype    = v,
    intercept   = summary(lm.r)$coefficient[1,1],
    slope       = summary(lm.r)$coefficient[2,1],
    correlation = cor.coeff,
    slope.p     = summary(lm.r)$coefficient[2,4]  
  ))
}

regressions = regressions[order(regressions$slope.p), ]
regressions$description=sprintf("%s (cor=%0.2f, p=%0.7f)", 
                                regressions$celltype, regressions$correlation, regressions$slope.p)
annotation = regressions[, c("celltype", "description")]
ExprClin = merge(ExprClin, annotation, by.x="celltype", by.y="celltype")
ExprClin$description = ordered(ExprClin$description, regressions$description)

##figure regression (min 1149 SNmet)
ggplot(data=subset(ExprClin2, type%in%c("LKmet", "Sn","SNmet")&variable=="TIGIT"), aes(log2(TC_CyTOF_percent_corrigated),value)) +
  geom_point(aes(color=TC_CyTOF, size=1),position = position_dodge(width=0.9),size=1) +
  ggtitle("Regression_TC% vs subfrac")+geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(celltype~description, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=12, vjust=0.5))+ theme_bw()

################################################################################################



################# overview all subtypes and all markers
ggplot(data=subset(ExprClin, type%in%c("LKmet", "Sn","SNmet")), aes(type,value)) +
  geom_boxplot(outlier.shape=NA,aes(fill=type),position = position_dodge(width=0.9)) +
  geom_jitter(aes(fill=type),position=position_jitter(width=0.1, height = .2),
              size=2) + ggtitle("overview all markers")+
  facet_grid(celltype~variable, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5))


###subpop Treg + expre
ggplot(data=subset(ExprClin, type%in%c("LKmet", "Sn","SNmet")&celltype=="memB"), aes(type,value)) +
  geom_boxplot(outlier.shape=NA,aes(fill=type),position = position_dodge(width=0.9))  + scale_fill_brewer(palette="Greens")+
  geom_jitter(aes(fill=type),position=position_jitter(width=0.1, height = .2),
              size=2) + ggtitle("overview Treg subPop")+
  facet_wrap(.~variable, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5))+ theme_bw()


#####
N.wilcox = length(unique(ExprClin$variable))
data.wilcox = subset(ExprClin, celltype=="CD4" &type!="ctr")[, c("variable", "type", "value")]
res.wilcox = data.frame(
  variable    = unique(data.wilcox$variable), 
  Sn.SNmet    = rep(NA, N.wilcox),
  Sn.LKmet    = rep(NA, N.wilcox),
  SNmet.LKmet = rep(NA, N.wilcox),
  stringsAsFactors = F)

for(i in 1:N.wilcox) {
  ct = res.wilcox$variable[i]
  d.Sn    = subset(data.wilcox, variable==ct & type=="Sn")$value
  d.SNmet = subset(data.wilcox, variable==ct & type=="SNmet")$value
  d.LKmet = subset(data.wilcox, variable==ct & type=="LKmet")$value
  res.wilcox$Sn.SNmet[i] = wilcox.test(d.Sn, d.SNmet)$p.value
  res.wilcox$Sn.LKmet[i] = wilcox.test(d.Sn, d.LKmet)$p.value
  res.wilcox$SNmet.LKmet[i] = wilcox.test(d.SNmet, d.LKmet)$p.value
  
}

####  expression_subtype

N.wilcox = length(unique(ExprClin$variable))
data.wilcox = subset(ExprClin, celltype=="CD4" &type!="ctr")[, c("variable", "TC_CyTOF_split", "value")]
res.wilcox = data.frame(
  variable    = unique(data.wilcox$variable), 
  TCneg.TClow    = rep(NA, N.wilcox),
  TCneg.TChigh    = rep(NA, N.wilcox),
  TClow.TChigh = rep(NA, N.wilcox),
  stringsAsFactors = F)

for(i in 1:N.wilcox) {
  ct = res.wilcox$variable[i]
  d.TCneg    = subset(data.wilcox, variable==ct & TC_CyTOF_split=="TCneg")$value
  d.TClow = subset(data.wilcox, variable==ct & TC_CyTOF_split=="TClow")$value
  d.TChigh = subset(data.wilcox, variable==ct & TC_CyTOF_split=="TChigh")$value
  res.wilcox$TCneg.TClow[i] = wilcox.test(d.TCneg, d.TClow)$p.value
  res.wilcox$TCneg.TChigh[i] = wilcox.test(d.TCneg, d.TChigh)$p.value
  res.wilcox$TClow.TChigh[i] = wilcox.test(d.TClow, d.TChigh)$p.value
  
}


############ clustermap
####### cluster of median/mean

setwd("C:/Users/ingah/Dropbox/clustermap")

library("gclus")
library("clusterGenomics")
library("foreach")
library("cluster")
library("doParallel")
library("iterators")
library("plyr")

# Set up everything example
source("clustermap.R") # Make package functions available

###abundance

AbunClin=merge(fracSub, clin, by.x="ID", by.y = "ID", all.x = T)

da=AbunClin[-c(17,19:20),] # remove 1149 Snmet
da=subset(da,type!="ctr")
#da=da[-c(1:16,19:57),] #removeall except
# da=subset(da, type%in%c("LKmet","Sn")&run!="run3")
#da=subset(da, type%in%c("LKmet","Sn", "SNmet"))
#da=subset(da, type2!="enz")
#clustData=da[,c (3:24,26:39)]


#normalization? scale between 0 and 1 

##the doit function take all the columns an rescale the numbers between 0 and 1. 0 corresponds to the lowest number and 
#and 1 correspond to the highest number.

dat=da[, c(2:28)] #all markers

doit <- function(dat) {(dat - min(dat, na.rm=TRUE))/(max(dat,na.rm=TRUE) -
                                                       min(dat, na.rm=TRUE))}

norm=as.data.frame(lapply(dat, doit)) # use lapply to apply doit() to every column in a data.frame

lapply(norm, range) #very that the range of al is [0,1]


clustData=norm[, c(1:5)] # min parental
clustData=norm[, c(8:16,18,20,22:24)] # min parental
clustData=norm[, c(1,4,5,6,7,9:14,16,17,22,24)] # min parental + double positive
#clustData=norm[, c(3,4,5,7,10,11,14:17,19:20,24:28,30:34,37)] #all markers min TC and TCR gd
#clustData=norm[, c(4,10,11,14,15,26,28,34)] #Tcell markers
#clustData=norm[, c(4,10,15,16,19,24,28,33,34)] #Tcell markers


type=da$type
subtype=as.factor(da$subtype)
ER=da$ER
TC=da$TC
run=da$run
enz=da$dig2
tc=da$TC_CyTOF
tc2=da$TC_CyTOF_split
rownames(clustData) = da$ID.int

####                   
#########################

#plot.hmap(clustData, colorscale="blue-white-red") 

# clustData = clustData - 0.5
#pdf("ClustHER2.pdf")
plot.init(tree=c(2,3), cbar=3, text=5)
hcluster(clustData, clust="col", distance="euclidean", linkage="complete")
hcluster(clustData, clust="row", distance="euclidean", linkage="complete")
subclust(3, clust="col", method="part") # make identified subcluster
subclust(3, clust="row", method="part") # make subclusteres
plot.hmap(clustData, as.image = F,  colorscale=c("blue", "white", "red"))
z1 = set.color(type, label="type", type="discrete", color=c("tomato1", "pink","blanchedalmond","orange"))
z2 = set.color(subtype, label="subtype", type="discrete", color=c("red", "green","blue"))
z3 = set.color(ER, label="ER", type="discrete", color=c("aquamarine4","aquamarine3"))
z4 = set.color(run, label="run", type="discrete", color=c("red","green","blue","cyan","pink","yellow","orange", "purple","brown"))
z5 = set.color(enz, label="enz",type="discrete",color=c("grey", "white", "pink"))
#z6 = set.color(TC,  label="TC", type="discrete", color=c("aquamarine4","aquamarine1"))
z7 = set.color(tc,  label="TC Cyt", type="discrete", color=c("aquamarine2","pink"))
z8 = set.color(tc2, label="TC cyt2", type="discrete",color=c("red","grey","cyan"))
plot.cbar(z1,z2,z3,z4,z5,z7,z8,  pvalue=T, pvalue.method="fisher", border="black", side=2, cex=1)
plot.tree(side=2)
plot.tree(side=3)
plot.text(colnames(clustData), cex=0.7, side=1)
plot.text(rownames(clustData), cex=0.7, side=4)
plot.hmap.key()
#plot.cbar.key()
###Fig3
#Comb=subset(ClinInfo, Herceptin.Pri=="yes")

da$clust= .CLUSTERMAP$rowgroup    #extract the cluster groups                 


#sett in cluster 1 og 2 i scatterplot
ggplot(da, aes(clust,TC_CyTOF_percent_corrigated))+geom_jitter()


fisher.test(da$clust, da$type)

t.test(clust~TC_CyTOF_percent_corrigated, data=subset(da, type!="ctr" ,clust!="2"))

wilcox.test(TC_CyTOF_percent_corrigated~clust, data=subset(da, type!="ctr"),alternative="two.sided")

####clustering of expression

CD4clin=merge(CD4, clin, by.x="ID", by.y = "ID", all.x=T)
CD4clin= subset(CD4clin, type!="ctr")
CD4clin=CD4clin[-c(17:19),] #remove 1149 SNmet+ one extra copy of 1149 LKmet

            
ClustCD4=CD4clin[, c(2:38)] #all markers

doit <- function(ClustCD4) {(ClustCD4 - min(ClustCD4, na.rm=TRUE))/(max(ClustCD4,na.rm=TRUE) -
                                                       min(ClustCD4, na.rm=TRUE))}
ClustCD4=as.data.frame(lapply(ClustCD4, doit)) # use lapply to apply doit() to every column in a data.frame

lapply(ClustCD4, range) #very that the range of al is [0,1]



clust=ClustCD4[, c(1:37)] # min parental
clustData=clust[, c(3,5,7,11,14:16,19:21,26:27,32:33)] # min breast markers 


type=CD4clin$type
subtype=as.factor(CD4clin$subtype)
ER=CD4clin$ER
TC=CD4clin$TC
run=CD4clin$run
enz=CD4clin$dig2
tc=CD4clin$TC_CyTOF
tc2=CD4clin$TC_CyTOF_split
rownames(clustData) = CD4clin$ID.int

####                   
#########################

plot.hmap(clustData, colorscale="blue-white-red") 

# clustData = clustData - 0.5
#pdf("ClustHER2.pdf")
plot.init(tree=c(2,3), cbar=2, text=5)
hcluster(clustData, clust="col", distance="euclidean", linkage="complete")
hcluster(clustData, clust="row", distance="euclidean", linkage="complete")
subclust(3, clust="col", method="part") # make identified subcluster
subclust(3, clust="row", method="part") # make subclusteres
plot.hmap(clustData,as.image=F,  colorscale=c("blue", "white", "red"))
z1 = set.color(type, label="type", type="discrete", color=c("tomato1", "pink","blanchedalmond","orange"))
z2 = set.color(subtype, label="subtype", type="discrete", color=c("red", "green","blue"))
z3 = set.color(ER, label="ER", type="discrete", color=c("aquamarine4","aquamarine3"))
z4 = set.color(run, label="run", type="discrete", color=c("red","green","blue","cyan","pink","yellow","orange", "purple","brown"))
z5 = set.color(enz, label="enz",type="discrete",color=c("grey", "white", "pink"))
#z6 = set.color(TC,  label="TC", type="discrete", color=c("aquamarine4","aquamarine1"))
z7 = set.color(tc,  label="TC Cyt", type="discrete", color=c("aquamarine2","pink"))
z8 = set.color(tc2, label="TC cyt2", type="discrete",color=c("red","grey","cyan"))
plot.cbar(z1,z2,z3,z4,z5,z7,z8,  pvalue=T, pvalue.method="fisher", border="black", side=2, cex=0.5)
plot.tree(side=2)
plot.tree(side=3)
plot.text(colnames(clustData), cex=1, side=1)
plot.text(rownames(clustData), cex=1, side=4)
plot.hmap.key()
#plot.cbar.key()
###Fig3
#Comb=subset(ClinInfo, Herceptin.Pri=="yes")



########clustering 1149 SNmet vs LKmet


##abundance
AbunClin=merge(fracSub, clin, by.x="ID", by.y = "ID", all.x = T)

daRed=AbunClin[c(17,20),] # keep only 1149 Snmet and LKmet


#normalization? scale between 0 and 1 

##the doit function take all the columns an rescale the numbers between 0 and 1. 0 corresponds to the lowest number and 
#and 1 correspond to the highest number.

daRed2=daRed[, c(2:28)] #all markers

#alternative:range of all numbers between 0 and 1
doit <- function(daRed2) {(daRed2 - min(daRed2, na.rm=TRUE))/(max(daRed2,na.rm=TRUE) -
                                                       min(daRed2, na.rm=TRUE))}

norm=as.data.frame(lapply(daRed2, doit)) # use lapply to apply doit() to every column in a data.frame

lapply(norm, range) #very that the range of al is [0,1]


clustData=norm[, c(1:5)] #  parental
clustData=norm[, c(8:16,18,20,22:24)] # min parental
clustData=norm[, c(1,4,5,6,7,9:14,16,17,22,24)] # min parental + double positive
#clustData=norm[, c(3,4,5,7,10,11,14:17,19:20,24:28,30:34,37)] #all markers min TC and TCR gd
#clustData=norm[, c(4,10,11,14,15,26,28,34)] #Tcell markers
#clustData=norm[, c(4,10,15,16,19,24,28,33,34)] #Tcell markers


subtype=as.factor(da$subtype)
TC=da$TC
run=da$run
enz=da$dig2
tc=da$TC_CyTOF
tc2=da$TC_CyTOF_split
rownames(clustData) = da$type

####                   
#########################

#plot.hmap(clustData, colorscale="blue-white-red") 

# clustData = clustData - 0.5
#pdf("ClustHER2.pdf")
plot.init(tree=c(2,3), cbar=3, text=5)
hcluster(clustData, clust="col", distance="euclidean", linkage="complete")
hcluster(clustData, clust="row", distance="euclidean", linkage="complete")
subclust(3, clust="col", method="part") # make identified subcluster
subclust(3, clust="row", method="part") # make subclusteres
plot.hmap(clustData, as.image = F,  colorscale=c("blue", "white", "red"))
z1 = set.color(type, label="type", type="discrete", color=c("tomato1", "pink","blanchedalmond","orange"))
z2 = set.color(subtype, label="subtype", type="discrete", color=c("red", "green","blue"))
z4 = set.color(run, label="run", type="discrete", color=c("red","green","blue","cyan","pink","yellow","orange", "purple","brown"))
z5 = set.color(enz, label="enz",type="discrete",color=c("grey", "white", "pink"))
#z6 = set.color(TC,  label="TC", type="discrete", color=c("aquamarine4","aquamarine1"))
z7 = set.color(tc,  label="TC Cyt", type="discrete", color=c("aquamarine2","pink"))
z8 = set.color(tc2, label="TC cyt2", type="discrete",color=c("red","grey","cyan"))
plot.cbar(z1,z2,z4,z5,z7,z8,  pvalue=T, pvalue.method="fisher", border="black", side=2, cex=1)
plot.tree(side=2)
plot.tree(side=3)
plot.text(colnames(clustData), cex=0.7, side=1)
plot.text(rownames(clustData), cex=0.7, side=4)
plot.hmap.key()

#############histogram of fig4. barplots of CD4 CD8
hist1=subset(dataClin2, ID%in% c("1026","1326","1363","1182"))
hist2=subset(hist1, variable %in% c("CD4","CD8","CD4CD8"))

ggplot(data=subset(hist2, type%in%c("LKmet", "Sn","SNmet")&subfrac=="fracSub"), aes(ID,value)) +
  geom_bar(outlier.shape=NA,aes(fill=type),stat="identity")+reorder(type) + scale_fill_brewer(palette="Reds")+
    facet_wrap(~variable, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5)) + theme_bw()

#NB if you want several bars sideways use position="dodge" after stat="identity"

###splitting on subtype
ggplot(data=subset(dataClin2, type%in%c("LKmet","SNmet")&subfrac=="fracSub"&variable%in%c("Plasmacells","memBcell")), aes(subtype,value)) +
  geom_boxplot(outlier.shape=NA,aes(fill=subtype),position = position_dodge(width=0.5)) + scale_fill_brewer(palette="Greens")+
  geom_jitter(aes(color=as.factor(type)),position=position_dodge(width=0.1),
              size=1) + ggtitle("amn 1149SNmet, split on subtype")+
  facet_wrap(~variable, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5)) + theme_bw()

N.wilcox = length(unique(dataClin2$variable))
data.wilcox = subset(dataClin2, subfrac=="fracSub" &type%in%c("LKmet","SNmet"))[, c("variable", "subtype", "value")]
res.wilcox = data.frame(
  variable    = unique(data.wilcox$variable), 
  ER.HER2    = rep(NA, N.wilcox),
  ER.TN    = rep(NA, N.wilcox),
  HER2.TN = rep(NA, N.wilcox),
  stringsAsFactors = F)

for(i in 1:N.wilcox) {
  ct = res.wilcox$variable[i]
  d.ER    = subset(data.wilcox, variable==ct & subtype=="ER")$value
  d.HER2 = subset(data.wilcox, variable==ct & subtype=="HER2")$value
  d.TN = subset(data.wilcox, variable==ct & subtype=="TN")$value
  res.wilcox$ER.HER2[i] = wilcox.test(d.ER, d.HER2)$p.value
  res.wilcox$ER.TN[i] = wilcox.test(d.ER, d.TN)$p.value
  res.wilcox$HER2.TN[i] = wilcox.test(d.HER2, d.TN)$p.value
  
}

#merge TN neg and HER files for comparison
dataClin2$subtype2=rep(NA, nrow(dataClin2))
dataClin2$subtype2[dataClin2$subtype=="HER2"]="HER2+TN"
dataClin2$subtype2[dataClin2$subtype=="TN"]="HER2+TN"
dataClin2$subtype2[dataClin2$subtype=="ER"]="ER"
dataClin2$subtype2<-factor(dataClin2$subtype2, levels=c("HER2+TN", "ER"))


#merge from two columns





#### comparison of mass cytomerty data and IHC measured by pathvision

CompIHCMC=read.table("../Rscripts_files/comparisonIHCMassCyt.txt", header=T, stringsAsFactors=F, sep="\t",na.strings = "NA")
CompIHCMC$ID= as.character(CompIHCMC$ID)
CompIHCMC$variable= as.factor(CompIHCMC$variable)
CompIHCMC$celltype<-character(CompIHCMC$celltype, levels=c("CD4", "CD8","CD4/CD8"))

#make comparison histogram between mass cytometry with whole section field fraction from immunopath
ggplot(data=subset(CompIHCMC, variable%in%c("fracSub","fieldFrac")), aes(ID,value)) +
  geom_bar(outlier.shape=NA,aes(fill=variable),stat="identity", position="dodge") + scale_fill_brewer(palette="Reds")+
  
  facet_wrap(~celltype, scales="free_y") +
  theme(axis.text.x=element_text(angle=-90, hjust=1, size=15, vjust=0.5)) + theme_bw()

##line plot

CompIHCMC$ID=as.factor(CompIHCMC$ID)
ggplot(data=subset(CompIHCMC, variable%in%c("fieldFrac","fracTot")),aes(x=ID, y=value, group=variable, color=variable)) +
  geom_line() +geom_point()+
  scale_color_brewer(palette="paired") +
  facet_wrap(~celltype, scales="free_y") +
  theme_minimal()

## 

################################ survival analysis###################################

clin2=clin[-25,] # remove the extra column for 1149 SNmet, we are only using the ALNmet sample

#clin$MonthsToLastUpdate=as.numeric(as.character(clin$MonthsToLastUpdate))
#clin$mestasticStatus=as.numeric(as.character(clin$mestasticStatus))

#plot(survfit(Surv(clin$ MonthsToLastUpdate, 
 #                 clin$mestasticStatus)~clin$subtype))

#fit=survfit(Surv(Clust$MontSurgFinal, 
 #                Clust$Survival.Status)~Clust$clustGen)
#plot(fit, col=c(1:3),xlab="time",mark.time = T, ylab="survival probability", main="surv Gen") 
#legend(30, 1, c("1", "2","3"), col= (1:3), lwd=0.5)
#summary(fit)
#survdiff(Surv(Clust$MontSurgFinal, Clust$Survival.Status)
#         ~Clust$clustGen, rho=0)


fit<-survfit(Surv(MonthsToLastUpdate, metastaticStatus)~subtype, data=clin2)
 ggsurvplot(fit, pval=TRUE,conf.int=TRUE, risk.table = TRUE, data=clin)

 res<-ggsurvplot(fit, pval=TRUE,
                 conf.int=TRUE,
                 risk.table = TRUE,
               
                 data=clin)

###########boxplot type vs TC cytof percentage
 
 clin2=clin[-25, ]# remove 1149 SNmet
 clin2=subset(clin2, type!="ctr")
 p<-ggplot(clin2, aes(x=type, y=TC_CyTOF_percent_corrigated))+ geom_jitter(position=position_jitter(width=0.1)) +scale_y_log10() + geom_hline(yintercept = 0.02) 
 
 #calculate Fisher.exact pathology vs mass cytometry
 
 table<-matrix(c(14,2,7,8,3,15),ncol=2, byrow=TRUE)
 fisher.test(table)
