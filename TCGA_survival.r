rm(list=ls())
options(stringsAsFactors = F)
read.table("GBM_TCGA/HT_HG-U133A",header=T,row.names=1,sep="\t") -> U133A_exp
U133A_exp <- U133A_exp[,-(grep("11$",colnames(U133A_exp)))]
load("sig_res/markers100.RData")

library(ssgsea.GBM.classification)
source("~/ssgseaMOD.r")
marker <- rep(1,length(markers100$mes))
names(marker) <- markers100$mes
mod.generate(marker,out="mes.mod")
marker <- rep(1,length(markers100$pn))
names(marker) <- markers100$pn
mod.generate(marker,out="pn.mod")
marker <- rep(1,length(markers100$oxphos))
names(marker) <- markers100$oxphos
mod.generate(marker,out="oxphos.mod")
score <- as.data.frame(mod.analyze2(U133A_exp,c('mes','pn','oxphos'),permN=1000))
save(score,file="TCGA_U133A_sig_score.RData")

## Determine the subtype
type_identifier <- function(score){  ## score was the output of ssGSEA, including norm_score and p_val
   pval <- score[,grep("pval",colnames(score))]
   pt_type <- rownames(pval)
   names(pt_type) <- rownames(pval)
   pt_type[apply(pval,1,function(x){min(x)>0.05})] <- "None"
   pval[intersect(rownames(pval),pt_type),] -> tmp
   apply(tmp,1,function(x){length(which(x==min(x)))>1}) -> tmp_index
   pt_type[rownames(tmp)[tmp_index]] <- "None"
   pval[intersect(rownames(pval),pt_type),] -> tmp
   apply(tmp,1,function(x){which(x==min(x))}) -> tmp_index
   pt_type[rownames(tmp)] <- colnames(tmp)[tmp_index]
   gsub("_pval","",pt_type) -> pt_type
   return(pt_type)
}
type_identifier(score) -> pt_type
score$type <- pt_type[rownames(score)]
## heatmap
library(pheatmap)
library(RColorBrewer)
library(dplyr)
dat <- U133A_exp[intersect(c(markers100$mes,markers100$pn,markers100$oxphos),rownames(U133A_exp)),rownames(subset(score,type%in%c("mes","pn","oxphos")))]
ann_col <- subset(score,select=type)
ann_colors <- list(type=c(mes="red",pn="purple",oxphos="mediumblue"))
tiff(file="TCGA_heatmap.tiff",width=8,height=8,units="in",compression="lzw",res=300)
#pdf("TCGA_heatmap.pdf",8,6)
pheatmap(dat,clustering_method="ward.D2",cutree_cols = 3,cutree_rows = 3,legend=F,annotation_legend = F,annotation_names_col=F,annotation_col=ann_col,annotation_colors = ann_colors,cluster_rows=T, show_rownames=F,show_colnames =F,cluster_cols=T,scale="row",clustering_distance_cols = "correlation",clustering_distance_rows = "correlation",col=colorRampPalette(c("blue", "white", "red"))(100),breaks=seq(-4,4,length=100))
dev.off()

rownames(score) <- gsub("[.]","-",rownames(score))
read.table("GBM_TCGA/GBM_survival.txt",sep="\t",header = T,row.names = 1) -> TCGA_survival
dat <- subset(TCGA_survival,rownames(TCGA_survival)%in%rownames(score))
dat$type <- score[rownames(dat),10]
dat <- dat %>% filter(type != "None")
dat$OS.time <- dat$OS.time/30
dat$PFI.time <- dat$PFI.time/30
library(survminer)
library(survival)
sfit <- survfit(Surv(OS.time,OS)~type, data=dat)
ggsurv <- ggsurvplot(sfit, data = dat,
                     pval = T,
                     surv.median.line = "hv",
                     risk.table = T,
                     xlim = c(0,24),
                     break.time.by = 5,risk.table.height = 0.2,surv.plot.height=0.5,
                     palette = c("red","mediumblue","purple"),
                     size = 0.6,fontsize=3)
pdf("TCGA_survival.pdf",10,7)
print(ggsurv, newpage = FALSE)
dev.off()