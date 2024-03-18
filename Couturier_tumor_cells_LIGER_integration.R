################################################### integrate data from different patients using LIGER
library(Seurat)
library(rliger)
library(tidyverse)
library(SeuratWrappers)
library(UCell)
library(RColorBrewer)
packageVersion("rliger")   ## version 1.0.0
setwd("Couturier_scRNA/GBM_res")
readRDS("GBM_liger.rds") -> GBM_liger
GBM_tumor <- subset(GBM_liger,idents = c(0:3,5:7,10,12,14))
obj.list <- SplitObject(GBM_tumor, split.by = "orig.ident")
tumor_list <- lapply(obj.list,function(x){x <- CreateSeuratObject(x@assays$RNA@counts)})
tumor_merge <- merge(tumor_list[[1]],tumor_list[2:length(tumor_list)])
tumor_merge <- NormalizeData(tumor_merge) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% ScaleData(split.by = "orig.ident", do.center = FALSE)
nFactors <- 20
tumor_merge <- RunOptimizeALS(tumor_merge,k=nFactors,split.by="orig.ident")
tumor_merge <- RunQuantileNorm(tumor_merge, split.by="orig.ident")
tumor_merge <- FindNeighbors(tumor_merge, reduction="iNMF", dims=1:nFactors) %>% FindClusters(resolution = 0.4)
tumor_merge <- RunUMAP(tumor_merge, dims=1:nFactors, reduction="iNMF")
saveRDS(tumor_merge,file="tumor_merge.RDS")

read.table("top150_markers.txt",header=T,sep="\t") -> top150_markers
top150_markers <- as.list(top150_markers)
load("cellular_state_score/cell_state_marker2.RData")
tumor_merge <- AddModuleScore_UCell(tumor_merge, features = cell_state_marker2)
colnames(tumor_merge@meta.data)[7] <- "MES_Neftel_UCell"
tumor_merge <- AddModuleScore_UCell(tumor_merge, features = top150_markers)

source("type_identifier.R")
load("cellular_state2_score.RData")
tumor_merge$cellular_state2 <- type_identifier(score)
load("three_subtype_score.RData")
tumor_merge$subtype <- type_identifier(score)

GSC <- subset(tumor_merge,subset = orig.ident%in%c("BT322-GSC","BT324-GSC","BT326-GSC","BT333-GSC","BT363-GSC","BT368-GSC"))
tumor <- subset(tumor_merge,subset = orig.ident%in%c("BT322-GSC","BT324-GSC","BT326-GSC","BT333-GSC","BT363-GSC","BT368-GSC"),invert=T)
mycolors <- c(brewer.pal(12,'Paired'),"lightgray")
pdf("tumor_umap.pdf")
DimPlot(tumor, reduction = "umap",label=T,pt.size = 0.8)
DimPlot(tumor, reduction = "umap",label=F,group.by="subtype",pt.size = 0.8,cols = scales::alpha(c("red","lightgray","mediumblue","purple"),0.4))
DimPlot(tumor, reduction = "umap",label=F,group.by="cellular_state2",pt.size = 0.8,cols = alpha(mycolors[c(5,6,13,1,2)],0.4))
FeaturePlot(tumor, reduction = "umap", features = c("MES_UCell", "PN_UCell","OXPHOS_UCell"),ncol = 2,pt.size = 0.8,min.cutoff = "q03",max.cutoff = "q99",order=T)
FeaturePlot(tumor, reduction = "umap", features = c("MES_Neftel_UCell","AC_UCell","NPC_UCell","OPC_UCell"),ncol = 2,pt.size = 0.8,min.cutoff = "q03",max.cutoff = "q99",order=T)
dev.off()
as.data.frame(table(tumor$orig.ident,tumor$cellular_state2)) -> dat1
pdf("state2_tumor_stat.pdf", width = 9, height = 8)
ggplot(data = dat1) + geom_bar(aes(x = Var1, y = Freq, fill = factor(Var2)),position = "fill",stat = "identity",alpha=0.7,width = 0.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=mycolors[c(5,6,13,1,2)])
dev.off()
as.data.frame(table(tumor$orig.ident,tumor$subtype)) -> dat2
pdf("subtype_tumor_stat.pdf", width = 9, height = 8)
ggplot(data = dat2) + geom_bar(aes(x = Var1, y = Freq, fill = factor(Var2)),position = "fill",stat = "identity",alpha=0.7,width = 0.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c("red","lightgray","mediumblue","purple"))
dev.off()

pdf("GSC_umap.pdf")
DimPlot(GSC, reduction = "umap",label=F,pt.size = 0.8)
DimPlot(GSC, reduction = "umap",label=F,group.by="subtype",pt.size = 0.8,cols = alpha(c("red","lightgray","mediumblue","purple"),0.4))
DimPlot(GSC, reduction = "umap",label=F,group.by="cellular_state2",pt.size = 0.8,cols = alpha(mycolors[c(5,6,13,1,2)],0.4))
FeaturePlot(GSC, reduction = "umap", features = c("MES_UCell", "PN_UCell","OXPHOS_UCell"),ncol = 2,pt.size = 0.8,min.cutoff = "q03",max.cutoff = "q99",order=T)
FeaturePlot(GSC, reduction = "umap", features = c("MES_Neftel_UCell","AC_UCell","NPC_UCell","OPC_UCell"),ncol = 2,pt.size = 0.8,min.cutoff = "q03",max.cutoff = "q99",order=T)
dev.off()
as.data.frame(table(GSC$orig.ident,GSC$cellular_state2)) -> dat3
pdf("state2_GSC_stat.pdf", width = 9, height = 8)
ggplot(data = dat3) + geom_bar(aes(x = Var1, y = Freq, fill = factor(Var2)),position = "fill",stat = "identity",alpha=0.7,width = 0.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=mycolors[c(5,6,13,1,2)])
dev.off()
as.data.frame(table(GSC$orig.ident,GSC$subtype)) -> dat4
pdf("subtype_GSC_stat.pdf", width = 9, height = 8)
ggplot(data = dat4) + geom_bar(aes(x = Var1, y = Freq, fill = factor(Var2)),position = "fill",stat = "identity",alpha=0.7,width = 0.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c("red","lightgray","mediumblue","purple"))
dev.off()

rbind(dat1,dat3) %>% filter(Var1%in%c("BT333","BT333-GSC","BT363","BT363-GSC","BT368","BT368-GSC")) -> dat5
dat5$Var1 <- factor(dat5$Var1,levels=c("BT333","BT333-GSC","BT363","BT363-GSC","BT368","BT368-GSC"))
pdf("state2_GSC_tumor_compare.pdf")
ggplot(data = dat5) + geom_bar(aes(x = Var1, y = Freq, fill = factor(Var2)),position = "fill",stat = "identity",alpha=0.7,width = 0.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=mycolors[c(5,6,13,1,2)])
dev.off()

rbind(dat2,dat4) %>% filter(Var1%in%c("BT333","BT333-GSC","BT363","BT363-GSC","BT368","BT368-GSC")) -> dat6
dat6$Var1 <- factor(dat6$Var1,levels=c("BT333","BT333-GSC","BT363","BT363-GSC","BT368","BT368-GSC"))
pdf("subtype_GSC_tumor_compare.pdf")
ggplot(data = dat6) + geom_bar(aes(x = Var1, y = Freq, fill = factor(Var2)),position = "fill",stat = "identity",alpha=0.7,width = 0.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c("red","lightgray","mediumblue","purple"))
dev.off()