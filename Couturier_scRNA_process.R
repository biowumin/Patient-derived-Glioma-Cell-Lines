rm(list=ls())
library(Seurat)
library(rliger)
library(SeuratWrappers)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
setwd("Couturier_scRNA/filtered")
dir <- dir("./")
samples_name <- gsub(".filtered_gene_matrices","",dir)
samples_name <- gsub("_.*","",samples_name)
samples_name <- gsub("BT322","BT322-GSC",samples_name)

## read the scRNA-seq data of tumor tissue
scRNA_GBM <- list()
for(i in 1:23){
    counts <- Read10X(data.dir=dir[i])
    scRNA_GBM[[i]] <- CreateSeuratObject(counts = counts, project = samples_name[i], min.cells = 3, min.features = 200)
    scRNA_GBM[[i]] <- RenameCells(scRNA_GBM[[i]],add.cell.id = samples_name[i])
    scRNA_GBM[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNA_GBM[[i]], pattern = "^MT-")
    scRNA_GBM[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNA_GBM[[i]], pattern = "^RP[SL]")
}
sum(sapply(scRNA_GBM,function(aa){dim(aa@assays$RNA)[2]}))
GBM_merge <- merge(scRNA_GBM[[1]],scRNA_GBM[2:length(scRNA_GBM)])
table(GBM_merge$orig.ident)
saveRDS(GBM_merge,file="Couturier_scRNA/GBM_res/GBM_merge.rds")

## quality control
readRDS("GBM_merge.rds") -> GBM_merge
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb")
group = "orig.ident"
plots = list()
theme.set2 <- theme(axis.title.x=element_blank())
for(i in seq_along(plot.features)){
    plots[[i]] <- VlnPlot(GBM_merge,group.by=group,pt.size = 0,features = plot.features[i]) + theme.set2 + NoLegend()
}
violin <- wrap_plots(plots=plots,nrow=2)
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf",plot=violin,width = 9, height = 8)

## set QC threshold
minGene <- 800
maxGene <- 7500
minUMI <- 1500
maxUMI <- 50000
pctMT <- 10
pctRB <- 30
GBM_QC <- subset(GBM_merge,subset = nCount_RNA < maxUMI & nCount_RNA > minUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.rb < pctRB)
saveRDS(GBM_QC,file="Couturier_scRNA/GBM_res/GBM_QC.rds")

plots = list()
for(i in seq_along(plot.features)){
    plots[[i]] <- VlnPlot(GBM_QC,group.by=group,pt.size = 0,features = plot.features[i]) + theme.set2 + NoLegend()
}
violin <- wrap_plots(plots=plots,nrow=2)
ggsave("QC/vlnplot_after_qc.pdf",plot=violin,width = 9, height = 8)

## integrate data from different patients
readRDS("GBM_QC.rds") -> GBM_QC
obj.list <- SplitObject(GBM_QC, split.by = "orig.ident")
tumor_list <- lapply(obj.list,function(x){x <- CreateSeuratObject(x@assays$RNA@counts)})
GBM_liger <- merge(tumor_list[[1]],tumor_list[2:length(tumor_list)])
GBM_liger <- NormalizeData(GBM_liger) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% ScaleData(split.by = "orig.ident", do.center = FALSE)
nFactors <- 20
GBM_liger <- RunOptimizeALS(GBM_liger,k=nFactors,split.by="orig.ident")
GBM_liger <- RunQuantileNorm(GBM_liger, split.by="orig.ident")
GBM_liger <- FindNeighbors(GBM_liger, reduction="iNMF", dims=1:nFactors) %>% FindClusters(resolution = 0.3)
GBM_liger <- RunUMAP(GBM_liger, dims=1:nFactors, reduction="iNMF")
saveRDS(GBM_liger,file="Couturier_scRNA/GBM_res/GBM_liger.rds")

GBM_liger$cell_anno <- as.character(GBM_liger$seurat_clusters)
GBM_liger$cell_anno[GBM_liger$cell_anno%in%c("4","8","11")] <- "Myeloid"
GBM_liger$cell_anno[GBM_liger$cell_anno%in%c("9")] <- "Oligodendrocyte"
GBM_liger$cell_anno[GBM_liger$cell_anno%in%c("13")] <- "Endothelial"
GBM_liger$cell_anno[GBM_liger$cell_anno%in%c("0","1","2","3","5","6","7","10","12","14")] <- "tumor"

mycolors <- c(brewer.pal(12,'Paired'),"lightgray")
pdf("all_cells_umap.pdf")
DimPlot(GBM_liger, reduction = "umap",label=T)
tiff(file="cell_anno.tiff",width=12,height=8,units="in",compression="lzw",res=300)
DimPlot(GBM_liger, reduction = "umap",label=F,group.by="cell_anno") + theme(plot.title=element_blank())
dev.off()
tiff(file="cell_marker.tiff",width=12,height=8,units="in",compression="lzw",res=300)
FeaturePlot(GBM_liger, features = c("MCAM","ESAM","CD53","CD68","MOG","MBP"))
dev.off()
VlnPlot(GBM_liger,features = c("MCAM","ESAM","CD53","CD68","MOG","MBP"),group.by="cell_anno",pt.size=0)
dev.off()