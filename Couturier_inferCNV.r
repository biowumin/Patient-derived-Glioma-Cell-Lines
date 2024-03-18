## infercnv
## get the genome position of each gene from gtf file
wget https://github.com/broadinstitute/infercnv/blob/master/scripts/gtf_to_position_file.py ./
###################################################### Run_infercnv.R begin
rm(list = ls())
library(Seurat)
library(infercnv)
setwd("Couturier_scRNA/GBM_res")
read.table("GRCh38.gtf_gene_pos.txt",sep="\t",stringsAsFactors = F) -> gene_order
gsub("^([0-9XY])","chr\\1",gene_order$V2) -> gene_order$V2
readRDS("GBM_liger.rds") -> GBM_liger
Idents(object = GBM_liger) -> dat
type <- ifelse(dat%in%c(4,8,9,11,13),"normal","malignant")
data.frame(cell=names(dat),type=type,stringsAsFactors = F) -> cell_anno
rownames(cell_anno) <- cell_anno$cell
as.matrix(GetAssayData(object = GBM_liger, slot = "counts")) -> expr_matrix
cell_anno <- cell_anno[colnames(expr_matrix),]
write.table(cell_anno,file="cell_anno.txt",sep="\t",quote = F,row.names = F,col.names=F)
write.table(gene_order,file="gene_order.txt",sep="\t",quote = F,row.names = F,col.names=F)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix=expr_matrix,
  annotations_file="cell_anno.txt",
  gene_order_file= "gene_order.txt",
  delim="\t",
  ref_group_names= c("normal"))

out_dir = "/media/usb3/Couturier_scRNA/GBM_res/infercnv/"
infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out_dir,
  cluster_by_groups=TRUE,
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=T,
  num_threads=60,
  plot_probabilities=F,
  no_plot=T,
  no_prelim_plot=T,
  output_format=NA
)
readRDS("infercnv/run.final.infercnv_obj") -> infercnv_obj
plot_cnv(infercnv_obj)    ## generate some files for CNV plots
###################################################### visualization the inferCNV results
setwd("Couturier_scRNA/GBM_res/infercnv")
readRDS("run.final.infercnv_obj") -> infercnv_obj
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices$normal
tumor_loc <- infercnv_obj@observation_grouped_cell_indices$malignant
norm_expr <- expr[,normal_loc]
tumor_expr <- expr[,tumor_loc]
used <- grep("GSC",colnames(tumor_expr),invert = T)
tumor_expr <- tumor_expr[,used]
readRDS("../GBM_liger.rds") -> GBM_liger
GBM_liger$orig.ident -> tmp
tmp <- tmp[tumor_loc]
cell_to_patient <- data.frame(patient=as.character(tmp),row.names = names(tmp),stringsAsFactors = F)
cell_to_patient <- subset(cell_to_patient,patient%in%grep("GSC",cell_to_patient$patient,value=T,invert=T))
## begin to draw the heatmap
library("ComplexHeatmap")
library("RColorBrewer")
library(circlize)
pos_control <- c("EGFR","PDGFRA","CDK4","CDKN2A","PTEN")
n <- t(norm_expr[c("MYC","MYCN",pos_control),])
t <- t(tumor_expr[c("MYC","MYCN",pos_control),])
col_fun = colorRamp2(c(0.8,0.938,1.018,1.124,1.21), c("#4169E1","white","white","white", "#DC143C"))
ht_normal = Heatmap(n,
                    col = col_fun,
                    border = "black",
                    cluster_rows = F,
                    cluster_columns = F,
                    show_heatmap_legend=F,
                    show_column_names = T,
                    show_row_names = F,
                    column_title = NULL,
                    row_title = "normal")
ht_tumor = Heatmap(t,
                    col = col_fun,
                    border = "black",
                    cluster_rows=F,
                    cluster_columns = F,
                    row_split = cell_to_patient$patient,
                    show_heatmap_legend=F,
                    show_column_names = T,
                    show_row_names = F,
                    column_title = NULL,row_title_gp = gpar(fontsize = 7),
                    row_title_rot=0)
ht_list <- ht_normal %v% ht_tumor
pdf("infercnv.pdf")
draw(ht_list)
dev.off()
###################################################### Compare the CNV level between GSC and tissue
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
read.table("Couturier_scRNA/GBM_res/infercnv/infercnv.references.txt",header = T,row.names = 1,sep = " ",stringsAsFactors = F) -> ref_CNV
save(ref_CNV,file="ref_CNV.RData")
ref_CNV <- ref_CNV[c("MYC","MYCN"),c(grep("BT333[-]*",colnames(ref_CNV)),grep("BT363[-]*",colnames(ref_CNV)),grep("BT368[-]*",colnames(ref_CNV)))]
ref_CNV <- ref_CNV -1
apply(ref_CNV,2,function(x){sum(x^2)}) -> ref_score
ref_anno <- rep("Normal",length(ref_score))
read.table("Couturier_scRNA/GBM_res/infercnv/infercnv.observations.txt",header = T,row.names = 1,sep = " ",stringsAsFactors = F) -> obs_CNV
save(obs_CNV,file="obs_CNV.RData")
obs_CNV <- obs_CNV[c("MYC","MYCN"),c(grep("BT333[-]*",colnames(obs_CNV)),grep("BT363[-]*",colnames(obs_CNV)),grep("BT368[-]*",colnames(obs_CNV)))]
obs_CNV <- obs_CNV -1
apply(obs_CNV,2,function(x){sum(x^2)}) -> obs_score
obs_anno <- gsub("_.*","",names(obs_score))
dat <- data.frame(CNV=c(ref_score,obs_score),anno=c(ref_anno,obs_anno))

obs_CNV <- as.data.frame(t(obs_CNV))
obs_CNV$anno <- gsub("_.*","",rownames(obs_CNV))
dat <- rbind(ref_CNV,obs_CNV)
df <- pivot_longer(dat,cols = MYC:PTEN, names_to = "gene",values_to = "CNV")
ggbarplot(df, x = "anno", y = "CNV", facet.by = "gene", add = "mean_se",color = "anno", position = position_dodge(0.9),width=0.3) + theme(legend.position = "none")




