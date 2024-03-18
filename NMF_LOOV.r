rm(list=ls())
library(data.table)
load("CPM_exp.RData")
CPM_exp -> PDGC_CPM_exp
############################################################# step 0
mads <- apply(PDGC_CPM_exp,1,mad)   ## calculate MAD value
exp_by_mad <- PDGC_CPM_exp[rev(order(mads)),]   ## sort MAD from high to low
library(NMF)
nmf(exp_by_mad[1:4500,], rank = 2:7, seed = 123456, nrun = 50,.opt='vp100') -> nmf_res
coph <- nmf_res$measures$cophenetic
pdf("NMF_res/select_k.pdf")
plot(2:7,coph,type="b",col="red",ylim=c(0.95,1))
dev.off()
nmf(exp_by_mad[1:4500,], rank = 4, seed = 123456, nrun = 100,.opt='vp100') -> nmf_res
pdf("NMF_res/k4_consensusmap.pdf")
consensusmap(nmf_res)
dev.off()
pdf("NMF_res/k4_coefmap.pdf")
coefmap(nmf_res)
dev.off()
pdf("NMF_res/k4_silhouette.pdf")
plot(silhouette(nmf_res))
dev.off()
data.frame(NMF_cluster = sort(predict(nmf_res,what="consensus")),stringsAsFactors=F) -> NMF_cluster
save(NMF_cluster,file="NMF_res/NMF_cluster.RData")
silhouette(nmf_res) -> silh
save(silh,file="NMF_res/silh.RData")
############################################################# step 1
runNMF2 <- function(exp_matrix, k=4, gene_num=4500){
    mads <- apply(exp_matrix,1,mad)   ## calculate MAD value
    exp_by_mad <- exp_matrix[rev(order(mads)),]   ## sort MAD from high to low
    exp_by_mad <- exp_by_mad[1:gene_num,]
    rm(exp_matrix)
    library(NMF)
    res <- lapply(gene_num, function(r,k){
        nmf(exp_by_mad[1:r,], rank = k, seed = 123456, nrun = 100,.opt='vp100') -> nmf_res
        list(coph = cophcor(nmf_res),cls=as.numeric(predict(nmf_res,,what="consensus")))
    },k=k)
    names(res) <- paste0("gene",gene_num)
    return(res)
}
############################################################# step 2
sapply(1:ncol(PDGC_CPM_exp),function(x,CPM_exp){
    tmp <- CPM_exp[,-x]
    runNMF2(tmp) -> aa
    sapply(aa,'[[','cls') -> bb  ## 聚类结果
    rownames(bb) <- colnames(tmp)
    sapply(aa,'[[','coph') -> cc
    rm(list=c("tmp","aa"))
    bb[,names(which(cc==max(cc)))] -> final_cls   ## 选出最佳聚类结果
    add <- "N"
    names(add) <- setdiff(colnames(CPM_exp),names(final_cls))
    final_cls <- c(final_cls,add)
    final_cls <- final_cls[colnames(CPM_exp)]
},CPM_exp=PDGC_CPM_exp) -> final_cNMF_res
save(final_cNMF_res,file="NMF_res/final_cNMF_res.RData")
############################################################# step 3
count_cls <- function(x){
    stat_matrix <- matrix(0,length(x),length(x),dimnames=list(names(x),names(x)))
    unique(x) -> tmp
    for(i in tmp){
        names(x)[which(x==i)] -> name
        stat_matrix[name,name] <- 1
    }
    return(stat_matrix)
}
final_stat <- matrix(0,nrow(final_cNMF_res),nrow(final_cNMF_res))
for(i in 1:ncol(final_cNMF_res)){
    count_cls(final_cNMF_res[,i]) -> tmp
    final_stat <- final_stat + tmp
}
write.table(final_stat,file="NMF_res/final_stat.txt",sep="\t",quote=F)
############################################################# step 4
read.table("NMF_res/final_stat.txt",sep="\t",header=T) -> final_stat
final_stat <- final_stat/ncol(PDGC_CPM_exp)
library(pheatmap)
pdf("NMF_res/LOOV_clustering.pdf",width=8,height=6)
aa <- pheatmap(final_stat,fontsize_row=8,fontsize_col=8,border="lightgray",cluster_rows=T, show_rownames=T,cutree_col = 4,show_colnames =T,cluster_cols=T,scale="none",clustering_distance_rows = "correlation",breaks=seq(0,1,length=100),color = colorRampPalette(c("steelblue", "white", "red"))(100))
dev.off()
############################################################# step 5
cutree(aa$tree_row,k=4) -> groups
groups[groups==1] <- "mes"
groups[groups==2] <- "pn"
groups[groups==3] <- "oxphos"
groups[groups==4] <- "other"
save(groups,file="NMF_res/NMF_LOOV_cluster.RData")