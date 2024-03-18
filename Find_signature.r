rm(list=ls())
library(data.table)
load("CPM_exp.RData")
CPM_exp -> PDGC_CPM_exp
load("NMF_res/NMF_LOOV_cluster.RData")
###########################################################################
DEAgroups <- function(ddata,
                      groups,
                      method = c("MWW", "t-test")
                      ){
  gg <- sort(unique(groups))
  ans <- vector("list", length(gg))
  names(ans) <- gg
  for(i in 1:length(gg)){
    whichOfInterest <- names(groups)[groups == gg[i]]
    theOthers <- setdiff(colnames(ddata), whichOfInterest)
    diffActivity <- apply(ddata, 1, function(x){
      if(method == "MWW") suppressWarnings(a <- wilcox.test(x[whichOfInterest], x[theOthers]))
      if(method == "t-test") suppressWarnings(a <- t.test(x[whichOfInterest], x[theOthers]))
      a <- c(as.numeric(a$statistic), a$p.value)
      return(a)
    })
    diffActivity <- t(diffActivity)
    colnames(diffActivity) <- c("wTest", "pValue")
    fc <- rowMeans(ddata[, whichOfInterest]) - rowMeans(ddata[, theOthers])
    qValue <- p.adjust(diffActivity[, "pValue"], method = "fdr")
    diffActivity <- data.frame(statistic = diffActivity[, "wTest"], dm = fc, p.value = diffActivity[, "pValue"], fdr = qValue)
    diffActivity <- diffActivity[, -1]
    colnames(diffActivity) <- c("logFC", "pValue", "qValue")
    ans[[i]] <- diffActivity
  }
  return(ans)
}
###########################################################################
findSig <- function(geData,groups){
	groups <- as.factor(groups)
	rankedLists <- vector("list", length(levels(groups)))
    names(rankedLists) <- levels(groups)
    for(gg in 1:length(rankedLists)){
    	whichOfInterest <- names(groups)[groups == names(rankedLists)[gg]]   ## 感兴趣样本的名字
    	theOthers <- setdiff(names(groups), whichOfInterest)	## 其他样本的名字
    	ans <- apply(geData, 1, function(x) wilcox.test(x[whichOfInterest], x[theOthers]))
    	rankedList <- unlist(sapply(ans, function(x) x$statistic))/length(whichOfInterest)/length(theOthers)
    	names(rankedList) <- gsub("\\.W", "", names(rankedList))
    	rankedList <- log2(rankedList/(1-rankedList))
    	rankedList <- sort(rankedList, decreasing = TRUE)
    	rankedLists[[gg]] <- rankedList
    	print(gg)
    }
    return(rankedLists)
}
###########################################################################
DEAgroups(PDGC_CPM_exp,groups,"MWW") -> diff_exp_of_type
findSig(PDGC_CPM_exp,groups) -> signature_of_type

mes <- diff_exp_of_type$mes
colnames(mes) <- paste0("mes_",colnames(mes))
mes$mes_NES <- signature_of_type$mes[rownames(mes)]

pn <- diff_exp_of_type$pn
colnames(pn) <- paste0("pn_",colnames(pn))
pn$pn_NES <- signature_of_type$pn[rownames(pn)]

oxphos <- diff_exp_of_type$oxphos
colnames(oxphos) <- paste0("oxphos_",colnames(oxphos))
oxphos$oxphos_NES <- signature_of_type$oxphos[rownames(oxphos)]

other <- diff_exp_of_type$other
colnames(other) <- paste0("other_",colnames(other))
other$other_NES <- signature_of_type$other[rownames(other)]

cbind(mes,pn,oxphos) -> diff_gene
write.table(diff_gene,file="sig_res/pdgc_diff_gene.txt",sep="\t",quote = F)
save(diff_gene,file="sig_res/pdgc_diff_gene.RData")
