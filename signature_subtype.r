rm(list=ls())
load("CPM_exp.RData")
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
signature_score <- as.data.frame(mod.analyze2(CPM_exp,c('mes','pn','oxphos'),permN=1000))
save(signature_score,file="subtype_res/signature_score.RData")
write.table(signature_score,file="subtype_res/signature_score.txt",sep="\t",quote = F)

## Determine the subtype
type_identifier <- function(score){  ## score was the output of ssGSEA, including norm_score and p_val
   pval <- score[,grep("pval",colnames(score))]
   norm <- score[,grep("norm",colnames(score))]
   apply(pval,1,function(x){which(x==min(x))}) -> tmp
   if(class(tmp)=="list"){
       unlist(lapply(tmp,function(x){if(length(x)==1){return(names(x))};if(length(x)>1){return("None")}})) -> tmp_type
       names(tmp_type)[tmp_type=="None"] -> indx
       tmp_type <- tmp_type[tmp_type!="None"]
       apply(norm[indx,],1,which.max) -> tmp_type2
       tmp_type2 <- colnames(norm)[tmp_type2]
       names(tmp_type2) <- indx
       tmp_type <- c(tmp_type,tmp_type2)
       return(tmp_type)
   }
   else{return(colnames(pval)[tmp])}
}
type_identifier(signature_score) -> subtype
signature_subtype <- gsub("_[a-z].*","",subtype)
save(signature_subtype,file="subtype_res/signature_subtype.RData")

load("NMF_res/NMF_LOOV_cluster.RData")
signature_subtype <- signature_subtype[names(groups)]
data.frame(subtype=signature_subtype,NMF=groups) -> NMF_sig_res
save(NMF_sig_res,file="subtype_res/NMF_sig_res.RData")
