rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggalluvial)
read.table("GBM_TCGA/HT_HG-U133A",header=T,row.names=1,sep="\t") -> U133A_exp
U133A_exp <- U133A_exp[,-(grep("11$",colnames(U133A_exp)))]
library(ssgsea.GBM.classification)
source("~/ssgseaMOD.r")
Wang_score <- as.data.frame(mod.analyze2(U133A_exp,c('MES','PN','CLS'),permN=1000))
save(Wang_score,file="TCGA_U133A_Wang_score.RData")

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

load("TCGA_U133A_Wang_score.RData")
type_identifier(Wang_score) -> Wang_subtype
load("TCGA_U133A_sig_score.RData")
type_identifier(score) -> Chen_subtype

dat <- data.frame(Chen=Chen_subtype,Wang=Wang_subtype[names(Chen_subtype)])
dat <- dat %>% filter(Chen!="None") %>% filter(Wang!="None")
dat %>% group_by(Chen,Wang) %>% summarise(count=n()) %>% as.data.frame() -> sankey
sankey$Wang <- paste("Wang",sankey$Wang,sep="-")
p <- ggplot(data = sankey,
  aes(axis1 = Chen, axis2 = Wang,y = count)) +
  scale_x_discrete(limits = c("Chen", "Wang"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = Chen)) +
  geom_stratum(fill=c("purple","mediumblue","red","purple","red","orange")) +
  geom_text(colour="white",stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("red","mediumblue","purple")) +
  theme_void() + theme(legend.position = 'none')
ggsave(p,filename = "sankey.pdf")


     Chen  Wang    count
1    mes Wang-CLS    17
2    mes Wang-MES    69
3    mes  Wang-PN     2
4 oxphos Wang-CLS    30
5 oxphos Wang-MES    32
6 oxphos  Wang-PN    30
7     pn Wang-CLS    80
8     pn Wang-MES     4
9     pn  Wang-PN    62