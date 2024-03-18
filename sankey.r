rm(list=ls())
library("ggalluvial")
library(dplyr)
load("../NMF_res/NMF_LOOV_cluster.RData")
load("Wang_subtype.RData")
load("Neftel_subtype.RData")
load("Garafano_subtype.RData")

dat <- data.frame(Wang=Wang_subtype[names(groups)],NMF=groups,stringsAsFactors=F,row.names=names(groups))
dat %>% filter(NMF != "other") -> dat
dat %>% group_by(NMF,Wang) %>% summarise(count=n()) %>% as.data.frame() -> sankey

     NMF Wang count
1    mes  CLS     1
2    mes  MES    15
3 oxphos  CLS     1
4 oxphos  MES     5
5 oxphos   PN     7
6     pn  CLS    10
7     pn   PN     6

p <- ggplot(data = sankey,
  aes(axis1 = NMF, axis2 = Wang,y = count)) +
  scale_x_discrete(limits = c("NMF", "Wang"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = NMF)) +
  geom_stratum(fill=c("purple","mediumblue","red","purple","red","orange")) +
  geom_text(colour="white",stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("red","mediumblue","purple")) +
  theme_void() + theme(legend.position = 'none')
ggsave(p,file="sankey_NMF_Wang.pdf")

dat <- data.frame(Neftel=Neftel_subtype[names(groups)],NMF=groups,stringsAsFactors=F,row.names=names(groups))
dat %>% filter(NMF != "other") -> dat
dat %>% group_by(NMF,Neftel) %>% summarise(count=n()) %>% as.data.frame() -> sankey

      NMF Neftel count
1     mes   scAC     1
2     mes  scMES    14
3     mes  scNPC     1
4  oxphos   scAC     2
5  oxphos  scMES     9
6  oxphos  scNPC     2
7      pn   scAC     8
8      pn  scMES     1
9      pn  scNPC     4
10     pn  scOPC     3

p <- ggplot(data = sankey,
  aes(axis1 = NMF, axis2 = Neftel,y = count)) +
  scale_x_discrete(limits = c("NMF", "Neftel"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = NMF)) +
  geom_stratum(fill=c("purple","mediumblue","red","blue","steelblue","red","orange")) +
  geom_text(colour="white",stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("red","mediumblue","purple")) +
  theme_void() + theme(legend.position = 'none')
ggsave(p,file="sankey_NMF_Neftel.pdf")

dat <- data.frame(Garafano=Garafano_subtype[names(groups)],NMF=groups,stringsAsFactors=F,row.names=names(groups))
dat %>% filter(NMF != "other") -> dat
dat %>% group_by(NMF,Garafano) %>% summarise(count=n()) %>% as.data.frame() -> sankey

      NMF Garafano count
1    mes      GPM     9
2    mes      MTC     1
3    mes      PPR     6
4 oxphos      GPM     1
5 oxphos      MTC    10
6 oxphos      PPR     2
7     pn      MTC     1
8     pn      NEU     5
9     pn      PPR    10

p <- ggplot(data = sankey,
  aes(axis1 = NMF, axis2 = Garafano,y = count)) +
  scale_x_discrete(limits = c("NMF", "Garafano"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = NMF)) +
  geom_stratum(fill=c("purple","mediumblue","red","red","green","blue","cyan")) +
  geom_text(colour="white",stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("red","mediumblue","purple")) +
  theme_void() + theme(legend.position = 'none')
ggsave(p,file="sankey_NMF_Garafano.pdf")