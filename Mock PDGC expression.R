load('PDGC_CPM_exp.RData')
source('~/Workshop/18-09-12_Multicentric/codes/rfunc/ssgseaMOD.r')

tableS1 = read.table('supplementary_table1.txt',header=T,sep='\t',row.names=1)
PDGC_CPM_exp.f = PDGC_CPM_exp[,which(tableS1$In.house.subtype!='other')]
rs = mod.analyze2(PDGC_CPM_exp.f,c('OXPHOS','PN','MES','Proneural','Classical','Mesenchymal','GPM','MTC','NEU','PPR'),permN=1000)
rs = as.data.frame(rs)
rs.inhouse = type_identifier(rs[,grep('(OXPHOS|PN|MES)',colnames(rs))])
rs.wang    = type_identifier(rs[,grep('(Proneural|Classical|Mesenchymal)',colnames(rs))])
rs.anto    = type_identifier(rs[,grep('(GPM|MTC|NEU|PPR)',colnames(rs))])

rmimic0.9 = as.data.frame(sapply(1:ncol(PDGC_CPM_exp.f),function(i){
  x = as.vector(unlist(PDGC_CPM_exp.f[,i]))
  set.seed(i*10)
  tmp = PDGC_CPM_exp.f[,-i]
  id = sample(ncol(tmp),1)
  y = unlist(tmp[,id])
  x*0.9 + 0.1*y
}))
rrs0.9 = mod.analyze2(rmimic0.9,c('OXPHOS','PN','MES','Proneural','Classical','Mesenchymal','GPM','MTC','NEU','PPR'),permN=1000)
rrs0.9 = as.data.frame(rrs0.9)
rrs0.9.inhouse = type_identifier(rrs0.9[,grep('(OXPHOS|PN|MES)',colnames(rrs0.9))])
rrs0.9.wang = type_identifier(rrs0.9[,grep('(Proneural|Classical|Mesenchymal)',colnames(rrs0.9))])
rrs0.9.anto = type_identifier(rrs0.9[,grep('(GPM|MTC|NEU|PPR)',colnames(rrs0.9))])

rmimic0.8 = as.data.frame(sapply(1:ncol(PDGC_CPM_exp.f),function(i){
  x = as.vector(unlist(PDGC_CPM_exp.f[,i]))
  set.seed(i*10)
  tmp = PDGC_CPM_exp.f[,-i]
  id = sample(ncol(tmp),1)
  y = unlist(tmp[,id])
  x*0.8 + 0.2*y
}))
rrs0.8 = mod.analyze2(rmimic0.8,c('OXPHOS','PN','MES','Proneural','Classical','Mesenchymal','GPM','MTC','NEU','PPR'),permN=1000)
rrs0.8 = as.data.frame(rrs0.8)
rrs0.8.inhouse = type_identifier(rrs0.8[,grep('(OXPHOS|PN|MES)',colnames(rrs0.8))])
rrs0.8.wang = type_identifier(rrs0.8[,grep('(Proneural|Classical|Mesenchymal)',colnames(rrs0.8))])
rrs0.8.anto = type_identifier(rrs0.8[,grep('(GPM|MTC|NEU|PPR)',colnames(rrs0.8))])

rmimic0.7 = as.data.frame(sapply(1:ncol(PDGC_CPM_exp.f),function(i){
  x = as.vector(unlist(PDGC_CPM_exp.f[,i]))
  set.seed(i*10)
  tmp = PDGC_CPM_exp.f[,-i]
  id = sample(ncol(tmp),1)
  y = unlist(tmp[,id])
  x*0.7 + 0.3*y
}))
rrs0.7 = mod.analyze2(rmimic0.7,c('OXPHOS','PN','MES','Proneural','Classical','Mesenchymal','GPM','MTC','NEU','PPR'),permN=1000)
rrs0.7 = as.data.frame(rrs0.7)
rrs0.7.inhouse = type_identifier(rrs0.7[,grep('(OXPHOS|PN|MES)',colnames(rrs0.7))])
rrs0.7.wang = type_identifier(rrs0.7[,grep('(Proneural|Classical|Mesenchymal)',colnames(rrs0.7))])
rrs0.7.anto = type_identifier(rrs0.7[,grep('(GPM|MTC|NEU|PPR)',colnames(rrs0.7))])

rmimic0.6 = as.data.frame(sapply(1:ncol(PDGC_CPM_exp.f),function(i){
  x = as.vector(unlist(PDGC_CPM_exp.f[,i]))
  set.seed(i*10)
  tmp = PDGC_CPM_exp.f[,-i]
  id = sample(ncol(tmp),1)
  y = unlist(tmp[,id])
  x*0.6 + 0.4*y
}))
rrs0.6 = mod.analyze2(rmimic0.6,c('OXPHOS','PN','MES','Proneural','Classical','Mesenchymal','GPM','MTC','NEU','PPR'),permN=1000)
rrs0.6 = as.data.frame(rrs0.6)
rrs0.6.inhouse = type_identifier(rrs0.6[,grep('(OXPHOS|PN|MES)',colnames(rrs0.6))])
rrs0.6.wang = type_identifier(rrs0.6[,grep('(Proneural|Classical|Mesenchymal)',colnames(rrs0.6))])
rrs0.6.anto = type_identifier(rrs0.6[,grep('(GPM|MTC|NEU|PPR)',colnames(rrs0.6))])

rmimic.stat = data.frame(
  proportion = c(0.6,0.7,0.8,0.9),
  inhouse = c((length(rs.inhouse) - sum(diag(table(gsub('_.*$','',rs.inhouse),gsub('_.*$','',rrs0.6.inhouse)))))/length(rs.inhouse),
              (length(rs.inhouse) - sum(diag(table(gsub('_.*$','',rs.inhouse),gsub('_.*$','',rrs0.7.inhouse)))))/length(rs.inhouse),
              (length(rs.inhouse) - sum(diag(table(gsub('_.*$','',rs.inhouse),gsub('_.*$','',rrs0.8.inhouse)))))/length(rs.inhouse),
              (length(rs.inhouse) - sum(diag(table(gsub('_.*$','',rs.inhouse),gsub('_.*$','',rrs0.9.inhouse)))))/length(rs.inhouse)
  ),
  wang = c((length(rs.wang) - sum(diag(table(gsub('_.*$','',rs.wang),gsub('_.*$','',rrs0.6.wang)))))/length(rs.wang),
           (length(rs.wang) - sum(diag(table(gsub('_.*$','',rs.wang),gsub('_.*$','',rrs0.7.wang)))))/length(rs.wang),
           (length(rs.wang) - sum(diag(table(gsub('_.*$','',rs.wang),gsub('_.*$','',rrs0.8.wang)))))/length(rs.wang),
           (length(rs.wang) - sum(diag(table(gsub('_.*$','',rs.wang),gsub('_.*$','',rrs0.9.wang)))))/length(rs.wang)
  ),
  anto = c((length(rs.anto) - sum(diag(table(gsub('_.*$','',rs.anto),gsub('_.*$','',rrs0.6.anto)))))/length(rs.anto),
           (length(rs.anto) - sum(diag(table(gsub('_.*$','',rs.anto),gsub('_.*$','',rrs0.7.anto)))))/length(rs.anto),
           (length(rs.anto) - sum(diag(table(gsub('_.*$','',rs.anto),gsub('_.*$','',rrs0.8.anto)))))/length(rs.anto),
           (length(rs.anto) - sum(diag(table(gsub('_.*$','',rs.anto),gsub('_.*$','',rrs0.9.anto)))))/length(rs.anto)
  )
)

plot(1-rmimic.stat$proportion,rmimic.stat$inhouse,type='o',col='red',ylim=c(0,1),ylab='inconsistence',xlab='mixture')
lines(1-rmimic.stat$proportion,rmimic.stat$wang,type='o',col='blue',ylim=c(0,1))
lines(1-rmimic.stat$proportion,rmimic.stat$anto,type='o',col='green',ylim=c(0,1))
