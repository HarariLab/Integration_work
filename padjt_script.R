library(openxlsx)
library(enrichR)
library(stringr)

setwd("~/integration_project/step_16/part_2/")


files = system('ls /home/a.garg/integration_project/step_15/microglia*xlsx',intern=T)
# read and save into excel
pages = lapply(files,function(f){
      name<- gsub(".*/ *", "\\1\\2",f)
})

# padj matrix --------------------------------------------------------------------------------------------

Hits = data.frame(Term = NA)

tmp = lapply(pages,function(p){
 d = read.xlsx('/home/a.garg/integration_project/step_15/Microglia_nebula_ouput.xlsx',sheet=p)
# d = d[d$logFC_cluster > 0,]
 d = d[d$logFC_cluster > 0.25,]

 if(nrow(d) >= 2){
  r = enrichr(d$gene,c('KEGG_2021_Human','GO_Biological_Process_2021'))
  r = rbind(r$KEGG_2021_Human, r$GO_Biological_Process_2021)
  r[,str_remove(p, ".xlsx")] = r$Adjusted.P.value
  Hits <<- merge(Hits, r[,c('Term',str_remove(p, ".xlsx"))], by='Term', all=T)
 }
})


outfile = '~/integration_project/step_16/part_2/EnricherConsensus_cellStates_upOnly_padj.csv'
Hits = Hits[!is.na(Hits$Term),]
rownames(Hits) = Hits$Term
Hits = Hits[,-1]
write.table(Hits,outfile,sep=',')
