rm(list = ls())

library(dplyr)
library(ggplot2)
library(devtools)
library(rlist)
library(ggord)
library(ggfortify)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(purrr)
library(enrichR)
library(UpSetR)
library(tidyr)
ptrn='conn_HEPG2_'; types=c('kd','oe')
rm( list = Filter( exists, c('dataset') ) )

hepg2_conndata=connectivitydatafile_compiler(ptrn = ptrn,types = types)
# hepg2_conndata[grepl(paste(c('prcp', 'ptpn','bnip3','vabp','isl', 'znf366','tfap2a','gale','ptgs2'), collapse = '|'),hepg2_conndata$Name_type,ignore.case = T),grepl(paste(c('mestranol', 'warfarin'), collapse = '|'), colnames(hepg2_conndata))]
hepg2_conndata[grepl(paste(k, collapse = '|'),hepg2_conndata$Name_type,ignore.case = T),grepl(paste(c('mestranol', 'pantoprazole'), collapse = '|'), colnames(hepg2_conndata))]


##EnrichR analyses
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

kegg="KEGG_2019_Human"
go_bio="GO_Biological_Process_2018"
rnaseq_disease_geo="RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"

alt_dbs = list(kegg=kegg,go_bio=go_bio, rnaseq_disease_geo=rnaseq_disease_geo)
drug_list=colnames(hepg2_conndata)[!colnames(hepg2_conndata)%in%'Name_type']

#pro-coagulant drugs only
drug_pros=drug_list[grepl(paste(c('estradiol','E2','mestranol','icariin','memantine','methylandrostenediol','procarbazine'), collapse = '|'),drug_list)]
#anti-coagulant drugs only
drug_antis=drug_list[grepl(paste(c('brucine','hyperoside','alprazolam','phenothiazine','pantoprazole','rimcazole','trifluoperazine','warfarin'), collapse = '|'),drug_list)]

rm( list = Filter( exists, c('go_output','k','d','pcutoff','dbs','result') ) )

kegg_data_summ=enrichr_output(drug_list = drug_list ,dataset = hepg2_conndata,dbs = alt_dbs[1],cmapranks = T,cmaprank_cutoff = 97)
go_data_summ=enrichr_output(drug_list = drug_list ,dataset = hepg2_conndata,dbs = alt_dbs[2],cmapranks = T,cmaprank_cutoff = 97)
rnaseq_disease_data_summ=enrichr_output(drug_list = drug_list ,dataset = summ_conndata,dbs = alt_dbs[3],cmapranks = T,cmaprank_cutoff = 97)

# y=rbind(kegg_data_summ[[1]][[3]], kegg_data_summ[[3]][[3]])
# y[grepl('p53 signaling pathway',y$Term),grepl(paste(c('mestranol','warfarin'),collapse = '|'),colnames(y))]

# x=rbind(go_data_summ[[2]][[3]], go_data_summ[[3]][[3]])
# x[grepl('negative regulation of intracellular estrogen receptor',x$Term),grepl(paste(c('mestranol','warfarin'),collapse = '|'),colnames(x))]

hepg2_conndata_mestwar=hepg2_conndata[,grepl(paste(c('estradiol','mestranol','warfarin'),collapse = '|'),colnames(hepg2_conndata))]
part_hepg2_conndata_mestwar=hepg2_conndata_mestwar[abs(hepg2_conndata_mestwar$`_estradiol`)>80|abs(hepg2_conndata_mestwar$`_mestranol`)>80,]
partgenes=unique(gsub('_kd|_oe','',rownames(part_hepg2_conndata_mestwar)))

part_KEGG_BOTHpathcounts_mestEstr=read.csv(file = '20210611_part_KEGG_BOTHpathcounts_mestEstr.csv')

KEGG_subpartUP=partgenes[grepl(paste(unlist(strsplit(part_KEGG_BOTHpathcounts_mestEstr$freqGenesOverlap_E2vsDMSO_UP,';')), collapse = '|'),partgenes)]
KEGG_subpartDOWN=partgenes[grepl(paste(unlist(strsplit(part_KEGG_BOTHpathcounts_mestEstr$freqGenesOverlap_E2vsDMSO_DOWN,';')), collapse = '|'),partgenes)]
KEGG_subparts=append(KEGG_subpartUP,KEGG_subpartDOWN)
KEGG_subparts[!duplicated(KEGG_subparts)]

KEGG_p=part_KEGG_BOTHpathcounts_mestEstr[grepl(paste(KEGG_subparts[!duplicated(KEGG_subparts)], collapse = '|'),paste0(part_KEGG_BOTHpathcounts_mestEstr$freqGenesOverlap_proUPDOWN,';',part_KEGG_BOTHpathcounts_mestEstr$freqGenesOverlap_proDOWNUP)),]

KEGG_CMapPros=hepg2_conndata[grepl(paste(KEGG_subparts,collapse = '|'),rownames(hepg2_conndata)),grepl(paste(c('estradiol','mestranol','warfarin'),collapse = '|'),colnames(hepg2_conndata))]
KEGG_CMapPros[KEGG_CMapPros==0]='---'

# write.csv(KEGG_p,file = paste0(gsub('-','',Sys.Date()),'_HepG2_CMapExpressionProfilesVsCMapSignatures_KEGG.csv'),row.names = F)
# write.csv(KEGG_CMapPros,file = paste0(gsub('-','',Sys.Date()),'_HepG2_CMapExpressionProfilesVsCMapSignatures_KEGG_ConnectivityScores.csv'))


part_GO_BOTHpathcounts_mestEstr=read.csv(file = '20210611_part_GO_BOTHpathcounts_mestEstr.csv')

GO_subpartUP=partgenes[grepl(paste(unlist(strsplit(part_GO_BOTHpathcounts_mestEstr$freqGenesOverlap_E2vsDMSO_UP,';')), collapse = '|'),partgenes)]
GO_subpartDOWN=partgenes[grepl(paste(unlist(strsplit(part_GO_BOTHpathcounts_mestEstr$freqGenesOverlap_E2vsDMSO_DOWN,';')), collapse = '|'),partgenes)]
GO_subparts=append(GO_subpartUP,GO_subpartDOWN)
GO_subparts[!duplicated(GO_subparts)]

GO_p=part_GO_BOTHpathcounts_mestEstr[grepl(paste(GO_subparts[!duplicated(GO_subparts)], collapse = '|'),paste0(part_GO_BOTHpathcounts_mestEstr$freqGenesOverlap_proUPDOWN,';',part_GO_BOTHpathcounts_mestEstr$freqGenesOverlap_proDOWNUP)),]

GO_CMapPros=hepg2_conndata[grepl(paste(GO_subparts,collapse = '|'),rownames(hepg2_conndata)),grepl(paste(c('estradiol','mestranol','warfarin'),collapse = '|'),colnames(hepg2_conndata))]
GO_CMapPros[GO_CMapPros==0]='---'

# write.csv(GO_p,file = paste0(gsub('-','',Sys.Date()),'_HepG2_CMapExpressionProfilesVsCMapSignatures_GO.csv'),row.names = F)
# write.csv(GO_CMapPros,file = paste0(gsub('-','',Sys.Date()),'_HepG2_CMapExpressionProfilesVsCMapSignatures_GO_ConnectivityScores.csv'))



kegg_mestwar_summ=slicenricher(enrichrlist = kegg_data_summ,combscore_cutoff = 50,cor_outputs = T,heatmap_outputs = T,cmapranks = T,refcompounds = c('mestranol','warfarin'),enrichmentdb = 'KEGG')
go_mestwar_summ=slicenricher(enrichrlist = go_data_summ,combscore_cutoff = 200,cor_outputs = T,heatmap_outputs = T,cmapranks = T,refcompounds = c('mestranol','warfarin'),enrichmentdb = 'GO_BIO')
rnaseq_disease_geo_mestwar_summ=slicenricher(enrichrlist = rnaseq_disease_data_summ,combscore_cutoff = 200,cor_outputs = T,heatmap_outputs = T,cmapranks = T,refcompounds = c('mestranol','warfarin'),enrichmentdb = 'RNAseq-disease-GEO')


upsetlist=list()

KEGG_upsetlist=upsetlister(slicenrichrlist = kegg_mestwar_summ,refcompounds =c('mestranol','warfarin') )
GO_upsetlist=upsetlister(slicenrichrlist = go_mestwar_summ,refcompounds =c('mestranol','warfarin') )
rnaseq_disease_geo_upsetlist=upsetlister(slicenrichrlist =rnaseq_disease_geo_mestwar_summ,refcompounds =c('mestranol','warfarin') )

library(devtools)
library(Vennerable)
KEGG_temp=Venn(KEGG_upsetlist)  #provide all your groups as list
GO_temp=Venn(GO_upsetlist) 
rnaseq_disease_geo_temp=Venn(rnaseq_disease_geo_upsetlist) 


names(KEGG_upsetlist)=factor(names(KEGG_upsetlist),levels = sort(names(KEGG_upsetlist)))
upset(fromList(KEGG_upsetlist), order.by = "freq",nsets = sum(lengths(KEGG_temp@IntersectionSets)>0))


names(GO_upsetlist)=factor(names(GO_upsetlist),levels = sort(names(GO_upsetlist)))
upset(fromList(GO_upsetlist), order.by = "freq",nsets = sum(lengths(GO_temp@IntersectionSets)>0))


names(rnaseq_disease_geo_upsetlist)=factor(names(rnaseq_disease_geo_upsetlist),levels = sort(names(rnaseq_disease_geo_upsetlist)))
upset(fromList(rnaseq_disease_geo_upsetlist), order.by = "freq",nsets = sum(lengths(rnaseq_disease_geo_temp@IntersectionSets)>0))

KEGG_querylist=list(list("_mestranol_complete", "_mestranol_up"),list("_mestranol_down", "_mestranol_complete"),list("_warfarin_complete","_mestranol_complete",'_warfarin_up','_mestranol_down'))
KEGG_querylist=list(list("_menadione_complete", "_menadione_up"),list("_menadione_down", "_menadione_complete"),list("_warfarin_complete","_menadione_complete",'_warfarin_up','_menadione_down'))

GO_querylist=list(list("_warfarin_complete","_warfarin_down"),list("_mestranol_complete", "_mestranol_up", '_warfarin_up'),list("_warfarin_complete",'_warfarin_up','_warfarin_down'),list("_mestranol_down", "_warfarin_up"))

head(termgenesmajor(enrichrlist = go_data_summ,refcompounds =  c("Mestranol","Warfarin")))

kegg_upsetlegend=upsetlistlegend(upsetlist = KEGG_upsetlist,querylist = KEGG_querylist,enrichrlist = kegg_data_summ,enrichmentdb = 'KEGG',cell_line = 'summ', refcompounds = c('mestranol', 'warfarin'),combscore_cutoff = 50, enrichmentTOOL = 'EnrichR')

go_upsetlegend=upsetlistlegend(upsetlist = GO_upsetlist,querylist = GO_querylist,enrichrlist = go_data_summ,enrichmentdb = 'GO_BIO',cell_line = 'summ', refcompounds = c('mestranol', 'warfarin'),combscore_cutoff = 200, enrichmentTOOL = 'EnrichR')


dataset$alprazolam.x=NULL
colnames(dataset)[colnames(dataset)=="alprazolam.y"]= "alprozalam"
dataset[is.na(dataset)]=0
dataset=dataset[!duplicated(dataset),]
rownames(dataset)=dataset$Name_type
dataset$Name_type=NULL
dataset=as.matrix.data.frame(dataset)

dataset_main=dataset
sim_lim=98
negs= vector()
for (i in 1:nrow(dataset)) {
  if (all(abs(dataset[i,])<sim_lim)) {
    negs=append(negs,i)  
  }
  
}

dataset=dataset[-negs,]
cols <- rep('black', ncol(dataset))
#turn red the specified rows in tf
proco=c("mestranol","memantine","icariin","linopirdine","anagrelide","donepezil","progesterone","testosterone","estradiol")
cols[colnames(dataset)%in%proco]= 'red'

empty.rows = unlist(lapply(row.names(dataset),function(x){a = " "}))

