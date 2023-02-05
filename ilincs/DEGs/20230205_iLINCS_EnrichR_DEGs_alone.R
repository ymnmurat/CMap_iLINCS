### EnrichR analyses for iLINCS data
#HA1E cell line
#Final (3rd) upgrade to MY   May 19, 2021

rm(list = ls())

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

if (exists('dataset')) {rm(dataset)}
files=list.files(pattern="iLINCS_complete_HA")
for (i in files) {  
  
  if (!exists("dataset")){
    dataset <- read.table(i,header = T,sep = "\t")[,c("Name_GeneSymbol", "Value_LogDiffExp","Significance_pvalue")]
    colnames(dataset)[colnames(dataset)=="Value_LogDiffExp"] = sub(".xls", "_logDE",sub(".*\\_HA1E_","",i));
    colnames(dataset)[colnames(dataset)=="Significance_pvalue"] = sub(".xls", "_pval",sub(".*\\_HA1E_","",i));
    colnames(dataset)[colnames(dataset)=='Name_GeneSymbol']='Genes'
  }
  if (exists("dataset")){
    temp_dataset <-read.table(i,header = T,sep = "\t")[,c("Name_GeneSymbol", "Value_LogDiffExp","Significance_pvalue")]
    colnames(temp_dataset)[colnames(temp_dataset)=="Value_LogDiffExp"] = sub(".xls", "_logDE",sub(".*\\_HA1E_","",i));
    colnames(temp_dataset)[colnames(temp_dataset)=="Significance_pvalue"] = sub(".xls", "_pval",sub(".*\\_HA1E_","",i));
    colnames(temp_dataset)[colnames(temp_dataset)=='Name_GeneSymbol']='Genes'
    if (!sum(colnames(temp_dataset)%in%colnames(dataset))==ncol(temp_dataset)) {
      dataset<-merge(dataset, temp_dataset, by="Genes",all= T)
    }
    rm(temp_dataset)
  }
  
}

# dataset=dataset[,-c(2:3)]
# colnames(dataset)[colnames(dataset)=="24_7-[4-(2-Hydroxyphenyl)-2-(trifluoromethyl)-1,3-dioxan-5-yl]hept-5-enoic acid_10uM_logDE"]= "24_Hydroxyphenyl, _10uM_logDE"
# colnames(dataset)[colnames(dataset)=="24_7-[4-(2-Hydroxyphenyl)-2-(trifluoromethyl)-1,3-dioxan-5-yl]hept-5-enoic acid_10uM_pval"]= "24_Hydroxyphenyl, _10uM_pval"
# colnames(dataset)[colnames(dataset)=="6_2-Propenoic acid, 3-[4-(1H-imidazol-1-ylmethyl)phenyl]-, (Z)-_10uM_logDE"]= "6_2-Propenoic acid, _10uM_logDE"
# colnames(dataset)[colnames(dataset)=="6_2-Propenoic acid, 3-[4-(1H-imidazol-1-ylmethyl)phenyl]-, (Z)-_10uM_pval"]= "6_2-Propenoic acid, _10uM_pval"
# colnames(dataset)[colnames(dataset)=="24_ETHYLENEDIAMINE TETRACETATE TETRA-ETHANOLAMINE_10uM_logDE"]= "24_ETHYLENEDIAMINE, _10uM_logDE"
# colnames(dataset)[colnames(dataset)=="24_ETHYLENEDIAMINE TETRACETATE TETRA-ETHANOLAMINE_10uM_pval"]= "24_ETHYLENEDIAMINE, _10uM_pval"
# colnames(dataset)[grepl('\\.y',colnames(dataset))]=sub('\\.y','',colnames(dataset)[grepl('\\.y',colnames(dataset))])

rownames(dataset)= dataset$Genes; dataset$Genes=NULL
dataset_pvals= dataset[,grepl('pval', colnames(dataset))]
dataset_logDE=dataset[which(apply(dataset_pvals, 1, function(x) any(x<=0.05))==T),grepl('logDE', colnames(dataset))]
dataset_logDE=as.matrix.data.frame(dataset_logDE)
dataset_main=dataset_logDE

##EnrichR analyses
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

kegg="KEGG_2019_Human"
go_bio="GO_Biological_Process_2018"

alt_dbs = list(kegg=kegg,go_bio=go_bio)

      # #EnrichR test
      # if (websiteLive) {
      #   data_mst_tmp=dataset_mest[dataset_mest$adj.P.Val<0.05,]
      #   enriched <- enrichr(data_mst_tmp[order(data_mst_tmp$adj.P.Val),1], alt_dbs[1])
      # }
      # library(openxlsx)
      # enriched[[kegg]]
      # write.xlsx(enriched[[go_bio]],file = '24h_mestranol_EnrichR_all.xlsx')
      # if (websiteLive) enriched[[kegg]]

drug_list=unique(sub('_logDE|_pval','',colnames(dataset)))

rm( list = Filter( exists, c('go_output','k','d','pcutoff','dbs','result') ) )
k=list();

#enrichr_output: an implemented function to retrieve EnrichR results. It requires a predetermined drug list, p-value cut off, name of the database used (KEGG, GO_BIO etc), and a specific direction (for the genes either up or down regulated, if inquired). If the direction is not specified, the calculations are done for all the DEGs + ups only + downs only). The output contains the lists for the EnrichR results in the specified directions separately. Each list contains three subcategories in data frame format: (1) List of pathways and associated combination scores, (2) Expanded list that contains the full range of EnrichR results, and (3) List of pathways and associated genes coming from the input gene expression data.

enrichr_output= function(drug_list,dataset,logDEcutoff=0,pcutoff,dbs,direction=NULL, cmapranks=F, cmaprank_cutoff){
  require(dplyr);require(purrr);require(enrichR)
  result=list()
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
    while (sum(grepl(ds, c('complete','up','down'),ignore.case = T))==0) { ds=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: complete/up/down?')))}
  }
  for (d in ds) {
    if (exists('go_output')) {rm(go_output)}
    if (exists('termgenes')) {rm(termgenes)}
    k=list()
    if (!cmapranks) {
     for (drugs in drug_list) {
      if (!exists("go_output")){
        go_input= dataset[,grepl(drugs, sub('_logDE|_pval','',colnames(dataset)))]
        go_input= go_input[order(go_input[,2]),]
        if (d=='complete') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&abs(go_input[,1])>abs(logDEcutoff),]), dbs)}
        else if (d=='up') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]>abs(logDEcutoff),]), dbs)}
        else {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]<(-abs(logDEcutoff)),]), dbs)}
        k[[drugs]]=enriched[[dbs[[1]]]]
        go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
        termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
        # rownames(go_output)=go_output$Term; go_output$Term=NULL
        colnames(go_output)=c('Term',drugs)
        colnames(termgenes)=c('Term',drugs)
        rm( list = Filter( exists, c('go_input','enriched') ) )}
      if (exists("go_output")){
        go_input <-dataset[,grepl(drugs, sub('_logDE|_pval','',colnames(dataset)))]
        go_input= go_input[order(go_input[,2]),]
        if (d=='complete') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&abs(go_input[,1])>abs(logDEcutoff),]), dbs)
        } else if (d=='up') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]>abs(logDEcutoff),]), dbs)
        } else {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]<(-abs(logDEcutoff)),]), dbs)}
        k[[drugs]]=enriched[[dbs[[1]]]]
        enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
        enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
        # rownames(enriched)=enriched$Term; enriched$Term=NULL
        colnames(enriched_go)=c('Term',drugs)
        colnames(enriched_go_termgenes)=c('Term',drugs)
        go_output= merge(go_output, enriched_go, by='Term',all= T)
        termgenes= merge(termgenes, enriched_go_termgenes, by='Term',all= T)
        rm( list = Filter( exists, c('go_input','enriched_go') ) )}}
    } else if(cmapranks) {
      for (drugs in drug_list) {
        if (!exists("go_output")){
          go_input= dataset[,c(1,grep(drugs,colnames(dataset)))]
          if (d=='complete') {go_input= go_input[grepl('_kd|_oe',go_input[,1]),]
            upranks=go_input[abs(go_input[,2])>cmaprank_cutoff,]
            upranks[order(abs(upranks[,2]),decreasing = T),]
            enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
            rm( list = Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {go_input= go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_oe',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_kd',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing = T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list = Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {go_input= go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_kd',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_oe',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing = T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list = Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          # rownames(go_output)=go_output$Term; go_output$Term=NULL
          colnames(go_output)=c('Term',drugs)
          colnames(termgenes)=c('Term',drugs)
          rm( list = Filter( exists, c('go_input','enriched') ) )}
        if (exists("go_output")){
          go_input= dataset[,c(1,grep(drugs,colnames(dataset)))]
          if (d=='complete') {upranks=go_input[abs(go_input[,2])>cmaprank_cutoff,]
          upranks[order(abs(upranks[,2]),decreasing = T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list = Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {go_input= go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_oe',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_kd',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing = T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list = Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {go_input= go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_kd',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_oe',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing = T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list = Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          # rownames(enriched)=enriched$Term; enriched$Term=NULL
          colnames(enriched_go)=c('Term',drugs)
          colnames(enriched_go_termgenes)=c('Term',drugs)
          go_output= merge(go_output, enriched_go, by='Term',all= T)
          termgenes= merge(termgenes, enriched_go_termgenes, by='Term',all= T)
          rm( list = Filter( exists, c('go_input','enriched_go') ) )}}
    }
    go_output=go_output[,-2]
    colnames(go_output)[grepl('\\.y',colnames(go_output))]=gsub('\\.y','',colnames(go_output)[grepl('\\.y',colnames(go_output))])
    colnames(go_output)[grepl('\\.',colnames(go_output))]=gsub('\\.','_',colnames(go_output)[grepl('\\.',colnames(go_output))])
    # rownames(go_output)=go_output$Term; go_output$Term=NULL
    go_output[is.na(go_output)]=0
    termgenes=termgenes[,-2]
    colnames(termgenes)[grepl('\\.y',colnames(termgenes))]=gsub('\\.y','',colnames(termgenes)[grepl('\\.y',colnames(termgenes))])
    colnames(termgenes)[grepl('\\.',colnames(termgenes))]=gsub('\\.','_',colnames(termgenes)[grepl('\\.',colnames(termgenes))])
    # rownames(go_output)=go_output$Term; go_output$Term=NULL
    termgenes[is.na(termgenes)]='---'
    result[[d]]=list(go_output, k, termgenes)
    names(result[[d]])=c(paste0(names(dbs),'_EnrichedTermsScores_',d),paste0(names(dbs),'_Expanded_',d),paste0(names(dbs),'_TermGenes_',d))
    rm(k)}
  for (d in ds) {
    for (i in 1:length(result[[d]])) {
      names(result[[d]][[i]])=gsub('\\.','_',names(result[[d]][[i]]))}}
return(result)}

keggtr=enrichr_output(drug_list = drug_list ,dataset = dataset,pcutoff = 0.05,dbs = alt_dbs[1])
gotr=enrichr_output(drug_list = drug_list ,dataset = dataset,pcutoff = 0.05,dbs = alt_dbs[2])


# enrichrlist=keggtr; direction="all"; combscore_cutoff=200;refcompounds = c("Mestranol","Warfarin");d='complete'

#slicenricher: an implemented function to slice and collect specified data from the enrichr_output function. Currently, it only utilizes the EnrichedTermsScores objects from the enrichr_output data. Specifications rely on the combscore cutoff values for each reference compound one by one (the inquiries can be provided directly, or the function will continuously prompt to the user to provide compounds, until the user indicates not to. The main outputfile contains a list for the pathway names passing the threshold and their combination scores for the respective compounds inquired by the user. Lists are given in separate sections for each direction (up/down/complete)). In addition, the correlation plots (by default 'spearman') and heatmaps are produced and saved in the current directory, unless indicated otherwise.  
slicenricher=function(enrichrlist, direction=NULL, combscore_cutoff=200, refcompounds=NULL, cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL){
  require(tidyverse); require(ggplot2); require(ggord);require(ggfortify);
  require(ComplexHeatmap);require(circlize);
  if (!cmapranks) {drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][1])))))  
  } else if (cmapranks) {drugnames=gsub('.*_','',names(enrichrlist[[1]][[1]]))}
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case = T))==0) { ds=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: complete/up/down?')))}
  }
  resp='y'
  while (is.null(refcompounds)) {refcompounds=as.character(readline(prompt = 'Please provide a valid compound name from the list: '))
  while (sum(grepl(paste(refcompounds,collapse = "|"), drugnames,ignore.case = T))==0) { refcompounds=as.character(readline(prompt = 'Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt = 'Name of the next compound: ')))
      while (sum(grepl(paste(nextcomp,collapse = "|"), drugnames,ignore.case = T))==0) { nextcomp=as.character(readline(prompt = 'Please recheck the compound name, and provide it once again: '))}
    refcompounds=append(refcompounds, nextcomp)
    resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))
    while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))}}
  }
  
  if (is.null(enrichmentdb)) {enrichmentdb=as.character(readline(prompt = 'Name of the gene enrichment database: '))
  while (length(enrichmentdb)==0) { enrichmentdb=as.character(readline(prompt = 'Please provide the name of the gene enrichment database that you have used: '))}
  }
  proco=c("mestranol","memantine","icariin","linopirdine","anagrelide","donepezil","progesterone","testosterone","Ethynyl Estradiol",'EthynylEstradiol','estradiol','17.betaestradiol','E2')
  antico_thr=c('2-Propenoic acid,','Hydroxyphenyl,','APIXABAN','Betrixaban','BRD-A49838158 ','Dabigatran Etexilate','ETHYLENEDIAMINE, ',
               'Gabexate','ICI 192,605','L 655240','Nafamostat Mesylate','Ozagrel','Picotamide','Skatole','SQ-29548','Torsemide','SQ.29548','L.655240')
  proco_thr=c('Desmopressin','Menadione','U-46619','U.46619')
  bloods=c('blood', 'thromb', 'agula', 'erythroc', 'leukoc', 'arteri','vein ', 'vascular')
  
  result=list() 
  corrplots=list() 
  heatmaps=list() 
  for(d in ds){tmp=as.data.frame(enrichrlist[names(enrichrlist)%in%d][[1]][1])
    colnames(tmp)=gsub(paste0(names(enrichrlist[names(enrichrlist)%in%d][[1]][1]),'.'),'',colnames(tmp))
    rownames(tmp)=tmp$Term; 
    ind=which(grepl(paste(refcompounds,collapse = "|"), colnames(tmp),ignore.case = T))
    pathways=vector()
    for (i in ind) {pathways= unique(append(pathways,rownames(tmp[which(tmp[,i]>combscore_cutoff),])))}
    tmp=tmp[rownames(tmp)%in%pathways,]
    tmp[,-1]=log2(tmp[,-1]) #Log scale the combination scores 
    tmp[tmp<0]=0
    result[[d]]=tmp
    tmp=tmp[,-1]
    tmp=tmp[,colSums(tmp)>0]
    if (length(tmp)>0) {
    
    #corplot
    cormat_go<-signif(cor(tmp,method = cormethod),2)  
    cols <- rep('black', ncol(tmp))
    cols[grepl(paste(proco, collapse = '|'), colnames(tmp),ignore.case = T)]= 'red4'
    cols[grepl(paste(antico_thr, collapse = '|'), colnames(tmp),ignore.case = T)]='dimgray'
    cols[grepl(paste(proco_thr, collapse = '|'), colnames(tmp),ignore.case = T)]='red'
    col_paths=rep('black',nrow(tmp))
    col_paths[grepl(paste(bloods, collapse = '|'), rownames(tmp),ignore.case = T)]= 'blue'
    if (cor_outputs) {
    par(cex.main=1)
    lgd_list=Legend(labels = c("Pro-coagulants", "Anti-coagulants", 'iLINCS: Pro-coagulants','iLINCS: Anti-coagulants'), title = "Compound types", type = "points", pch = 18, legend_gp = gpar(col = c("red4","black",'red', 'dimgray'), fontsize=24)) 
    f1 = colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length = 3), c("#3f7fff", "#EEEEEE", "red"))
    jpeg(paste0(deparse(substitute(enrichrlist)),'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substring(paste(refcompounds,collapse = '-'),1,15),toupper(cormethod),'corplot_',toupper(d),".jpeg"), units="in", width=18, height=12, res=300)
    ht1=Heatmap(cormat_go,clustering_method_rows = "ward.D",clustering_method_columns  = "ward.D", col = f1, name= "Pathway combination score correlations (iLINCS - EnrichR)", row_names_gp = gpar(fontsize = 8, fontface = "bold", col=cols),column_names_gp = gpar(fontsize = 8, fontface = "bold", col=cols),
                heatmap_legend_param=list(legend_direction="horizontal", title_position= "topcenter", grid_height = unit(0.6, "cm"),legend_width = unit(8, "cm"),labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold")),
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                  grid.text(cormat_go[i, j], x, y,gp = gpar(fontsize=4.5, fontface='bold'))})
    tmpfig= draw(ht1, heatmap_legend_side = "bottom", column_title= paste0("Corplot - ",d,'\n',cormethod," correlations on iLINCS pathways \n (EnrichR ",paste(refcompounds,collapse = '-'), " Combination Scores >", combscore_cutoff,", n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line = 5), annotation_legend_list = lgd_list,annotation_legend_side = "right")
    k=paste0('Corplot','_',d)
    corrplots[[k]]<-tmpfig
    dev.off()}
    
    #heatmaps
    if (heatmap_outputs) {
    f2 = colorRamp2(seq(max(tmp), min(tmp), length=3),  c("#3f7fff", "#EEEEEE", "red"))
    jpeg(paste0(deparse(substitute(enrichrlist)),'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substring(paste(refcompounds,collapse = '-'),1,15),"_Heatmap_",toupper(d),".jpeg"), units="in", width=18, height=12, res=300)
    
    ht2=Heatmap(as.matrix.data.frame(tmp),clustering_method_rows = "ward.D",clustering_method_columns  = "ward.D", col = f2, name = paste0(toupper(d), " - Pathway combinations scores (EnrichR ",paste(refcompounds,collapse = '-'), " Combination Scores >", combscore_cutoff,", n=", nrow(tmp),")"), row_names_gp = gpar(fontsize = 10-nrow(tmp)*0.05, fontface = "bold", col=col_paths),column_names_gp = gpar(fontsize = 12, fontface = "bold",col=cols),
                       heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height = unit(0.8, "cm"),legend_width = unit(10, "cm"),labels_gp = gpar(fontsize = 12),title_gp = gpar(fontsize = 18, fontface = "bold"),title_position = "topcenter"),use_raster = F)
    tmpheat=draw(ht2, heatmap_legend_side = "top",  column_title_gp=gpar(fontface="bold", line = 5),annotation_legend_list = lgd_list,annotation_legend_side = "bottom")
    k=paste0('Heatmap','_',d)
    heatmaps[[k]]<-tmpheat
    dev.off()
    rm(ind); rm(pathways); rm(tmp)
    }}}
  rm(drugnames)
  return(result)
}  

slicenricher_updown=function(enrichrlist, combscore_cutoff=200, contr1,contr2,classContr1=NULL,classContr2=NULL,cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL, repositoryName=NULL){
  require(tidyverse); require(ggplot2); require(ggord);require(ggfortify);
  require(ComplexHeatmap);require(circlize);
  if (!cmapranks) {drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][1])))))  
  } else if (cmapranks) {drugnames=gsub('.*_','',names(enrichrlist[[1]][[1]]))}
  if (is.null(enrichmentdb)) {enrichmentdb=as.character(readline(prompt = 'Name of the gene enrichment database: '))
  while (length(enrichmentdb)==0) { enrichmentdb=as.character(readline(prompt = 'Please provide the name of the gene enrichment database that you have used: '))}
  }
  if (is.null(classContr1)) {classContr1=as.character(readline(prompt = 'What is the class of the first contrasts (Pro-agulant / Anti-coagulant etc): '))
  while (length(classContr1)==0) { classContr1=as.character(readline(prompt = 'What is the class of the first contrasts (Pro-agulant / Anti-coagulant etc): '))}}
  if (is.null(classContr2)) {classContr2=as.character(readline(prompt = 'What is the class of the second contrasts (Pro-agulant / Anti-coagulant etc): '))
  while (length(classContr2)==0) { classContr2=as.character(readline(prompt = 'What is the class of the second contrasts (Pro-agulant / Anti-coagulant etc): '))}}
  if (is.null(repositoryName)) {repositoryName=as.character(readline(prompt = 'What is the source of the data (CMap, iLINCS etc): '))
  while (length(repositoryName)==0) { repositoryName=as.character(readline(prompt = 'What is the source of the data (CMap, iLINCS etc): '))}}
  proco=c("mestranol","memantine","icariin","linopirdine","anagrelide","donepezil","progesterone","testosterone","Ethynyl Estradiol",'EthynylEstradiol','estradiol','17.betaestradiol','E2')
  antico_thr=c('2-Propenoic acid,','Hydroxyphenyl,','APIXABAN','Betrixaban','BRD-A49838158 ','Dabigatran Etexilate','ETHYLENEDIAMINE, ',
               'Gabexate','ICI 192,605','L 655240','Nafamostat Mesylate','Ozagrel','Picotamide','Skatole','SQ-29548','Torsemide','SQ.29548','L.655240')
  proco_thr=c('Desmopressin','Menadione','U-46619','U.46619')
  bloods=c('blood', 'thromb', 'agula', 'erythroc', 'leukoc', 'arteri','vein ', 'vascular')
  
  ds=list(c('up','down'),c('down','up'))
  result=list() 
  corrplots=list() 
  heatmaps=list() 
  for(d in 1:length(ds)){tmp1=as.data.frame(enrichrlist[names(enrichrlist)%in%ds[[d]][1]][[1]][1])
  colnames(tmp1)=paste0(gsub(paste0(names(enrichrlist[names(enrichrlist)%in%ds[[d]][1]][[1]][1]),'.'),'',colnames(tmp1)),'_',ds[[d]][1])
  colnames(tmp1)[1]='Term'
  rownames(tmp1)=tmp1$Term;
  tmp1=tmp1[,grepl(paste(c('Term',contr1),collapse = '|'),colnames(tmp1))]
  tmp2=as.data.frame(enrichrlist[names(enrichrlist)%in%ds[[d]][2]][[1]][1])
  colnames(tmp2)=paste0(gsub(paste0(names(enrichrlist[names(enrichrlist)%in%ds[[d]][2]][[1]][1]),'.'),'',colnames(tmp2)),'_',ds[[d]][2])
  colnames(tmp2)[1]='Term'
  rownames(tmp2)=tmp2$Term;
  tmp2=tmp2[,grepl(paste(c('Term',contr2),collapse = '|'),colnames(tmp2))]
  tmp=merge(tmp1,tmp2,by='Term',all = T)
  rownames(tmp)=tmp$Term
  ind=which(grepl(paste(c(contr1,contr2),collapse = "|"), colnames(tmp),ignore.case = T))
  pathways=vector()
  for (i in ind) {pathways= unique(append(pathways,rownames(tmp[which(tmp[,i]>combscore_cutoff),])))}
  tmp=tmp[rownames(tmp)%in%pathways,]
  tmp[,-1]=log2(tmp[,-1]) #Log scale the combination scores 
  tmp[tmp<0]=0
  result[[paste(c(ds[[d]][1],ds[[d]][2]),collapse = '_')]]=tmp
  tmp=tmp[,-1]
  tmp=tmp[,colSums(tmp)>0]
  if (length(tmp)>0) {
    
    #corplot
    cormat_go<-signif(cor(tmp,method = cormethod),2)  
    cols <- rep('black', ncol(tmp))
    cols[grepl(paste(proco, collapse = '|'), colnames(tmp),ignore.case = T)]= 'red4'
    cols[grepl(paste(antico_thr, collapse = '|'), colnames(tmp),ignore.case = T)]='dimgray'
    cols[grepl(paste(proco_thr, collapse = '|'), colnames(tmp),ignore.case = T)]='red'
    col_paths=rep('black',nrow(tmp))
    col_paths[grepl(paste(bloods, collapse = '|'), rownames(tmp),ignore.case = T)]= 'blue'
    if (cor_outputs) {
      par(cex.main=1)
      lgd_list=Legend(labels = c("Pro-coagulants", "Anti-coagulants", 'iLINCS: Pro-coagulants','iLINCS: Anti-coagulants'), title = "Compound types", type = "points", pch = 18, legend_gp = gpar(col = c("red4","black",'red', 'dimgray'), fontsize=24)) 
      f1 = colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length = 3), c("#3f7fff", "#EEEEEE", "red"))
      jpeg(paste0(deparse(substitute(enrichrlist)),'_',paste(c(repositoryName,enrichmentdb,enrichmentTOOL),collapse = '-'),'_',substr(paste(c(contr1,contr2),collapse = '-'),1,30),'_',toupper(cormethod),'corplot_',toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse = '_')),".jpeg"), units="in", width=18, height=12, res=300)
      ht1=Heatmap(cormat_go,clustering_method_rows = "ward.D",clustering_method_columns  = "ward.D", col = f1, name= paste0("Pathway combination score correlations (",repositoryName," - ", enrichmentTOOL,")"), row_names_gp = gpar(fontsize = 8, fontface = "bold", col=cols),column_names_gp = gpar(fontsize = 8, fontface = "bold", col=cols),
                  heatmap_legend_param=list(legend_direction="horizontal", title_position= "topcenter", grid_height = unit(0.6, "cm"),legend_width = unit(8, "cm"),labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold")),
                  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                    grid.text(cormat_go[i, j], x, y,gp = gpar(fontsize=4.5, fontface='bold'))})
      tmpfig= draw(ht1, heatmap_legend_side = "bottom", column_title= paste0("Corplot - ",toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse = 'vs')),'\n',cormethod," correlations on ",repositoryName,'-',enrichmentdb," pathways \n (", enrichmentTOOL," ",substr(paste(c(contr1,contr2),collapse = '-'),1,30), " Combination Scores >", combscore_cutoff,", n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line = 5), annotation_legend_list = lgd_list,annotation_legend_side = "right")
      k=paste0('Corplot','_',d)
      corrplots[[k]]<-tmpfig
      dev.off()}
    
    #heatmaps
    if (heatmap_outputs) {
      f2 = colorRamp2(seq(max(tmp), min(tmp), length=3),  c("#3f7fff", "#EEEEEE", "red"))
      jpeg(paste0(deparse(substitute(enrichrlist)),'_',repositoryName,'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substr(paste(c(contr1,contr2),collapse = '-'),1,30),'____Heatmap_',toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse = '_')),".jpeg"), units="in", width=18, height=12, res=300)
      
      ht2=Heatmap(as.matrix.data.frame(tmp),clustering_method_rows = "ward.D",clustering_method_columns  = "ward.D", col = f2, name = paste0(toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse = 'vs')), " - Pathway combinations scores (",paste(c(repositoryName,enrichmentdb,enrichmentTOOL),collapse = '-')," ",substr(paste(c(contr1,contr2),collapse = '-'),1,30), "... Combination Scores >", combscore_cutoff,", n=", nrow(tmp),")"), row_names_gp = gpar(fontsize = 10-nrow(tmp)*0.05, fontface = "bold", col=col_paths),column_names_gp = gpar(fontsize = 12, fontface = "bold",col=cols), heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height = unit(0.8, "cm"),legend_width = unit(10, "cm"),labels_gp = gpar(fontsize = 12),title_gp = gpar(fontsize = 10, fontface = "bold"),title_position = "topcenter"),use_raster = F)
      tmpheat=draw(ht2, heatmap_legend_side = "top",  column_title_gp=gpar(fontface="bold", line = 5),annotation_legend_list = lgd_list,annotation_legend_side = "bottom")
      k=paste0('Heatmap','_',d)
      heatmaps[[k]]<-tmpheat
      dev.off()
      rm(ind); rm(pathways); rm(tmp)
    }}}
  rm(drugnames)
  return(result)
}  

kegg_mestwar2=slicenricher(enrichrlist = keggtr,combscore_cutoff = 20,refcompounds = drug_list,enrichmentdb = 'KEGG',enrichmentTOOL = 'EnrichR')
go_mestwar2=slicenricher(enrichrlist = gotr,combscore_cutoff = 200,refcompounds = drug_list,enrichmentdb = 'GO_BIO',enrichmentTOOL = 'EnrichR')
#Complete evaluations

# KEGG_combscore_cutoff=100
# GO_combscore_cutoff=200
# mestwar=c('24_mestranol_10uM','6_mestranol_10uM','24_warfarin_004uM')

upsetlist=list()


##upsetlister, an implemented function to provide the input for the upset plots and venn interaction terms. Slicenricher outputs are preferred in producing the upsetlists
upsetlister=function(slicenrichrlist,refcompounds=NULL,direction=NULL, combscore_cutoff=200,updowncontrasts=F){
  drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(slicenrichrlist[[1]])))))
  if (!updowncontrasts) {
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case = T))==0) { ds=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: complete/up/down?')))}
  }}else if (updowncontrasts) {
    if (is.null(direction)) {ds=c('up_down','down_up')}
    else {ds=tolower(direction)
    while (sum(grepl(ds, c('up_down','down_up'),ignore.case = T))==0) { ds=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: up_down/down_up?')))}
    }}
  resp='y'
  if (is.null(refcompounds)) {refcompounds=as.character(readline(prompt = 'Please provide a valid compound name from the list: '))
  while (sum(grepl(paste(refcompounds,collapse = "|"), drugnames,ignore.case = T))==0) { refcompounds=as.character(readline(prompt = 'Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt = 'Name of the next compound: ')))
  while (sum(grepl(paste(nextcomp,collapse = "|"), drugnames,ignore.case = T))==0) { nextcomp=as.character(readline(prompt = 'Please recheck the compound name, and provide it once again: '))}
  refcompounds=append(refcompounds, nextcomp)
  resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))}}
  }
  upsetlist=list();  
  for(d in ds){tmp=as.data.frame(slicenrichrlist[names(slicenrichrlist)%in%d][[1]])
  ind=which(grepl(paste(refcompounds,collapse = "|"), colnames(tmp),ignore.case = T))  
  if (!updowncontrasts) {colnames(tmp)[ind]=paste0(colnames(tmp[,ind]),'_',d)}
  for (i in ind) {comp=list(rownames(tmp[tmp[,i]>log2(combscore_cutoff),]))
      names(comp)=colnames(tmp)[i]
      upsetlist=append(upsetlist, comp)
      }
    }  
    return(upsetlist)}

KEGG_upsetlist=upsetlister(slicenrichrlist = kegg_mestwar2, combscore_cutoff = 20,refcompounds = c('menadione','warfarin', 'skatole'))
GO_upsetlist=upsetlister(slicenrichrlist = go_mestwar2,refcompounds = c('mestranol','warfarin'))

library(devtools)
# install_github("js229/Vennerable"); BiocManager::install('RBGL')
library(Vennerable)
KEGG_temp=Venn(KEGG_upsetlist)  #provide all your groups as list
GO_temp=Venn(GO_upsetlist) 


names(KEGG_upsetlist)=factor(names(KEGG_upsetlist),levels = sort(names(KEGG_upsetlist)))
upset(fromList(KEGG_upsetlist), order.by = "freq",nsets = sum(lengths(KEGG_temp@IntersectionSets)>0))


names(GO_upsetlist)=factor(names(GO_upsetlist),levels = sort(names(GO_upsetlist)))
upset(fromList(GO_upsetlist), order.by = "freq",nsets = sum(lengths(GO_temp@IntersectionSets)>0))


# BiocManager::install(c("RBGL","graph"), force=T, dependencies=T)
# install.packages("Vennerable", repos="http://R-Forge.R-project.org")
 #provide all your groups as list
names(KEGG_upsetlist)
names(GO_upsetlist)
#access the elements of a single intersection list
KEGG_temp@IntersectionSets[paste(as.integer(names(KEGG_upsetlist) %in% c("24_warfarin_004uM_up", "6_mestranol_10uM_down")), collapse = '')]
KEGG_temp@IntersectionSets[paste(as.integer(names(KEGG_upsetlist) %in% c("24_Menadione_3uM_up", '24_warfarin_004uM_down')), collapse = '')]
GO_temp@IntersectionSets[paste(as.integer(names(KEGG_upsetlist) %in% c("24_warfarin_004uM_down", "6_mestranol_10uM_up")), collapse = '')]

KEGG_querylist=list(list("6_mestranol_10uM_up", '6_mestranol_10uM_complete'),list("6_mestranol_10uM_down", '6_mestranol_10uM_complete'),list("24_warfarin_004uM_down","24_warfarin_004uM_complete"),list("24_mestranol_10uM_up","24_warfarin_004uM_up",'6_mestranol_10uM_down'),list("6_mestranol_10uM_complete","24_mestranol_10uM_down"),list("6_mestranol_10uM_up", '6_mestranol_10uM_complete','6_mestranol_10uM_down'))
KEGG_querylist=list(list("24_Menadione_3uM_up", '24_warfarin_004uM_down'),list("6_mestranol_10uM_down", '6_mestranol_10uM_complete'),list("24_warfarin_004uM_down","24_warfarin_004uM_complete"),list("24_mestranol_10uM_up","24_warfarin_004uM_up",'6_mestranol_10uM_down'),list("6_mestranol_10uM_complete","24_mestranol_10uM_down"),list("6_mestranol_10uM_up", '6_mestranol_10uM_complete','6_mestranol_10uM_down'))

GO_querylist=list(list("24_warfarin_004uM_down","6_mestranol_10uM_up"),list("24_warfarin_004uM_down","6_mestranol_10uM_down"),list("24_warfarin_004uM_down","24_warfarin_004uM_complete","24_mestranol_10uM_down"), list("24_warfarin_004uM_down", "6_mestranol_10uM_up", '6_mestranol_10uM_complete','24_warfarin_004uM_complete'),list("24_warfarin_004uM_down", "6_mestranol_10uM_up",'24_mestranol_10uM_up'))


#paths_list, an implemented function retrieving intersections of the pathway names across the user defined queries and upsetlists. Further implemented in the upsetlistlegend function
paths_list= function(myqueries, allqueries,intersectionlist){
  require(devtools); require(Vennerable)
  unlist(intersectionlist[paste(as.numeric(allqueries %in% myqueries),collapse = "")])}

#termgenesmajor, an implemented function retrieving the names of the DEGs relating with the enrichr pathway outputs. Further utilized in the upsetlistlegend function
termgenesmajor=function(enrichrlist,direction=NULL, refcompounds=NULL){
  drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][3])))))
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case = T))==0) { ds=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: complete/up/down?')))}
  }
  resp='y'
  if (is.null(refcompounds)) {refcompounds=as.character(readline(prompt = 'Please provide a valid compound name from the list: '))
  while (sum(grepl(paste(refcompounds,collapse = "|"), drugnames,ignore.case = T))==0) { refcompounds=as.character(readline(prompt = 'Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt = 'Name of the next compound:')))
  while (sum(grepl(paste(nextcomp,collapse = "|"), drugnames,ignore.case = T))==0) { nextcomp=as.character(readline(prompt = 'Please recheck the compound name, and provide it once again: '))}
  refcompounds=append(refcompounds, nextcomp)
  resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more compounds: yes (Y) / no (N)   ')))}}
  }
  majorlist=setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("Term"))
  for(d in ds){tmp=enrichrlist[[d]][[3]]
  ind=which(grepl(paste(refcompounds,collapse = "|"), colnames(tmp),ignore.case = T))  
  colnames(tmp)[ind]=paste0(colnames(tmp[,ind]),'_',d)
  majorlist=merge(majorlist, as.data.frame(tmp[,c(1,ind)]), by='Term', all=T)
  }
  return(majorlist)
}

head(termgenesmajor(enrichrlist = keggtr,refcompounds =  c("Mestranol","Warfarin")))

#upsetlistlegend, an implemented function for producing the upsetplots with appropriate legend annotations that are defined by the users' querylists. In addition, pathways of interests and their respective DEGs are extracted as excel sheets to the working directory.
upsetlistlegend=function(upsetlist, querylist, enrichrlist,outputlist=NULL,colorpal='Dark2', ylabel='No. of Pathways', cell_line='NA',combscore_cutoff='NA',enrichmentdb=NULL,refcompounds=NULL,upsetFigure=T,enrichmentTOOL,repositoryName='iLINCS'){
  require(stringr)
  require(Vennerable)
  require(RColorBrewer)
  require(UpSetR)
  require(openxlsx)
  
  temp=Venn(upsetlist)  #provide all your groups as list
  
  names(upsetlist)=factor(names(upsetlist),levels = sort(names(upsetlist)))
  
  # Define the number of colors you want
  nb.cols <- length(querylist)
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
  
  pathtables=termgenesmajor(enrichrlist, direction = NULL,refcompounds =  refcompounds)
  gridtable=setNames(data.frame(matrix(ncol = ncol(pathtables), nrow = 0)),colnames(pathtables))
  queries=list()
  
  for (i in 1:length(querylist)) {
    queries[[i]]=list(query=intersects,params=querylist[[i]], color=mycolors[i], active=T)}
  
  cnts=0
  if (upsetFigure) {
    
  jpeg(paste0(gsub('-','',Sys.Date()),"_",as.character(repositoryName),'_',as.character(enrichmentTOOL),cell_line,'_',paste(substr(str_to_title(refcompounds),1,4),collapse = ''),"_CombScore",combscore_cutoff,"_",enrichmentdb,"_UpsetPlot.jpeg"), units="in", width=14, height=10, res=300)
  
upsetpl=upset(fromList(upsetlist), nsets = sum(lengths(temp@IntersectionSets)>0), queries =queries, order.by = "freq", mainbar.y.label = ylabel)  
print(upsetpl)
grid.text(paste0(enrichmentdb,' Upset plot - ', paste(str_to_title(refcompounds),collapse = ' & '),' Pathways (Combination Score >', combscore_cutoff,')'),x = 0.4, y=0.98, gp=gpar(fontsize=10, fontface='bold'),just = "left")}

  for (i in 1:length(querylist)) {
      plist=paths_list(myqueries = unlist(querylist[[i]]), allqueries = names(upsetlist), temp@IntersectionSets)
      tmptable=pathtables[pathtables$Term %in% unlist(plist),]
      gridtable=rbind(gridtable,tmptable)
      
      if (upsetFigure) {for (j in 1:length(plist)) {
        grid.text(plist[[j]],x = 0.78, y=(0.96-cnts*0.011), gp=gpar(fontsize=6.5, col=mycolors[i], fontface='bold'),just = "left")
        cnts=cnts+1}}}
  if (upsetFigure) {grid.polygon(x=c(0.75,1,1,0.75), y=c((0.95-cnts*0.011),(0.95-cnts*0.011),0.97,0.97), gp = gpar(fill="gray", alpha=0.2))
  dev.off()}
  write.xlsx(x = gridtable,file = paste0(gsub('-','',Sys.Date()),"_",as.character(repositoryName),'_',as.character(enrichmentTOOL),'_',cell_line,'_',paste(substr(str_to_title(refcompounds),1,4),collapse = ''),"_CombScore",combscore_cutoff,"_",enrichmentdb,"_UpsetPlotLegend.xlsx"), row.names = F)
  return(gridtable)
}

querylister=function(intersectionSets,upsetlist, lowerlimit=0, upperlimit=max(lengths(intersectionSets))){qs=list()
ques=intersectionSets[lengths(intersectionSets)>lowerlimit & lengths(intersectionSets)<=max(lengths(intersectionSets))]
quenames=names(ques[order(lengths(ques),decreasing = T)])
for (ints in quenames) {
  ints_membs=strsplit(ints, '')[[1]]
  tmplist=list()
  for (j in 1:length(ints_membs)) {
    if(as.numeric(ints_membs[j])){tmplist=append(tmplist, names(upsetlist[j])) }}
  qs[[ints]]=tmplist}
rm(ques);rm(quenames);names(qs)=NULL;
return(qs)}

KEGG_querylist=querylister(upsetlist = KEGG_upsetlist,intersectionSets = KEGG_temp@IntersectionSets)

kegg_upsetlegend=upsetlistlegend(upsetlist = KEGG_upsetlist,querylist = KEGG_querylist,enrichrlist = keggtr,enrichmentdb = 'KEGG',cell_line = 'HA1E', refcompounds = c('menadione', 'warfarin'),enrichmentTOOL = 'EnrichR',combscore_cutoff = 20)

go_upsetlegend=upsetlistlegend(upsetlist = GO_upsetlist,querylist = GO_querylist,enrichrlist = gotr,enrichmentdb = 'GO_BIO',cell_line = 'HA1E', refcompounds = c('mestranol', 'warfarin'),enrichmentTOOL = 'EnrichR',combscore_cutoff = 200)

freqfunc <- function(x, n=30){
  tail(sort(table(unlist(strsplit(as.character(x), ";")))), n)
}
# library(RColorBrewer)
# # Define the number of colors you want
# nb.cols_kegg <- length(KEGG_querylist)
# nb.cols_go <- length(GO_querylist)
# mycolors_kegg <- colorRampPalette(brewer.pal(nb.cols_kegg, "Dark2"))(nb.cols_kegg)
# mycolors_go <- colorRampPalette(brewer.pal(nb.cols_go, "Dark2"))(nb.cols_go)
# 
# kegg_queries=list()
# go_queries=list()
# 
# for (i in 1:length(KEGG_querylist)) {
#   kegg_queries[[i]]=list(query=intersects,params=KEGG_querylist[[i]], color=mycolors_kegg[i], active=T)}
# 
# for (i in 1:length(GO_querylist)) {
#   go_queries[[i]]=list(query=intersects,params=GO_querylist[[i]], color=mycolors_go[i], active=T)}
# 
# paths_list= function(myqueries, allqueries,intersectionlist){
#   unlist(intersectionlist[paste(as.numeric(allqueries %in% myqueries),collapse = "")])}
# 
# as.vector(paths_list(myqueries = unlist(KEGG_querylist[[1]]), allqueries = names(KEGG_upsetlist),KEGG_temp@IntersectionSets))[1:2]
# 
# KEGG_temp@IntersectionSets$`000001010`  #params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up")
# KEGG_temp@IntersectionSets$``  #params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up", '6_mestranol_10uM_together','24_warfarin_004uM_together')
# 
# #number of columns in the upsetplot
# sum(lengths(KEGG_temp@IntersectionSets)>0)
# 
# 
# #UPDATED Upset plot with pathway annotations for GO_BIO
# jpeg("20210523_iLINCS_EnrichR_HA1E_mestWarSign_GO_Bio_UpsetPlot.jpeg", units="in", width=14, height=10, res=300)
# upset(fromList(GO_upsetlist), sets = names(GO_upsetlist),  queries = go_queries, order.by = "freq", mainbar.y.label = 'No. of Pathways')
# cnts=0
# grid.text(paste0('GO_BIO Upset plot - Mestranol & Warfarin Pathways (Combination Score >', combscore_cutoff,')'),x = 0.4, y=0.98, gp=gpar(fontsize=12, fontface='bold'),just = "left")
# for (i in 1:length(GO_querylist)) {
#   plist=paths_list(myqueries = unlist(GO_querylist[[i]]), allqueries = names(GO_upsetlist),GO_temp@IntersectionSets)
#   for (j in 1:length(plist)) {
#     grid.text(plist[[j]],x = 0.65, y=(0.9-cnts*0.02), gp=gpar(fontsize=10, col=mycolors_go[i], fontface='bold'),just = "left")
#     cnts=cnts+1
#   }
# }
# grid.polygon(x=c(0.63,1,1,0.63), y=c((0.9-cnts*0.02),(0.9-cnts*0.02),0.92,0.92), gp = gpar(fill="gray", alpha=0.2))
# dev.off()
# 
# #UPDATED Upset plot with pathway annotations for KEGG
# jpeg("20210524_iLINCS_EnrichR_HA1E_mestWarSign_KEGG_UpsetPlot.jpeg", units="in", width=14, height=10, res=300)
# upset(fromList(KEGG_upsetlist), nsets = sum(lengths(KEGG_temp@IntersectionSets)>0), queries =kegg_queries, order.by = "freq", mainbar.y.label = 'No. of Pathways')
# 
# 
# 
# uplistlegendary=function(upsetlist, querylist,enrichmentdb=NULL,refcompounds=NULL,combscore_cutoff=NULL){
# grid.text(paste0('KEGG - Upset plot - Mestranol & Warfarin Pathways (Combination Score >', combscore_cutoff,')'),x = 0.4, y=0.98, gp=gpar(fontsize=12, fontface='bold'),just = "left")
# cnts=0
# for (i in 1:length(KEGG_querylist)) {
#   plist=paths_list(myqueries = unlist(KEGG_querylist[[i]]), allqueries = names(KEGG_upsetlist),KEGG_temp@IntersectionSets)
#   for (j in 1:length(plist)) {
#     grid.text(plist[[j]],x = 0.65, y=(0.9-cnts*0.02), gp=gpar(fontsize=10, col=mycolors_kegg[i], fontface='bold'),just = "left")
#     cnts=cnts+1
#   }
# }
# grid.polygon(x=c(0.63,1,1,0.63), y=c((0.9-cnts*0.02),(0.9-cnts*0.02),0.92,0.92), gp = gpar(fill="gray", alpha=0.2))
# dev.off()
# }
# 
# 
# 
# # queries = list(list(query = intersects, params = list("24_warfarin_004uM_down","6_mestranol_10uM_up"), color = mycolors[1], active = T), list(query = intersects, params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up", '6_mestranol_10uM_complete','24_warfarin_004uM_complete'), color = mycolors[2], active = T), list(query = intersects, params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up",'24_mestranol_10uM_up'),color=mycolors[3], active = T),list(query = intersects, params = list("24_warfarin_004uM_down","6_mestranol_10uM_down"),color=mycolors[4], active = T),list(query = intersects, params = list("24_warfarin_004uM_down","24_warfarin_004uM_complete","24_mestranol_10uM_down"),color=mycolors[5], active = T))
# 
# # upset(fromList(upsetlist), sets = names(upsetlist),  queries = list(list(query = intersects, params = list("24_warfarin_004uM_down","6_mestranol_10uM_up"), color = "blue", active = T),
# #         list(query = intersects, params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up", '6_mestranol_10uM_together','24_warfarin_004uM_together'), color = "red", active = T), 
# #         list(query = intersects, params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up",'24_mestranol_10uM_up'),color='orange', active = T)),
# #         order.by = "freq", mainbar.y.label = 'No. of Pathways')
# 
# 
# temp=Venn(upsetlist)  #provide all your groups as list
# names(upsetlist)
# 
# #access the elements of a single intersection list
# temp@IntersectionSets[paste(as.integer(names(upsetlist) %in% c("24_warfarin_004uM_down", "6_mestranol_10uM_up")), collapse = '')]
# 
# 
# querylist=list(list("24_warfarin_004uM_down","6_mestranol_10uM_up"),list("24_warfarin_004uM_down", "6_mestranol_10uM_up", '6_mestranol_10uM_complete','24_warfarin_004uM_complete'),
#                list("24_warfarin_004uM_down", "6_mestranol_10uM_up",'24_mestranol_10uM_up'),
#                list("24_warfarin_004uM_down","6_mestranol_10uM_down"),
#                list("24_warfarin_004uM_down","24_warfarin_004uM_complete","24_mestranol_10uM_down"))
# 
# library(RColorBrewer)
# # Define the number of colors you want
# nb.cols <- length(querylist)
# mycolors <- colorRampPalette(brewer.pal(nb.cols, "Set1"))(nb.cols)
# queries=list()
# 
# for (i in 1:length(querylist)) {
#   queries[[i]]=list(query=intersects,params=querylist[[i]], color=mycolors[i], active=T)}
# 
# paths_list= function(myqueries, allqueries,intersectionlist){
#   unlist(intersectionlist[paste(as.numeric(allqueries %in% myqueries),collapse = "")])}
# 
# paths_list(myqueries = unlist(querylist[[1]]), allqueries = names(upsetlist),temp@IntersectionSets)[[1]]
# 
# temp@IntersectionSets$`000001010`  #params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up")
# temp@IntersectionSets$``  #params = list("24_warfarin_004uM_down", "6_mestranol_10uM_up", '6_mestranol_10uM_together','24_warfarin_004uM_together')
# 
# #number of columns in the upsetplot
# sum(lengths(temp@IntersectionSets)>0)
# 
# 
# #UPDATED Upset plot with pathway annotations
# jpeg("20210523_iLINCS_EnrichR_HA1E_mestWarSign_GO_Bio_UpsetPlot.jpeg", units="in", width=14, height=10, res=300)
# upset(fromList(upsetlist), sets = names(upsetlist),  queries = queries, order.by = "freq", mainbar.y.label = 'No. of Pathways')
# cnts=0
# grid.text(paste0('Upset plot - Mestranol & Warfarin Pathways (Combination Score >', combscore_cutoff,')'),x = 0.4, y=0.98, gp=gpar(fontsize=12, fontface='bold'),just = "left")
# for (i in 1:length(querylist)) {
#   plist=paths_list(myqueries = unlist(querylist[[i]]), allqueries = names(upsetlist),temp@IntersectionSets)
#   for (j in 1:length(plist)) {
#     grid.text(plist[[j]],x = 0.65, y=(0.9-cnts*0.02), gp=gpar(fontsize=10, col=mycolors[i], fontface='bold'),just = "left")
#     cnts=cnts+1
#   }
# }
# grid.polygon(x=c(0.63,1,1,0.63), y=c((0.9-cnts*0.02),(0.9-cnts*0.02),0.92,0.92), gp = gpar(fill="gray", alpha=0.2))
# dev.off()
# 
# 
