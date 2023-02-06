#' Data analysis workflow functions for the xls/txt documents that are retrieved from CMap/iLINCS/SwissTargetPrediction platforms
#'
#' @param ptrn Defaults to _conn_HA1E_, just change it the way you kept your data, preferably in seperate folder to avoid confusion
#' @param types Defaults to kd and oe
#' @keywords mylib
#' connectivitydatafile_compiler()

##connectivitydatafile_compiler: Compiles gene connectivity data retrieved from CMap.Currently, it works locally and utilizes names of the downloaded files to parse through. 
##Future formats can include a logical variable whether the data should be accessed locally or from CMap directly. 
##Accordingly, curl commands can be posted to and retrieved from CMap Clue.
connectivitydatafile_compiler=function(ptrn=NULL, types=NULL){
  if(is.null(ptrn)){pattern='_conn_HA1E_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt = 'Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}

  if (is.null(types)) {typewise=c('kd','oe')
  }  else {typewise=tolower(types)
  while (sum(grepl(paste(typewise, collapse = '|'), c('kd','oe','cp','cc'),ignore.case = T))==0) { typewise=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: knockdown (kd)/overexpression (oe) / compound (cp) / CMap Class (cc)?  ')))}
  }
  files=list.files(pattern = pattern)
  if (sum(typewise%in%'cp')==0) {for (i in files) {  if (!exists("dataset")){
    dataset <- read.table(i,header = T,sep = "\t")[,c("Name", "Type","Score")]
    colnames(dataset)[colnames(dataset)=="Score"] = sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
    dataset=dataset[dataset$Type%in%typewise,]
    dataset$Name_type=paste0(dataset$Name,"_",dataset$Type);
    dataset$Name_type=gsub(' ','_',dataset$Name_type)
    dataset$Name=NULL; dataset$Type=NULL
  }
    if (exists("dataset")){
      temp_dataset <-read.table(i,header = T,sep = "\t")[,c("Name", "Type","Score")]
      colnames(temp_dataset)[colnames(temp_dataset)=="Score"] = sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
      temp_dataset$Name_type=paste0(temp_dataset$Name,"_",temp_dataset$Type);
      temp_dataset=temp_dataset[temp_dataset$Type%in%typewise,]
      temp_dataset$Name_type=gsub(' ','_',temp_dataset$Name_type)
      temp_dataset$Name=NULL; temp_dataset$Type=NULL
      dataset<-merge(dataset, temp_dataset, by="Name_type",all= T)
      rm(temp_dataset)
    }}} else {for (i in files) {

      if (!exists("dataset")){
        dataset <- read.table(i,header = T,sep = "\t")[,c("Name", "Type","Score",'Description')]
        colnames(dataset)[colnames(dataset)=="Score"] = sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
        dataset=dataset[dataset$Type=='cp',]
        dataset$Name_type=paste0(dataset$Name,"_",dataset$Description);
        dataset$Name_type=gsub(' ','_',dataset$Name_type)
        dataset$Name=NULL; dataset$Type=NULL; dataset$Description=NULL
      }
      if (exists("dataset")){
        temp_dataset <-read.table(i,header = T,sep = "\t")[,c("Name", "Type","Score")]
        colnames(temp_dataset)[colnames(temp_dataset)=="Score"] = sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
        temp_dataset=temp_dataset[temp_dataset$Type=='cp']
        temp_dataset$Name_type=paste0(temp_dataset$Name,"_",temp_dataset$Description);
        temp_dataset$Name_type=gsub(' ','_',temp_dataset$Name_type)
        temp_dataset$Name=NULL; temp_dataset$Type=NULL;dataset$Desctiption=NULL
        dataset<-merge(dataset, temp_dataset, by="Name_type",all= T)
        rm(temp_dataset)
      }}}
  dataset=dataset[,-2]
  colnames(dataset)[grepl('\\.y',colnames(dataset))]=gsub('\\.y','',colnames(dataset)[grepl('\\.y',colnames(dataset))])
  colnames(dataset)[grepl('\\.',colnames(dataset))]=gsub('\\.','_',colnames(dataset)[grepl('\\.',colnames(dataset))])
  dataset=dataset[!duplicated(dataset),]
  rownames(dataset)=dataset$Name_type
  dataset[is.na(dataset)]=0
  return(dataset)
}

##iLINCS_kd_connectivitydatafile_compiler: Similar to connectivitydatafile_compiler function, and compiles data from iLINCS instead. 
##Further developments by online queries are needed. So that, connectivity information would be collected for specified compounds.
iLINCS_kd_connectivitydatafile_compiler=function(ptrn=NULL){
  if(is.null(ptrn)){pattern='_conn_HA1E_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt = 'Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}
  files=list.files(pattern = pattern)
  for (i in files) {  if (!exists("dataset")){
    dataset <- read.table(i,header = T,sep = "\t")[,c("GeneTarget", "Correlation","zScore")]
    dataset[dataset=='+']=1; dataset[dataset=='-']=-1
    dataset$zScore=dataset$zScore*as.numeric(dataset$Correlation)
    colnames(dataset)[colnames(dataset)=="zScore"] = sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
    dataset$Correlation=NULL
  }
    if (exists("dataset")){
      temp_dataset <-read.table(i,header = T,sep = "\t")[,c("GeneTarget", "Correlation","zScore")]
      temp_dataset[temp_dataset=='+']=1; temp_dataset[temp_dataset=='-']=-1
      temp_dataset$zScore=temp_dataset$zScore*as.numeric(temp_dataset$Correlation)
      colnames(temp_dataset)[colnames(temp_dataset)=="zScore"] = sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
      temp_dataset$Correlation=NULL
      dataset<-merge(dataset, temp_dataset, by="GeneTarget",all= T)
      rm(temp_dataset)
    }}
  dataset=dataset[,-2]
  colnames(dataset)[grepl('\\.y',colnames(dataset))]=gsub('\\.y','',colnames(dataset)[grepl('\\.y',colnames(dataset))])
  colnames(dataset)[grepl('\\.',colnames(dataset))]=gsub('\\.','_',colnames(dataset)[grepl('\\.',colnames(dataset))])
  dataset=dataset[!duplicated(dataset),]
  rownames(dataset)=dataset$Name
  # dataset[is.na(dataset)]=0
  return(dataset)
}

#" enrichr_output: Retrieves EnrichR results. It requires a vector of drug list, p-value cut off, name of the database used (kegg="KEGG_2019_Human","GO_Biological_Process_2018" (see EnrichR sources)), 
##and a specific direction (for the genes either up or down regulated, if inquired). If the direction is not specified, the calculations are done for all the DEGs (up and downregulated ones separately). 
##The output file contains the lists for the EnrichR results in the specified directions separately. Each list contains three subcategories in data frame format: 
##(1) List of pathways and associated combination scores, (2) Expanded list that contains the full range of EnrichR results, and (3) List of pathways and associated genes coming from the input gene expression data.

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

# enrichrlist=keggtr; direction="all"; combscore_cutoff=200;refcompounds = c("Mestranol","Warfarin");d='complete'

#slicenricher: Slices and collects specified data (based on Combined score cutoff) from the enrichr_output function. 
##Currently, it only utilizes the EnrichedTermsScores objects from the enrichr_output data. 
##Specifications rely on the combscore cutoff values for each reference compound one by one (the inquiries can be provided 
##directly, or the function will continuously prompt to the user to provide compounds, until the user indicates not to. 
##Main outputfile contains a list for the pathway names passing the threshold and their combination scores for the respective compounds inquired by the user. 
##Lists are given in separate sections for each direction (up/down/complete)). In addition, the correlation plots (by default 'spearman') and 
##heatmaps are produced and saved in the current directory, unless indicated otherwise.
##Heatmap label aesthetics are also improved to indicate whether the drugs of interest are coagulant in the literature and/or in our in-house phenotypic screenings. 
##Pathways onthologies that matches with thrombosis keyword patterns are also manually included. A warning function can be included to change the cutoff value, in case there are not many enriched pathways, hence errors based on non-finite calculations
slicenricher=function(enrichrlist, direction=NULL, combscore_cutoff=200, refcompounds=NULL, cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL='EnrichR'){
  require(tidyverse); require(ggplot2); require(ggord);require(ggfortify);
  require(ComplexHeatmap);require(circlize);
  if (!cmapranks) {drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][1])))))
  } else if (cmapranks) {drugnames=gsub('.*_','',names(enrichrlist[[1]][[1]]))}
  if (is.null(direction)) {ds=c('complete','up','down')}
    else {ds=tolower(direction)
    while (sum(grepl(ds, c('complete','up','down'),ignore.case = T))==0) { 
      ds=tolower(as.character(readline(prompt = 'Please check out misspelling, if any: complete/up/down?')))}
  }
  resp='y'
  while (is.null(refcompounds)) {refcompounds=as.character(readline(prompt = 'Please provide a valid compound name from the list: '))
    while (sum(grepl(paste(refcompounds,collapse = "|"), drugnames,ignore.case = T))==0) { 
      refcompounds=as.character(readline(prompt = 'Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
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
  proco=c("heparin","thrombin",'phylloquinone', 'menadione')
  antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
  proco_thr=c("heparin","thrombin",'phylloquinone', 'menadione')
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

#slicenricher_updown: Similar to slicenricher. It allows comparing different groups of drugs (pro-coagulants vs anti-coagulants) to see whether there are any pathways that move in different directions
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
  proco=c("heparin","thrombin",'phylloquinone', 'menadione')
  antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
  proco_thr=c("heparin","thrombin",'phylloquinone', 'menadione')
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
        cols[grepl(pattern = paste0(proco, collapse = '|'), x= colnames(tmp),ignore.case = T)]= 'red4'
        cols[grepl(pattern = paste0(antico_thr, collapse = '|'), x=colnames(tmp),ignore.case = T)]='dimgray'
        cols[grepl(pattern = paste0(proco_thr, collapse = '|'), x = colnames(tmp),ignore.case = T)]='red'
        col_paths=rep('black',nrow(tmp))
        col_paths[grepl(pattern = paste0(bloods, collapse = '|'), x = rownames(tmp),ignore.case = T)]= 'blue'
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



##upsetlister, an implemented function to provide the input for the upset plots and venn interaction terms. So that, from many pathways in the list, users can see how many of them are shared among the reference compounds and across different directions. If the lists do not seem complicated (meaning that they are not enriched with multiple pathways and directions) upsetlist and further plots will not be needed.Slicenricher outputs are preferred and transformed into upset plot convenient list forms while collecting pathways that are above a certain threshold. Reference compounds and regulation direction (up, down or complete) can be further specified.
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


#paths_list: Retrieves intersections of the pathway names across the user defined queries and upsetlists. Further implemented in the upsetlistlegend function
paths_list= function(myqueries, allqueries,intersectionlist){
  require(devtools); require(Vennerable)
  unlist(intersectionlist[paste(as.numeric(allqueries %in% myqueries),collapse = "")])}


#termgenesmajor, retrieves the names of the differential genes relating with the previous enrichr pathway outputs. 
##Hence, it allows seeing overexpression or knockdown of which genes (from CMap/iLINCS) can be further important. 
##The function is also nested inside the upsetlistlegend function
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

#querylister allows conversion of queries to compatible formats to be processed along with upsetlists and upset plots. Requires intersectionSets from Venn function (see vennerable)
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


#upsetlistlegend creates the upsetplots with appropriate legend annotations that are defined by the users' querylists. 
##In addition, pathways of interests and their respective genes from CMap/iLINCS studies are extracted in excel sheets to the working directory.
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

#Adopted function (reference needed)
freqfunc2 <- function(x, n=30){
  tail(sort(table( unlist(strsplit(c(na.omit(as.character(unlist(x, use.names=FALSE)))), ";",fixed = T)))), n)
}

#Adopted function (reference needed)
strsplit2 <- function(x,
                     split,
                     type = "remove",
                     perl = FALSE,
                     ...) {
  if (type == "remove") {
    # use base::strsplit
    out <- base::strsplit(x = x, split = split, perl = perl, ...)
  } else if (type == "before") {
    # split before the delimiter and keep it
    out <- base::strsplit(x = x,
                          split = paste0("(?<=.)(?=", split, ")"),
                          perl = TRUE,
                          ...)
  } else if (type == "after") {
    # split after the delimiter and keep it
    out <- base::strsplit(x = x,
                          split = paste0("(?<=", split, ")"),
                          perl = TRUE,
                          ...)
  } else {
    # wrong type input
    stop("type must be remove, after or before!")
  }
  return(out)
}

###pairs with p-val (imported function (Michael Love)). Online link is needed for reference
cols_BR = brewer.pal(11, "RdBu")   # goes from red to white to blue
pal = colorRampPalette(cols_BR)
cor_colors = data.frame(correlation = seq(-1,1,0.01),
                        correlation_color = pal(201)[1:201])  # assigns a color for each r correlation value
cor_colors$correlation_color = as.character(cor_colors$correlation_color)

panel.cor <- function(x, y, digits=2, cex.cor)
{
  par(usr = c(0, 1, 0, 1))
  u <- par('usr')
  names(u) <- c("xleft", "xright", "ybottom", "ytop")
  r <- cor(x, y,method="pearson",use="complete.obs")
  test <- cor.test(x,y)
  bgcolor = cor_colors[2+(-r+1)*100,2]    # converts correlation into a specific color
  do.call(rect, c(col = bgcolor, as.list(u))) # colors the correlation box

  if (test$p.value> 0.05){
    text(0.5,0.5,"Insignificant",cex=1.5)
  } else{
    text(0.5, 0.75, paste("r=",round(r,2)),cex=2.5) # prints correlatoin coefficient
    text(.5, .25, paste("p=",formatC(test$p.value, format = "e", digits = 1)),cex=2)
    abline(h = 0.5, lty = 2) # draws a line between correlatoin coefficient and p value
  }

}
panel.smooth<-function (x, y, col = "black", bg = NA, pch = 19, cex = 1.2,
                        col.smooth = "blue", span = 2/3, iter = 3, ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    abline(lm(y~x), lwd=2.5, col = col.smooth, ...)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# bind Ensembl ID to results and name the columns (imported function (Michael Love)). Online link is needed for reference
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


