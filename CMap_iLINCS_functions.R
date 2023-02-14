#' Data analysis workflow functions for the xls/txt documents that are retrieved from CMap/iLINCS/SwissTargetPrediction platforms
#'
#' @param ptrn Defaults to _conn_HA1E_, just change it the way you kept your data, preferably in seperate folder to avoid confusion
#' @param types Defaults to kd and oe
#' @keywords mylib
#' connectivitydatafile_compiler()

##connectivitydatafile_compiler: Compiles CMap connectivity data.Currently, it works locally and utilizes names of the downloaded files to parse through. Main output file contains gene connectivity data from knock-down/over-expression studies. 
##In addition, connectivity score thresholds can be further specified and heatmaps/correlation plots can be produced, if called. Users can also query for compound-compound connectivity and compound-compound class connectivity data. 
##Future formats can include a logical variable whether the data should be accessed locally or from CMap directly.
##Accordingly, curl commands can be posted to and retrieved from CMap Clue.

connectivitydatafile_compiler=function(ptrn=NULL, types=NULL,cormethod='spearman', cmaprank_cutoff=98, cor_outputs=F, heatmap_outputs=F, width_heat=18, height_heat=12, res_heat=300){
  fullds=c('kd','oe','cp','cc')
  #Confirm input parameters and update when needed
  if(is.null(ptrn)){pattern='_conn_HA1E_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt='Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}
  
  resp='y'
  if (is.null(types)) {typewise=c('kd','oe')
  }  else {typewise=unique(tolower(types))
  if (sum(tolower(typewise) %in% fullds)!=length(typewise)) {typewise=fullds[which(fullds%in%typewise)]
  print('Warning: One or more types arguement is missing. Calculation will be continued with the correct arguements available.')
  }}
  while (sum(tolower(typewise) %in% fullds)==0) { 
    typewise=tolower(as.character(readline(prompt='Please check out misspelling, if any: knockdown (kd)/overexpression (oe) / compound (cp) / CMap Class (cc)?  ')))
    resp=tolower(as.character(readline(prompt='Wanna add more connectivity types: yes (Y) / no (N)   ')))
    while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more connectivity types: yes (Y) / no (N)   ')))}
    while (resp!='n') {nextype=tolower(as.character(readline(prompt='Name of the next connectivity types: ')))
    while (sum(tolower(nextype) %in% fullds)==0) { nextype=as.character(readline(prompt='Please check out misspelling, if any: knockdown (kd)/overexpression (oe) / compound (cp) / CMap Class (cc)?  '))}
    typewise=tolower(unique(append(typewise, nextype)))
    resp=tolower(as.character(readline(prompt='Wanna add more connectivity types: yes (Y) / no (N)   ')))
    while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more connectivity types: yes (Y) / no (N)   ')))}
    }}
  
  #compile files in the cmap directory (data is derived from CMap connections in the .txt format)
  files=list.files(pattern=pattern)
  for (i in files) {if (!exists("dataset")){
    dataset=read.csv(i,header=T,sep="\t")[,c("Name", "Type","Score")]
    colnames(dataset)[colnames(dataset)=="Score"]=sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
    dataset$Name=gsub(x=dataset$Name, pattern="[[:punct:]]", replacement="_")
    dataset$Name_type=paste0(dataset$Name,"_",dataset$Type);
    dataset$Name=NULL; dataset$Type=NULL
    dataset=dataset %>% group_by(Name_type) %>% summarise_all(funs(mean))
  }
    if (exists("dataset")){
      temp_dataset=read.csv(i,header=T,sep="\t")[,c("Name", "Type","Score")]
      colnames(temp_dataset)[colnames(temp_dataset)=="Score"]=sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
      temp_dataset$Name=gsub(x=temp_dataset$Name, pattern="[[:punct:]]", replacement="_")
      temp_dataset$Name_type=paste0(temp_dataset$Name,"_",temp_dataset$Type);
      temp_dataset$Name=NULL; temp_dataset$Type=NULL
      temp_dataset=temp_dataset %>% group_by(Name_type) %>% summarise_all(funs(mean))
      if (sum(colnames(temp_dataset)%in%colnames(dataset))!=ncol(temp_dataset)) {
        dataset=merge(dataset, temp_dataset, by="Name_type",all=T)}
      rm(temp_dataset)
    }}
  dataset=dataset[!duplicated(dataset),]
  rownames(dataset)=dataset$Name_type
  #replace special characters with underscore to avoid further issues
  colnames(dataset)=gsub(x=colnames(dataset), pattern="[[:punct:]]", replacement="_")
  dataset[is.na(dataset)]=0
  
  cpcc=which(typewise%in%c('cp', 'cc'))
  kdoe=which(typewise%in%c('kd','oe'))
  
  #Iterate through connectivity types (kd, oe, cp, cc) and produce correlation plots/heatmaps for each separately
  for (ds in typewise) {
    require(strex)
    tmp=dataset[grepl(x=rownames(dataset),paste0("^.+_", ds,"$")),]
    tmp$Name_type=NULL
    tmp=as.matrix.data.frame(tmp)
    if (ds %in% typewise[cpcc]) {
      rownames(tmp)=strex::str_before_last(pattern = paste0('_',ds),string = rownames(tmp)) 
    }
    plotlab=''
    if (ds=='cp') {plotlab='ConnectivityMap matching compounds'
    } else if (ds=='cc') {plotlab='ConnectivityMap matching compound classes'
    } else if (ds=='kd') {plotlab='ConnectivityMap gene knock-down matches'
    } else if (ds=='oe') {plotlab='ConnectivityMap gene over-expression matches'
    }
    
    #corplots. Correlation plots take full list of connections into account when comparing different samples
    cormat_go=signif(cor(tmp,method=cormethod),2)
    proco=c("thrombin",'phylloquinone', 'menadione','desmopressin')
    antico_thr=c("heparin",'Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
    proco_thr=c("thrombin",'phylloquinone', 'menadione','desmopressin')
    bloods=c('blood', 'thromb', 'agula', 'erythroc', 'leukoc', 'arteri','vein ', 'vascular')
    
    corrplots=list()
    heatmaps=list()
    
    cols=rep('black', ncol(tmp))
    cols[grepl(paste(proco, collapse='|'), colnames(tmp),ignore.case=T)]='red4'
    cols[grepl(paste(antico_thr, collapse='|'), colnames(tmp),ignore.case=T)]='dimgray'
    cols[grepl(paste(proco_thr, collapse='|'), colnames(tmp),ignore.case=T)]='red'
    col_paths=rep('black',nrow(tmp))
    col_paths[grepl(paste(bloods, collapse='|'), rownames(tmp),ignore.case=T)]='blue'
    if (cor_outputs) {
      par(cex.main=1)
      lgd_list=Legend(labels=c("Pro-coagulants", "Anti-coagulants", 'iLINCS: Pro-coagulants','iLINCS: Anti-coagulants'), title="Compound types", type="points", pch=18, legend_gp=gpar(col=c("red4","black",'red', 'dimgray'), fontsize=24))
      f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
      jpeg(paste0(gsub('-','',Sys.Date()),'_',toupper(ds),'_','threshold',cmaprank_cutoff,'_',toupper(cormethod),'corplot',".jpeg"), units="in", width=18, height=12, res=300)
      ht1=Heatmap(cormat_go,clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f1, name=paste0(plotlab,"\n |CMap threshold| > ", cmaprank_cutoff), row_names_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold", col=cols),column_names_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold", col=cols),
                  heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(8, "cm"),labels_gp=gpar(fontsize=nrow(cormat_go)*2),title_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold")),
                  cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
                    grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=nrow(cormat_go)*2, fontface='bold'))})
      tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - ",plotlab,'\n',cormethod," correlations \n |Connectivity Scores| > ", cmaprank_cutoff,", n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
      k=paste0('Corplot','_',ds)
      corrplots[[k]]=tmpfig
      dev.off()}
    
    #heatmaps. Heatmaps evaluate the connections that are above an absolute threshold in at least one of the samples
    if (heatmap_outputs) {
      tmp2=as.data.frame(tmp)
      tmp2=as.data.frame.matrix(tmp[rowSums(abs(tmp)>cmaprank_cutoff)>=1,])
      col_paths2=rep('black',nrow(tmp2))
      col_paths2[grepl(paste(bloods, collapse='|'), rownames(tmp2),ignore.case=T)]='blue'
      f2=colorRamp2(seq(max(tmp2), min(tmp2), length=3),  c("#3f7fff", "#EEEEEE", "red"))
      jpeg(paste0(gsub('-','',Sys.Date()),'_',toupper(ds),'_','threshold',cmaprank_cutoff,'_','Heatmap','.jpeg'), units="in", width=width_heat, height=height_heat, res=res_heat)
      
      ht2=Heatmap(as.matrix.data.frame(tmp2),clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f2, name=paste0(plotlab,'\n',cormethod," correlations \n |Connectivity Scores| > ", cmaprank_cutoff,", n=", nrow(tmp2),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp2)*0.05,2), fontface="bold", col=col_paths), show_row_names = T,column_names_gp=gpar(fontsize=12, fontface="bold",col=cols),
                  heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"), legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold"), title_position="topcenter"),use_raster=F)
      tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
      k=paste0('Heatmap','_',ds)
      heatmaps[[k]]=tmpheat
      dev.off()
      # rm(ind); rm(pathways); rm(tmp)
    }
    rm(tmp)}
  
  #Retrieve the data for kd/oe experiments only to process further with gene enrichment analyses
  dataset=dataset[grepl(x=rownames(dataset),paste0("^.+_(kd|oe)$")),]
  return(dataset)
  rm( list=Filter( exists, c('dataset','typewise', 'pattern','ds') ) )
}


##iLINCS_kd_connectivitydatafile_compiler: Similar to connectivitydatafile_compiler function, and compiles data from iLINCS instead.
##Further developments by online queries are needed. So that, connectivity information would be collected for specified compounds.
iLINCS_kd_connectivitydatafile_compiler=function(ptrn=NULL){
  if(is.null(ptrn)){pattern='_conn_HA1E_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt='Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}
  files=list.files(pattern=pattern)
  for (i in files) {  if (!exists("dataset")){
    dataset=read.table(i,header=T,sep="\t")[,c("GeneTarget", "Correlation","zScore")]
    dataset[dataset=='+']=1; dataset[dataset=='-']=-1
    dataset$zScore=dataset$zScore*as.numeric(dataset$Correlation)
    colnames(dataset)[colnames(dataset)=="zScore"]=sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
    dataset$Correlation=NULL
  }
    if (exists("dataset")){
      temp_dataset=read.table(i,header=T,sep="\t")[,c("GeneTarget", "Correlation","zScore")]
      temp_dataset[temp_dataset=='+']=1; temp_dataset[temp_dataset=='-']=-1
      temp_dataset$zScore=temp_dataset$zScore*as.numeric(temp_dataset$Correlation)
      colnames(temp_dataset)[colnames(temp_dataset)=="zScore"]=sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
      temp_dataset$Correlation=NULL
      if (sum(colnames(temp_dataset)%in%colnames(dataset))!=ncol(temp_dataset)) {
        dataset=merge(dataset, temp_dataset, by="GeneTarget",all=T)
      }
      rm(temp_dataset)
    }}

  dataset=dataset[!duplicated(dataset),]
  dataset[is.na(dataset)]=0
  rownames(dataset)=dataset$Name
  #replace special characters with underscore to avoid further issues
  colnames(dataset)=gsub(x=colnames(dataset), pattern="[[:punct:]]", replacement="_")
  return(dataset)
  rm( list=Filter(exists, c('dataset','pattern')))
}

##iLINCS_kd_connectivitydatafile_compiler: Similar to connectivitydatafile_compiler function, and compiles data from iLINCS instead.
##Further developments by online queries are needed. So that, connectivity information would be collected for specified compounds.
iLINCS_DEGdata_compiler=function(ptrn=NULL){
  if(is.null(ptrn)){pattern='iLINCS_complete_HA1E_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt='Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}
  files=list.files(pattern=pattern)
  for (i in files) {  if (!exists("dataset")){
    dataset=read.table(i,header=T,sep="\t")[,c("Name_GeneSymbol", "Value_LogDiffExp","Significance_pvalue")]
    colnames(dataset)[colnames(dataset)=="Value_LogDiffExp"]=sub(".xls", "_logDE",sub(paste0(".*\\",pattern),"",i));
    colnames(dataset)[colnames(dataset)=="Significance_pvalue"]=sub(".xls", "_pval",sub(paste0(".*\\",pattern),"",i));
    colnames(dataset)[colnames(dataset)=='Name_GeneSymbol']='Genes'
  }
    if (exists("dataset")){
      temp_dataset=read.table(i,header=T,sep="\t")[,c("Name_GeneSymbol", "Value_LogDiffExp","Significance_pvalue")]
      colnames(temp_dataset)[colnames(temp_dataset)=="Value_LogDiffExp"]=sub(".xls", "_logDE",sub(paste0(".*\\",pattern),"",i));
      colnames(temp_dataset)[colnames(temp_dataset)=="Significance_pvalue"]=sub(".xls", "_pval",sub(paste0(".*\\",pattern),"",i));
      colnames(temp_dataset)[colnames(temp_dataset)=='Name_GeneSymbol']='Genes'
      if (sum(colnames(temp_dataset)%in%colnames(dataset))!=ncol(temp_dataset)) {
        dataset=merge(dataset, temp_dataset, by="Genes",all=T)
      }
      rm(temp_dataset)
    }}
  
  rownames(dataset)=dataset$Genes; dataset$Genes=NULL
  #replace special characters with underscore to avoid further issues
  colnames(dataset)=gsub(x=colnames(dataset), pattern="[[:punct:]]", replacement="_")
  return(dataset)
  rm(list=Filter(exists, c('dataset', 'pattern')))
}

##swisstargetdatafile_compiler: Compiles SwissTargetPrediction derived data (.xlsx). Currently, it works locally and utilizes names of the downloaded files to parse through. Each input file has a first row which is omitted by the startRow arguement. 
##If the input file is manually modified, you may need to change this arguement. Main output file contains gene connectivity data from knock-down/over-expression studies. 
##In addition, targetscore_cutoff argument inherits a threshold for filtering data where each target and for any compound should meet a certain minimum value (it is 0.5 by default here)
##If you have a group of compounds to focus and see how the data from rest aligns with them, you can change the focused argument to TRUE, specify a set of compounds (by focuscandies arguement), and further check on the focuscandiesThreshold value to change the resolution based on the focused compounds
##Future formats can include a logical variable whether the data should be accessed locally or from SwissTargetPrediction directly.
##Although slightly slow for a long range of compounds, SwissSimilarity tool from the same platform has a curl and bash based accession pipeline. It would be worthwhile to check

swisstargetdatafile_compiler=function(ptrn=NULL,startRow=2, targetscore_cutoff=0.5, focused=F, focuscandies=NULL, focuscandiesThreshold=0.2,heatmap_outputs=F, clustering_method='ward.D',width_heat=18, height_heat=18, res_heat=300){
  require(openxlsx)
  require(ComplexHeatmap)
  require(dplyr)
  require(RColorBrewer); require(circlize)
  #Confirm input parameters and update when needed
  if(is.null(ptrn)){pattern='SwissTargetPrediction_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt='Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}
  
  dataset=data.frame();rm(dataset)
  #compile files in the swisstarget directory (data is derived from SwissTargetPrediction in the .xlsx format)
  files=list.files(pattern=ptrn)
  
  for (i in files) {  
    
    if (!exists("dataset")){
      dataset=read.xlsx(i,startRow = startRow)[,c("Target", "Probability*")]
      colnames(dataset)[colnames(dataset)=="Probability*"] = sub(".xlsx", "",sub(pattern,"",i));
      dataset_info=read.xlsx(i,startRow = startRow)[,c('Target',"Common.name")]
    }
    if (exists("dataset")){
      temp_dataset=read.xlsx(i,startRow = startRow)[,c("Target", "Probability*")]
      colnames(temp_dataset)[colnames(temp_dataset)=="Probability*"] = sub(".xlsx", "",sub(pattern,"",i));
      if (sum(colnames(temp_dataset)%in%colnames(dataset))!=ncol(temp_dataset)) {
        dataset=merge(dataset, temp_dataset, by="Target",all=T)
      }
      temp_dataset_info=read.xlsx(i,startRow = startRow)[,c('Target',"Common.name")]
      dataset_info=rbind(dataset_info,temp_dataset_info)
      rm(temp_dataset);rm(temp_dataset_info)
    }
    
  }
  
  dataset_info=dataset_info[!duplicated(dataset_info),]
  dataset=dataset[,!duplicated(colnames(dataset))]
  dataset[is.na(dataset)]=0
  dataset=dataset[!duplicated(dataset),]
  rownames(dataset)=dataset$Target
  dataset$Target=NULL
  ##get the targets that have at least some considerable score for at least one of the compounds
  dataset=dataset[rowSums(dataset>targetscore_cutoff)>=1,]
  dataset=as.matrix.data.frame(dataset)
  colnames(dataset)=gsub('_iso','',colnames(dataset))
  colnames(dataset)=tolower(colnames(dataset))
  
  if (focused) {
    if (is.null(focuscandies)) {focuscandies=colnames(dataset)
    }  else {focuscandies=unique(tolower(focuscandies))
    if (sum(tolower(focuscandies) %in% colnames(dataset))!=length(focuscandies)) {focuscandies=colnames(dataset)[which(colnames(dataset)%in%focuscandies)]
    print('Warning: One or more focuscandies arguement is missing. Calculation will be continued with the correct arguements available.')
    }}
    if (sum(tolower(focuscandies) %in% colnames(dataset))==0) {
      print('Warning: Focuscandies arguement does not match with the query compounds. Calculations will be continued for the whole data')
      focuscandies=colnames(dataset)}
    
    dataset=as.data.frame(dataset)
    dataset[rowSums(dataset[which(colnames(dataset)%in%focuscandies)])>focuscandiesThreshold,]}
  
  dataset=as.data.frame(dataset)
  dataset$Target=rownames(dataset)
  dataset=merge(dataset, dataset_info, by='Target', all.x=T)
  
  #remove homologous genes to avoid matching row names in the next steps
  dataset=dataset[!grepl('by homology', dataset$Target),]
  dataset$Target=NULL
  rownames(dataset)=dataset$Common.name; dataset$Common.name=NULL
  dataset=as.matrix.data.frame(dataset)
  
  if (heatmap_outputs) {
    heatmaps=list()
    cols <- rep('black', ncol(dataset))
    #turn red the specified rows in tf
    proco=c("thrombin",'phylloquinone', 'menadione','desmopressin')
    cols[colnames(dataset)%in%proco]= 'red4'
    
    lgd_list=Legend(labels = c("Pro-coagulants", "Anti-coagulants"), title = "Compound types", type = "points", pch = 18, legend_gp = gpar(col = c("red4","black"), fontsize=24))
    f1 = colorRamp2(seq(max(dataset), min(dataset), length=5), c("dark red","red3","red2","red","white"))
    
    jpeg(paste0(gsub('-','',Sys.Date()),'_SwissTarget_ProAntiCoagulants_threshold',targetscore_cutoff,'_Heatmap_',clustering_method,'.jpeg'), units="in", width=width_heat, height=height_heat, res=res_heat)
    
    ht=Heatmap(dataset,clustering_method_rows = clustering_method,clustering_method_columns  = clustering_method, col = f1, name = "SwissTarget similarity scores", row_names_gp = gpar(fontsize = max(10-nrow(dataset)*0.05,2), fontface = "bold"),column_names_gp = gpar(fontsize = 12, fontface = "bold",col=cols),
               heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height = unit(0.8, "cm"),legend_width = unit(10, "cm"),labels_gp = gpar(fontsize = 15),title_gp = gpar(fontsize = 15, fontface = "bold"),title_position = "topcenter"), border = T)
    tmpheat=draw(ht, heatmap_legend_side = "top",  column_title_gp=gpar(fontface="bold", line = 5),annotation_legend_list = lgd_list,annotation_legend_side = "bottom")
    heatmaps[['Heatmap']]=tmpheat
    dev.off()
  }
  
  #Retrieve the data for kd/oe experiments only to process further with gene enrichment analyses
  return(dataset)
  rm( list=Filter( exists, c('dataset','temp_dataset','focuscandies', 'pattern') ) )
}

#" enrichr_output: Retrieves EnrichR results. It requires a vector of drug list, p-value cut off, name of the database used (kegg="KEGG_2019_Human","GO_Biological_Process_2018" (see EnrichR sources)),
##and a specific direction (for the genes either up or down regulated, if inquired). If the direction is not specified, the calculations are done for all the DEGs (up and downregulated ones separately).
##The output file contains the lists for the EnrichR results in the specified directions separately. Each list contains three subcategories in data frame format:
##(1) List of pathways and associated combination scores, (2) Expanded list that contains the full range of EnrichR results, and (3) List of pathways and associated genes coming from the input gene expression data.

enrichr_output=function(drug_list,dataset,logDEcutoff=0,pcutoff,dbs,direction=NULL, cmapranks=F, cmaprank_cutoff){
  require(dplyr);require(purrr);require(enrichR)
  result=list()
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case=T))==0) { ds=tolower(as.character(readline(prompt='Please check out misspelling, if any: complete/up/down?')))}
  }
  for (d in ds) {
    if (exists('go_output')) {rm(go_output)}
    if (exists('termgenes')) {rm(termgenes)}
    k=list()
    if (!cmapranks) {
      for (drugs in drug_list) {
        if (!exists("go_output")){
          go_input=dataset[,grepl(drugs, sub('_logDE|_pval','',colnames(dataset)))]
          go_input=go_input[order(go_input[,2]),]
          if (d=='complete') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&abs(go_input[,1])>abs(logDEcutoff),]), dbs)}
          else if (d=='up') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]>abs(logDEcutoff),]), dbs)}
          else {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]<(-abs(logDEcutoff)),]), dbs)}
          k[[drugs]]=enriched[[dbs[[1]]]]
          go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          # rownames(go_output)=go_output$Term; go_output$Term=NULL
          colnames(go_output)=c('Term',drugs)
          colnames(termgenes)=c('Term',drugs)
          rm( list=Filter( exists, c('go_input','enriched') ) )}
        if (exists("go_output")){
          go_input=dataset[,grepl(drugs, sub('_logDE|_pval','',colnames(dataset)))]
          go_input=go_input[order(go_input[,2]),]
          if (d=='complete') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&abs(go_input[,1])>abs(logDEcutoff),]), dbs)
          } else if (d=='up') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]>abs(logDEcutoff),]), dbs)
          } else {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]<(-abs(logDEcutoff)),]), dbs)}
          k[[drugs]]=enriched[[dbs[[1]]]]
          enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          # rownames(enriched)=enriched$Term; enriched$Term=NULL
          colnames(enriched_go)=c('Term',drugs)
          colnames(enriched_go_termgenes)=c('Term',drugs)
          go_output=merge(go_output, enriched_go, by='Term',all=T)
          termgenes=merge(termgenes, enriched_go_termgenes, by='Term',all=T)
          rm( list=Filter( exists, c('go_input','enriched_go') ) )}}
    } else if(cmapranks) {
      for (drugs in drug_list) {
        if (!exists("go_output")){
          go_input=dataset[,c(1,grep(drugs,colnames(dataset)))]
          if (d=='complete') {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks=go_input[abs(go_input[,2])>cmaprank_cutoff,]
          upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_oe',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_kd',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_kd',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_oe',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          # rownames(go_output)=go_output$Term; go_output$Term=NULL
          colnames(go_output)=c('Term',drugs)
          colnames(termgenes)=c('Term',drugs)
          rm( list=Filter( exists, c('go_input','enriched') ) )}
        if (exists("go_output")){
          go_input=dataset[,c(1,grep(drugs,colnames(dataset)))]
          if (d=='complete') {upranks=go_input[abs(go_input[,2])>cmaprank_cutoff,]
          upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_oe',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_kd',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>cmaprank_cutoff,]
          upranks1=upranks1[grepl('_kd',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-cmaprank_cutoff),]
          upranks2=upranks2[grepl('_oe',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          # rownames(enriched)=enriched$Term; enriched$Term=NULL
          colnames(enriched_go)=c('Term',drugs)
          colnames(enriched_go_termgenes)=c('Term',drugs)
          go_output=merge(go_output, enriched_go, by='Term',all=T)
          termgenes=merge(termgenes, enriched_go_termgenes, by='Term',all=T)
          rm( list=Filter( exists, c('go_input','enriched_go') ) )}}
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

# enrichrlist=keggtr; direction="all"; combscore_cutoff=200;refcompounds=c("Mestranol","Warfarin");d='complete'

#slicenricher: CMap ranks are default to FALSE, unless users provide a rank list of kd/oe studies from CMap connectivity data. The functions, slices and collects specified data (based on Combined score cutoff) from the enrichr_output function. Accordingly, it retrieves pathway enrichment scores from up or downregulated genes ('complete' option accounts both up and downregulated genes together). Based on the enrichment scores, users can visually inspect pathway enrichment profiles across the treatments either by heatmaps or corplots
##Currently, it only utilizes the EnrichedTermsScores objects from the enrichr_output data.
##Specifications rely on the combscore cutoff values for each reference compound. List of the drugs can be fetched into data as a drug_list, or the function prompts users to specify the names of the compounds to be compared.
##Combination scores are presented on log2 scale. Accordingly, minimum combination score value is capped to 1 (hence log2(1)=0 across the representations)
##Main outputfile contains a list for the pathway names passing the threshold and their combination scores for the respective compounds inquired by the user.
##Lists are given in separate sections for each direction (up/down/complete)). In addition, the correlation plots (by default 'spearman') and
##heatmaps are produced and saved in the current directory, unless indicated otherwise.
##Heatmap label aesthetics are also improved to indicate whether the drugs of interest are coagulant in the literature and/or in our in-house phenotypic screenings.
##Pathways onthologies that matches with thrombosis keyword patterns are also manually included. A warning function can be included to change the cutoff value, in case there are not many enriched pathways, hence errors based on non-finite calculations
##Heatmap files' width, height and resolution parameters can be user specified, as well. For the future updates, the fucntion can be incorporated into a dashboard or a shinyapp to give better control over the aesthetics and dynamic features
slicenricher=function(enrichrlist, direction=NULL, combscore_cutoff=200, refcompounds=NULL, cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL='EnrichR', width_heat=18, height_heat=12, res_heat=300){
  require(tidyverse); require(ggplot2); require(ggord);require(ggfortify);
  require(ComplexHeatmap);require(circlize);
  tryCatch({if (!cmapranks) {drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][1])))))
  } else if (cmapranks) {drugnames=gsub('.*_','',names(enrichrlist[[1]][[1]]))}
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case=T))==0) {
    ds=tolower(as.character(readline(prompt='Please check out misspelling, if any: complete/up/down?')))}
  }
  resp='y'
  while (is.null(refcompounds)) {refcompounds=as.character(readline(prompt='Please provide a valid compound name from the list: '))
  while (sum(grepl(paste(refcompounds,collapse="|"), drugnames,ignore.case=T))==0) {
    refcompounds=as.character(readline(prompt='Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt='Name of the next compound: ')))
  while (sum(grepl(paste(nextcomp,collapse="|"), drugnames,ignore.case=T))==0) { nextcomp=as.character(readline(prompt='Please recheck the compound name, and provide it once again: '))}
  refcompounds=append(refcompounds, nextcomp)
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}}
  }

  if (is.null(enrichmentdb)) {enrichmentdb=as.character(readline(prompt='Name of the gene enrichment database (KEGG/GO/REACT etc.): '))
  while (length(enrichmentdb)==0) { enrichmentdb=as.character(readline(prompt='Please provide the name of the gene enrichment database (KEGG/GO/REACT etc) that you have used: '))}
  }
  proco=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin')
  antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
  proco_thr=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin')
  bloods=c('blood', 'thromb', 'agula', 'erythroc', 'leukoc', 'arteri','vein ', 'vascular')

  result=list()
  corrplots=list()
  heatmaps=list()
  for(d in ds){tmp=as.data.frame(enrichrlist[names(enrichrlist)%in%d][[1]][1])
  colnames(tmp)=gsub(paste0(names(enrichrlist[names(enrichrlist)%in%d][[1]][1]),'.'),'',colnames(tmp))
  rownames(tmp)=tmp$Term;
  ind=which(grepl(paste(refcompounds,collapse="|"), colnames(tmp),ignore.case=T))
  pathways=vector()
  for (i in ind) {pathways=unique(append(pathways,rownames(tmp[which(tmp[,i]>combscore_cutoff),])))}
  tmp=tmp[rownames(tmp)%in%pathways,ind]
  tmp[tmp<=1]=1
  tmp=log2(tmp) #Log scale the combination scores
  result[[d]]=tmp
  tmp$Term=NULL
  tmp=tmp[,colSums(tmp)>0]
  if (length(tmp)>0) {

    #corplot
    cormat_go=signif(cor(tmp,method=cormethod),2)
    cols=rep('black', ncol(tmp))
    cols[grepl(paste(proco, collapse='|'), colnames(tmp),ignore.case=T)]='red4'
      cols[grepl(paste(antico_thr, collapse='|'), colnames(tmp),ignore.case=T)]='dimgray'
        cols[grepl(paste(proco_thr, collapse='|'), colnames(tmp),ignore.case=T)]='red'
          col_paths=rep('black',nrow(tmp))
          col_paths[grepl(paste(bloods, collapse='|'), rownames(tmp),ignore.case=T)]='blue'
            if (cor_outputs) {
              par(cex.main=1)
              lgd_list=Legend(labels=c("Pro-coagulants", "Anti-coagulants", 'iLINCS: Pro-coagulants','iLINCS: Anti-coagulants'), title="Compound types", type="points", pch=18, legend_gp=gpar(col=c("red4","black",'red', 'dimgray'), fontsize=24))
              f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
              jpeg(paste0(deparse(substitute(enrichrlist)),'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substring(paste(refcompounds,collapse='-'),1,15),toupper(cormethod),'corplot_',toupper(d),".jpeg"), units="in", width=18, height=12, res=300)
              ht1=Heatmap(cormat_go,clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f1, name="Pathway combination score correlations (iLINCS - EnrichR)", row_names_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold", col=cols),column_names_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold", col=cols),
                          heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(8, "cm"),labels_gp=gpar(fontsize=nrow(cormat_go)*2),title_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold")),
                          cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
                            grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=nrow(cormat_go)*2, fontface='bold'))})
              tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - ",d,'\n',cormethod," correlations on iLINCS pathways \n (EnrichR ",paste(refcompounds,collapse='-'), "\n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
              k=paste0('Corplot','_',d)
              corrplots[[k]]=tmpfig
              dev.off()}

          #heatmaps
          if (heatmap_outputs) {
            f2=colorRamp2(seq(max(tmp), min(tmp), length=3),  c("#3f7fff", "#EEEEEE", "red"))
            jpeg(paste0(deparse(substitute(enrichrlist)),'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substring(paste(refcompounds,collapse='-'),1,15),"_Heatmap_",toupper(d),".jpeg"), units="in", width=width_heat, height=height_heat, res=res_heat)

            ht2=Heatmap(as.matrix.data.frame(tmp),clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f2, name=paste0(toupper(d), " - Pathway combinations scores \n","(EnrichR ",paste(refcompounds,collapse='-'), "\n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp2)*0.05,2), fontface="bold", col=col_paths),column_names_gp=gpar(fontsize=12, fontface="bold",col=cols),
                        heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"),legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold"),title_position="topcenter"),use_raster=F)
            tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
            k=paste0('Heatmap','_',d)
            heatmaps[[k]]=tmpheat
            dev.off()
            rm(ind); rm(pathways); rm(tmp)
          }}}
  rm(drugnames)
  return(result)},error=function(e){print('Warning: There was an error. Please, check if the parameters are correct. Combination score cut-off might be too stringent for some of the compounds. Please, check the outputs or you can also try lower combination score thresholds')})
}

#slicenricher_updown: Similar to slicenricher. It allows comparing different groups of drugs (pro-coagulants vs anti-coagulants) to see whether there are any pathways that move in same/opposite directions.
slicenricher_updown=function(enrichrlist, combscore_cutoff=200, contr1=NULL,contr2=NULL,classContr1=NULL,classContr2=NULL,cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL='EnrichR', repositoryName=NULL, noContrastwarning=T){
  require(tidyverse); require(ggplot2); require(ggord);require(ggfortify);
  require(ComplexHeatmap);require(circlize);
  tryCatch({if (!cmapranks) {drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][1])))))
  } else if (cmapranks) {drugnames=gsub('.*_','',names(enrichrlist[[1]][[1]]))}
  if (is.null(enrichmentdb)) {enrichmentdb=as.character(readline(prompt='Name of the gene enrichment database (KEGG/GO/REACT): '))
  while (length(enrichmentdb)==0) { enrichmentdb=as.character(readline(prompt='Please provide the name of the gene enrichment database that you have used (KEGG/GO/REACT etc): '))}
  }
  resp='y'
  while (is.null(contr1)) {contr1=as.character(readline(prompt='Please provide a valid compound name from the list to contrast for: '))
  while (sum(grepl(paste(contr1,collapse="|"), drugnames,ignore.case=T))==0) {
    contr1=as.character(readline(prompt='Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt='Wanna add more compounds for the same direction: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds for the same direction: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt='Name of the next compound: ')))
  while (sum(grepl(paste(nextcomp,collapse="|"), drugnames,ignore.case=T))==0) { nextcomp=as.character(readline(prompt='Please recheck the compound name, and provide it once again: '))}
  contr1=append(contr1, nextcomp)
  resp=tolower(as.character(readline(prompt='Wanna add more compounds from the same direction: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds from the same direction: yes (Y) / no (N)   ')))}}
  }
  while (is.null(contr2)) {contr2=as.character(readline(prompt='Please provide a valid compound name from the list to contrast against the first compound(s): '))
  while (sum(grepl(paste(contr2,collapse="|"), drugnames,ignore.case=T))==0) {
    contr2=as.character(readline(prompt='Please recheck the compound names, and provide the name of the second compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt='Name of the next compound: ')))
  while (sum(grepl(paste(nextcomp,collapse="|"), drugnames,ignore.case=T))==0) { nextcomp=as.character(readline(prompt='Please recheck the compound name, and provide it once again: '))}
  contr2=append(contr2, nextcomp)
  resp=tolower(as.character(readline(prompt='Wanna add more compounds from the same direction: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds from the same direction: yes (Y) / no (N)   ')))}}
  }
  if (is.null(classContr1)) {classContr1=as.character(readline(prompt='What is the class of the first contrasts (Pro-coagulant / Anti-coagulant etc): '))
  while (length(classContr1)==0) { classContr1=as.character(readline(prompt='What is the class of the first contrasts (Pro-coagulant / Anti-coagulant etc): '))}}
  if (is.null(classContr2)) {classContr2=as.character(readline(prompt='What is the class of the second contrasts (Pro-coagulant / Anti-coagulant etc): '))
  while (length(classContr2)==0) { classContr2=as.character(readline(prompt='What is the class of the second contrasts (Pro-coagulant / Anti-coagulant etc): '))}}
  if (is.null(repositoryName)) {repositoryName=as.character(readline(prompt='What is the source of the data (CMap, iLINCS etc): '))
  while (length(repositoryName)==0) { repositoryName=as.character(readline(prompt='What is the source of the data (CMap, iLINCS etc): '))}}

  proco=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin')
  antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
  proco_thr=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin')
  bloods=c('blood', 'thromb', 'agula', 'erythroc', 'leukoc', 'arteri','vein ', 'vascular')

  ds=list(c('up','down'),c('down','up'))
  result=list()
  corrplots=list()
  heatmaps=list()
  for(d in 1:length(ds)){tmp1=as.data.frame(enrichrlist[names(enrichrlist)%in%ds[[d]][1]][[1]][1])
  colnames(tmp1)=paste0(gsub(paste0(names(enrichrlist[names(enrichrlist)%in%ds[[d]][1]][[1]][1]),'.'),'',colnames(tmp1)),'_',ds[[d]][1])
  colnames(tmp1)[1]='Term'
  rownames(tmp1)=tmp1$Term;
  tmp1=tmp1[,grepl(paste(c('Term',contr1),collapse='|'),colnames(tmp1),ignore.case = T)]
  tmp2=as.data.frame(enrichrlist[names(enrichrlist)%in%ds[[d]][2]][[1]][1])
  colnames(tmp2)=paste0(gsub(paste0(names(enrichrlist[names(enrichrlist)%in%ds[[d]][2]][[1]][1]),'.'),'',colnames(tmp2)),'_',ds[[d]][2])
  colnames(tmp2)[1]='Term'
  rownames(tmp2)=tmp2$Term;
  tmp2=tmp2[,grepl(paste(c('Term',contr2),collapse='|'),colnames(tmp2),ignore.case = T)]
  tmp=merge(tmp1,tmp2,by='Term',all=T)
  rownames(tmp)=tmp$Term
  #prevent any possible issues with generated NAs from the merges; and capping the minimum value to 1 (hence log2(1)=0 in the follow up representations)
  tmp[is.na(tmp)]=1
  tmp[tmp<1]=1
  ind=which(grepl(paste(c(contr1,contr2),collapse="|"), colnames(tmp),ignore.case=T))
  pathways=vector()
  for (i in ind) {pathways=unique(append(pathways,rownames(tmp[which(tmp[,i]>combscore_cutoff),])))}
  tmp=tmp[rownames(tmp)%in%pathways,]
  tmp[,-1]=log2(tmp[,-1]) #Log scale the combination scores
  # tmp[tmp<0]=0
  result[[paste(c(ds[[d]][1],ds[[d]][2]),collapse='_')]]=tmp
  tmp$Term=NULL
  tmp=tmp[,colSums(tmp)>0]
  if (length(tmp)>0) {

    #corplot
    cormat_go=signif(cor(tmp,method=cormethod),2)
    cols=rep('black', ncol(tmp))
    cols[grepl(pattern=paste0(proco, collapse='|'), x=colnames(tmp),ignore.case=T)]='red4'
      cols[grepl(pattern=paste0(antico_thr, collapse='|'), x=colnames(tmp),ignore.case=T)]='dimgray'
        cols[grepl(pattern=paste0(proco_thr, collapse='|'), x=colnames(tmp),ignore.case=T)]='red'
          col_paths=rep('black',nrow(tmp))
          col_paths[grepl(pattern=paste0(bloods, collapse='|'), x=rownames(tmp),ignore.case=T)]='blue'
            if (cor_outputs) {
              par(cex.main=1)
              lgd_list=Legend(labels=c("Pro-coagulants", "Anti-coagulants", 'iLINCS: Pro-coagulants','iLINCS: Anti-coagulants'), title="Compound types", type="points", pch=18, legend_gp=gpar(col=c("red4","black",'red', 'dimgray'), fontsize=24))
              f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
              jpeg(paste0(deparse(substitute(enrichrlist)),'_',paste(c(repositoryName,enrichmentdb,enrichmentTOOL),collapse='-'),'_',substr(paste(c(contr1,contr2),collapse='-'),1,30),'_',toupper(cormethod),'corplot_',toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='_')),".jpeg"), units="in", width=18, height=12, res=300)
              ht1=Heatmap(cormat_go,clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f1, name=paste0("Pathway combination score correlations (",repositoryName," - ", enrichmentTOOL,")"), row_names_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold", col=cols),column_names_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold", col=cols),
                          heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(nrow(cormat_go), "cm"),labels_gp=gpar(fontsize=nrow(cormat_go)*2),title_gp=gpar(fontsize=nrow(cormat_go)*2, fontface="bold")),
                          cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
                            grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=nrow(cormat_go)*2, fontface='bold'))})
              tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - ",toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='vs')),'\n',cormethod," correlations on ",repositoryName,'-',enrichmentdb," pathways \n (", enrichmentTOOL," ",substr(paste(c(contr1,contr2),collapse='-'),1,30), "\n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
              k=paste0('Corplot','_',d)
              corrplots[[k]]=tmpfig
              dev.off()}
          #heatmaps
          if (heatmap_outputs) {
            f2=colorRamp2(seq(max(tmp), min(tmp), length=3),  c("#3f7fff", "#EEEEEE", "red"))
            jpeg(paste0(deparse(substitute(enrichrlist)),'_',repositoryName,'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substr(paste(c(contr1,contr2),collapse='-'),1,30),'__Heatmap_',toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='_')),".jpeg"), units="in", width=18, height=12, res=300)

            ht2=Heatmap(as.matrix.data.frame(tmp),clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f2, name=paste0(toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='vs')), " - Pathway combinations scores \n (",paste(c(repositoryName,enrichmentdb,enrichmentTOOL),collapse='-'),"\n",substr(paste(c(contr1,contr2),collapse='-'),1,30), "... \n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp2)*0.05,2), fontface="bold", col=col_paths),column_names_gp=gpar(fontsize=12, fontface="bold",col=cols), heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"),legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=10, fontface="bold"),title_position="topcenter"),use_raster=F)
            tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
            k=paste0('Heatmap','_',d)
            heatmaps[[k]]=tmpheat
            dev.off()
            rm(ind); rm(pathways); rm(tmp)
          }}}
  rm(drugnames)
  return(result)
}, error=function(e){print('Warning: There was an error. Please, check if the parameters are correct. Otherwise, some contrasts may not have been processed (due to lack of overlaps between the query compounds). Please, check the outputs or you can also try lower combination score thresholds')})}



##upsetlister, an implemented function to provide the input for the upset plots and venn interaction terms. So that, from many pathways in the list, users can see how many of them are shared among the reference compounds and across different directions. If the lists do not seem complicated (meaning that they are not enriched with multiple pathways and directions) upsetlist and further plots will not be needed.Slicenricher outputs are preferred and transformed into upset plot convenient list forms while collecting pathways that are above a certain threshold. Reference compounds and regulation direction (up, down or complete) can be further specified.
upsetlister=function(slicenrichrlist,refcompounds=NULL,direction=NULL, combscore_cutoff=200){
  drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(slicenrichrlist[[1]])))))
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case=T))==0) { ds=tolower(as.character(readline(prompt='Please check out misspelling, if any: complete/up/down?')))}}
  resp='y'
  if (is.null(refcompounds)) {refcompounds=as.character(readline(prompt='Please provide a valid compound name from the list: '))
  while (sum(grepl(paste(refcompounds,collapse="|"), drugnames,ignore.case=T))==0) { refcompounds=as.character(readline(prompt='Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt='Name of the next compound: ')))
  while (sum(grepl(paste(nextcomp,collapse="|"), drugnames,ignore.case=T))==0) { nextcomp=as.character(readline(prompt='Please recheck the compound name, and provide it once again: '))}
  refcompounds=append(refcompounds, nextcomp)
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}}
  }
  upsetlist=list();
  for(d in ds){tmp=as.data.frame(slicenrichrlist[names(slicenrichrlist)%in%d][[1]])
  ind=which(grepl(paste(refcompounds,collapse="|"), colnames(tmp),ignore.case=T))
  colnames(tmp)[ind]=paste0(colnames(tmp[,ind]),'_',d)
  for (i in ind) {comp=list(rownames(tmp[tmp[,i]>log2(combscore_cutoff),]))
  names(comp)=colnames(tmp)[i]
  upsetlist=append(upsetlist, comp)
  }
  }
  return(upsetlist)}


#paths_list: Retrieves intersections of the pathway names across the user defined queries and upsetlists. Further implemented in the upsetlistlegend function
paths_list=function(myqueries, allqueries,intersectionlist){
  require(devtools); require(Vennerable)
  unlist(intersectionlist[paste(as.numeric(allqueries %in% myqueries),collapse="")])}


#termgenesmajor, retrieves the names of the differential genes relating with the previous enrichr pathway outputs.
##Hence, it allows seeing overexpression or knockdown of which genes (from CMap/iLINCS) can be further important.
##The function is also nested inside the upsetlistlegend function
termgenesmajor=function(enrichrlist,direction=NULL, refcompounds=NULL){
  drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][3])))))
  if (is.null(direction)) {ds=c('complete','up','down')}
  else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case=T))==0) { ds=tolower(as.character(readline(prompt='Please check out misspelling, if any: complete/up/down?')))}
  }
  resp='y'
  if (is.null(refcompounds)) {refcompounds=as.character(readline(prompt='Please provide a valid compound name from the list: '))
  while (sum(grepl(paste(refcompounds,collapse="|"), drugnames,ignore.case=T))==0) { refcompounds=as.character(readline(prompt='Please recheck the compound names, and provide the name of the first compound \n Do not worry about capital letters, time/conc of exposures: '))}
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}
  while (resp!='n') {nextcomp=tolower(as.character(readline(prompt='Name of the next compound:')))
  while (sum(grepl(paste(nextcomp,collapse="|"), drugnames,ignore.case=T))==0) { nextcomp=as.character(readline(prompt='Please recheck the compound name, and provide it once again: '))}
  refcompounds=append(refcompounds, nextcomp)
  resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))
  while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt='Wanna add more compounds: yes (Y) / no (N)   ')))}}
  }
  majorlist=setNames(data.frame(matrix(ncol=1, nrow=0)), c("Term"))
  for(d in ds){tmp=enrichrlist[[d]][[3]]
  ind=which(grepl(paste(refcompounds,collapse="|"), colnames(tmp),ignore.case=T))
  colnames(tmp)[ind]=paste0(colnames(tmp[,ind]),'_',d)
  majorlist=merge(majorlist, as.data.frame(tmp[,c(1,ind)]), by='Term', all=T)
  }
  return(majorlist)
}

#querylister allows conversion of queries to compatible formats to be processed along with upsetlists and upset plots. Requires intersectionSets from Venn function (see vennerable). Numbers of queries for each comparisons can be adjusted with the lowerlimit and upperlimit options
querylister=function(intersectionSets,upsetlist, lowerlimit=0, upperlimit=max(lengths(intersectionSets))){qs=list()
ques=intersectionSets[lengths(intersectionSets)>lowerlimit & lengths(intersectionSets)<=max(lengths(intersectionSets))]
quenames=names(ques[order(lengths(ques),decreasing=T)])
for (ints in quenames) {
  ints_membs=strsplit(ints, '')[[1]]
  tmplist=list()
  for (j in 1:length(ints_membs)) {
    if(as.numeric(ints_membs[j])){tmplist=append(tmplist, names(upsetlist[j])) }}
  qs[[ints]]=tmplist}
rm(ques);rm(quenames);names(qs)=NULL;
return(qs)}


#upsetlistlegend creates the upsetplots with appropriate legend annotations that are defined by the users' querylists.
##In addition, pathways of interests and their respective genes from CMap/iLINCS studies are extracted in excel sheets to the working directory. Maximum total numbers of pathway terms in the querylists (npathways arguement) are defaulted to 50 to prevent overcrowding the upsetplot figure while the full range of pathways between the queries are still accessible within the output excel sheet. 
upsetlistlegend=function(upsetlist, querylist, enrichrlist,outputlist=NULL,colorpal=NULL, ylabel='No. of Pathways', cell_line='NA',combscore_cutoff=200,enrichmentdb=NULL,refcompounds=NULL,upsetFigure=T,enrichmentTOOL,repositoryName='iLINCS', npathways=50){
  require(stringr)
  require(Vennerable)
  require(RColorBrewer)
  require(UpSetR)
  require(openxlsx)
  if (is.null(cell_line)) {cell_line=as.character(readline(prompt='What was the name of the cell line used in the study (please double-check, this is for representation purposes only): '))
  while (length(cell_line)==0) { cell_line=as.character(readline(prompt='What was the name of the cell line used in the study (please double-check, this is for representation purposes only): '))}}
  if (is.null(combscore_cutoff)) {combscore_cutoff=as.character(readline(prompt='What was the ncombination score threshold that you applied (please double-check, this is for representation purposes only): '))
  while (length(combscore_cutoff)==0) { combscore_cutoff=as.character(readline(prompt='What was the ncombination score threshold that you applied (please double-check, this is for representation purposes only): '))}}

  temp=Venn(upsetlist)  #provide all your groups as list

  names(upsetlist)=factor(names(upsetlist),levels=sort(names(upsetlist)))

  # Define the number of colors you want
  if (is.null(colorpal)) {
    colorpal='Dark2'
  }
  nb.cols=min(npathways, length(querylist))
  mycolors=colorRampPalette(brewer.pal(8, colorpal))(nb.cols)
  pathtables=termgenesmajor(enrichrlist, direction=NULL,refcompounds= refcompounds)
  gridtable=setNames(data.frame(matrix(ncol=ncol(pathtables), nrow=0)),colnames(pathtables))
  queries=list()

  for (i in 1:length(querylist)) {
    queries[[i]]=list(query=intersects,params=querylist[[i]], color=mycolors[i], active=T)}

 
  if (upsetFigure) {

    jpeg(paste0(gsub('-','',Sys.Date()),"_",as.character(repositoryName),'_',as.character(enrichmentTOOL),cell_line,'_',paste(substr(str_to_title(refcompounds),1,4),collapse=''),"_CombScore",combscore_cutoff,"_",enrichmentdb,"_UpsetPlot.jpeg"), units="in", width=14, height=10, res=300)

    upsetpl=upset(fromList(upsetlist), nsets=sum(lengths(temp@IntersectionSets)>0), queries=queries, order.by="freq", mainbar.y.label=ylabel)
    print(upsetpl)
    grid.text(paste0(enrichmentdb,' Upset plot - Pathways \n (log2(Combination Score) >', combscore_cutoff,')'),x=0.5, y=0.94, gp=gpar(fontsize=10, fontface='bold'),just="centre")}
  
  cnts=0
  
  for (i in 1:length(querylist)) {
    plist=paths_list(myqueries=unlist(querylist[[i]]), allqueries=names(upsetlist), temp@IntersectionSets)
    tmptable=pathtables[pathtables$Term %in% unlist(plist),]
    gridtable=rbind(gridtable,tmptable)
    if(nrow(gridtable)<=npathways & upsetFigure){
    for (j in 1:length(plist)) {
      grid.text(plist[[j]],x=0.78, y=(0.96-cnts*0.011), gp=gpar(fontsize=6.5, col=mycolors[i], fontface='bold'),just="left")
      cnts=cnts+1}}}
  if (upsetFigure) {grid.polygon(x=c(0.75,1,1,0.75), y=c((0.95-cnts*0.011),(0.95-cnts*0.011),0.97,0.97), gp=gpar(fill="gray", alpha=0.2))
    dev.off()}
  write.xlsx(x=gridtable,file=paste0(gsub('-','',Sys.Date()),"_",as.character(repositoryName),'_',as.character(enrichmentTOOL),'_',cell_line,'_',paste(substr(str_to_title(refcompounds),1,4),collapse=''),"_CombScore",combscore_cutoff,"_",enrichmentdb,"_UpsetPlotLegend.xlsx"), rowNames=F)
  return(gridtable)
}

#Adopted function (reference needed)
freqfunc2=function(x, n=30){
  tail(sort(table( unlist(strsplit(c(na.omit(as.character(unlist(x, use.names=FALSE)))), ";",fixed=T)))), n)
}

#Adopted function (reference needed)
strsplit2=function(x,
                      split,
                      type="remove",
                      perl=FALSE,
                      ...) {
  if (type=="remove") {
    # use base::strsplit
    out=base::strsplit(x=x, split=split, perl=perl, ...)
  } else if (type=="before") {
    # split before the delimiter and keep it
    out=base::strsplit(x=x,
                          split=paste0("(?<=.)(?=", split, ")"),
                          perl=TRUE,
                          ...)
  } else if (type=="after") {
    # split after the delimiter and keep it
    out=base::strsplit(x=x,
                          split=paste0("(?<=", split, ")"),
                          perl=TRUE,
                          ...)
  } else {
    # wrong type input
    stop("type must be remove, after or before!")
  }
  return(out)
}

###pairs with p-val (imported function (Michael Love)). Online link is needed for reference
library(RColorBrewer)
cols_BR=brewer.pal(11, "RdBu")   # goes from red to white to blue
pal=colorRampPalette(cols_BR)
cor_colors=data.frame(correlation=seq(-1,1,0.01),
                        correlation_color=pal(201)[1:201])  # assigns a color for each r correlation value
cor_colors$correlation_color=as.character(cor_colors$correlation_color)

panel.cor=function(x, y, digits=2, cex.cor)
{
  par(usr=c(0, 1, 0, 1))
  u=par('usr')
  names(u)=c("xleft", "xright", "ybottom", "ytop")
  r=cor(x, y,method="pearson",use="complete.obs")
  test=cor.test(x,y)
  bgcolor=cor_colors[2+(-r+1)*100,2]    # converts correlation into a specific color
  do.call(rect, c(col=bgcolor, as.list(u))) # colors the correlation box

  if (test$p.value> 0.05){
    text(0.5,0.5,"Insignificant",cex=1.5)
  } else{
    text(0.5, 0.75, paste("r=",round(r,2)),cex=2.5) # prints correlatoin coefficient
    text(.5, .25, paste("p=",formatC(test$p.value, format="e", digits=1)),cex=2)
    abline(h=0.5, lty=2) # draws a line between correlatoin coefficient and p value
  }

}
panel.smooth=function (x, y, col="black", bg=NA, pch=19, cex=1.2,
                        col.smooth="blue", span=2/3, iter=3, ...) {
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
  ok=is.finite(x) & is.finite(y)
  if (any(ok))
    abline(lm(y~x), lwd=2.5, col=col.smooth, ...)
}
panel.hist=function(x, ...)
{
  usr=par("usr"); on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5) )
  h=hist(x, plot=FALSE)
  breaks=h$breaks; nB=length(breaks)
  y=h$counts; y=y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# bind Ensembl ID to results and name the columns (imported function (Michael Love)). Online link is needed for reference
convertIDs=function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple=match.arg( ifMultiple )
  suppressWarnings( selRes=AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple=="putNA" ) {
    duplicatedIds=selRes[ duplicated( selRes[,1] ), 1 ]
    selRes=selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}



## Modified plotPCA from DESeq2 package. Shows the Names of the Samples (the first col of SampleTable), and uses ggrepel pkg to plot them conveniently.
# @SA 10.02.2017
library(genefilter)
library(ggplot2)
library(ggrepel)

plotPCA.san=function (object, intgroup="condition", ntop=500, returnData=FALSE, PCno1=1, PCno2=2)
{
  rv=rowVars(assay(object))
  select=order(rv, decreasing=TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca=prcomp(t(assay(object)[select, ]))
  percentVar=pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df=as.data.frame(colData(object)[, intgroup, drop=FALSE])
  group=if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse=" : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d=data.frame(PC1=pca$x[, PCno1], PC2=pca$x[, PCno2], group=group,
                  intgroup.df, name=colData(rld)[,1])
  if (returnData) {
    attr(d, "percentVar")=percentVar
    return(d)
  }
  ggplot(data=d, aes_string(x=paste0("PC", as.character(PCno1)), y=paste0("PC", as.character(PCno2)), color="group", label="name")) + geom_point(size=3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3)

}
