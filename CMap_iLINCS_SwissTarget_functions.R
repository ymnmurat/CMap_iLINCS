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
  require(ComplexHeatmap)
  require(circlize); require(RColorBrewer)
  require(dplyr)
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
    dataset=dataset %>% group_by(Name_type) %>% summarise_all(list(~mean(., na.rm=T)))
    dataset=as.data.frame(dataset)
  }
    if (exists("dataset")){
      temp_dataset=read.csv(i,header=T,sep="\t")[,c("Name", "Type","Score")]
      colnames(temp_dataset)[colnames(temp_dataset)=="Score"]=sub(".[^.]*$", "",sub(paste0(".*\\",pattern),"",i));
      temp_dataset$Name=gsub(x=temp_dataset$Name, pattern="[[:punct:]]", replacement="_")
      temp_dataset$Name_type=paste0(temp_dataset$Name,"_",temp_dataset$Type);
      temp_dataset$Name=NULL; temp_dataset$Type=NULL
      temp_dataset=temp_dataset %>% group_by(Name_type) %>% summarise_all(list(~mean(., na.rm=T)))
      temp_dataset=as.data.frame(temp_dataset)
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
  
  #Iterate through connectivity types (kd, oe, cp, cc) and produce correlation plots/heatmaps for each of them separately
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
    } else if (ds=='kd') {plotlab='ConnectivityMap gene knock-down matches \n (gene-based)'
    } else if (ds=='oe') {plotlab='ConnectivityMap gene over-expression matches \n (gene-based)'
    }
    
    #some color coding is also included. So that, labels for each compound will show whether they are pro-coagulant (red), anti-coagulant (dimgray) or black if they belong to none. 
    
    proco=c("thrombin",'phylloquinone', 'menadione','desmopressin', 'U46619')
    antico_thr=c("heparin",'Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
    cols=rep('black', ncol(tmp))
    cols[grepl(paste(proco, collapse='|'), colnames(tmp),ignore.case=T)]='red4'
      cols[grepl(paste(antico_thr, collapse='|'), colnames(tmp),ignore.case=T)]='dimgray'
        
      #Plot legened corresponding to color coding of the compounds
      lgd_list=Legend(labels=c("Pro-coagulants", "Anti-coagulants", "Other compounds"), title="Compound types", type="points", pch=18, legend_gp=gpar(col=c("red4",'dimgray',"black"), fontsize=24))
      
      #corplots. Correlation plots take full list of connections into account when comparing different samples.
      cormat_go=signif(cor(tmp,method=cormethod),2)
      corrplots=list()
      
      if (cor_outputs) {
        par(cex.main=1)
        
        f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
        jpeg(paste0(gsub('-','',Sys.Date()),'_',toupper(ds),'_CMapConnectivityScores_',toupper(cormethod),'corplot',".jpeg"), units="in", width=width_heat, height=height_heat, res=res_heat)
        ht1=Heatmap(cormat_go,clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f1, name=paste0(plotlab,"\n |CMap threshold| > ", cmaprank_cutoff), row_names_gp=gpar(fontsize=max(10-nrow(cormat_go)*0.05,2), fontface="bold", col=cols),column_names_gp=gpar(fontsize=max(12-ncol(cormat_go)*0.05,2), fontface="bold", col=cols), heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(8, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold")),cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
       grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=max(12-nrow(cormat_go)*0.05,2), fontface='bold'))})
        #Draw the figure and store it in a list
        tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - ",plotlab,'\n',cormethod," correlations \n |Connectivity Scores| > ", cmaprank_cutoff,", n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
        k=paste0('Corplot','_',ds)
        corrplots[[k]]=tmpfig
        dev.off()}
      
      #heatmaps. Heatmaps evaluate the connections that are above an absolute threshold in at least one of the samples. Blue color is also assigned for pathways that contain coagulation related terms
      heatmaps=list()
      
      if (heatmap_outputs) {
        tmp2=as.data.frame(tmp)
        tmp2=as.data.frame.matrix(tmp[rowSums(abs(tmp)>cmaprank_cutoff)>=1,])

        f2=colorRamp2(seq(max(tmp2), min(tmp2), length=3),  c("#3f7fff", "#EEEEEE", "red"))
        
        jpeg(paste0(gsub('-','',Sys.Date()),'_',toupper(ds),'_CMapConnectivityScores_threshold',cmaprank_cutoff,'_','Heatmap','.jpeg'), units="in", width=width_heat, height=height_heat, res=res_heat)
        
        ht2=Heatmap(as.matrix.data.frame(tmp2),clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f2, name=paste0(plotlab,"\n |Connectivity Scores| > ", cmaprank_cutoff,", n=", nrow(tmp2),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp2)*0.05,2), fontface="bold", col='black'), show_row_names = T,column_names_gp=gpar(fontsize=max(12-ncol(tmp2)*0.05,2), fontface="bold",col=cols), heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"), legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold"), title_position="topcenter"),use_raster=F)
        #Draw the figure and store it in a list
        tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
        k=paste0('Heatmap','_',ds)
        heatmaps[[k]]=tmpheat
        dev.off()
      }
      rm(tmp)}
  
  #Retrieve the data for kd/oe experiments only to process further i.e. gene enrichment analyses
  dataset=dataset[grepl(x=rownames(dataset),paste0("^.+_(kd|oe)$")),]
  return(dataset)
  rm( list=Filter( exists, c('dataset','typewise', 'pattern','ds','ptrn') ) )
}

##iLINCS_kd_connectivitydatafile_compiler: Similar to connectivitydatafile_compiler function, and compiles data from iLINCS instead.
#Main output file will contain gene knockdown based connectivity data and z-scores. You can also check out additional arguments of this function below, e.g. to produce heatmaps/correlation plots for the connectivity data based on compound-gene knock-down z-scores
##ptrn arguement is by default '_conn_HA1E_'. You can apply similar naming convention for your own files. 
##Further developments by online queries are needed. So that, connectivity information would be collected for user-selected compounds.
iLINCS_kd_connectivitydatafile_compiler=function(ptrn=NULL, zscores_cutoff=2, clusterMethodRow="ward.D", clusterMethodColumn="ward.D",cormethod='spearman', cor_outputs=F, heatmap_outputs=F, width_heat=18, height_heat=12, res_heat=300){
  if(is.null(ptrn)){pattern='_conn_HA1E_'; ptrn=pattern
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
  
  #replace special characters with underscore to avoid further issues
  colnames(dataset)=gsub(x=colnames(dataset), pattern="[[:punct:]]", replacement="_")
  tmp=dataset; rownames(tmp)=tmp$GeneTarget; tmp$GeneTarget=NULL; tmp=as.data.frame.matrix(tmp)
  
  #some color coding is also included. So that, labels for each compound will show whether they are pro-coagulant (red), anti-coagulant (dimgray) or black if they belong to none. 
  proco=c("thrombin",'phylloquinone', 'menadione','desmopressin', 'U46619')
  antico_thr=c("heparin",'Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
  cols=rep('black', ncol(tmp))
  cols[grepl(paste(proco, collapse='|'), colnames(tmp),ignore.case=T)]='red4'
    cols[grepl(paste(antico_thr, collapse='|'), colnames(tmp),ignore.case=T)]='dimgray'
      #   
    # col_paths=rep('black',nrow(tmp))
    # col_paths[grepl(paste(bloods, collapse='|'), rownames(tmp),ignore.case=T)]='blue'
    
    #Plot legened corresponding to color coding of the compounds
    lgd_list=Legend(labels=c("Pro-coagulants", "Anti-coagulants", "Other compounds"), title="Compound types", type="points", pch=18, legend_gp=gpar(col=c("red4",'dimgray',"black"), fontsize=24))
    
    
    #corplots. Correlation plots take full list of connections into account when comparing different samples
    corrplots=list()
    cormat_go=signif(cor(tmp,method=cormethod),2)
    # proco_thr=c("thrombin",'phylloquinone', 'menadione','desmopressin','U46619')
    # bloods=c('blood', 'thromb', 'agula', 'erythroc', 'leukoc', 'arteri','vein ', 'vascular')
    if (cor_outputs) {
      par(cex.main=1)
      
      f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
      jpeg(paste0(gsub('-','',Sys.Date()),'_iLINCSkd_GeneConnectivityScores_',toupper(cormethod),'corplot',".jpeg"), units="in", width=18, height=12, res=300)
      ht1=Heatmap(cormat_go,clustering_method_rows=clusterMethodRow,clustering_method_columns =clusterMethodColumn, col=f1, name=paste0("iLINCS Gene Connectivity \n |zScores| > ", zscores_cutoff), row_names_gp=gpar(fontsize=max(10-nrow(cormat_go)*0.05,2), fontface="bold", col=cols),column_names_gp=gpar(fontsize=max(12-ncol(cormat_go)*0.05,2), fontface="bold", col=cols), heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(8, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold")),cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
        grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=max(12-nrow(cormat_go)*0.05,2), fontface='bold'))})
      
      #Draw the figure and store it in a list
      tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - iLINCS Gene connectivity\n",cormethod," correlations \n |zScores| > ", zscores_cutoff,", n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
      k=paste0('Corplot')
      corrplots[[k]]=tmpfig
      dev.off()}
    
    #heatmaps. Heatmaps evaluate the connections that are above an absolute threshold in at least one of the samples
    if (heatmap_outputs) {
      heatmaps=list()
      tmp2=as.data.frame(tmp)
      tmp2=as.data.frame.matrix(tmp[rowSums(abs(tmp)>zscores_cutoff)>=1,])
      # col_paths2=rep('black',nrow(tmp2))
      # col_paths2[grepl(paste(bloods, collapse='|'), rownames(tmp2),ignore.case=T)]='blue'
      f2=colorRamp2(seq(max(tmp2), min(tmp2), length=3),  c("#3f7fff", "#EEEEEE", "red"))
      jpeg(paste0(gsub('-','',Sys.Date()),'_iLINCSkd_GeneConnectivityScores_thresholdAbs',zscores_cutoff,'_','Heatmap','.jpeg'), units="in", width=width_heat, height=height_heat, res=res_heat)
      
      ht2=Heatmap(as.matrix.data.frame(tmp2),clustering_method_rows=clusterMethodRow,clustering_method_columns =clusterMethodColumn, col=f2, name=paste0("iLINCS Gene Connectivity \n |zScores| > ", zscores_cutoff,", n=", nrow(tmp2),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp2)*0.05,2), fontface="bold", col='black'), show_row_names = T,column_names_gp=gpar(fontsize=max(12-ncol(tmp2)*0.05,2), fontface="bold",col=cols),heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"), legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold"), title_position="topcenter"),use_raster=F)
      
      #Draw the figure and store it in a list
      tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
      k=paste0('Heatmap')
      heatmaps[[k]]=tmpheat
      dev.off()
    }
    #Retrieve the data for kd experiments to process further i.e. gene enrichment analyses
    return(dataset)
    rm( list=Filter(exists, c('dataset','pattern','ptrn','tmp','tmp2')))
}


##iLINCS_DEGdata_compiler: Allows compiling differential gene expression data from iLINCS complete lists (.xls). Returned file contains differential expression values (log scaled) and pvalues which can be further utilized to do follow-up processes i.e. gene enrichment analyses

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
  rm(list=Filter(exists, c('dataset', 'pattern','ptrn')))
}


##swisstargetdatafile_compiler: Compiles SwissTargetPrediction derived data (.xlsx). Currently, it works locally and utilizes names of the downloaded files to parse through. Each input file has a first row which is omitted by the startRow argument. 
##If the input file is manually modified, you may need to change this argument. Main output file contains SwissTargetPrediction combined scores(2D/3D). Please, visit http://www.swisstargetprediction.ch/faq.php#Why_is_it_useful_to_combine_different_similarity_measures for details.
##In addition, targetscore_cutoff argument inherits a threshold for filtering data where each possible target for any of the compounds should meet a certain minimum value (it is 0.5 by default here)
##If you have a group of compounds to focus and see how the rest of the data aligns with them, you can change the focused argument to TRUE, specify a set of compounds (by focuscandies argument), and further check on the focuscandiesThreshold value to change the resolution based on the focused compounds
##Future formats can include a logical variable whether the data should be accessed locally or from SwissTargetPrediction directly.
##Although slightly slow for a long range of compounds, SwissSimilarity tool from the same platform has a curl and bash based accession pipeline. It would be worthwhile to check

swisstargetdatafile_compiler=function(ptrn=NULL,startRow=2, clippedSwissData=F, targetscore_cutoff=0.5, focused=F, focuscandies=NULL, focuscandiesThreshold=0.2,heatmap_outputs=F, clustering_method='ward.D',width_heat=18, height_heat=18, res_heat=300){
  require(openxlsx)
  require(dplyr)
  #Confirm input parameters and update when needed
  if(is.null(ptrn)){pattern='SwissTargetPrediction_'
  }else {pattern=ptrn
  while (length(list.files(pattern=pattern))==0) {pattern=as.character(readline(prompt='Please check out the current directory, and  if any, misspelling, capital letters etc.  '))}}
  if (clippedSwissData) {startRow=1}
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
  
  ## Some targets contain multiple common names delimited by a space. We can produce new rows by splitting them.
  dataset_info=dataset_info %>% tidyr::separate_rows(Common.name,sep = ' ')
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
  
  #If interested you can additionally include a filter for some reference compounds of your choice, and see how the rest of the data aligns with their targets
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
  dataset=dataset %>% group_by(Common.name) %>% summarise_all(list(~mean(., na.rm=T)))
  dataset=as.data.frame(dataset)
  rownames(dataset)=dataset$Common.name; dataset$Common.name=NULL
  dataset=as.matrix.data.frame(dataset)
  
  
  if (heatmap_outputs) {
    require(ComplexHeatmap)
    require(RColorBrewer); require(circlize)
    heatmaps=list()
    
    #Color code pro-coagulant (red) and anti-coagulant (gray) for the representations
    proco=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin','U46619')
    antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
    cols=rep('black', ncol(dataset))
    cols[grepl(paste(proco, collapse='|'), colnames(dataset),ignore.case=T)]='red4'
      cols[grepl(paste(antico_thr, collapse='|'), colnames(dataset),ignore.case=T)]='dimgray'
        # cols[grepl(paste(proco_thr, collapse='|'), colnames(tmp),ignore.case=T)]='red'
      
      lgd_list=Legend(labels=c("Pro-coagulants", "Anti-coagulants", "Other compounds"), title="Compound types", type="points", pch=18, legend_gp=gpar(col=c("red4",'dimgray',"black"), fontsize=24))
      
      #Make a color ramp of reds where low scores have fader tones of red and white in the lowest end
      f1 = colorRamp2(seq(max(dataset), min(dataset), length=5), c("dark red","red3","red2","red","white"))
      
      jpeg(paste0(gsub('-','',Sys.Date()),'_SwissTarget_ProAntiCoagulants_threshold',targetscore_cutoff,'_Heatmap_',clustering_method,'.jpeg'), units="in", width=width_heat, height=height_heat, res=res_heat)
      
      ht=Heatmap(dataset,clustering_method_rows = clustering_method,clustering_method_columns  = clustering_method, col = f1, name = "SwissTarget similarity scores", row_names_gp = gpar(fontsize = max(12-nrow(dataset)*0.05,2), fontface = "bold"),column_names_gp = gpar(fontsize = max(15-ncol(dataset)*0.05,2), fontface = "bold",col=cols),
                 heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height = unit(0.8, "cm"),legend_width = unit(10, "cm"),labels_gp = gpar(fontsize =15),title_gp = gpar(fontsize = 15, fontface = "bold"),title_position = "topcenter"), border = T)
      tmpheat=draw(ht, heatmap_legend_side = "top",  column_title_gp=gpar(fontface="bold", line = 5),annotation_legend_list = lgd_list,annotation_legend_side = "bottom")
      heatmaps[['Heatmap']]=tmpheat
      dev.off()
  }
  
  #Retrieve the data for kd/oe experiments only to process further with gene enrichment analyses
  return(dataset)
  rm( list=Filter( exists, c('dataset','temp_dataset','focuscandies', 'pattern','ptrn') ) )
}

##multi_heatmaps_SwissCMapILINCS compiles SwissTargetPrediction, iLINCS and CMap derived data, and produces a combined heatmap figure from them and from a query list of genes. Optionally, users can provide a reference data (RNAseq based logFC, gene-phenotype correlation scores etc) as an annotation for the figure. Please, use similar file name conventions and directory allocations (see folders; swisstarget, ilincs, cmap under github/cmap_ilincs directory). #CellLine is by default HA1E, you can additionally include secondCellLine argument. Optionally, you can provide a refdata to annotate your heatmaps. Final heatmaps are saved inside the main  github/cmap_ilincs directory automatically. Here you can also change the name for your single column reference data annotation on the heatmap, colors of the NA data (default 'gray') and width/height/resolution of your heatmap figure
multi_heatmaps_SwissCMapILINCS=function(swisstargetDir='swisstarget', iLINCSconnectivityDir='../ilincs/', cmapConnectDir='../cmap/',singleAnnotationTag='Measurements', annotationLegendTitle='user input logFC data',refdata=NULL, candygenes, swissFocus=T, CellLine='HA1E', secondCellLine=NULL, nacols=NULL,width_heat=20, height_heat=14, res_heat=600){
  
  #set and change current directory to main github/cmap_ilincs folder 
  currentdir=getwd()
  # setwd(gsub(x=currentdir,pattern = '(github_stuff).*','\\1'))
  setwd(gsub(x=currentdir,pattern = '(github/cmap_ilincs).*','\\1'))
  #SwissTargetPrediction collections
  require(ComplexHeatmap)
  require(dplyr)
  require(circlize); require(RColorBrewer)
  #Capitalize letters of each gene for consistency and exact matching  
  candygenes=toupper(candygenes)
  candygenes=candygenes[order(candygenes)]
  if (length(candygenes)>=30) {print('Warning: Length of the candidate genes is too long. You may need to refine your query for finer images')}
  
  #If provided, modify the reference data based on the list of query genes (candygenes). Add aditional rows for the missing genes, and arrange the data in a proper shape. so that, they can be annotated on the heatmaps. Otherwise, it will be removed from the analysis pipeline
  if(!is.null(refdata)&is.data.frame(refdata)){
    miss_genes_ref=candygenes[!candygenes%in%toupper(rownames(refdata))]
    miss_mat_ref=matrix(data=as.numeric(NA),nrow = length(miss_genes_ref),ncol = ncol(refdata))
    rownames(miss_mat_ref)=miss_genes_ref; colnames(miss_mat_ref)=colnames(refdata)
    refdata=rbind(refdata,miss_mat_ref)
    if (ncol(refdata)==1) {
      refdata=t(refdata); refdata=refdata[,toupper(colnames(refdata))%in%candygenes]
      refdata=as.data.frame(refdata[order(names(refdata))]); colnames(refdata)=singleAnnotationTag
    }else{refdata=refdata[toupper(rownames(refdata))%in%candygenes,]
    refdata=refdata[order(rownames(refdata)),]}
  } else if(!is.null(refdata)&!is.data.frame(refdata)){refdata=NULL
  print('Warning: refdata is not a data frame and removed from the analysis')}
  
  #Compile SwissTargetPrediction data for the canidate genes. If numbers of compounds are less than two, the function will not proceed. Since the main heatmap is from SwissTargePrediction data, users will also have the option to make the observations based on the available compounds from SwissTarget and other platforms. If the overlap of compound is limited, or the users are not willing to focus on SwissTarget solely, this approach will be cancelled and the comparisons will proceed with no bias for SwissTargetPrediction compounds.
  setwd(swisstargetDir)
  stp=swisstargetdatafile_compiler(ptrn = 'SwissTargetPrediction_',targetscore_cutoff = 0)
  if(ncol(stp)<2){print('Number of samples are limited for the SwissTarget derived data. Please, retrieve additional data from SwissTarget to proceed')
    quit(save = 'no')}
  miss_genes_stp=candygenes[!candygenes%in%toupper(rownames(stp))]
  miss_mat_stp=matrix(data=as.numeric(NA),nrow = length(miss_genes_stp),ncol = ncol(stp))
  rownames(miss_mat_stp)=miss_genes_stp; colnames(miss_mat_stp)=colnames(stp)
  stp=rbind(stp,miss_mat_stp)
  stp=stp[toupper(rownames(stp))%in%candygenes,]
  stp=stp[order(rownames(stp)),]
  
  #Retrieve and compile iLINCS connectivity collections for the query genes.
  setwd(iLINCSconnectivityDir)
  
  ptrn1=paste0('_conn_',CellLine,'_')
  conndata=iLINCS_kd_connectivitydatafile_compiler(ptrn = ptrn1)
  if(ncol(conndata)<3){print(paste0('Number of samples are limited for the iLINCS derived ', CellLine,' data. Please, retrieve additional data from iLINCS to proceed'))
    quit(save = 'no')}else if(ncol(conndata)>9){print('Number of samples might be large for the iLINCS derived ', CellLine,' data. You may need to consider reducing amount of data you are providing')}
  if(swissFocus & sum(grepl(pattern = paste0(toupper(colnames(stp)), collapse = '|'),toupper(colnames(conndata))))>1) {
    conndata=conndata[,grepl(pattern = toupper(paste(c('GeneTarget',colnames(stp)), collapse = '|')), x = toupper(colnames(conndata)),ignore.case = T)]
  }else{swissFocus=F; print(paste0('Numbers of mutual compounds between SwissTarget and iLINCS ', CellLine,' data are limited (<2). Representations will ignore SwissTarget focused refinements further on'))} 
  miss_genes_ilincsConn=rownames(stp)[!toupper(rownames(stp))%in%toupper(conndata$GeneTarget)]
  rownames(conndata)=conndata$GeneTarget; conndata$GeneTarget=NULL
  miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_ilincsConn),ncol = ncol(conndata))
  rownames(miss_mat)=miss_genes_ilincsConn; colnames(miss_mat)=colnames(conndata)
  conndata=rbind(conndata,miss_mat)
  conndata_stp=as.matrix(conndata[toupper(rownames(conndata))%in%toupper(rownames(stp)),])
  conndata_stp=conndata_stp[order(rownames(conndata_stp)),]
  
  #Collect data from the second cell line, if provided
  if (!is.null(secondCellLine)) {
    ptrn2=paste0('_conn_',secondCellLine,'_')
    conndata=iLINCS_kd_connectivitydatafile_compiler(ptrn = ptrn2)
    if(ncol(conndata)<3){print(paste0('Number of samples are limited for the iLINCS derived ', secondCellLine,' data. Please, retrieve additional data from iLINCS to proceed'))
      quit(save = 'no')}else if(ncol(conndata)>9){print('Number of samples might be large for the iLINCS derived ', secondCellLine,' data. You may need to consider reducing amount of data you are providing')}
    if(swissFocus & sum(grepl(pattern = paste0(toupper(colnames(stp)), collapse = '|'),toupper(colnames(conndata))))>1) {
      conndata=conndata[,grepl(pattern = toupper(paste(c('GeneTarget',colnames(stp)), collapse = '|')), x = toupper(colnames(conndata)),ignore.case = T)]
    }else{swissFocus=F; print(paste0('Numbers of mutual compounds between SwissTarget and iLINCS ', secondCellLine,' data are limited (<2). Representations will ignore SwissTarget focused refinements further on'))}
    miss_genes_ilincsConn=rownames(stp)[!toupper(rownames(stp))%in%toupper(conndata$GeneTarget)]
    rownames(conndata)=conndata$GeneTarget; conndata$GeneTarget=NULL
    miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_ilincsConn),ncol = ncol(conndata))
    rownames(miss_mat)=miss_genes_ilincsConn; colnames(miss_mat)=colnames(conndata)
    conndata=rbind(conndata,miss_mat)
    conndata_stp2=as.matrix(conndata[toupper(rownames(conndata))%in%toupper(rownames(stp)),])
    conndata_stp2=conndata_stp[order(rownames(conndata_stp)),]
  }else{conndata_stp2=NULL}
  
  #Collect CMap connectivity data for the query genes
  setwd(cmapConnectDir)
  types=c('kd','oe')
  rm( list = Filter( exists, c('dataset') ) )
  ptrn1=paste0('_conn_',CellLine,'_')
  conndata_cmap=connectivitydatafile_compiler(ptrn = ptrn1,types = types)
  if(ncol(conndata_cmap)<3){print(paste0('Number of samples are limited for the CMap Connectivity derived ', CellLine,' data. Please, retrieve additional data from CMap Connections to proceed'))
    quit(save = 'no')}else if(ncol(conndata_cmap)>9){print('Number of samples might be large for the CMap Connectivity derived ', CellLine,' data. You may need to consider reducing amount of data you are providing')}
  if(swissFocus==T & sum(grepl(pattern = paste0(toupper(colnames(stp)), collapse = '|'),toupper(colnames(conndata_cmap))))>1) {
    conndata_cmap=conndata_cmap[,grepl(pattern = paste(gsub('6_|24_','',colnames(stp)), collapse = '|'),colnames(conndata_cmap), ignore.case = T)]
    conndata_cmap=conndata_cmap[,grepl(pattern = paste(colnames(stp), collapse = '|'), x = colnames(conndata_cmap),ignore.case = T)]
  }else{swissFocus=F; print(paste0('Numbers of mutual compounds between SwissTarget and CMap Connectivity ', CellLine,' data are limited (<2). Representations will ignore SwissTarget focused refinements further on'))} 
  conndata_cmap$Name_type=NULL
  
  #Subset CMap connectivity data for knockdown (kd) genes only 
  conndata_cmap_kd=conndata_cmap[grepl('_kd', rownames(conndata_cmap)),]; rownames(conndata_cmap_kd)=gsub('_kd','',rownames(conndata_cmap_kd))
  miss_genes_cMapConn_kd=candygenes[!candygenes%in%toupper(rownames(conndata_cmap_kd))]
  miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_kd),ncol = ncol(conndata_cmap_kd))
  rownames(miss_mat)=miss_genes_cMapConn_kd; colnames(miss_mat)=colnames(conndata_cmap_kd)
  conndata_cmap_kd=rbind(conndata_cmap_kd,miss_mat)
  conndata_cmap_kd=conndata_cmap_kd[order(rownames(conndata_cmap_kd)),]
  conndata_cmap_kd_stp=as.matrix(conndata_cmap_kd[toupper(rownames(conndata_cmap_kd))%in%candygenes,])
  conndata_cmap_kd_stp=conndata_cmap_kd_stp[order(rownames(conndata_cmap_kd_stp)),]
  rm(miss_mat)
  
  #Subset CMap connectivity data for overexpression (oe) genes only 
  conndata_cmap_oe=conndata_cmap[grepl('_oe', rownames(conndata_cmap)),]; rownames(conndata_cmap_oe)=gsub('_oe','',rownames(conndata_cmap_oe))
  miss_genes_cMapConn_oe=candygenes[!candygenes%in%toupper(rownames(conndata_cmap_oe))]
  miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_oe),ncol = ncol(conndata_cmap_oe))
  rownames(miss_mat)=miss_genes_cMapConn_oe; colnames(miss_mat)=colnames(conndata_cmap_oe)
  conndata_cmap_oe=rbind(conndata_cmap_oe,miss_mat)
  conndata_cmap_oe=conndata_cmap_oe[order(rownames(conndata_cmap_oe)),]
  conndata_cmap_oe_stp=conndata_cmap_oe[toupper(rownames(conndata_cmap_oe))%in%candygenes,]
  conndata_cmap_oe_stp=conndata_cmap_oe_stp[order(rownames(conndata_cmap_oe_stp)),]
  
  #If provided, do the same compilations for the second cell line as well
  if (!is.null(secondCellLine)) {
    rm( list = Filter( exists, c('dataset') ) )
    
    conndata_cmap=connectivitydatafile_compiler(ptrn = ptrn2,types = types)
    if(ncol(conndata_cmap)<3){print(paste0('Number of samples are limited for the CMap Connectivity derived ', secondCellLine,' data. Please, retrieve additional data from CMap Connections to proceed'))
      quit(save = 'no')}else if(ncol(conndata_cmap)>9){print('Number of samples might be large for the CMap Connectivity derived ', secondCellLine,' data. You may need to consider reducing amount of data you are providing')}
    if(swissFocus==T & sum(grepl(pattern = paste0(toupper(colnames(stp)), collapse = '|'),toupper(colnames(conndata_cmap))))>1) {
      conndata_cmap=conndata_cmap[,grepl(pattern = paste(gsub('6_|24_','',colnames(stp)), collapse = '|'),colnames(conndata_cmap), ignore.case = T)]
      conndata_cmap=conndata_cmap[,grepl(pattern = paste(colnames(stp), collapse = '|'), x = colnames(conndata_cmap),ignore.case = T)]
    }else{swissFocus=F; print(paste0('Numbers of mutual compounds between SwissTarget and CMap Connectivity ', secondCellLine,' data are limited (<2). Representations will ignore SwissTarget focused refinements further on'))} 
    conndata_cmap$Name_type=NULL
    
    #CMap kd 
    conndata_cmap_kd=conndata_cmap[grepl('_kd', rownames(conndata_cmap)),]; rownames(conndata_cmap_kd)=gsub('_kd','',rownames(conndata_cmap_kd))
    miss_genes_cMapConn_kd=candygenes[!candygenes%in%toupper(rownames(conndata_cmap_kd))]
    miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_kd),ncol = ncol(conndata_cmap_kd))
    rownames(miss_mat)=miss_genes_cMapConn_kd; colnames(miss_mat)=colnames(conndata_cmap_kd)
    conndata_cmap_kd=rbind(conndata_cmap_kd,miss_mat)
    conndata_cmap_kd=conndata_cmap_kd[order(rownames(conndata_cmap_kd)),]
    conndata_cmap_kd_stp2=as.matrix(conndata_cmap_kd[toupper(rownames(conndata_cmap_kd))%in%candygenes,])
    conndata_cmap_kd_stp2=conndata_cmap_kd_stp[order(rownames(conndata_cmap_kd_stp)),]
    rm(miss_mat)
    #CMap oe
    conndata_cmap_oe=conndata_cmap[grepl('_oe', rownames(conndata_cmap)),]; rownames(conndata_cmap_oe)=gsub('_oe','',rownames(conndata_cmap_oe))
    miss_genes_cMapConn_oe=candygenes[!candygenes%in%toupper(rownames(conndata_cmap_oe))]
    miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_oe),ncol = ncol(conndata_cmap_oe))
    rownames(miss_mat)=miss_genes_cMapConn_oe; colnames(miss_mat)=colnames(conndata_cmap_oe)
    conndata_cmap_oe=rbind(conndata_cmap_oe,miss_mat)
    conndata_cmap_oe=conndata_cmap_oe[order(rownames(conndata_cmap_oe)),]
    conndata_cmap_oe_stp2=conndata_cmap_oe[toupper(rownames(conndata_cmap_oe))%in%candygenes,]
    conndata_cmap_oe_stp2=conndata_cmap_oe_stp[order(rownames(conndata_cmap_oe_stp)),]
  }else{conndata_cmap_kd_stp2=NULL;conndata_cmap_oe_stp2=NULL}
  
  #Confirm the colors for the representations (default: 'gray' for NAs, reds for SwissTargetPrediction, blues for iLINCS and CMap data)
  if(is.null(nacols)){nacols='gray'}
  reds=c("dark red","red3","red2","red","white")
  blues=c("#3f7fff", "#EEEEEE", "red")
  
  ##heatmaptags function: Setting the color scales and adjusting heatmap legends can become redundant and tricky. Especially, if any of the compiled data does not contain viable numbers but all NAs, we can run into some issues. heatmaptags function mainly handles the redundancy for setting the color scales and produces reference lists for us to handle potentially problematic cases. For example, if whole set of the compiled data is NA, it will convert them to 0s, assign a color scale of grays between a minimal range c(0,0.0001) as labels, include a flag tag to update the heatmap legends and store them in a list. Otherwise, it will keep the data as it is, assign a color range of choice ('reds', 'blues') and keep the flag silent (F) to not to update the legend with NAs. 
  heatmaptags=function(dataset,colorscales=NULL,nacols,ranges=T){
    if(is.null(colorscales)){colorscales=c("#3f7fff", "#EEEEEE", "red")}
    details=list()
    if (sum(is.na(dataset))<length(as.matrix(dataset))) {
      if(ranges){f1 = colorRamp2(seq(max(dataset,  na.rm = TRUE)*-1,max(dataset, na.rm = TRUE), length=length(colorscales)), colorscales)
      labels=round(seq(max(dataset,  na.rm = TRUE)*-1,max(dataset, na.rm = TRUE), length=length(colorscales)),2)}else{f1 = colorRamp2(seq(max(dataset,na.rm = T), min(dataset,na.rm = T), length=length(colorscales)), colorscales);labels=rev(round(seq(max(dataset,na.rm = T), min(dataset,na.rm = T), length=length(colorscales)),digits = 2))}
      f1tag=F; at=labels} else {dataset[is.na(dataset)]=0
      f1 = colorRamp2(breaks = c(0,0.0001),c(nacols,nacols)); f1tag=T; at=seq(0,0.0001,length=length(colorscales))
      labels=rep('',length(colorscales)); labels[median(1:length(labels))]='NAs'}
    details[['f1']]=f1; details[['dataset']]=dataset; details[['f1tag']]=f1tag; details[['at']]=at; details[['labels']]=labels;
    return(details)
  }
  
  #Generate tag data for the heatmaps
  stp_tags=heatmaptags(stp, reds, nacols=nacols,ranges = F)
  conndata_stp_tags=heatmaptags(conndata_stp, nacols=nacols)
  
  if(!is.null(refdata)) {refdata_tags=heatmaptags(refdata, nacols=nacols)}
  conndata_cmap_kd_stp_tags=heatmaptags(conndata_cmap_kd_stp, nacols=nacols)
  conndata_cmap_oe_stp_tags=heatmaptags(conndata_cmap_oe_stp, nacols=nacols)
  
  if (!is.null(conndata_stp2)|!is.null(conndata_cmap_kd_stp2)|!is.null(conndata_cmap_oe_stp2)) {
    conndata_stp2_tags=heatmaptags(conndata_stp2, nacols=nacols)
    conndata_cmap_kd_stp2_tags=heatmaptags(conndata_cmap_kd_stp2, nacols=nacols)
    conndata_cmap_oe_stp2_tags=heatmaptags(conndata_cmap_oe_stp2, nacols=nacols)
  }
  
  #assign color codes for the compound names whether they are pro-coagulant (red) or anti-coagulant (dimgray). Assigning colors for multiple data can become redundant. Here is a function that will make the writing process cleaner
  
  colgenaid=function(x, proco=NULL, antico=NULL, proco_col='red4', antico_col='dimgray', base_col='black'){
    if (is.null(proco)) {proco=c("thrombin",'phylloquinone', 'menadione','desmopressin', 'U46619')}
    if (is.null(antico)) {antico=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')}
    cols=rep(base_col, ncol(x));
    cols[grepl(paste(toupper(proco), collapse = '|'),toupper(colnames(x)),ignore.case = T)]= proco_col;
    cols[grepl(paste(toupper(antico), collapse = '|'),toupper(colnames(x)),ignore.case = T)]= antico_col;
    return(cols)
  }
  
  cols=colgenaid(stp); cols_ilincs_Con=colgenaid(conndata_stp); cols_cmap_kd_Con=colgenaid(conndata_cmap_kd_stp); cols_cmap_oe_Con=colgenaid(conndata_cmap_oe_stp)
   if (!is.null(conndata_stp2)|!is.null(conndata_cmap_kd_stp2)|!is.null(conndata_cmap_oe_stp2)) {cols_ilincs_Con2=colgenaid(conndata_stp2); cols_cmap_kd_Con2=colgenaid(conndata_cmap_kd_stp2); cols_cmap_oe_Con2=colgenaid(conndata_cmap_oe_stp2)}
 
  # Make a legend for the color coded names. Additionally, another legend for the missing data points is also included (lgd2), and a final legend list is generated by using both
  lgd1=Legend(labels = c("Pro-coagulants", "Anti-coagulants", "Other drugs"), title = "Compound types", type = "points", pch = 15, legend_gp = gpar(col = c("red4","dimgray","black"), fontsize=48),grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp =gpar(fontsize = 15, fontface = "bold"),title_gp = gpar(fontsize = 15, fontface = "bold") )
  lgd2 = Legend(labels = '', legend_gp = gpar(col=nacols, fill=nacols, fontsize=48), title = "NA data\n color",  pch=15, grid_height = unit(0.5, "cm"),grid_width = unit(1, "cm"),labels_gp =gpar(fontsize = 15, fontface = "bold"),title_gp = gpar(fontsize = 15, fontface = "bold"))
  
  lgd_list = packLegend(lgd1, lgd2,  direction = "horizontal", column_gap = unit(5, "mm"), row_gap = unit(1, "cm"))
  
  heatmaps=list()        
  #Generate a heatmap annotation with our reference data, of provided. Since heatmap annotations can get tricky, here we are handling them based on the column numbers of the reference data while assigning appropriate legends and colors that were generated by heatmaptags function previously
  if (!is.null(refdata)) {if (ncol(refdata)==1) {
    ha = HeatmapAnnotation(E2_JAS=as.matrix(refdata_tags$dataset), col = list(E2_JAS=as.function(refdata_tags$f1)), annotation_legend_param = list(E2_JAS= list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7, fontface = "bold"),title_gp = gpar(fontsize = 7, fontface = "bold"), title=annotationLegendTitle)), na_col = nacols, annotation_name_side = 'right'); if(refdata_tags$f1tag){ha@anno_list$E2_JAS@legend_param$at=refdata_tags$at; ht1@anno_list$E2_JAS@legend_param$labels=refdata_tags$labels}
  }else{tmp=t(apply(refdata_tags$dataset, 1, cbind)); colnames(tmp)=colnames(refdata_tags$dataset)
  ha = HeatmapAnnotation(E2_JAS=tmp, col = list(E2_JAS=as.function(refdata_tags$f1)), annotation_legend_param = list(E2_JAS= list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7, fontface = "bold"),title_gp = gpar(fontsize = 7, fontface = "bold"), title=annotationLegendTitle)), na_col = nacols, annotation_name_side = 'right'); if(refdata_tags$f1tag){ha@anno_list$E2_JAS@legend_param$at=refdata_tags$at; ht1@anno_list$E2_JAS@legend_param$labels=refdata_tags$labels}
  }} else{ha=NULL}
  
  # Change directory to keep the heatmap figures which will not overwrite one another
  # setwd(gsub(x=currentdir,pattern = '(github_stuff).*','\\1'))
  setwd(gsub(x=currentdir,pattern = '(github/cmap_ilincs).*','\\1')) 
  filenames=paste0(gsub('-','',Sys.Date()),'_',paste0(c(CellLine,secondCellLine), collapse = '_'),ifelse(is.null(ha),'','_Annotated'),'_CombinedHeatmaps_1.jpeg')
  while (file.exists(filenames)) {
    filenames=paste0(gsub('_[^_]*$', '', filenames),'_', as.numeric(gsub('\\.jpeg','',gsub("^.+_", "", filenames)))+1,'.jpeg')
  }
  
  #Generate heatmaps
  jpeg(filenames, units="in", width=width_heat, height=height_heat, res=res_heat)  
  
  ht1=Heatmap(t(stp_tags$dataset),cluster_rows = F,  col = as.function(stp_tags$f1), name= "SwissTarget similarity scores",top_annotation=ha,row_names_gp = gpar(fontsize = 15, fontface = "bold", col=cols),column_names_gp = gpar(fontsize = 15, fontface = "bold"),
              heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")),cluster_columns=F,row_title = "SwissTargetPrediction", row_title_rot = 0, height = 3,row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"), na_col = nacols,use_raster=F); if(stp_tags$f1tag){ht1@matrix_legend_param$at=stp_tags$at; ht1@matrix_legend_param$labels=stp_tags$labels}
  
  ht2=Heatmap(t(conndata_stp_tags$dataset),cluster_rows = F, cluster_columns = F, col = as.function(conndata_stp_tags$f1), name= paste0(CellLine," z-scores \niLINCS Connectivity (kd)"),row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_ilincs_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = paste0(CellLine," iLINCS Connectivity (kd)"), row_title_rot = 0, row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"), heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = nacols,height = 1.5,use_raster=F); if(conndata_stp_tags$f1tag){ht2@matrix_legend_param$at=conndata_stp_tags$at;ht2@matrix_legend_param$labels=conndata_stp_tags$labels}
  
  ht3=Heatmap(t(conndata_cmap_kd_stp_tags$dataset),cluster_rows = F, cluster_columns = F, col = as.function(conndata_cmap_kd_stp_tags$f1), name= paste0(CellLine, " CMap Connectivity (kd)"),row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_kd_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = paste0(CellLine," CMap Connectivity (kd)"), row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = nacols,height = 1.5,use_raster=F); if(conndata_cmap_kd_stp_tags$f1tag){ht3@matrix_legend_param$at=conndata_cmap_kd_stp_tags$at;ht3@matrix_legend_param$labels=conndata_cmap_kd_stp_tags$labels}
  
  ht4=Heatmap(t(conndata_cmap_oe_stp_tags$dataset),cluster_rows = F, cluster_columns = F, col = as.function(conndata_cmap_oe_stp_tags$f1), name= paste0(CellLine, " CMap Connectivity (oe)"),row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_oe_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = paste0(CellLine, " CMap Connectivity (oe)"), row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = nacols,height = 1.5,use_raster=F); if(conndata_cmap_oe_stp_tags$f1tag){ht4@matrix_legend_param$at=conndata_cmap_oe_stp_tags$at;ht4@matrix_legend_param$labels=conndata_cmap_oe_stp_tags$labels}
  
  if (!is.null(conndata_stp2)|!is.null(conndata_cmap_kd_stp2)|!is.null(conndata_cmap_oe_stp2)) {
    ht12=Heatmap(t(conndata_stp2_tags$dataset),cluster_rows = F, cluster_columns = F,col = as.function(conndata_cmap_kd_stp2_tags$f1), name= paste0(secondCellLine," z-scores \niLINCS Connectivity (kd)"),row_names_gp = gpar(fontsize = 8, fontface = "bold", col=cols_ilincs_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = paste0(secondCellLine," iLINCS Connectivity (kd)"), row_title_rot = 0,row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"), heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = nacols,height = 2,use_raster=F); if(conndata_stp2_tags$f1tag){ht12@matrix_legend_param$at=conndata_stp2_tags$at;ht12@matrix_legend_param$labels=conndata_stp2_tags$labels}
    
    ht13=Heatmap(t(conndata_cmap_kd_stp2_tags$dataset),cluster_rows = F, cluster_columns = F, col = as.function(conndata_cmap_kd_stp2_tags$f1), name= paste0(secondCellLine, " CMap Connectivity (kd)"),row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_kd_Con2),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = paste0(secondCellLine, "  CMap Connectivity (kd)"), row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = nacols,height = 1.5,use_raster=F); if(conndata_cmap_kd_stp2_tags$f1tag){ht13@matrix_legend_param$at=conndata_cmap_kd_stp2_tags$at;ht13@matrix_legend_param$labels=conndata_cmap_kd_stp2_tags$labels}
    
    ht14=Heatmap(t(conndata_cmap_oe_stp2_tags$dataset),cluster_rows = F, cluster_columns = F, col = as.function(conndata_cmap_oe_stp2_tags$f1), name= paste0(secondCellLine, " CMap Connectivity (oe)"),row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_oe_Con2),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = paste0(secondCellLine, " CMap Connectivity (oe)"), row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = nacols,height = 1.5,use_raster=F); if(conndata_cmap_oe_stp2_tags$f1tag){ht14@matrix_legend_param$at=conndata_cmap_oe_stp2_tags$at;ht14@matrix_legend_param$labels=conndata_cmap_oe_stp2_tags$labels}
    
    #Store the heatmaps as a list in vertical order, then draw them
    ht_list=ht1 %v% ht2  %v% ht3  %v% ht4  %v% ht12  %v% ht13  %v% ht14
  } else{ht_list=ht1 %v% ht2  %v% ht3  %v% ht4}
  
  tmpheat=draw(ht_list, heatmap_legend_side = "bottom", main_heatmap="SwissTarget similarity scores",annotation_legend_list = lgd_list,annotation_legend_side = "bottom",merge_legend = TRUE, column_title = paste0(paste0(c(CellLine,secondCellLine),collapse = '&')," - SwissTarget and Connectivity data compilation"), column_title_gp = gpar(fontsize = 16, fontface='bold'))
  heatmaps[['Multimap']]=tmpheat
  dev.off()
}

#enrichr_output: Retrieves EnrichR results. It requires a vector of drug list, p-value cut off, name of the database used (kegg="KEGG_2019_Human","GO_Biological_Process_2018" (see EnrichR sources)),
##and a specific direction (for the genes either up or down regulated, if inquired). If the direction is not specified, the calculations are done for all the DEGs (up and downregulated ones separately). The workflow can handle data from iLINCS DEGs, iLINCS connectivity and iLINCS CMap results. Accordingly, please be careful while assigning cmapranks (TRUE for CMap connectivity only) and ilincsranks (TRUE for iLINCS gene knock-down connectivity data only)
##The output file contains the lists for the EnrichR results in the specified directions. Each list has three subcategories in data frame format:
##(1) List of pathways and associated combination scores, (2) Expanded list that contains the full range of EnrichR results, and (3) List of pathways and associated genes coming from the input gene expression data. 
## Note: there is some redundancy in the code. In addition, if there are replicates of the same compound, the enrichment analysis only takes into the first occurring one. It currently runs fine for the limited set of data that we have, but an update is needed for cleaner and more comprehensive analyses 

enrichr_output=function(drug_list,dataset,logDEcutoff=0,pcutoff,dbs,direction=NULL, cmapranks=F, cmaprank_cutoff,ilincsranks=F,zscores_cutoff=5){
  require(dplyr);require(purrr);require(enrichR)
  result=list()
  if (is.null(direction)) {ds=c('complete','up','down')
  }else {ds=tolower(direction)
  while (sum(grepl(ds, c('complete','up','down'),ignore.case=T))==0) { ds=tolower(as.character(readline(prompt='Please check out misspelling, if any: complete/up/down?')))}
  }
  for (d in ds) {
    if (exists('go_output')) {rm(go_output)}
    if (exists('termgenes')) {rm(termgenes)}
    k=list()
    #Analyses for the differentially expressed genes data retrieved from iLINCS
    if (!cmapranks&!ilincsranks) {
      for (drugs in drug_list) {
        if (!exists("go_output")){
          #Filter the data based on log fold changes and p-values. Store name of the terms, enrichment scores, relevant genes from the associated pathways and full list of results
          go_input=dataset[,grepl(drugs, sub('_logDE|_pval','',colnames(dataset),ignore.case = T))]
          go_input=go_input[order(go_input[,2]),]
          #Store the results separately for (i) all DEGs ('complete'), (ii) only upregulated genes ('up') and only downregulated genes ('down')
          if (d=='complete') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&abs(go_input[,1])>abs(logDEcutoff),]), dbs)}
          else if (d=='up') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]>abs(logDEcutoff),]), dbs)}
          else {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]<(-abs(logDEcutoff)),]), dbs)}
          k[[drugs]]=enriched[[dbs[[1]]]]
          go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          colnames(go_output)=c('Term',drugs)
          colnames(termgenes)=c('Term',drugs)
          rm( list=Filter( exists, c('go_input','enriched') ) )}
        if (exists("go_output")){
          go_input=dataset[,grepl(drugs, sub('_logDE|_pval','',colnames(dataset)),ignore.case = T)]
          go_input=go_input[order(go_input[,2]),]
          if (d=='complete') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&abs(go_input[,1])>abs(logDEcutoff),]), dbs)
          } else if (d=='up') {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]>abs(logDEcutoff),]), dbs)
          } else {enriched=enrichr(rownames(go_input[go_input[,2]<pcutoff&go_input[,1]<(-abs(logDEcutoff)),]), dbs)}
          k[[drugs]]=enriched[[dbs[[1]]]]
          enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          colnames(enriched_go)=c('Term',drugs)
          colnames(enriched_go_termgenes)=c('Term',drugs)
          go_output=merge(go_output, enriched_go, by='Term',all=T)
          termgenes=merge(termgenes, enriched_go_termgenes, by='Term',all=T)
          rm( list=Filter( exists, c('go_input','enriched_go') ) )}}
    } else if(cmapranks&!ilincsranks) {
      #Analyses for the kd/oe based gene signature data from ConnectivityMap Connections
      for (drugs in drug_list) {
        if (!exists("go_output")){
          #Similarly, collect the results for the pathways associating with (i) complete: all gene signatures above ABSOLUTE cmapranks threshold, (ii) up: all overexpression signatures above the positive ABSOLUTE threshold plus all knockdown signatures below the negative ABSOLUTE threshold, and (iii) down: all overexpression signatures below the negative ABSOLUTE threshold and all knockdown signatures above the positive ABSOLUTE threshold
          go_input=dataset[,c(1,grep(drugs,colnames(dataset),ignore.case = T))]
          if (d=='complete') {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks=go_input[abs(go_input[,2])>abs(cmaprank_cutoff),]
          upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>abs(cmaprank_cutoff),]
          upranks1=upranks1[grepl('_oe',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-abs(cmaprank_cutoff)),]
          upranks2=upranks2[grepl('_kd',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>abs(cmaprank_cutoff),]
          upranks1=upranks1[grepl('_kd',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-abs(cmaprank_cutoff)),]
          upranks2=upranks2[grepl('_oe',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          colnames(go_output)=c('Term',drugs)
          colnames(termgenes)=c('Term',drugs)
          rm( list=Filter( exists, c('go_input','enriched') ) )}
        if (exists("go_output")){
          go_input=dataset[,c(1,grep(drugs,colnames(dataset),ignore.case = T))]
          if (d=='complete') {upranks=go_input[abs(go_input[,2])>abs(cmaprank_cutoff),]
          upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>abs(cmaprank_cutoff),]
          upranks1=upranks1[grepl('_oe',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-abs(cmaprank_cutoff)),]
          upranks2=upranks2[grepl('_kd',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {go_input=go_input[grepl('_kd|_oe',go_input[,1]),]
          upranks1=go_input[go_input[,2]>abs(cmaprank_cutoff),]
          upranks1=upranks1[grepl('_kd',row.names(upranks1)),]
          upranks2=go_input[go_input[,2]<(-abs(cmaprank_cutoff)),]
          upranks2=upranks2[grepl('_oe',rownames(upranks2)),]
          upranks=rbind(upranks1,upranks2)
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(gsub('_kd|_oe','',rownames(upranks)), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          colnames(enriched_go)=c('Term',drugs)
          colnames(enriched_go_termgenes)=c('Term',drugs)
          go_output=merge(go_output, enriched_go, by='Term',all=T)
          termgenes=merge(termgenes, enriched_go_termgenes, by='Term',all=T)
          rm( list=Filter( exists, c('go_input','enriched_go') ) )}}
    }else if(!cmapranks&ilincsranks) {
      #Analyses for the kd based gene signature data from iLINCS perturbations data
      for (drugs in drug_list) {
        if (!exists("go_output")){
          #Similar to the previous steps, collect the results for the pathways associating with (i) complete: all gene signatures above ABSOLUTE ilincsranks threshold (zscores_cutoff), (ii) up:  all knockdown signatures below the negative ABSOLUTE threshold, and (iii) down: all knockdown signatures above the positive ABSOLUTE threshold
          go_input=dataset[,c(1,grep(pattern = drugs,x = colnames(dataset),ignore.case = T))]
          rownames(go_input)=go_input[,1]
          if (d=='complete') {upranks=go_input[abs(go_input[,2])>abs(zscores_cutoff),]
          upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(rownames(upranks), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') { upranks=go_input[go_input[,2]<(-abs(zscores_cutoff)),]
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(rownames(upranks), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          else {upranks=go_input[go_input[,2]>abs(zscores_cutoff),]
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(rownames(upranks), dbs)
          rm( list=Filter( exists, c('go_input','upranks','upranks1','upranks2') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          go_output=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          colnames(go_output)=c('Term',drugs)
          colnames(termgenes)=c('Term',drugs)
          rm( list=Filter( exists, c('go_input','enriched') ) )}
        if (exists("go_output")){
          go_input=dataset[,c(1,grep(drugs,colnames(dataset),ignore.case = T))]
          rownames(go_input)=go_input[,1]
          if (d=='complete') {upranks=go_input[abs(go_input[,2])>abs(zscores_cutoff),]
          upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(rownames(upranks), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else if (d=='up') {upranks=go_input[go_input[,2]<(-abs(zscores_cutoff)),]
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(rownames(upranks), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          else {upranks=go_input[go_input[,2]>zscores_cutoff,]
          upranks=upranks[order(abs(upranks[,2]),decreasing=T),]
          enriched=enrichr(rownames(upranks), dbs)
          rm( list=Filter( exists, c('go_input','upranks') ) )}
          k[[drugs]]=enriched[[dbs[[1]]]]
          enriched_go=enriched[[dbs[[1]]]][,c('Term', 'Combined.Score')]
          enriched_go_termgenes=enriched[[dbs[[1]]]][,c('Term', 'Genes')]
          colnames(enriched_go)=c('Term',drugs)
          colnames(enriched_go_termgenes)=c('Term',drugs)
          go_output=merge(go_output, enriched_go, by='Term',all=T)
          termgenes=merge(termgenes, enriched_go_termgenes, by='Term',all=T)
          rm( list=Filter( exists, c('go_input','enriched_go') ) )}}
    }
    go_output=go_output[,-2]
    colnames(go_output)[grepl('\\.y',colnames(go_output))]=gsub('\\.y','',colnames(go_output)[grepl('\\.y',colnames(go_output))])
    colnames(go_output)[grepl('\\.',colnames(go_output))]=gsub('\\.','_',colnames(go_output)[grepl('\\.',colnames(go_output))])
    go_output[is.na(go_output)]=0
    termgenes=termgenes[,-2]
    colnames(termgenes)[grepl('\\.y',colnames(termgenes))]=gsub('\\.y','',colnames(termgenes)[grepl('\\.y',colnames(termgenes))])
    colnames(termgenes)[grepl('\\.',colnames(termgenes))]=gsub('\\.','_',colnames(termgenes)[grepl('\\.',colnames(termgenes))])
    termgenes[is.na(termgenes)]='---'
    result[[d]]=list(go_output, k, termgenes)
    names(result[[d]])=c(paste0(names(dbs),'_EnrichedTermsScores_',d),paste0(names(dbs),'_Expanded_',d),paste0(names(dbs),'_TermGenes_',d))
    rm(k)}
  for (d in ds) {
    for (i in 1:length(result[[d]])) {
      names(result[[d]][[i]])=gsub('\\.','_',names(result[[d]][[i]]))}}
  return(result)}


#slicenricher: CMap ranks are default to FALSE, unless users provide a rank list of kd/oe studies from CMap connectivity data. The functions, slices and collects specified data (based on Combined score cutoff) from the enrichr_output function. Accordingly, it retrieves pathway enrichment scores from up or downregulated genes ('complete' option accounts both up and downregulated genes together). Based on the enrichment scores, users can visually inspect pathway enrichment profiles across the treatments either by heatmaps or corplots
##Currently, it only utilizes the EnrichedTermsScores objects from the enrichr_output data.
##Specifications rely on the combscore cutoff values (see https://maayanlab.cloud/Enrichr/help#background for more information of this metric) for each reference compound. List of the drugs can be fetched into data as a drug_list, or the function prompts users to specify the names of the compounds to be compared.
##Combined scores are presented on log2 scale. Accordingly, minimum combination score value is capped to 1 (hence log2(1)=0 across the representations)
##Main outputfile contains a list for the pathway names passing the threshold and their combination scores for the respective compounds inquired by the user.
##Lists are given in separate sections for each direction (up/down/complete)). In addition, the correlation plots (by default 'spearman') and
##heatmaps are produced and saved in the current directory, unless indicated otherwise.
##Heatmap label aesthetics are also improved to indicate whether the drugs of interest are coagulant in the literature and/or in our in-house phenotypic screenings.
##Pathways onthologies that matches with thrombosis keyword patterns are also manually included. A warning can be included to change the cutoff value, in case there are not many enriched pathways, hence errors based on non-finite calculations
##Heatmap files' width, height and resolution parameters can be user specified, as well. For the future updates, the fucntion can be incorporated into a dashboard or a shinyapp to give better control over the aesthetics and dynamic features
slicenricher=function(enrichrlist, direction=NULL, combscore_cutoff=200, refcompounds=NULL, cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL='EnrichR', width_heat=18, height_heat=12, res_heat=300){
  require(tidyverse); require(ggplot2); require(ggord);require(ggfortify);
  require(ComplexHeatmap);require(circlize);
  tryCatch({if (!cmapranks) {drugnames=gsub('\\.','_',gsub('.*_','',gsub('_[^_]*$','',colnames(as.data.frame(enrichrlist[[1]][1])))))
  } else {drugnames=gsub('.*_','',names(enrichrlist[[1]][[1]]))}
    if (is.null(direction)) {ds=c('complete','up','down')
    }else {ds=tolower(direction)
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
    proco=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin','U46619')
    antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
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
      col_paths=rep('black',nrow(tmp))
      col_paths[grepl(paste(bloods, collapse='|'), rownames(tmp),ignore.case=T)]='blue'
        
      lgd1=Legend(labels = c("Pro-coagulants", "Anti-coagulants", "Other drugs"), title = "Compound types", type = "points", pch = 18, legend_gp = gpar(col = c("red4","dimgray","black"), fontsize=24))
      lgd2 = Legend(labels ='', title = "Compound types", pch = 15, legend_gp = gpar(col = 'blue', fill='blue', fontsize=24))
          
      lgd_list = packLegend(lgd1, lgd2,  direction = "horizontal", column_gap = unit(5, "mm"), row_gap = unit(1, "cm"))
          
      if (cor_outputs) {
          par(cex.main=1)
            
          f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
          jpeg(paste0(deparse(substitute(enrichrlist)),'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substring(paste(refcompounds,collapse='-'),1,15),toupper(cormethod),'corplot_',toupper(d),".jpeg"), units="in", width=18, height=12, res=300)
            ht1=Heatmap(cormat_go,clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f1, name="Pathway combination score correlations (EnrichR)", row_names_gp=gpar(fontsize=max(10-nrow(cormat_go)*0.05,2), fontface="bold", col=cols),column_names_gp=gpar(fontsize=max(12-ncol(cormat_go)*0.05,2), fontface="bold", col=cols),heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(8, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold")), cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
              grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=max(12-nrow(cormat_go)*0.05,2), fontface='bold'))})
            tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - ",d,'\n',cormethod," correlations on gene pathways \n (EnrichR ",paste(refcompounds,collapse='-'), "\n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
            k=paste0('Corplot','_',d)
            corrplots[[k]]=tmpfig
            dev.off()}
          
          #heatmaps
          if (heatmap_outputs) {
            f2=colorRamp2(seq(max(tmp), min(tmp), length=3),  c("#3f7fff", "#EEEEEE", "red"))
            jpeg(paste0(deparse(substitute(enrichrlist)),'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substring(paste(refcompounds,collapse='-'),1,15),"_Heatmap_",toupper(d),".jpeg"), units="in", width=width_heat, height=height_heat, res=res_heat)
            
            ht2=Heatmap(as.matrix.data.frame(tmp),clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f2, name=paste0(toupper(d), " - Pathway combinations scores \n","(EnrichR ",paste(refcompounds,collapse='-'), "\n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp)*0.05,2), fontface="bold", col=col_paths),column_names_gp=gpar(fontsize=max(12-nrow(tmp)*0.05,2), fontface="bold",col=cols), heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"),legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=18, fontface="bold"),title_position="topcenter"),use_raster=F)
            tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
            k=paste0('Heatmap','_',d)
            heatmaps[[k]]=tmpheat
            dev.off()
            rm(ind); rm(pathways); rm(tmp)
          }}}
    rm(drugnames)
    return(result)},error=function(e){print('Warning: There was an error. Please, check if the parameters are correct. Combined score cut-off might be too stringent for some of the compounds. Please, check the outputs or you can also try lower combination score thresholds')})
}

#slicenricher_updown: Similar to slicenricher. It allows comparing different groups of drugs (such as pro-coagulants vs anti-coagulants) to see whether there are any pathways that move in same/opposite directions.
slicenricher_updown=function(enrichrlist, combscore_cutoff=200, contr1=NULL,contr2=NULL,classContr1=NULL,classContr2=NULL,cmapranks=F,cormethod='spearman',enrichmentdb=NULL, cor_outputs=T, heatmap_outputs=T, enrichmentTOOL='EnrichR', repositoryName=NULL, width_heat=18, height_heat=12, res_heat=300){
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
    
    proco=c("heparin","thrombin",'phylloquinone', 'menadione','desmopressin','U46619')
    antico_thr=c('Warfarin','APIXABAN','Dabigatran', 'argatroban', 'skatole', 'phenindione')
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
    result[[paste(c(ds[[d]][1],ds[[d]][2]),collapse='_')]]=tmp
    tmp$Term=NULL
    tmp=tmp[,colSums(tmp)>0]
    if (length(tmp)>0) {
      
      #corplot
      cormat_go=signif(cor(tmp,method=cormethod),2)
      cols=rep('black', ncol(tmp))
      cols[grepl(paste(proco, collapse='|'), colnames(tmp),ignore.case=T)]='red4'
      cols[grepl(paste(antico_thr, collapse='|'), colnames(tmp),ignore.case=T)]='dimgray'
      col_paths=rep('black',nrow(tmp))
      col_paths[grepl(paste(bloods, collapse='|'), rownames(tmp),ignore.case=T)]='blue'
      lgd1=Legend(labels = c("Pro-coagulants", "Anti-coagulants", "Other drugs"), title = "Compound types", type = "points", pch = 18, legend_gp = gpar(col = c("red4","dimgray","black"), fontsize=24))
      lgd2 = Legend(labels ='', title = "Compound types", pch = 15, legend_gp = gpar(col = 'blue', fill='blue', fontsize=24))
          
      lgd_list = packLegend(lgd1, lgd2,  direction = "horizontal", column_gap = unit(5, "mm"), row_gap = unit(1, "cm"))
          
      if (cor_outputs) {
          par(cex.main=1)
          f1=colorRamp2(seq(-abs(max(cormat_go)), abs(max(cormat_go)), length=3), c("#3f7fff", "#EEEEEE", "red"))
          jpeg(paste0(deparse(substitute(enrichrlist)),'_',paste(c(repositoryName,enrichmentdb,enrichmentTOOL),collapse='-'),'_',substr(paste(c(contr1,contr2),collapse='-'),1,30),'_',toupper(cormethod),'corplot_',toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='_')),".jpeg"), units="in", width=width_heat, height=height_heat, res=res_heat)
          ht1=Heatmap(cormat_go,clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f1, name=paste0("Pathway combination score correlations (",repositoryName," - ", enrichmentTOOL,")"), row_names_gp=gpar(fontsize=max(10-nrow(cormat_go)*0.05,2), fontface="bold", col=cols),column_names_gp=gpar(fontsize=max(12-ncol(cormat_go)*0.05,2), fontface="bold", col=cols), heatmap_legend_param=list(legend_direction="horizontal", title_position="topcenter", grid_height=unit(0.6, "cm"),legend_width=unit(max(10-nrow(cormat_go)*0.05,2), "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=15, fontface="bold")), cell_fun=function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(cormat_go[i, j], x, y,gp=gpar(fontsize=max(12-nrow(cormat_go)*0.05,2), fontface='bold'))})
          tmpfig=draw(ht1, heatmap_legend_side="bottom", column_title=paste0("Corplot - ",toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='vs')),'\n',cormethod," correlations on ",repositoryName,'-',enrichmentdb," pathways \n (", enrichmentTOOL," ",substr(paste(c(contr1,contr2),collapse='-'),1,30), "\n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), column_title_gp=gpar(fontface="bold", line=5), annotation_legend_list=lgd_list,annotation_legend_side="right")
            k=paste0('Corplot','_',d)
            corrplots[[k]]=tmpfig
            dev.off()}
          #heatmaps
          if (heatmap_outputs) {
            f2=colorRamp2(seq(max(tmp), min(tmp), length=3),  c("#3f7fff", "#EEEEEE", "red"))
            jpeg(paste0(deparse(substitute(enrichrlist)),'_',repositoryName,'_',enrichmentdb,'_',as.character(enrichmentTOOL),'_',substr(paste(c(contr1,contr2),collapse='-'),1,30),'__Heatmap_',toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='_')),".jpeg"), units="in", width=width_heat, height=height_heat, res=res_heat)
            
            ht2=Heatmap(as.matrix.data.frame(tmp),clustering_method_rows="ward.D",clustering_method_columns ="ward.D", col=f2, name=paste0(toupper(paste(c(ds[[d]][1],ds[[d]][2]),collapse='vs')), " - Pathway combinations scores \n (",paste(c(repositoryName,enrichmentdb,enrichmentTOOL),collapse='-'),"\n",substr(paste(c(contr1,contr2),collapse='-'),1,30), "... \n log2(Combination Scores >", combscore_cutoff,"), n=", nrow(tmp),")"), row_names_gp=gpar(fontsize=max(10-nrow(tmp)*0.05,2), fontface="bold", col=col_paths),column_names_gp=gpar(fontsize=max(12-ncol(tmp)*0.05,2), fontface="bold",col=cols), heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height=unit(0.8, "cm"),legend_width=unit(10, "cm"),labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize=15, fontface="bold"),title_position="topcenter"),use_raster=F)
            tmpheat=draw(ht2, heatmap_legend_side="top",  column_title_gp=gpar(fontface="bold", line=5),annotation_legend_list=lgd_list,annotation_legend_side="bottom")
            k=paste0('Heatmap','_',d)
            heatmaps[[k]]=tmpheat
            dev.off()
            rm(ind); rm(pathways); rm(tmp)
          }}}
    rm(drugnames)
    return(result)
  }, error=function(e){print('Warning: There was an error. Please, check if the parameters are correct. Otherwise, some contrasts may not have been processed (due to lack of overlaps between the query compounds). Please, check the outputs or you can also try lower combination score thresholds')})}



##upsetlister, an implemented function to provide the input for the upset plots and venn interaction terms. So that, from many pathways in the list, users can see how many of them are shared among the reference compounds and across different directions. If the lists do not seem complicated (meaning that they are not enriched with multiple pathways and directions) upsetlist and further plots will not be needed.Slicenricher outputs are preferred and transformed into upset plot compatible lists while collecting pathways that are above a certain threshold. Reference compounds and regulation direction (up, down or complete) can be further specified.
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


#termgenesmajor, retrieves the names of the significant genes relating with the previous enrichr pathway outputs.
##Hence, it allows seeing overexpression or downregulation of which genes (from CMap/iLINCS) can be further important.
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

#querylister allows conversion of queries to upsetplot and upsetlist compatible formats. It further requires intersectionSets from Venn function (see vennerable). Numbers of queries for each comparisons can be adjusted with the lowerlimit and upperlimit options
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
##In addition, pathways of interests and their respective genes from CMap/iLINCS studies are imported into excel sheets inside the working directory, in the same order that they appear on the upsetplot. Maximum total numbers of pathway terms in the querylists (npathways argument) are defaulted to 50 to prevent overcrowding the upsetplot figure while the full range of pathways between the queries are still accessible within the output excel sheet. As the outputs can be enriched with multiple pathways, it may become hard to navigate the plot and figure out which pathway is on the intersection of the relevant exposures. Future developments can include interactive upsetplots to make the navigation easier - upsetplots and interactive plots might have some incompatibilities as I checked. It may take a while to implement this. For now, you can focus on the pathway enriched intersections (size>2), check the excel sheet outputs for smaller intersections or in some cases making additional upsetplots for only selected numbers of compounds (you can do this refinement during the upsetlister() function step) 

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
  
  #Get your color palette for the upsetplot representations
  nb.cols=min(npathways, length(querylist))
  mycolors=colorRampPalette(brewer.pal(8, colorpal))(nb.cols)
  
  #Generate tables to store information from the upsetplot intersections. So that, you can see the relevant information (pathways, genes) inside the excel sheets
  pathtables=termgenesmajor(enrichrlist, direction=NULL,refcompounds= refcompounds)
  gridtable=setNames(data.frame(matrix(ncol=ncol(pathtables), nrow=0)),colnames(pathtables))
  queries=list()
  
  #Make a list of your query compounds and see how their signatures overlap with one another
  for (i in 1:length(querylist)) {
    queries[[i]]=list(query=intersects,params=querylist[[i]], color=mycolors[i], active=T)}
  
  
  if (upsetFigure) {
    
    jpeg(paste0(gsub('-','',Sys.Date()),"_",as.character(repositoryName),'_',as.character(enrichmentTOOL),cell_line,'_',paste(substr(str_to_title(refcompounds),1,4),collapse=''),"_CombScore",combscore_cutoff,"_",enrichmentdb,"_UpsetPlot.jpeg"), units="in", width=14, height=10, res=300)
    
    upsetpl=upset(fromList(upsetlist), nsets=sum(lengths(temp@IntersectionSets)>0), queries=queries, order.by="freq",nintersects=length(querylist), mainbar.y.label=ylabel)
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

#Implemented and modified function (https://stackoverflow.com/questions/25892912/matching-words-in-strings-with-variables-in-r) - sometimes your data can contain some cells that has multiple strings attached side by side with a delimiter. With this function, you can count how many times each specific string occurs. It was useful for me to retrieve gene names from groups of genes and get the most frequently observed genes across these groups.  
freqfunc2=function(x, n=30){
  tail(sort(table( unlist(strsplit(c(na.omit(as.character(unlist(x, use.names=FALSE)))), ";",fixed=T)))), n)
}


###pairs with p-val (https://stackoverflow.com/questions/31851537/set-backgroud-color-of-a-panel-in-pairs-function-call). It allows making nice looking correlation plots. 
pairs.pal=function(cols_BR=NULL){require(RColorBrewer)
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
    } else{ text(0.5, 0.75, paste("r=",round(r,2)),cex=2.5) # prints correlatoin coefficient
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
  }}

# bind Ensembl ID to results and name the columns (http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html) identifiers of a list of genes from one type to another (such as gene id -> gene name) 
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

#Adopted function (https://www.r-bloggers.com/2018/04/strsplit-but-keeping-the-delimiter/) allows fine splitting of long appended string into pieces by using specific delimiters
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

##biomart_orthologs: Retrieves ortholog genes for two organisms. Zebrafish and human are the default species here. You can change try alternatives like 'mmusculus' for mice. See BioMaRT for details.
##In addition, for some reason getLDS() does not work with recent versions of Ensembl BioMart requests. We can refer to one of the older repositories. 
## You can also access the data directly from the Ensembl BioMaRT webpage (http://useast.ensembl.org/biomart/martview/), and handle the remaining parts manually
biomart_orthologs=function(species1='drerio', species2='hsapiens', host="https://dec2021.archive.ensembl.org", repoOLDorNEW=NULL){
  require(biomaRt)
  if(isTRUE(repoOLDorNEW)&&repoOLDorNEW=='new') {
    zebra.mart=useEnsembl(biomart="ensembl", dataset=paste0(species1,"_gene_ensembl"))
    human.mart=useEnsembl(biomart="ensembl", dataset=paste0(species2,"_gene_ensembl"))
  } else { zebra.mart=biomaRt::useMart(host=host, "ENSEMBL_MART_ENSEMBL", dataset=paste0(species1,"_gene_ensembl"))
  human.mart=biomaRt::useMart(host=host, "ENSEMBL_MART_ENSEMBL", dataset=paste0(species2,"_gene_ensembl"))}
  annot_table=getLDS(mart = zebra.mart, attributes = c('ensembl_gene_id','external_gene_name'),
                     martL = human.mart, attributesL = c('ensembl_gene_id','external_gene_name','chromosome_name','gene_biotype'))
  annot_table=annot_table[,c('Gene.name','Gene.name.1')]
  annot_table=annot_table[!duplicated(annot_table),]; annot_table=annot_table[!rowSums(annot_table=='')>0,];
  return(annot_table)
  rm( list=Filter(exists, c('annot_table')))
}

## Modified plotPCA from DESeq2 package (https://www.biostars.org/p/243695/). Shows the Names of the Samples (the first col of SampleTable), and uses ggrepel pkg to plot them conveniently.
# @Santosh Annand 10.02.2017
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
