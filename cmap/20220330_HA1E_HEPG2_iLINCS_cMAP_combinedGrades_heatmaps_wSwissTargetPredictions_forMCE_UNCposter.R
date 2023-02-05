rm(list=ls())
gc()

library(cmapR)
library(tidyverse)
library(dplyr)
# BiocManager::install("org.Hs.eg.db")
# suppressMessages(library(org.Hs.eg.db))
# install.packages('reshape')
# install.packages("reshape2")
library(reshape)
library(reshape2)
library(ggord)
#devtools::install_github("sinhrks/ggfortify")
# install.packages('ggfortify')
library(ggfortify)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(rlist)

setwd('../SwissTargetPrediction/')
stp=readRDS('20220330_SwissTarPred_heatmap_matrix_genenames_forMCE_UNCposter.RDS')
stp=stp[order(rownames(stp)),]
setwd('../20210830_HA1E_HEPG2_iLINCS_cMAP_combinedGrades/')
#HEPG2

zf_e2=read.csv('E2vsDMSO_hsa_formeta_ALL.csv')
zf_e2=zf_e2[,c(1:3)];zf_e2=zf_e2[!duplicated(zf_e2),];colnames(zf_e2)[colnames(zf_e2)=='human_gene_symbol']='symbol'
zf_e2$ppow=log10(zf_e2$E2vsDMSO_pval)
zf_e2=zf_e2 %>% group_by(zf_e2$symbol) %>% mutate_each(funs(mean),-1) %>% distinct
zf_e2$ppow=10^zf_e2$ppow;zf_e2$E2vsDMSO_pval=zf_e2$ppow;zf_e2$ppow=NULL;zf_e2$`zf_e2$symbol`=NULL
zf_e2=as.data.frame(zf_e2)
rownames(zf_e2)=zf_e2$symbol; zf_e2$symbol=NULL;zf_e2=as.data.frame.matrix(zf_e2)
zf_e2=zf_e2[order(rownames(zf_e2)),]
miss_genes_zf=rownames(stp)[!rownames(stp)%in%rownames(zf_e2)]
miss_mat_zf=matrix(data=as.numeric(NA),nrow = length(miss_genes_zf),ncol = ncol(zf_e2))
rownames(miss_mat_zf)=miss_genes_zf; colnames(miss_mat_zf)=colnames(zf_e2)
zf_e2=rbind(zf_e2,miss_mat_zf)
zf_e2_moltar=zf_e2[rownames(zf_e2)%in%rownames(stp),]
zf_e2_moltar=zf_e2_moltar[order(rownames(zf_e2_moltar)),]

#iLINCS connectivity data

##HEPG2
setwd('../iLINCS_conn_kdoe/')
ptrn='iLINCS_complete_conn_HEPG2_'
conndata=iLINCS_kd_connectivitydatafile_compiler(ptrn = ptrn)
conndata=conndata[,grepl(toupper(paste(c('GeneTarget',colnames(stp)), collapse = '|')), toupper(colnames(conndata)),ignore.case = T)]
miss_genes_ilincsConn=rownames(stp)[!rownames(stp)%in%conndata$GeneTarget]
rownames(conndata)=conndata$GeneTarget; conndata$GeneTarget=NULL
miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_ilincsConn),ncol = ncol(conndata))
rownames(miss_mat)=miss_genes_ilincsConn; colnames(miss_mat)=colnames(conndata)
conndata=rbind(conndata,miss_mat)

conndata=conndata[order(rownames(conndata)),]
# conndata=conndata[,c(which(grepl('estradiol|mestranol',colnames(conndata),ignore.case = T)),which(grepl('prazol',colnames(conndata))),which(grepl('calcipotriol|tranilast',colnames(conndata),ignore.case = T)))]
conndata_stp=as.matrix(conndata[rownames(conndata)%in%rownames(stp),])
conndata_stp=conndata_stp[order(rownames(conndata_stp)),]

#CMap connectivity data
setwd('../GeneConnectivityScores/')
ptrn='conn_HEPG2_'; types=c('kd','oe')
rm( list = Filter( exists, c('dataset') ) )

conndata_cmap=connectivitydatafile_compiler(ptrn = ptrn,types = types)
conndata_cmap=conndata_cmap[,grepl(toupper(paste(gsub('6_|24_','',colnames(stp)), collapse = '|')),toupper(colnames(conndata_cmap)), ignore.case = T)]
conndata_cmap=conndata_cmap[,grepl(paste(colnames(stp), collapse = '|'), colnames(conndata_cmap),ignore.case = T)]
# conndata_cmap=conndata_cmap[,c(which(grepl('estradiol|mestranol',colnames(conndata_cmap),ignore.case = T)),which(grepl('prazol',colnames(conndata_cmap),ignore.case = T)),which(grepl('calcipotriol|tranilast',colnames(conndata_cmap),ignore.case = T)))]

#CMap kd 
conndata_cmap_kd=conndata_cmap[grepl('_kd', rownames(conndata_cmap)),]; rownames(conndata_cmap_kd)=gsub('_kd','',rownames(conndata_cmap_kd))
miss_genes_cMapConn_kd=rownames(stp)[!rownames(stp)%in%rownames(conndata_cmap_kd)]
miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_kd),ncol = ncol(conndata_cmap_kd))
rownames(miss_mat)=miss_genes_cMapConn_kd; colnames(miss_mat)=colnames(conndata_cmap_kd)
conndata_cmap_kd=rbind(conndata_cmap_kd,miss_mat)
# for (i in 1:length(miss_genes_hepg2_ilincsConn)) {hepg2_conndata=rbind(hepg2_conndata,c(miss_genes_hepg2_ilincsConn[i],rep(as.numeric(NA),ncol(hepg2_conndata)-1)))}
conndata_cmap_kd=conndata_cmap_kd[order(rownames(conndata_cmap_kd)),]
conndata_cmap_kd_stp=as.matrix(conndata_cmap_kd[rownames(conndata_cmap_kd)%in%rownames(stp),])
conndata_cmap_kd_stp=conndata_cmap_kd_stp[order(rownames(conndata_cmap_kd_stp)),]
rm(miss_mat)
#CMap oe
conndata_cmap_oe=conndata_cmap[grepl('_oe', rownames(conndata_cmap)),]; rownames(conndata_cmap_oe)=gsub('_oe','',rownames(conndata_cmap_oe))
miss_genes_cMapConn_oe=rownames(stp)[!rownames(stp)%in%rownames(conndata_cmap_oe)]
miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_oe),ncol = ncol(conndata_cmap_oe))
rownames(miss_mat)=miss_genes_cMapConn_oe; colnames(miss_mat)=colnames(conndata_cmap_oe)
conndata_cmap_oe=rbind(conndata_cmap_oe,miss_mat)

conndata_cmap_oe=conndata_cmap_oe[order(rownames(conndata_cmap_oe)),]
conndata_cmap_oe_stp=conndata_cmap_oe[rownames(conndata_cmap_oe)%in%rownames(stp),]
conndata_cmap_oe_stp=conndata_cmap_oe_stp[order(rownames(conndata_cmap_oe_stp)),]


f1 = colorRamp2(seq(max(stp), min(stp), length=5), c("dark red","red3","red2","red","white"))
f2 = colorRamp2(seq(max(conndata_stp,  na.rm = TRUE)*-1,max(conndata_stp,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))
f3= colorRamp2(seq(max(zf_e2_moltar$E2vsDMSO_logDE,  na.rm = TRUE)*-1,max(zf_e2_moltar$E2vsDMSO_logDE,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))
f4= colorRamp2(seq(max(conndata_cmap_kd_stp,  na.rm = TRUE)*-1,max(conndata_cmap_kd_stp,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))
f5= colorRamp2(seq(max(conndata_cmap_oe_stp,  na.rm = TRUE)*-1,max(conndata_cmap_oe_stp,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))

# ha_mix_top = HeatmapAnnotation(heatmap_hepg2 = anno_density(t(rseqVal), type = 'heatmap'),heatmap = anno_density(t(rseqVal_ha1e), type = 'heatmap'))
cols<- rep('black', ncol(stp))
cols_ilincs_Con=rep('black', ncol(conndata_stp))
cols_cmap_kd_Con=rep('black', ncol(conndata_cmap_kd_stp))
cols_cmap_oe_Con=rep('black', ncol(conndata_cmap_oe_stp))
# cols_ha1e=rep('black', ncol(rseqVal_ha1e))
#turn red the specified rows in tf
proco=c("mestranol","Ethynyl Estradiol",'EthynylEstradiol','estradiol','17.betaestradiol','17-beta-estradiol','methylandr', 'methandr')
# prev=c('xemestane','irarubicin', 'ivozanib', 'intedanib', 'elamanid')

cols[grepl(paste(proco, collapse = '|'),colnames(stp),ignore.case = T)]= 'red4'
# cols[grepl(paste(prev, collapse = '|'),colnames(stp),ignore.case = T)]= 'grey41'
cols_ilincs_Con[grepl(paste(proco, collapse = '|'),colnames(conndata_stp),ignore.case = T)]= 'red4'
# cols_ilincs_Con[grepl(paste(prev, collapse = '|'),colnames(conndata_stp),ignore.case = T)]= 'grey41'
cols_cmap_kd_Con[grepl(paste(proco, collapse = '|'),colnames(conndata_cmap_kd_stp),ignore.case = T)]= 'red4'
# cols_cmap_kd_Con[grepl(paste(prev, collapse = '|'),colnames(conndata_cmap_kd_stp),ignore.case = T)]= 'grey41'
cols_cmap_oe_Con[grepl(paste(proco, collapse = '|'),colnames(conndata_cmap_oe_stp),ignore.case = T)]= 'red4'
# cols_cmap_oe_Con[grepl(paste(prev, collapse = '|'),colnames(conndata_cmap_oe_stp),ignore.case = T)]= 'grey41'

# cols_ha1e[grepl(paste(proco, collapse = '|'),colnames(rseqVal_ha1e),ignore.case = T)]= 'red4'
lgd_list=Legend(labels = c("E2 derivatives", "Candidate drugs"), title = "Compound types", type = "points", pch = 15, legend_gp = gpar(col = c("red4","black"), fontsize=48),grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp =gpar(fontsize = 15, fontface = "bold"),title_gp = gpar(fontsize = 15, fontface = "bold") )

ha = HeatmapAnnotation(E2_JAS=zf_e2_moltar$E2vsDMSO_logDE, col = list(E2_JAS=f3),
                       annotation_legend_param = list(E2_JAS= list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7, fontface = "bold"),title_gp = gpar(fontsize = 7, fontface = "bold"), title='logFC E2_JAS')), na_col = "gray", annotation_name_side = 'right')


ht1=Heatmap(t(stp),cluster_rows = F, clustering_method_columns = 'average', col = f1, name= "SwissTarget similarity scores",top_annotation = ha,row_names_gp = gpar(fontsize = 15, fontface = "bold", col=cols),column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")),cluster_columns=TRUE,row_title = "SwissTargetPrediction", row_title_rot = 0, height = 3,row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"))
# ht2=Heatmap(t(hepg2_m_cast_top25),cluster_rows = F, col = f2, name= "z-scores iLINCS HepG2",row_names_gp = gpar(fontsize = 15, fontface = "bold", col=cols_hepg2),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HepG2",
#             heatmap_legend_param=list(legend_direction="vertical",grid_height = unit(0.8, "cm"),grid_width = unit(0.8, "cm"),labels_gp = gpar(fontsize = 15),title_gp = gpar(fontsize = 15, fontface = "bold")), height = 3)

ht2=Heatmap(t(conndata_stp),cluster_rows = F, cluster_columns = F, col = f2, name= "HEPG2 z-scores \niLINCS Connectivity (kd)",row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_ilincs_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HEPG2 iLINCS Connectivity (kd)", row_title_rot = 0, row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"), heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = "gray",height = 1.5)

ht3=Heatmap(t(conndata_cmap_kd_stp),cluster_rows = F, cluster_columns = F, col = f4, name= "HEPG2 CMap Connectivity (kd)",row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_kd_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HEPG2 CMap Connectivity (kd)", row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = "gray",height = 1.5)

ht4=Heatmap(t(conndata_cmap_oe_stp),cluster_rows = F, cluster_columns = F, col = f5, name= "HEPG2 CMap Connectivity (oe)",row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_oe_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HEPG2 CMap Connectivity (oe)", row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = "gray",height = 1.5)


#HA1E
setwd('../20210830_HA1E_HEPG2_iLINCS_cMAP_combinedGrades/')
#iLINCS connectivity data

##HA1E
setwd('../iLINCS_conn_kdoe/')
ptrn='_conn_HA1E_'
conndata=iLINCS_kd_connectivitydatafile_compiler(ptrn = ptrn)
conndata=conndata[,grepl(toupper(paste(c('GeneTarget',colnames(stp)), collapse = '|')), toupper(colnames(conndata)),ignore.case = T)]
miss_genes_ilincsConn=rownames(stp)[!rownames(stp)%in%conndata$GeneTarget]
rownames(conndata)=conndata$GeneTarget; conndata$GeneTarget=NULL
miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_ilincsConn),ncol = ncol(conndata))
rownames(miss_mat)=miss_genes_ilincsConn; colnames(miss_mat)=colnames(conndata)
conndata=rbind(conndata,miss_mat)
# for (i in 1:length(miss_genes_hepg2_ilincsConn)) {hepg2_conndata=rbind(hepg2_conndata,c(miss_genes_hepg2_ilincsConn[i],rep(as.numeric(NA),ncol(hepg2_conndata)-1)))}
conndata=conndata[order(rownames(conndata)),]
# conndata=conndata[,c(which(grepl('estradiol|mestranol',colnames(conndata),ignore.case = T)),which(grepl('prazol',colnames(conndata))),which(grepl('calcipotriol|tranilast',colnames(conndata),ignore.case = T)))]
conndata_stp=conndata[rownames(conndata)%in%rownames(stp),]

#CMap connectivity data
setwd('../GeneConnectivityScores/')
ptrn='conn_HA1E_'; types=c('kd','oe')
rm( list = Filter( exists, c('dataset') ) )

conndata_cmap=connectivitydatafile_compiler(ptrn = ptrn,types = types)
conndata_cmap=conndata_cmap[,grepl(paste(gsub('6_|24_','',colnames(stp)), collapse = '|'),colnames(conndata_cmap), ignore.case = T)]
conndata_cmap=conndata_cmap[,grepl(paste(colnames(stp), collapse = '|'), colnames(conndata_cmap),ignore.case = T)]
# conndata_cmap=conndata_cmap[,c(which(grepl('estradiol|mestranol',colnames(conndata_cmap),ignore.case = T)),which(grepl('prazol',colnames(conndata_cmap),ignore.case = T)),which(grepl('calcipotriol|tranilast',colnames(conndata_cmap),ignore.case = T)))]

#CMap kd 
conndata_cmap_kd=conndata_cmap[grepl('_kd', rownames(conndata_cmap)),]; rownames(conndata_cmap_kd)=gsub('_kd','',rownames(conndata_cmap_kd))
miss_genes_cMapConn_kd=rownames(stp)[!rownames(stp)%in%rownames(conndata_cmap_kd)]
miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_kd),ncol = ncol(conndata_cmap_kd))
rownames(miss_mat)=miss_genes_cMapConn_kd; colnames(miss_mat)=colnames(conndata_cmap_kd)
conndata_cmap_kd=rbind(conndata_cmap_kd,miss_mat)
# for (i in 1:length(miss_genes_hepg2_ilincsConn)) {hepg2_conndata=rbind(hepg2_conndata,c(miss_genes_hepg2_ilincsConn[i],rep(as.numeric(NA),ncol(hepg2_conndata)-1)))}
conndata_cmap_kd=conndata_cmap_kd[order(rownames(conndata_cmap_kd)),]
conndata_cmap_kd_stp=conndata_cmap_kd[rownames(conndata_cmap_kd)%in%rownames(stp),]
rm(miss_mat)
#CMap oe
conndata_cmap_oe=conndata_cmap[grepl('_oe', rownames(conndata_cmap)),]; rownames(conndata_cmap_oe)=gsub('_oe','',rownames(conndata_cmap_oe))
miss_genes_cMapConn_oe=rownames(stp)[!rownames(stp)%in%rownames(conndata_cmap_oe)]
miss_mat=matrix(data=as.numeric(NA),nrow = length(miss_genes_cMapConn_oe),ncol = ncol(conndata_cmap_oe))
rownames(miss_mat)=miss_genes_cMapConn_oe; colnames(miss_mat)=colnames(conndata_cmap_oe)
conndata_cmap_oe=rbind(conndata_cmap_oe,miss_mat)

conndata_cmap_oe=conndata_cmap_oe[order(rownames(conndata_cmap_oe)),]
conndata_cmap_oe_stp=conndata_cmap_oe[rownames(conndata_cmap_oe)%in%rownames(stp),]


# f11 = colorRamp2(seq(min(stp), max(stp), length=3), c("#3f7fff", "#EEEEEE", "red"))
f12 = colorRamp2(seq(max(conndata_stp,  na.rm = TRUE)*-1,max(conndata_stp,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))
# f13= colorRamp2(seq(max(zf_e2_top25$E2vsDMSO_logDE,  na.rm = TRUE)*-1,max(zf_e2_top25$E2vsDMSO_logDE,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))
f14= colorRamp2(seq(max(conndata_cmap_kd_stp,  na.rm = TRUE)*-1,max(conndata_cmap_kd_stp,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))
f15= colorRamp2(seq(max(conndata_cmap_oe_stp,  na.rm = TRUE)*-1,max(conndata_cmap_oe_stp,  na.rm = TRUE), length = 3), c("#3f7fff", "#EEEEEE", "red"))

# ha_mix_top = HeatmapAnnotation(heatmap_hepg2 = anno_density(t(rseqVal), type = 'heatmap'),heatmap = anno_density(t(rseqVal_ha1e), type = 'heatmap'))
cols<- rep('black', ncol(stp))
cols_ilincs_Con=rep('black', ncol(conndata_stp))
cols_cmap_kd_Con=rep('black', ncol(conndata_cmap_kd_stp))
cols_cmap_oe_Con=rep('black', ncol(conndata_cmap_oe_stp))
# cols_ha1e=rep('black', ncol(rseqVal_ha1e))
#turn red the specified rows in tf
proco=c("mestranol","Ethynyl Estradiol",'EthynylEstradiol','estradiol','17.betaestradiol','17-beta-estradiol','methylandr', 'methandr')
# prev=c('calcipotriol','tranilast')

cols[grepl(paste(proco, collapse = '|'),colnames(stp),ignore.case = T)]= 'red4'
# cols[grepl(paste(prev, collapse = '|'),colnames(stp),ignore.case = T)]= 'grey41'
cols_ilincs_Con[grepl(paste(proco, collapse = '|'),colnames(conndata_stp),ignore.case = T)]= 'red4'
# cols_ilincs_Con[grepl(paste(prev, collapse = '|'),colnames(conndata_stp),ignore.case = T)]= 'grey41'
cols_cmap_kd_Con[grepl(paste(proco, collapse = '|'),colnames(conndata_cmap_kd_stp),ignore.case = T)]= 'red4'
# cols_cmap_kd_Con[grepl(paste(prev, collapse = '|'),colnames(conndata_cmap_kd_stp),ignore.case = T)]= 'grey41'
cols_cmap_oe_Con[grepl(paste(proco, collapse = '|'),colnames(conndata_cmap_oe_stp),ignore.case = T)]= 'red4'
# cols_cmap_oe_Con[grepl(paste(prev, collapse = '|'),colnames(conndata_cmap_oe_stp),ignore.case = T)]= 'grey41'

setwd('../20210830_HA1E_HEPG2_iLINCS_cMAP_combinedGrades/')
jpeg("20220330_SwissTarPred_HEPG2andHA1E_iLINCS_cMAP_combinedGrades_forMCE_UNCposter.jpeg", units="in", width=20, height=14, res=600)
# ht11=Heatmap(t(m_cast_stp),cluster_rows = F, clustering_method_columns = 'average', col = f1, name= "z-scores HA1E gene expression",top_annotation = ha,row_names_gp = gpar(fontsize = 15, fontface = "bold", col=cols),column_names_gp = gpar(fontsize = 15, fontface = "bold"),
#             heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")),cluster_columns=TRUE,row_title = "HA1E gene expression", row_title_rot = 0, height = 3,row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"))
# ht2=Heatmap(t(hepg2_m_cast_top25),cluster_rows = F, col = f2, name= "z-scores iLINCS HepG2",row_names_gp = gpar(fontsize = 15, fontface = "bold", col=cols_hepg2),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HepG2",
#             heatmap_legend_param=list(legend_direction="vertical",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 15),title_gp = gpar(fontsize = 15, fontface = "bold")), height = 3)

ht12=Heatmap(t(conndata_stp),cluster_rows = F, cluster_columns = F,col = f2, name= "HA1E z-scores \niLINCS Connectivity (kd)",row_names_gp = gpar(fontsize = 8, fontface = "bold", col=cols_ilincs_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HA1E iLINCS Connectivity (kd)", row_title_rot = 0,row_title_gp = gpar(fontsize=20,fontface='bold', fill = "red", col = "white", border = "blue"), heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = "gray",height = 2)

ht13=Heatmap(t(conndata_cmap_kd_stp),cluster_rows = F, cluster_columns = F, col = f4, name= "HA1E CMap Connectivity (kd)",row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_kd_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HA1E CMap Connectivity (kd)", row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = "gray",height = 1.5)

ht14=Heatmap(t(conndata_cmap_oe_stp),cluster_rows = F, cluster_columns = F, col = f5, name= "HA1E CMap Connectivity (oe)",row_names_gp = gpar(fontsize = 10, fontface = "bold", col=cols_cmap_oe_Con),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HA1E CMap Connectivity (oe)", row_title_rot = 0, heatmap_legend_param=list(legend_direction="horizontal",grid_height = unit(0.3, "cm"),grid_width = unit(0.3, "cm"),labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")), na_col = "gray",height = 1.5)
# ht5=Heatmap(t(rseqVal_ha1e),clustering_method_rows = "average", col = f4, name= "z-scores iLINCS HA1E",row_names_gp = gpar(fontsize = 15, fontface = "bold", col=cols_ha1e),column_names_gp = gpar(fontsize = 15, fontface = "bold"),row_title = "HA1E",
#             heatmap_legend_param=list(legend_direction="vertical",grid_height = unit(0.8, "cm"),grid_width = unit(0.8, "cm"),labels_gp = gpar(fontsize = 15),title_gp = gpar(fontsize = 15, fontface = "bold")), height = 1.5)


ht_list=ht1 %v% ht2  %v% ht3  %v% ht4  %v% ht12  %v% ht13  %v% ht14
draw(ht_list, heatmap_legend_side = "bottom", main_heatmap="SwissTarget similarity scores",annotation_legend_list = lgd_list,annotation_legend_side = "bottom",merge_legend = TRUE, column_title = "HEPG2 & HA1E - SwissTarget and Connectivity data compilation", column_title_gp = gpar(fontsize = 16, fontface='bold'), ht_gap = unit(c(5, 5, 5, 7, 5, 5), "mm"))
dev.off()


#MolTarPred
k=readRDS('20211006_MolTarPred_matix.RDS')


lgd_list=Legend(labels = c("Pro-coagulants", "Anti-coagulants"), title = "Compound types", type = "points", pch = 18, legend_gp = gpar(col = c("red4","black"), fontsize=24))
f1 = colorRamp2(seq(max(dataset), min(dataset), length=5), c("dark red","red3","red2","red","white"))

jpeg("20211006_MolTarPred_AntiProCoagulants_PPIs_MePS_cutoff2.jpeg", units="in", width=18, height=12, res=300)

ht=Heatmap(dataset,clustering_method_rows = "ward.D",clustering_method_columns  = "ward.D", col = f1, name = "MolTarPred - Reliability of Prediction", row_names_gp = gpar(fontsize = 7.5, fontface = "bold"),column_names_gp = gpar(fontsize = 18, fontface = "bold",col=cols),
           heatmap_legend_param=list(color_bar="continuous",legend_direction="horizontal",grid_height = unit(0.8, "cm"),legend_width = unit(10, "cm"),labels_gp = gpar(fontsize = 18),title_gp = gpar(fontsize = 18, fontface = "bold"),title_position = "topcenter"))
draw(ht, heatmap_legend_side = "top",  column_title_gp=gpar(fontface="bold", line = 5),annotation_legend_list = lgd_list,annotation_legend_side = "bottom")
dev.off()


