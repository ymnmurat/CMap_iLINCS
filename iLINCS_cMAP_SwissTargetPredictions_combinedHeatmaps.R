##Combining heatmap graphs from multiple sources (ConnectivityMap, iLINCS and SwissTargetPrediction)
##You can check the folders named with 'cmap', 'ilincs', 'DEGs' and 'swisstarget' to access their data, R scripts and see how their heatmaps from original data and pathway enriched data would look like. We have not demonstrated SwissTargetPrediction based data heatmap previously. You can also test it here.
##Main objective is to combine data from multiple sources into a single heatmap figure to facilitate data-driven decision processes to follow-up with. Accordingly, we believe that multiple sources of information can complement one another, such as possible target annotations from SwissTargetPrediction and knock-down/over-expression signatures from ConnectivityMap and iLINCS, as well as user queried data (gene expression fold changes, or gene-phenotype correlation scores)
#The users should also approach these representations cautiously, as their respective data are derived from different sources (separate databases, cell lines, exposure times and even the algorithms that their respective data was pre-exposed to). Therefore, the data presented here should not be taken as a definitive answer, but as plausible hints for the follow-up studies. If possible additional literature on these hints from alternative sources would further complement your approach.
# Feb 21, 2023 - Murat Yaman

rm(list=ls())
gc()

suppressPackageStartupMessages(source('CMap_iLINCS_SwissTarget_functions.R'))

#Let's test our swisstargetdatafile_compiler() function. It mainly: 
# (i) gathers the SwissTargetPrediction derived data (.xlsx) together, 
# (ii) clips the first row due to SwissTargetPrediction derived data conventions, unless modified by the user
# (iii) applies a minimum vakue threshold to refine the data, and further finetunes for a specific group of compounds (if specified)
# (iv) generates a heatmap of these comparisons
setwd('swisstarget/')
swisstarget=swisstargetdatafile_compiler(ptrn ='SwissTargetPrediction_', targetscore_cutoff=0.2, heatmap_outputs = T)

# We can now continue to evaluate our function (multi_heatmaps_SwissCMapILINCS()) to genrate and gather heatmap figures together. For that, we will make use of user queried data (let's say gene expression data) as the heatmap annotation, SwissTargetPrediction heatmap as the primary output, and accordingly subsidiary heatmaps from gene knock-down/over-expression studies in ConnectivityMap Connections and iLINCS platforms. For the latter, we will utilize data from HA1E and HEPG2 cell lines. If needed, you can provide alternative lines. Please, just make sure of the naming conventions and the data formats.


## (1) Input data from zebrafish RNAseq. 
## We will create a dummy expression data here. In addition, you will need a list of ortholog and paralog genes for zebrafish and human, or your species of interest vise versa. Using human as the target species is strongly recommended in this pipeline for ortholog/paralog assessments. This is mainly because the heatmaps and the functions used here are more oriented towards the SwissTargetPrediction outputs and human gene names associated with them.
## Here BioMaRT package is utilized to get a fine grade ortholog/paralog data. Strangely, getLDS() function of this package is recently running into some issues, hence our custom function is currently referring to one of the older versions (Dec 2021 release). Note: It may take a while to get the data from BioMaRT.
## Alternatively, you can get this list by visiting Ensemble's BioMaRT services, or install a specific database distribution for your model organism (for zebrafish "org.Dr.eg.db"), or make use of alternative packages for ease of demonstration.
## orthogene library is one of them, and you can follow the code below, to retrieve ortholog genes between zebrafish and human

# library(orthogene)
# orth_zeb=orthogene::report_orthologs(target_species = "zebrafish",
#                                         reference_species = "human",
#                                         method_all_genes = 'gprofiler',
#                                         method_convert_orthologs = 'gprofiler')
# orth_zeb=orth_zeb$map[c('Gene.Symbol', 'ortholog_gene')]
# colnames(orth_zeb)=c('Gene.name', 'Gene.name.1')
# orth_zeb=orth_zeb[!duplicated(orth_zeb),]

## Mainly, it would be a better practice to get your data first, check your top genes, double-check their orthologs/paralogs to make informed decisions and follow up with these heatmap representations
library(biomaRt); library(tidyverse)
orth_zeb=biomart_orthologs()
zfs=unique(orth_zeb$Gene.name); zfs=as.data.frame(zfs); colnames(zfs)='Gene.name'; 

## Getting the top most interesting genes. Here, we generated a normally distributed gene expression data. Arcsin or log transformation on the gene counts or fold changes (vs controls) can provide you a similar kind of distribution.
set.seed(1994)
zfs$logFC=rnorm(nrow(zfs), mean=0, sd = .5)

orth_zeb=merge(orth_zeb,zfs, by='Gene.name')
##Since orthologs and paralogs can have one-to-many/many-to-one/many-to-many type of relationships, at some point you may need to summarize multiple occurrences of the same gene names across your tables. For that, you can group them by for a variable and take their average/geometric mean/median
orth_zeb$Gene.name=NULL; orth_zeb=orth_zeb%>%group_by(Gene.name.1)%>% summarise_all(list(~mean(., na.rm=T)))%>%as.data.frame()
orth_zeb=orth_zeb[order(orth_zeb$Gene.name.1),]
rownames(orth_zeb)=orth_zeb$Gene.name.1; 

##We can make a copy from our original data for alternative uses. Here, we can get the gene fold change data from one study alone and from two separate studies (orth_zeb1 and orth_zeb2, respectively). Please, be careful with keeping your data in a data frame format. This will be crucial for annotating your heatmaps in the follow-up step
orth_zeb1=orth_zeb;orth_zeb1$Gene.name.1=NULL; orth_zeb1=as.data.frame(orth_zeb1) 

orth_zeb2=orth_zeb; orth_zeb2$logFC2= rnorm(nrow(orth_zeb2), mean=0, sd = .5); orth_zeb2$Gene.name.1=NULL; orth_zeb2=as.data.frame(orth_zeb2) 

##Let's get a list of our candidate genes. Here, we collected the names of the top 20 genes with highest fold change differences. In the second list, 10 top genes from one group and 10 from the other one. Please, be careful here we referred absolute log fold differences. 
candygenes1=head(rownames(orth_zeb1)[order(abs(orth_zeb1[,1]),decreasing = T)],20)
candygenes2=c(head(rownames(orth_zeb2)[order(abs(orth_zeb2[,1]),decreasing = T)],10),head(rownames(orth_zeb2)[order(abs(orth_zeb2[,2]),decreasing = T)],10))

##Since our data was randomly generated previously, you may want to provide a list of genes that you are most interested in. Between 10-20 genes is good for representation and visual inspections. You can provide more, but that would also overcrowd your figures. If you do have a huge list to check for, you may consider ranking your genes and making batches of the gene lists, then checking each batch individually.  please also be extra careful about the naming convention/gene aliases for your candidate genes and their orthologs. Interestingly, instead of differential expression levels, genes that would correlate most with the observed phenotype and their correlation scores can be used here, too. Here, some made up gene names are also included to show that the heatmaps can also represent genes that have no matches with our reference. 
testcandygenes = c('ptgs2','nr3c1', 'kdr','flt4','esr1','drd1','cyp19a1','cyp39a1','lss','tsr1','kdrl','ar','pgr','f13','f2','kardo','PamPa','babayaro')
# candygenes=testcandygenes

## (2) Genenrating and combining heatmaps together. multi_heatmaps_SwissCMapILINCS is a custom made function that compiles data SwissTargetPredcition, iLINCS and CMap derived data. Please, use similar file name conventions and directory allocations (see folders; swisstarget, ilincs, cmap under github/cmap_ilincs directory). CellLine is by default HA1E, you can additionally include secondCellLine argument. Optionally, you can provide a refdata to annotate your heatmaps. Final heatmaps are saved in main github/cmap_ilincs directory automatically. Here you can also change the name for your single column reference top annotation ('Measurements' by default) on the heatmap, and change the colors of the NA data (default 'gray') and width/height/resolution of your heatmap figure

a=multi_heatmaps_SwissCMapILINCS(candygenes = candygenes1,CellLine = 'HA1E',secondCellLine = 'HEPG2',refdata = orth_zeb1)
b=multi_heatmaps_SwissCMapILINCS(candygenes = candygenes2,CellLine = 'HA1E',secondCellLine = 'HEPG2',refdata = orth_zeb2)
c=multi_heatmaps_SwissCMapILINCS(candygenes = testcandygenes,CellLine = 'HA1E',secondCellLine = 'HEPG2',refdata = orth_zeb2)
d=multi_heatmaps_SwissCMapILINCS(candygenes = testcandygenes,CellLine = 'HEPG2',refdata = orth_zeb1, singleAnnotationTag = 'CompoundX RNAseq')
e=multi_heatmaps_SwissCMapILINCS(candygenes = testcandygenes,CellLine = 'HEPG2',refdata = orth_zeb2, nacols = 'orange',width_heat = 20 ,height_heat = 18,res_heat = 600)
f=multi_heatmaps_SwissCMapILINCS(candygenes = testcandygenes,CellLine = 'HA1E',secondCellLine = 'HEPG2', nacols = 'orange',width_heat = 20 ,height_heat = 18,res_heat = 600)
