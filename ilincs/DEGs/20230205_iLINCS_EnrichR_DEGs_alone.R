### EnrichR analyses for iLINCS gene expression data
#The approach allows: 
# (i)compiling gene expression data (.xls) retrieved from iLINCS, 
# (ii) employ EnrichR based gene enrichment analyses on differentially expressed genes (DEGs), 
# (iii) retrieve pathways associated with these DEGs in either direction (complete: for complete set of DEGs, up: for the upregulated DEGs only, and down: for the downregulated genes only)
# (iv) compare and contrast different exposures based on the enriched pathways and their direction changes
# (v) obtain heatmaps, correlation plots and upset plots for visual inspections
# Feb 11, 2023 - Murat Yaman

rm(list = ls())

#Please be careful with the pattern and how you keep the names of your original files that you previously obtained from iLINCS gene expression data. Preferably, keep the names in a similar format as they are in the names of the sample files eg. iLINCS_complete_HA1E_24_warfarin_004uM.xls (iLINCS: data source, complete: complete list of genes, HA1E: name of the cell line, 24: 24h of exposure, warfarin: compound of interest, 004uM: 0.04 uM exposure concentration - becareful with using the underscores, too) - if the compound or your file name contains special characters like []{}#.,/ the compiler function converts them to underscores.
pattern="iLINCS_complete_HA1E_"
dataset=iLINCS_DEGdata_compiler(ptrn = pattern)

#in case you may need the list of genes with their p-values or only differentially expressed genes with their fold changes, you can follow the below comments
dataset_pvals= dataset[,grepl('pval', colnames(dataset))]
dataset_logDE=dataset[which(apply(dataset_pvals, 1, function(x) any(x<=0.05))==T),grepl('logDE', colnames(dataset))]
dataset_logDE=as.matrix.data.frame(dataset_logDE)

##EnrichR analyses
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

kegg="KEGG_2019_Human"
go_bio="GO_Biological_Process_2018"

alt_dbs = list(kegg=kegg,go_bio=go_bio)
#please visit enrichR manual for available enrichment tools. KEGG will be used as an example from here onwards

drug_list=unique(sub('_logDE|_pval','',colnames(dataset)))

rm( list = Filter( exists, c('go_output','k','d','pcutoff','dbs','result') ) )

##If you like you can refine our list to a couple of compounds only to allow better sample-wise contrasts and retrieve a list of a panel of pathways from differentially expressed genes
##Or, if you have several compounds only, you can run them along each other
drug_list2=drug_list[grepl(x = drug_list, pattern = paste0(c('menadione','warfarin','dabigatran'), collapse = '|'),ignore.case = T)]

keggtr=enrichr_output(drug_list = drug_list2 ,dataset = dataset,pcutoff = 0.05,dbs = alt_dbs[1])

##Let's obtain the list of the pathways with combination scores above a threshold for each group and direction, besides heatmaps and cormaps for visual inspection. Note: results are provided in a log2 scale for the respective combination score values
kegg_mestwar2=slicenricher(enrichrlist = keggtr,combscore_cutoff = 100,refcompounds = drug_list2)

##Let's see how our treatments compare to one another, based on pathway enrichments - but in opposite directions. This can give us some clues about drugs with opposing effects on the phenotype. Users will be prompted to provide some details about the analyses. If you like to make it faster, you can assign these values in the arguments. Note: you can feed vectors of compounds for contr1 and contr2 to compare groups of compounds eg contr1=c('warfarin', 'phenindione'), contr2=c('dabigatran', 'apixaban')
updownfunk=slicenricher_updown(enrichrlist = keggtr,combscore_cutoff = 100)
# updownfunk=slicenricher_updown(enrichrlist = keggtr,combscore_cutoff = 100,contr1 = 'menadione', contr2 = c('warfarin', 'dabigatran'), classContr1 = 'Pro-coagulants',classContr2 = 'Anti-coagulants',enrichmentdb = 'KEGG',repositoryName = 'iLINCS')

                                  
#Prepare data for upsetplot representations to see overlaps and contrasts between compounds' pathway enrichment profiles. Upsetplots allow nice representations when comparing more than 4 samples. It can still get overcrowded though. Hence, it would still be a good practice to keep the numbers of samples minimal (5-8)
KEGG_upsetlist=upsetlister(slicenrichrlist = kegg_mestwar2, combscore_cutoff = 100,refcompounds = drug_list2)

library(Vennerable)
KEGG_temp=Venn(KEGG_upsetlist)  #provide all your groups as list

#Retrieve all the points that contain enriched pathways 
KEGG_querylist=querylister(intersectionSets = KEGG_temp@IntersectionSets, upsetlist = KEGG_upsetlist)

#Make the plot and see how the profiles compare. List of pathways and DEGs for each comparison should be saved in an xlsx document in your directory along with the upsetplot. You can also provide a color palette of your choice. Please see, RColorBrewer for available palettes
kegg_upsetlegend=upsetlistlegend(upsetlist = KEGG_upsetlist,querylist = KEGG_querylist,enrichrlist = keggtr,enrichmentdb = 'KEGG',cell_line = 'HA1E', refcompounds = drug_list2,enrichmentTOOL = 'EnrichR',combscore_cutoff = 100)
