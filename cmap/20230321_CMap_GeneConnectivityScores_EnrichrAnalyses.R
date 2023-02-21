### EnrichR analyses for CMap drug-gene knockdown/overexpression connectivity
#The approach allows: 
# (i) compiling connectivity data (.txt) retrieved from CLUE Connections (clue.io /conn compounds of interest followed-up with accessing 'See all connections' button for ConnectivityMap data and representations. 
# (ii) employ EnrichR based gene enrichment analyses on genes that have high absolute connectivity scores with the compounds of interest based on knockdows.overexpression studies, 
# (iii) retrieve pathways associated with these genes moving in either direction (complete: for complete set of genes with high absolute connectivity scores, up: for the genes that positively correlate with the connectivity score 
## (if knock-down study; genes with negative connectivity score (<(-98), or if overexpression; genes with positive connectivity score above the threshold), 
## and down: for the genes that negatively correlate with the compounds connectivity signatures (if knock-down study; genes with positive connectivity score above the threshold, or if overexpression; genes with negative connectivity score below the threshold (<(-97)), 
# (iv) compare and contrast different exposures based on the enriched pathways and their direction changes
# (v) obtain heatmaps, correlation plots and upset plots for visual inspections
# Feb 21, 2023 - Murat Yaman
rm(list = ls())
suppressPackageStartupMessages(source('../20230221_mylib_functions.R'))
setwd('../cmap/')

ptrn='_conn_HA1E_'; types=c('kd','oe','cp','cc')
rm( list = Filter( exists, c('dataset') ) )

##Let's compile CMap connectivity .txt data together. Main output file will contain knockdown/overexpression based connectivity data and will be used for further gene enrichment analyses. In the meantime, you can also check out additional functionalities of this function below, e.g. to produce heatmaps/correlation plots for the connectivity data and explore compound/compound class based similarity scores
##ptrn arguement is by default '_conn_HA1E_'. You can apply similar naming convention for your own files. types arguement considers kd(knockdown) and oe(overexpression) signatures by default. You can test additional input types such as cp (compound) or cc (compound class) 
ha1e_conndata=connectivitydatafile_compiler(ptrn = ptrn,types = types, cmaprank_cutoff = 98,cor_outputs = T, heatmap_outputs = T)

# You can take quick look to see how knockdown/overexpression of a list of genes may relate with exposures to the compounds of interest. Some gene names may partially match with yours. You can try dplyr's %in% pipe and use toupper() or tolower() transformations to allow exact matches
candygenes=c('prcp', 'ptpn','bnip3','vabp','isl', 'znf366','tfap2a','gale','ptgs2')
drugs=c('menadione', 'warfarin', 'phylloquinone')
ha1e_conndata[grepl(paste(candygenes, collapse = '|'),ha1e_conndata$Name_type,ignore.case = T),grepl(paste(drugs, collapse = '|'), colnames(ha1e_conndata), ignore.case = T)]
drug_list=colnames(ha1e_conndata)[!colnames(ha1e_conndata)%in%'Name_type']

##EnrichR analyses
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive = TRUE
dbs = listEnrichrDbs()
if (is.null(dbs)) websiteLive = FALSE
if (websiteLive) head(dbs)

kegg="KEGG_2019_Human"
go_bio="GO_Biological_Process_2018"

alt_dbs = list(kegg=kegg,go_bio=go_bio)
#please visit enrichR manual for available enrichment tools. KEGG will be used as an example from here onwards

rm( list = Filter( exists, c('go_output','k','d','pcutoff','dbs','result') ) )

##If you like you can refine your list to a couple of compounds only to allow better sample-wise contrasts and retrieve a list of a panel of pathways from the candidate genes with absolute connectivity scores above a threshold
##Or, if you have several compounds only, you can run them along each other. Please, do not forget change the cmapranks arguement to TRUE to analyze this type of data
kegg_data_summ=enrichr_output(drug_list = drug_list ,dataset = ha1e_conndata,dbs = alt_dbs[1],cmapranks = T,cmaprank_cutoff = 98)

##In case you want to investigate only a limited numbers of compounds, you can refine your list and use it for the follow up
# drug_list2=drug_list[grepl(paste(c('menadione', 'warfarin', 'phylloquinone'), collapse = '|'), drug_list,ignore.case = T)]

##Let's obtain the list of the pathways with combination scores above a threshold for each group and direction, besides heatmaps and cormaps for visual inspection. Note: results are provided in a log2 scale for the respective combination score values
kegg_mestwar_summ=slicenricher(enrichrlist = kegg_data_summ,combscore_cutoff = 70,cor_outputs = T,heatmap_outputs = T,cmapranks = T,refcompounds = drug_list,enrichmentdb = 'KEGG')

##Let's see how our treatments compare to one another, based on pathway enrichments - but in opposite directions. This can give us some clues about drugs with opposing effects on the phenotype. Users will be prompted to provide some details about the analyses. If you like to make it faster, you can assign these values in the arguments. Note: you can feed vectors of compounds for contr1 and contr2 to compare groups of compounds eg contr1=c('menadione', 'phylloquinone'), contr2=c('warfarin', 'skatole', 'argatroban')
updownfunk=slicenricher_updown(enrichrlist = kegg_data_summ,combscore_cutoff = 70)
updownfunk=slicenricher_updown(enrichrlist = kegg_data_summ,combscore_cutoff = 70,contr1 = c('menadione', 'phylloquinone'), contr2 = c('warfarin', 'skatole', 'argatroban'), classContr1 = 'Pro-coagulants',classContr2 = 'Anti-coagulants',enrichmentdb = 'KEGG',repositoryName = 'CMap', cmapranks = T)

#Prepare data for upsetplot representations to see overlaps and contrasts between compounds' pathway enrichment profiles. Upsetplots allow nice representations when comparing more than 4 samples. It can still get overcrowded though. Hence, it would still be a good practice to keep the numbers of samples minimal (5-8)
KEGG_upsetlist=upsetlister(slicenrichrlist = kegg_mestwar_summ, combscore_cutoff = 70,refcompounds = drug_list)

library(Vennerable)
KEGG_temp=Venn(KEGG_upsetlist)  #provide all your groups as list

#Retrieve all the points that contain enriched pathways 
KEGG_querylist=querylister(intersectionSets = KEGG_temp@IntersectionSets, upsetlist = KEGG_upsetlist)

#Make the plot and see how the profiles compare. List of pathways and DEGs for each comparison should be saved in an xlsx document in your directory along with the upsetplot. Order of the pathways in the upsetplot from top to buttom are the same in the excel sheet. You can also provide a color palette of your choice. Please see, RColorBrewer for available palettes
kegg_upsetlegend=upsetlistlegend(upsetlist = KEGG_upsetlist,querylist = KEGG_querylist,enrichrlist = kegg_data_summ,enrichmentdb = 'KEGG',cell_line = 'HA1E', refcompounds = drug_list,enrichmentTOOL = 'EnrichR',combscore_cutoff = 70)

#Upon inspecting the upsetplot for the demo data here, Prolactin signaling pathway with changes in ESR1 gene might be an interesting one that moves different directions between menadione and skatole. Estrogen based therapies are known to relate with hypercoagulability/deep vein thrombosis 


