#CMap_iLINCS

Evaluating gene-compound relationships across multiple platforms can be hard. Here, we aimed to incorporate data from public sources together and facilitate analysis of your experimental data along with them.

The workflow includes analyses of data derived from ConnectivityMap, iLINCS and SwissTargetPrediction. From these platforms, (i) connectivity data with knock-down/over-expression/compound/compound class signatures, (ii) gene expression profiles, (iii) target predictions for the compounds of interests, were retrieved and analyzed. On the main folder, you can access major functions that are used throughout the workflow (see CMap_iLINCS_SwissTarget_functions.R). Inside each folder, you can also access data from the relevant databases and try out the R scripts to see how signatures from different compounds align or contrast with one another. This can become quite a useful approach if you have classes of compounds that have opposing actions. Pro-coagulant and anti-coagulant compounds have been tested for this purpose. Accordingly, heatmaps, correlation plots and upset plots have been frequently utilized on the gene level data as well as on the gene enrichment profiles. Hence, the users can see which genes or pathways move in different directions when comparing two groups of compounds.

In addition to that, in this main folder you can access 'iLINCS_cMAP_SwissTargetPredictions_combinedHeatmaps.R' script allowing you to combine data from the aforementioned platforms and test how your own experimental findings (like RNAseq experiments, or gene-phenotype association scores) would align with them. This is a sort of a meta-analysis approach. Accordingly, one can get insights from multiple sources at the same time, see which candidate genes or drugs have some potential to investigate further, make data-driven decisions, develop new hypotheses and test them further. Users should also be cautious that the figures and representations do not provide definitive answers, but just give you some angles from multiple studies to look at and derive new ideas to test for.

Initial pipeline is finalized by using local files. Potential follow ups:
1) Implementing the codes as R Markdown documents/Shiny apps for ease of demonstration and conveying the approaches better
2) Allowing html/curl queries across the platforms, so that users can simultaneously choose which data they want access to.
