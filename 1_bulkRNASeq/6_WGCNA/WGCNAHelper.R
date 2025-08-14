

### 1. === FUNCTION: MODULE MEMBERSHIP BOXPLOT


boxMM <- function(GeneMM=GeneMM, module=module, title=NULL) {
  
  # Arguments:
  # I. data: data frame specifying for each gene the module membership for each module (as MM+modulename, e.g. MMblue), as well as a 'Module' column with module assignment. 
  # Additional columns can be present.
  # II. module: string indicating the module to be represented
  # III. title: plot title.If null, a title with the specified module is generated
  # Dependencies: ggplot2
  
  # 1. Setting of plot title
  if (is.null(title)) {
    title <- paste0('Module Membership for ', toupper(module), ' Module')
  }  
  
  # 2. Boxplot of module membership stratified by module assignment
  ggplot(GeneMM, aes(x=Module, y=!!sym(paste0('MM', module)))) + 
    geom_boxplot(fill=module, alpha=0.5) + 
    ggtitle(title) + 
    theme_bw() +
    theme(plot.title = element_text(face='bold', colour='darkred', size=13, hjust=0.5), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}


### 2. === FUNCTION: MODULE EIGENGENE RIBBONPLOT

ribbonME <- function(METraits=METraits, module=module, ascisse='Day', title=NULL, Shape=18, DotSize=6.5) {
  
  # Arguments:
  # I. METraits: data frame specifying the module eigengene and the continuous trait of interest (Day); a Rep column must specify the sample grouping. Additional columns can be present.
  # II. Module: string indicating the module to be represented
  # III. title: plot title.If null, a title with the specified module is generated
  # Dependencies: ggplot2
  
  # 1. Setting of plot title
  if (is.null(title)) {
    title <- paste0('Module Eigengene for ', toupper(module), ' Module')
  }  
  
  # 2. Ribbonplot representing module eigengene behaviour through day
  ggplot(data=METraits, aes(x=!!sym(ascisse), y=!!sym(paste0('ME', module)), group=Rep, col=Rep)) +
    geom_point(stat='identity', size=DotSize, shape=Shape) +
    geom_line(size=3.5, alpha=0.55) +
    ggtitle(title) + 
    #geom_text(aes(label=Sample), size=2) + 
    scale_color_hue(h=c(120, 360)) + 
    theme_light()  +
    theme(plot.title = element_text(face='bold', colour='darkred', size=13, hjust=0.5)) 
}

# -----------



### 3. === FUNCTION: MODULE EIGENGENE BOXPLOT


boxME <- function(METraits=METraits, module=module, ascisse='Type', title=NULL, 
                  lines='No', reference=NA) {
  
  # Arguments:
  # I. METraits: data frame specifying the module eigengene and the continuous trait of interest (Day); a Rep column must specify the sample grouping. Additional columns can be present.
  # II. Module: string indicating the module to be represented
  # III. title: plot title.If null, a title with the specified module is generated
  # IV. lines: if lines should be drawn reporting the maximum and minimum level of a reference condition
  # V. reference: reference condition for the lines
  
  # Dependencies: ggplot2
  
  # 1. Setting of plot title
  if (is.null(title)) {
    title <- paste0('Module Eigengene for ', toupper(module), ' Module')
  }  
  
  if (lines == 'Yes'){
    ColSel= METraits[, names(METraits)==ascisse]
    Upline = max(METraits[ColSel==reference, paste0('ME', module)])
    Lowline = min(METraits[ColSel==reference, paste0('ME', module)])
  }
  
  
  # 2. Boxplot representing module eigengene behaviour through day
  Box <- ggplot(data=METraits, aes(x=!!sym(ascisse), y=!!sym(paste0('ME', module)), col=!!sym(ascisse))) +
    geom_boxplot(col=module, alpha=0.75) +
    geom_point(stat='identity', col=module, size=4)+ # size=DotSize, shape=Shape
    #geom_line(size=3.5, alpha=0.55) +
    ggtitle(title) + 
    #geom_text(aes(label=Sample), size=2) + 
    theme_light()  +
    theme(plot.title = element_text(face='bold', colour='darkred', size=13, hjust=0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  if(lines=='Yes'){
    Box <- Box +
      geom_hline(yintercept = Upline, linewidth=3, alpha=0.45) +
      geom_hline(yintercept = Lowline, linewidth=3, alpha=0.45)
  }
  
  return(Box)
  
}



# -----------


### 4. === FUNCTION: MEMBERSHIP-CONNECTIVITY SCATTERPLOT


plotGeneMetrics <- function(GeneMetrics=GeneMetrics, module=module, title=NULL, top=NULL) {
  
  # Arguments: 
  
  # I. GeneMetrics: dataframe containing gene annotation, module assignment and metrics for each gene. 
  #    Columns that are used: gene_biotype, external_gene_name, kWithin, module membership, Module
  # II. module: string specifying the module to be examined
  # III. title: string that will be incorporated in plot title
  # VI. top: integer specifying the maximum number of genes to be visualized, ranked according to kWithin
  
  # Dependencies: pheatmap, viridis
  
  
  # 1. Data selection
  # First a selection is implemented to identify only the genes of the specified module
  Metrics <- dplyr::filter(GeneMetrics, Module==module)
  # Then, if top is not null, only the top-specified genes are selected
  if (!is.null(top)) {
    Metrics <- dplyr::top_n(Metrics, top, kWithin)
  }
  
  # 2. Setting of plot title
  if (is.null(title)) {
    title <- paste('Membership-Connectivity Scatterplot for', dim(Metrics)[1],  toupper(module), 'Genes')
  }  
  
  
  # 3. Plotting 
  Scatter <- ggplot(Metrics, aes(x=log10(!!sym(paste0('MM', module))), y=log10(kWithin), data_id=Gene, color=-log10(!!sym(paste0('p_MM', module))))) +
    geom_point(size=3.5, shape=19) + 
    scale_color_continuous(low='blue', high='springgreen') +
    xlab('Log10 Module Membership') + ylab('Log10 Intramodular Connectivity') +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(face='bold', colour='darkred', size=20, hjust=0.5), 
          legend.text=element_text(size=13), legend.position='bottom', axis.title=element_text(size=16)) 
  
  Interactive <- plotly::ggplotly(Scatter)
  return(Interactive)
}


### 5. === FUNCTION InteractiveMetrics: generation of interactive table for main gene metrics

interactiveMetrics <- function(GeneMetrics, Module, searchURL='https://www.genecards.org/cgi-bin/carddisp.pl?gene='){
  
  # Arguments: 
  
  # I. GeneMetrics: dataframe containing gene annotation, module assignment and metrics for each gene. 
  # II. module: string specifying the module for which the table should be generated
  # III. searchURL: for the interactive link. Default to gene cards
  
  # Dependencies: dplyr, DT 
  
  GeneMetrics %>% 
    dplyr::filter(Module==ModName) %>%
    dplyr::top_n(TopGenes, kWithin) %>%
    # generation of the link
    dplyr::mutate(GeneLink=paste0('<a href="', searchURL, Gene, '">', Gene, '</a>')) %>% 
    # selection of columns to be shown
    dplyr::select(GeneLink, Module, Gene, GeneBiotype, kWithin, !!sym(paste0('MM', ModName)), !!sym(paste0('p_MM', ModName))) %>% 
    DT::datatable(class='hover', rownames=FALSE, caption=paste(ModName, 'genes'), filter='top', options=list(pageLength=10, autoWidth=TRUE), escape=FALSE) %>%
    DT::formatRound(5, 1) %>% DT::formatRound(6, 2) %>% DT::formatSignif(7, 2)
}



### 6. === FUNCTION topGOGeneVectors: generation of gene vectors for topGO analysis

#  Function modified from DEGCharacterization in order to be suited to take as input GeneMetrics data frame.


topGOGeneVectors <- function(GeneMetrics, modules){
  
  # Arguments: 
  
  # I. GeneMetrics: dataframe containing gene annotation, module assignment and metrics for each gene. 
  #    Columns that are used: gene_biotype, external_gene_name, kWithin, module membership, Module
  # II. modules: vector of strings specifying the modules for which the vector of genes should be generated
  
  # Dependencies: dplyr
  
  GeneVectors <- list()
  
  # 1. Generation of vectors for all tested genes, all modulated genes, up-reg and down-reg
  UniverseGenes <- GeneMetrics %>% dplyr::pull(Gene)
  
  # 2. For each gene module, generation of named vectors in the format required by TopGO
  for(i in 1:length(modules)){
    ModuleGenes <- dplyr::filter(GeneMetrics, Module==modules[i]) %>% dplyr::pull(Gene)
    GeneVectors[[i]] <- factor(as.integer(UniverseGenes%in%ModuleGenes))
    names(GeneVectors[[i]]) <- UniverseGenes
    names(GeneVectors)[i] <- modules[i]
  }
  
  return(GeneVectors)
}



### 7.  === FUNCTION topGOResults: topGO enrichment analysis

#  Function from DEGCharacterizationHSapiensGenecode. It generates the TopGO object, perform the TopGO enrichment analysis for a single statistical test and algorithm and return the two result table: ResAll reports the complete results: ResSel reports the terms after selection on the enrichment thresholdd for all the terms with pval < of PvalTh or for at least the minimun number of terms specified. 
# N.B. If less than minTerms terms respect the enrichment threshold, NA are returned (and could create issue in the generation pf the barplot) 

topGOResults <- function(Genes, gene2GO, ontology='BP', description=NULL, nodeSize=8, algorithm='weight01', statistic='fisher', PvalTh=0.01, EnTh=1.5, minTerms=15) {
  
  # Arguments:
  # Genes: named vector where names identify genes (e.g. gene symbol) and value gives information about differential expression (0: not differentially expressed; 1: differentially expressed). The vector gives this information about all the genes in the universe. 
  # gene2GO: 
  # Ontology: ontology domain
  # PvalTh: threshold on the pvalue for the selection of enriched terms
  # EnTh: threshold on the enrichment, calculated as significant/expected
  # minTerms reported if the significant threshold is not reached
  
  
  # 1. Setting description if null
  if (is.null(description)) {
    desc <- paste(ontology, deparse(substitute(Genes)))
  }
  
  # 2. Create list in which useful results are stored
  
  TopGORes <- list()
  
  # 3. TopGO object 
  TopGORes$GOdata <- new('topGOdata', ontology=ontology, allGenes=Genes, annotat=annFUN.gene2GO, gene2GO=gene2GO, description=desc, nodeSize=nodeSize)
  
  # 4. Statistical test
  TopGORes$Test <- runTest(TopGORes$GOdata, algorithm=algorithm, statistic=statistic)
  # message(Test@geneData)
  
  # 5. Generation of result table
  TopGORes$ResAll <- GenTable(TopGORes$GOdata, Statistics=TopGORes$Test, topNodes=length(TopGORes$Test@score)) 
  # Selection based on enrichment threshold
  ESel <- TopGORes$ResAll %>% dplyr::mutate(ER=round(Significant/Expected,2)) %>% dplyr::filter(ER > EnTh)
  n <- max(dim(dplyr::filter(ESel, as.numeric(Statistics) <= PvalTh))[1], minTerms)
  TopGORes$ResSel <- ESel[1:n,]
  # n for the selection of the results is set to the minimal number of terms if the selection
  # on pval gives a number of terms lower than the minimum number
  return(TopGORes)
}



### 8. === FUNCTION GOAnnotation: annotation for GO terms

#  Function from DEGCharacterizationHSapiensGenecode

GOAnnotation <- function(TopGOResults, GOdata, OutputFolder=getwd(), keytype='ENSEMBL') {
  # GenTableResults: oggetto risultante da GenTable con risultati statistica
  # GOdata: oggetto iniziale con ontology di riferimenti
  # OutputFolder: directory per storare tabelle
  # Keytype: key used to retrieve information about the gene. Default is Ensembl, to be changed to symbol for RefSeq
  GOvector <- as.vector(TopGOResults$GO.ID)
  for (i in 1:length(GOvector)) 
  {
    Genes <- genesInTerm(GOdata, as.character(GOvector[i]))
    Scores <- scoresInTerm(GOdata, as.character(GOvector[i]))
    GO <- as.data.frame(cbind(Genes[[1]], Scores[[1]]))
    names(GO) <- c("Gene", "Score")
    GOdifferentially <- GO[GO$Score==2,]
    GOannotation <- AnnotationDbi::select(org.Hs.eg.db,keys = as.character(GOdifferentially[,1]) ,columns=c("ENTREZID","SYMBOL","GENENAME","ENSEMBL"),keytype=keytype)
    GOTerm <- gsub( ':', '_', GOvector[i])
    write.table(GOannotation, file=paste0(OutputFolder, "/",GOTerm, "_annotation.txt"),quote=F, sep="\t" , row.names=F)
  }
}



### 8. === FUNCTION topGOBarplot: visualization of topGO results as barplot

# Function from DEGCharacterizationHSapiensGenecode for the visualization of topGO results as barplot

topGOBarplot <- function(TopGORes, terms=15, pvalTh=0.01, title=NULL, palette=NULL) {
  
  # TopGORes: object resulting from topGO analysis (produced by function above or from ...)
  # terms: number of terms to represent in the barplot
  # pvalTh: thresholf of pvalue to be plotted as line
  # title: plot title
  # palette: color for the barplot. If null, yellow-green is used
  # the function works considering that the results of the enrichment (pvalue) are in a column named Statistics
  
  # 1. Definition of the color palette
  if (is.null(palette)) {
    palette <- colorRampPalette(c('gold', 'forestgreen'))(terms)
  }
  
  # 2. Title definition
  if (is.null(title)) {
    title = deparse(substitute(TopGORes))
  }
  
  # Add full GO term
  GO <- as.list(GO.db::GOTERM)
  
  # add the complete GO term name
  TopGORes['extTerm'] <- stringr::str_wrap(sapply(as.character(TopGORes[['GO.ID']]),
                                                  function(x) X = GO[[x]]@Term))
  
  # 3. Dataframe reorder
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGORes$Statistics <- ifelse(grepl('<', TopGORes$Statistics), 1e-30, TopGORes$Statistics)
  # then I order the results
  ResOrdered <- transform(TopGORes, GO.ID=reorder(GO.ID, -as.numeric(Statistics)))[1:terms,]
  
  # 4. x-axis limit definition
  MaxVal <- round(max(-log10(as.numeric(ResOrdered$Statistics))), 0) +1
  
  # 4. BarPlot
  TopGOBarplot <- ggplot(data=ResOrdered, aes(x=GO.ID, y=-log10(as.numeric(Statistics)), fill=GO.ID)) + 
    geom_bar(stat='identity', aes(alpha=0.75)) +
    geom_text(aes(y=0), label=ResOrdered$extTerm, hjust=0) + 
    scale_y_continuous(breaks=seq(0,MaxVal,2), labels=abs(seq(0, MaxVal, 2)), limits=c(0,MaxVal), expand=c(0.025, 0.025)) +
    geom_hline(yintercept=-log10(pvalTh), col='darkred', lty='longdash') +
    coord_flip() + 
    scale_fill_manual(values=palette) +
    ylab('-log10 PValue') + xlab('') +
    ggtitle(paste('TopGO Enrichment results:\n', title)) +
    theme_bw() +
    theme(legend.position='none', axis.title.x = element_text(face = 'bold', colour = 'grey30', size=12), 
          plot.title= element_text(face='bold', colour='darkred', size=12))
}


### 9. === FUNCTION topGOBubbleplot: visualization of topGO results as bubbleplot

# Function for the visualization of topGO results as bubbleplot

# Dependencies: ggplots, ggrepel, 


topGOBubbleplot <- function(TopGORes, terms=15, Title=NULL, Col='blue', LabelSize=2) {
  
  # TopGORes: object resulting from topGO analysis (produced by function above or from ...)
  # terms: number of terms to represent in the barplot
  # pvalTh: thresholf of pvalue to be plotted as line
  # title: plot title
  # Color for the barplot.
  # the function works considering that the results of the enrichment (pvalue) are in a column named Statistics
  
  # 1. == Title definition
  if (is.null(Title)) {
    Title = deparse(substitute(TopGORes))
  }
  
  # 2. == Dataframe reorder
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGORes$Statistics <- as.numeric(ifelse(grepl('<', TopGORes$Statistics), 1e-30, TopGORes$Statistics))
  # then I order the results
  ResOrdered <- TopGORes %>% dplyr::top_n(terms, -Statistics) %>% dplyr::arrange(Statistics)
  # ResOrdered <- transform(TopGORes, GO.ID=reorder(GO.ID, -as.numeric(Statistics)))[1:terms,] 
  # checked that is equivalent to the one above
  
  # 3. == x-axis limit definition
  #MaxVal <- round(max(-log10(as.numeric(ResOrdered$Statistics))), 0) +1
  
  # 4. == BubblePlot
  TopGOBubble <- ggplot(data=ResOrdered, aes(x=ER, y=-log10(as.numeric(Statistics)))) + 
    geom_point(aes(size=ResOrdered$Significant), shape=1, color='gray35') + 
    geom_point(aes(size=ResOrdered$Significant), shape=19, color=Col, alpha=0.5) +
    ggrepel::geom_text_repel(label=ResOrdered$Term, cex=2, box.padding=0.5) +
    ylab('-log10 PValue') + xlab('Enrichment') + labs(size='Genes') +
    ggtitle(paste('TopGO Enrichment results: \n', Title)) +
    theme_bw() +
    theme(axis.title.x = element_text(face = 'bold', colour = 'grey30', size=12), 
          plot.title= element_text(face='bold', colour='darkred', size=12))
}




### 10. === FUNCTION: ENTREZ TEXT MINING FOR A VECTOR OF GENES


entrezGeneMining <- function(GeneVector, key='ENSEMBL'){
  
  # Arguments: 
  # I. GeneVector: vector of gene in Ensemble or other keys recognized by mapIds functions. 
  # II. key: specifies the key type, to be used by mapIds functions
  
  # Dependencies: rentrez, 
  
  # 1. Identifier change
  # duplicated values are discarded, as well as missing
  if(!is.null(key)) {
    EntrezGenes <-  mapIds(org.Hs.eg.db, keys=GeneVector, column="ENTREZID", keytype=key, multiVals="first") %>% unique()
    EntrezGenes <- EntrezGenes[!is.na(EntrezGenes)]}
  else{
    EntrezGenes <- GeneVector
  }
  
  # 2. Retrieval of information from entrez for each gene  
  EntrezSummary <- rentrez::entrez_summary(db="gene", id=EntrezGenes)
  
  # 3. All gene summary information is collapsed in a single string
  String <- ''
  for(i in 1:length(EntrezSummary)){
    String <- paste(String, EntrezSummary[[i]]$summary)
  }
  Res <- list()
  Res$GeneVector <- GeneVector
  Res$EntrezGenes <- EntrezGenes
  Res$EntrezSummary <- EntrezSummary
  Res$String <- String
  return(Res)
}




entrezMIMMining <- function(GeneVector, key='ENSEMBL'){
  
  # Arguments: 
  # I. GeneVector: vector of gene in Ensemble or other keys recognized by mapIds functions. 
  # II. key: specifies the key type, to be used by mapIds functions
  
  # Dependencies: rentrez, AnnotationDbi, org.Hs.eg.db;
  
  # 1. Identifier change
  # duplicated values are discarded, as well as missing
  if(!is.null(key)) {
    EntrezGenes <-  mapIds(org.Hs.eg.db, keys=GeneVector, column="ENTREZID", keytype=key, multiVals="first") %>% unique()
    EntrezGenes <- EntrezGenes[!is.na(EntrezGenes)]}
  else{
    EntrezGenes <- GeneVector
  }
  
  # 2. Retrieval of information from entrez for each gene  
  EntrezSummary <- rentrez::entrez_summary(db="gene", id=EntrezGenes)
  
  # 3. All gene summary information is collapsed in a single string
  String <- ''
  for(i in 1:length(EntrezSummary)){
    String <- paste(String, EntrezSummary[[i]]$summary)
  }
  Res <- list()
  Res$GeneVector <- GeneVector
  Res$EntrezGenes <- EntrezGenes
  Res$EntrezSummary <- EntrezSummary
  Res$String <- String
  return(Res)
}




moduleEntrezGene <- function(GeneMetrics, module, top=500, metric='kWithin'){
  
  # Arguments: 
  
  # I. GeneMetrics: dataframe containing gene annotation, module assignment and metrics for each gene. 
  #    Columns that are used: gene_biotype, external_gene_name, kWithin, module membership, Module
  # II. module: string specifying the module to be examined
  # III. title: string that will be incorporated in plot title
  # VI. top: integer specifying the maximum number of genes to be visualized, ranked according to kWithin
  
  # Dependencies: AnnotationDbi, org.Hs.eg.db;  relies also on function entrezGeneMining
  
  
  # 1. Data selection
  # First a selection is implemented to identify only the genes of the specified module
  Metrics <- dplyr::filter(GeneMetrics, Module==module)
  # Then, if top is not null, only the top-specified genes are selected
  if (!is.null(top)) {
    EnsGenes <- dplyr::top_n(Metrics, top, !!sym(metric)) %>% dplyr::pull(EnsGene) 
  } else {
    EnsGenes <- dplyr::pull(Metrics, EnsGene)   
  }
  
  entrezGeneMining(GeneVector=EnsGenes, key='ENSEMBL')
  
}




