plotGenesInTerm_v2 <- function(TopGOResults, GOdata, SE, nterms=12, ngenes=8,
                            plotTitle=NULL, Interactive=FALSE, fillCol='forestgreen'){
  
  #Prepare the DF for plotting
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGOResults$Statistics <- ifelse(grepl('<', TopGOResults$Statistics), 1e-30,
                                    TopGOResults$Statistics)
  # then I order the results
  ResOrdered <- transform(TopGOResults,
                          GO.ID=reorder(GO.ID,
                                        -as.numeric(Statistics)))[1:nterms,]
  #Create the data frame
  finalDF <- data.frame()
  
  for (i in seq(nterms)){
    if (!is.na(ResOrdered$GO.ID[i])) {
      
      #Create the basic dataframe
      GOid <- as.character(ResOrdered$GO.ID[i])
      GOterm <- as.character(ResOrdered$Term[i])
      Genes <- topGO::genesInTerm(GOdata, GOid)
      Scores <- topGO::scoresInTerm(GOdata, GOid)
      
      GOdf <- as.data.frame(cbind(Gene=unlist(Genes), Score=unlist(Scores)),
                            row.names = FALSE) %>%
        filter(Score == 2) %>%  #select genes that are present
        mutate(
          GOid = case_when(TRUE ~ GOid),   #repeat GOid
          GOterm = case_when(TRUE ~ GOterm)) %>%
        mutate(across(c(GOterm, Gene), as.factor)) %>%  #change types
        mutate(across(c(Score), as.integer))
      
      
      #Add FDR for DEG significance
      GOdf['FDR'] <- rowData(SE)[which(rownames(SE) %in% GOdf$Gene), c('FDR')]
      GOdf['logFC'] <- rowData(SE)[which(rownames(SE) %in% GOdf$Gene), c('logFC')]
      
      #Sort and add rank for x axis in plot
      GOdf <- GOdf %>%
        arrange(FDR) %>%
        slice_head(n = ngenes)
      
      #Bind to dataframe
      finalDF <- as.data.frame(rbind(finalDF, GOdf))
      
    } else{next}
  }
  
  ## Plot
  #Graphical stuff
  if (is.null(plotTitle)) { #Title
    plotTitle = "Genes in Term"
  }
  sub <- paste("Top", nterms, "GO terms | Top", ngenes, "genes by FDR")
  #palette <- colorRampPalette(c('#dfe5ef', '#0FBBE6', '#0635EE'))(nterms)
  
  BP <- ggplot(finalDF, aes(y=GOterm, x=.data$Score, fill=.data$FDR, label=.data$Gene, label2=GOid)) +
    geom_col(col = 'white', alpha = 0.6, width = 1) +
    #viridis::scale_fill_viridis(begin = 0.3) +
    scale_fill_gradient(low = fillCol, high = "grey90") + #mid = "white", #midpoint = .02
    scale_y_discrete(limits = rev(levels(finalDF$GOterm))) + #otherwise most significant below
    labs(x='Genes in Term', y='GO Terms',
         title=plotTitle, subtitle=sub) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
          plot.subtitle = element_text(size=12, hjust=0.5), #colour='darkred',
          text = element_text(size = 15),
          axis.text = element_text(size = 14),
          axis.text.x = element_blank(), #axis.text.y = element_blank(),
          axis.ticks.x =element_blank(), axis.ticks.y = element_blank())
  
  #4. Return interactive or static plot
  if (Interactive == TRUE){  #With ggplotly
    #ggplotly doesn't show subtitles: we use this workaround
    BP <- plotly::ggplotly(BP) %>%
      plotly::layout(title= list(text=paste0(plotTitle, '<br>', '<sup>', sub),
                                 tooltip=c("Gene", "GOid", "FDR"))) # selects the aesthetics to include in the tooltip
    return(BP)
    
  } else {  #Static plot
    BP <- BP + geom_text(size = 3.5, position = position_stack(vjust = 0.5), )
    return(BP)
  }
}
