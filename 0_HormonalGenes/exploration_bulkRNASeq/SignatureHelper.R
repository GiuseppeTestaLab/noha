options(stringsAsFactors = FALSE)

LogFpkmSignature <- function(MetaData, LogFpkmData, GeneSig){
  if (identical(MetaData$names, colnames(LogFpkmData))==FALSE){
    stop()
  }
  if(!is.character(GeneSig)){
    stop('Gene Signature is not a vector character')
  }
  if(!is.character(row.names(LogFpkmData))){
    stop('RowNames in LogFpKmData is not a vector character')
  }
  SigLogFpkm <- MetaData
  for (i in 1:length(GeneSig)){
    if(GeneSig[i] %in% row.names(LogFpkmData)){
      SigLogFpkm <- cbind(SigLogFpkm, as.vector(as.matrix(LogFpkmData[GeneSig[i],])))
      names(SigLogFpkm)[length(names(SigLogFpkm))] <- GeneSig[i]
    } else {
      SigLogFpkm <- cbind(SigLogFpkm, rep(0, dim(LogFpkmData)[2]))
      names(SigLogFpkm)[length(names(SigLogFpkm))] <- GeneSig[i]
    }
  }
  return(SigLogFpkm)
}

LogFpkmModule <- function(LogFpkmData, MetaData, Module, NSamples){
  if (identical(MetaData$names, colnames(LogFpkmData)[1:NSamples])==FALSE){
    stop()
  }
  if(!is.character(Module$EnsGene)){
    stop('EnsGene in Module is not a vector character')
  }
  if(!is.character(row.names(LogFpkmData))){
    stop('EnsGene in LogFpKmData is not a vector character')
  }
  ModuleLogFpkm <- dplyr::inner_join(LogFpkmData, Module, by='Gene')
  return(ModuleLogFpkm)
}

SigLogFpkmMean <- function(SigLogFpkm, StartCol=5, GroupCol='Stage', MetaInfo, Arrange='Population', Levels=NULL){ 
  if(is.null(Levels)){
    Levels <- levels(factor(dplyr::pull(MetaInfo, !!sym(Arrange))))
  }
  SigLogFpkmMean <- SigLogFpkm %>% tidyr::gather('GeneName', 'Exp', StartCol:dim(SigLogFpkm)[2]) %>%
    dplyr::group_by(GeneName, !!sym(GroupCol)) %>% dplyr::summarize(Count=n(), ExpMean=mean(Exp))
  if(!is.character(MetaInfo$GeneName) | !is.character(SigLogFpkmMean$GeneName)){
    stop('Gene Name in MetaInfo is not a vector character')
  }
  SigLogFpkmMean <- dplyr::left_join(SigLogFpkmMean, MetaInfo, by='GeneName') %>%
    dplyr::arrange(factor(!!sym(Arrange), levels=Levels))
  return(SigLogFpkmMean)
}

DeltaLogFpkmMean <- function(LogFpkmMean, GroupCol='Stage', DiscardStage=NULL, RefStage='S25', Arrange='Population', 
                             Levels=NULL, MetaCols=2, Th=-2, Substitute=NA){
  if(is.null(DiscardStage)){
    LogFpkmDelta <- LogFpkmMean 
  } else {
    LogFpkmDelta <- LogFpkmMean %>% filter(Stage != DiscardStage)
  }
  if(is.null(Levels)){
    Levels <- levels(factor(dplyr::pull(LogFpkmDelta, !!sym(Arrange))))
  }
  LogFpkmDelta <- LogFpkmDelta %>% dplyr::select(-Count) %>% spread(Stage, ExpMean) 
  StageNum <- dim(LogFpkmDelta)[2] - MetaCols 
  names(LogFpkmDelta)[(MetaCols+1):(StageNum+MetaCols)] <- paste0('S', names(LogFpkmDelta)[(MetaCols+1):(StageNum+MetaCols)]) 
  for (i in 1:StageNum) {
    if (names(LogFpkmDelta)[i+MetaCols] != RefStage){
      LogFpkmDelta <- LogFpkmDelta %>% mutate(!!sym(names(LogFpkmDelta)[i+MetaCols]) - !!sym(RefStage))
      names(LogFpkmDelta)[length(names(LogFpkmDelta))] <- paste0('D', i, '_', names(LogFpkmDelta)[i+MetaCols],'_', RefStage)
      if(is.null(Th)==FALSE){
        for(row in 1:dim(LogFpkmDelta)[1]){  
          LogFpkmDelta[row, names(LogFpkmDelta)[length(names(LogFpkmDelta))]] <- ifelse((LogFpkmDelta[row,RefStage] < Th & 
                                                                                           LogFpkmDelta[row,names(LogFpkmDelta)[i+MetaCols]] < Th), 
                                                                                        Substitute, 
                                                                                        LogFpkmDelta[row, names(LogFpkmDelta)[length(names(LogFpkmDelta))]])
        }
      }
    }
  }
  LogFpkmDelta <- LogFpkmDelta %>% dplyr::select(-c((MetaCols+1):(StageNum+MetaCols))) %>% 
    tidyr::gather('Stage', 'ExpMean', (MetaCols+1):(StageNum+MetaCols-1)) %>%
    dplyr::arrange(factor(!!sym(Arrange), levels=Levels))
  return(LogFpkmDelta)
}

LolliPopCols <- function(Data, PointSize=2, SegSize=1, cols=c('palegreen', 'turquoise', 'blue', 'purple'), yLow=0, 
                         yHigh=9, unit='LogFpkm', CeilUp=NULL, CeilDown=NULL, Shape=16, GeneSize=4, filterZero=FALSE){
  if(!is.null(CeilUp)){
    Data$ExpMean <- ifelse(Data$ExpMean > CeilUp, CeilUp, Data$ExpMean)
  }
  if(!is.null(CeilDown)){
    Data$ExpMean <- ifelse(Data$ExpMean < CeilDown, CeilDown, Data$ExpMean)
  }
  if(is.null(cols)){
    cols <- scales::hue_pal()(length(unique(Data$Population)))
  }
  if(filterZero==TRUE){
    Data <- dplyr::filter(Data, ExpMean!=0)
  }
  Plot <- ggplot(data=Data, aes(y=ExpMean, 
                                x=factor(Data$GeneName, levels=unique(Data$GeneName)), col=Population)) +
    geom_point(stat='identity', size=PointSize, shape=Shape)  +
    geom_segment(aes(x=GeneName, xend=GeneName, y=0, yend=ExpMean), size=SegSize) + 
    scale_color_manual(values=cols, limits=unique(Data$Population)) +
    xlab('') + ylab(paste('Expression', unit)) + ylim(yLow, yHigh) +
    theme_light()  +
    theme(panel.grid.major.x=element_blank(),
          axis.ticks.y=element_blank(), legend.position='bottom', axis.text=element_text(size=GeneSize, angle=45, hjust=1),
          axis.title.y=element_text(size=14, face='bold', colour='#35A2FF'))
  PlotGrid <- Plot + facet_grid(Data$Stage~.) + theme(strip.text.x = element_text(size=16, face='bold'), 
                                                      strip.background = element_rect(colour='grey', fill="#35A2FF"))
  return(PlotGrid)
}

LolliGene <- function(Data,  PointSize=2, SegSize=1, cols=c(100, 360), yLow=0, yHigh=9, unit='LogFpkm', Ceil=NULL, Shape=16, GeneSize=4){
  if(!is.null(Ceil)){
    Data$ExpMean <- ifelse(Data$ExpMean <= abs(Ceil), Data$ExpMean,
                           ifelse(Data$ExpMean > Ceil, Ceil, -Ceil))
  }
  Plot <- ggplot(data=Data, aes(y=ExpMean, 
                                x=factor(Data$GeneName, levels=unique(Data$GeneName)), col=GeneName)) +
    geom_point(stat='identity', size=PointSize, shape=Shape)  +
    geom_segment(aes(x=GeneName, xend=GeneName, y=0, yend=ExpMean), size=SegSize) + 
    scale_color_hue(h=cols, limits=Data$GeneName) + 
    xlab('') + ylab(paste('Expression', unit)) + ylim(yLow, yHigh) +
    theme_light()  +
    theme(panel.grid.major.x=element_blank(),
          axis.ticks.y=element_blank(), legend.position='bottom', axis.text=element_text(size=GeneSize, angle=45, hjust=1), 
          axis.title.y=element_text(size=16, face='bold', colour='#35A2FF')) 
  PlotGrid <- Plot + facet_grid(Data$Stage~.) + theme(strip.text.x = element_text(size=16, face='bold'), 
                                                      strip.background = element_rect(colour='grey', fill="#35A2FF"))
  return(PlotGrid)
}

boxStage <- function(Data, Signature, Title=NULL, Fill='turquoise', yHigh=NULL, yLow=NULL, 
                     filterZero=FALSE, sampleOrd=NULL) {
  if (is.null(Title)) {
    Title <- paste0('Expression levels for ', Signature, ' signature')
  }  
  Sig <- dplyr::filter(Data, Population==Signature)
  if(filterZero==TRUE){
    Sig <- dplyr::filter(Sig, ExpMean!=0)
  }
  if (is.null(yHigh)) {
    yHigh <- max(Sig$ExpMean) + 0.5
  } 
  if (is.null(yLow)) {
    yLow <- min(Sig$ExpMean) - 0.5
  } 
  if (is.null(sampleOrd)) {
    sampleOrd <- as.character(unique(Sig$Stage))
  } 
  Plot <- ggplot(Sig, aes(x=as.factor(Stage), y=ExpMean)) + 
    geom_boxplot(fill=Fill, alpha=0.5, outlier.size = 1E-10) + 
    geom_jitter(width=0.1, height=0, col='gray40', shape=16, size=0.75) +
    ggtitle(Title) + 
    ylim(yLow, yHigh) +
    scale_x_discrete(limits = sampleOrd) +
    theme_bw() +
    theme(plot.title = element_text(face='bold', colour='darkred', size=14, hjust=0.5))
  return(Plot)
}

scatterME <- function(MetaInfo, Predict, title=NULL, Shape=16, PointSize=5, Col='blue', yLow=-10, yHigh=10, Sign=1) {
  if (is.null(title)) {
    title <- 'Module Component 1'
  }  
  Pr <- data.frame(names=row.names(Predict), PC1=(Predict[,1]*Sign))
  Data <- dplyr::inner_join(MetaInfo, Pr, by='names')
  ggplot(data=Data, aes(x=Stage, y=PC1)) +
    geom_point(stat='identity', size=PointSize, shape=Shape, color=Col) +
    stat_summary(aes(y=PC1), fun.y=mean, colour=Col, alpha=0.45, size=3.4, geom='line') +
    ggtitle(title) + ylim(yLow, yHigh) + xlim(10, 210) +
    theme_light()  + 
    theme(plot.title = element_text(face='bold', colour='darkred', size=13, hjust=0.5)) 
}
