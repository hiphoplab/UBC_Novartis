
#library(XML)
library(dplyr)
library(igraph)
#library(xml2)
library(visNetwork)


########
# myhref1 = function(link,nameoflink){
#   link = paste("<a href=",link,">",nameoflink,"</a>",sep="")
#   #w = which(names(df)==links)
#   #df = df[,-w,drop=F]
#   link
# }
########

########putting the two tables together
# mygotable = function(lst){
#   myhref1 = function(link,nameoflink){
#     link = paste("<a href=",link,">",nameoflink,"</a>",sep="")
#     link
#   }
#   links = sapply(lst,function(x) x = unlist(x[2],use.names = F))
# 
#   names(links) = names(lst)
# 
#   links = myhref1(link = links,nameoflink = names(links))
# 
#   dinks = data.frame(condition = names(lst), GO_enrichment = links,stringsAsFactors = F)
# 
#   dinks
# }

#########
myreadSigs = function (sigFileName) 
{
  sigs <- scan(sigFileName, what = character(), sep = "\n")
  sigs <- strsplit(sigs, split = "\t")
  names(sigs) <- sapply(sigs, function(vec) {
    vec[1]
  })
  sigs <- lapply(sigs, function(vec) {
    setdiff(vec[-c(1:2)], "NULL")
  })
  lens <- sapply(sigs, length)
  sigs[lens > 0]
}


########
#bp_file = myreadSigs("BP_geneNames_current.gmt")
########
myrun_go_enrich1 = function (fdrThresh = 0.1, curr_exp, score, bp_path){
    
    fdrThresh = as.numeric(fdrThresh)
   
    bp_file = file.path(bp_path)
    scoreMat = score
    queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index >= 1)]))
    uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))
    bp <- myreadSigs(bp_file)
    enrichMat.mn <- myhyperG(querySet = queryGenes.mn, geneSets = bp,
      uni = uniGenes.mn, scoreMat = score, minSetSize = 5,
      maxSetSize = 300, uniSize = NA)
    queryGeneSets = list()
    queryGeneSets[[curr_exp]] = queryGenes.mn
    enrichMat.mn$filename <- curr_exp
    enrichMat_Ordered = enrichMat.mn[with(enrichMat.mn, order(FDR,
      -foldEnrichment)), ]
    scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),
      ]
    scoreMat <- scoreMat[match(uniGenes.mn, scoreMat$gene), "score",
      drop = F]
    rownames(scoreMat) <- uniGenes.mn
    colnames(scoreMat) <- curr_exp
    #head(scoreMat, 3)

    nonEnrichMat.mn <- mygenesNotInEnrichedTerm(queryGeneSets,
      enrichMat.mn, scoreMat, bp$NONSPECIFIC.TERMS, fdrThresh)

    bp <- lapply(bp, intersect, uniGenes.mn)
    lens <- sapply(bp, length)
    bp <- bp[lens >= 5 & lens <= 300]
   
    q = myclusterEnrich(enrichInfo = enrichMat.mn, geneSets = bp,
      outFile = curr_exp, fdrThresh = fdrThresh, overlapThresh = 0.5,
      nonEnrichInfo = nonEnrichMat.mn, barModGenes = NULL,
      scoreName = "score", plotForEachEnrichedTerm = T, goTable = NULL)
    return(list(enrichInfo = q$enrichInfo , edgeMat = q$edgeMat))
}

#####
# generates an enrichment map in xgmml format using hypergeometric test statistics
# enrichInfo - dataframe with enrichment stats for gene sets (one per row), with the following columns:
#              term, geneSetFraction, querySetFraction, FDR, overlapGenes, maxOverlapGeneScore
#            - see documentation for the output of hyperG() for descriptions of these columns
# geneSets - named list of gene sets tested for significant overlap w/ the query set,
#              restricted to genes in the universe
# outFile - the output xgmml file will be saved to this location
# fdrThresh - FDR threshold; only show gene sets that pass this significance threshold
# overlapThresh - an edge between a pair of enriched gene sets will only be shown if the overlap coefficient
#              is >= overlapThresh
# nonEnrichInfo - dataframe with info on gene sets that are *not* significantly enriched (one per row)
#             with the following columns:
#             term, overlapGenes, maxOverlapGeneScore, geneSetFraction, unenrichedGenes
#            - see documentation for the output of genesNotInEnrichedTerm() for descriptions of these columns
#            - can be NULL
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly in the
#              barplots, should they be in the top overlap genes
# scoreName - score label to use in top overlap gene barplots
# plotForEachEnrichedTerm - if TRUE, a top overlap gene barplot will be created for each enriched term;
#              if FALSE, a barplot will be created for each enriched term cluster
# goTable - a dataframe with the following columns describing GO terms:
#         - "term" (GO term), "id" (GOID)
#         - if provided (i.e. not NULL), the GO ID numbers of the enriched GO terms will be saved in
#           the output xgmml file as a node attribute called "GOID".
#         - the GOID allows for an easy link to the GO website page for the associated GO term
#########
#####
myclusterEnrich = function (enrichInfo, geneSets, outFile, fdrThresh = 0.1, overlapThresh = 0.5,
  nonEnrichInfo = NULL, barModGenes = NULL, scoreName = "Fitness defect score",
  plotForEachEnrichedTerm = F, goTable = NULL){
  nodeSizeRange <- c(10, 40)
  prunedCol <- "#BEBEBE"
  labelWidth <- 20
  edgeWidthRange <- c(1, 5)
  overlapCoeffRange <- c(overlapThresh, 1)
  
  if (!is.null(nonEnrichInfo)) {
    nonEnrichInfo$maxOverlapGeneScore <- round(nonEnrichInfo$maxOverlapGeneScore,
      2)
    nonEnrichInfo$geneSetFraction <- round(nonEnrichInfo$geneSetFraction *
        100, 1)
    if (is.null(nonEnrichInfo$nGenes)) {
      lens <- sapply(geneSets, length)
      i <- match(nonEnrichInfo$term, names(lens))
      nonEnrichInfo$nGenes <- lens[i]
    }
    tmp <- strsplit(nonEnrichInfo$overlapGenes, "\\|")
    w <- which(is.na(tmp))
    if(length(w)>0) tmp = tmp[-w]
    if (is.null(nonEnrichInfo$unenrichedGenes)) {
      nonEnrichInfo$overlapGenes <- sapply(tmp, paste,
        collapse = "| ")
    }
    else {
      unEnriched <- strsplit(nonEnrichInfo$unenrichedGenes,
        "\\|")
      tmp.mod <- sapply(1:length(tmp), function(termI) {
        vec <- tmp[[termI]]
        i <- match(unEnriched[[termI]], tmp[[termI]])
        vec[i] <- paste("<b>", vec[i], "</b>", sep = "")
        paste(vec, collapse = "| ")
      })

      nonEnrichInfo$overlapGenes <- tmp.mod
    }
    
    if (is.null(enrichInfo)) {
      
      return()
    }
  }
  ###### key step to exit
  enrichInfo <- enrichInfo[enrichInfo$FDR <= fdrThresh, , drop = F]
  #print(enrichInfo[1:5,1:5])
  #print(dim(enrichInfo))
  enrich = enrichInfo
  if (nrow(enrichInfo) == 0) {
    print("No enriched terms to cluster.")
    return()
  }
  enrichInfo$formattedLabel <- sapply(enrichInfo$term, function(curLabel) {
    curLabel <- strwrap(curLabel, labelWidth)
    paste(curLabel, collapse = "\n")
  })
  i <- match(enrichInfo$term, names(geneSets))
  if (any(is.na(i))) {
    stop("Could not find gene sets for ", sum(is.na(i)),
      " enriched terms.")
  }
  geneSets <- geneSets[i]
  if (is.null(enrichInfo$nGenes)) {
    enrichInfo$nGenes <- sapply(geneSets, length)
  }
  tmpSize <- -log10(enrichInfo$FDR)
  maxVal <- max(tmpSize[!is.infinite(tmpSize)])
  tmpSize[is.infinite(tmpSize)] <- maxVal + 2
  gsSizeRange <- range(tmpSize)
  if (gsSizeRange[1] == gsSizeRange[2]) {
      gsSizeRange[1] <- -log10(fdrThresh)
      gsSizeRange[2] <- gsSizeRange[2] + 1
  }
  tmpSize <- (tmpSize - gsSizeRange[1])/(gsSizeRange[2] - gsSizeRange[1])
  tmpSize <- nodeSizeRange[1] + tmpSize * (nodeSizeRange[2] -nodeSizeRange[1])
  enrichInfo$size <- round(tmpSize, 2)
  if (nrow(enrichInfo) == 1) {
    enrichInfo$cluster <- CLUST.COL[1]
    edgeMat <- NULL
  }
  else {
    pairI <- getUniquePairs(length(geneSets))
    distVal <- apply(pairI, 1, function(onePair) {
      myoverlapCoeff(geneSets[onePair])
    })
    distVal[distVal < overlapThresh] <- 0
    edgeMat <- data.frame(nodeA = pairI[, 1], nodeB = pairI[,
      2], coeff = distVal)
    enrichInfo$cluster <- prunedCol
    if (is.null(enrichInfo$pruneOutcome)) {
      termI <- 1:nrow(enrichInfo)
    }
    else {
      termI <- which(enrichInfo$pruneOutcome == enrichInfo$term)
    }
    if (length(termI) == 1) {
      enrichInfo$cluster[termI] <- CLUST.COL[1]
    }
    else {
      i <- which((edgeMat$nodeA %in% termI) & (edgeMat$nodeB %in%
          termI))
      enrichInfo$id = termI
      g=igraph::graph_from_data_frame(edgeMat[which(edgeMat$coeff!=0),],directed = F,vertices = enrichInfo$id)
      adj = igraph::as_adjacency_matrix(g)
      clusters = igraph::clusters(g)
      clusters = split(names(clusters$membership),clusters$membership)
      #clusters <- mclWrapper(edgeMat[i, , drop = F], dirname(outFile))
      mcl.in = edgeMat[i, , drop = F]
      if(is.null(mcl.in)) print("edge NULL")
      mcl.out = clusters
      clusters <- lapply(clusters, as.numeric)
      # if (length(clusters) > length(CLUST.COL)) {
      #   stop("Need more cluster colours!")
      # }
      lens <- sapply(clusters, length)
      clusters <- data.frame(id = unlist(clusters), cluster = CLUST.COL[rep(1:length(clusters),
        lens)], stringsAsFactors = F)
      enrichInfo$cluster[clusters$id] <- clusters$cluster
    }
    edgeMat <- edgeMat[edgeMat$coeff > 0, , drop = F]
    if (nrow(edgeMat) > 0) {
      edgeMat$size <- (edgeMat$coeff - overlapCoeffRange[1])/(overlapCoeffRange[2] -
          overlapCoeffRange[1])
      edgeMat$size <- edgeWidthRange[1] + edgeMat$size *
        (edgeWidthRange[2] - edgeWidthRange[1])
      edgeMat$coeff <- round(edgeMat$coeff, 2)
      edgeMat$size <- round(edgeMat$size, 2)
    }
    else {
      edgeMat <- NULL
    }
  }
  otherI <- order(enrichInfo$cluster)
  otherI <- otherI[order(enrichInfo$FDR[otherI])]
  termI <- which(enrichInfo$cluster[otherI] != prunedCol)
  if (length(termI) < length(otherI)) {
    otherI <- c(otherI[termI], otherI[-termI])
  }
  enrichInfo$id <- 1:nrow(enrichInfo)
  enrichInfo <- enrichInfo[otherI, , drop = F]
  enrichInfo$geneSetFraction <- round(enrichInfo$geneSetFraction *
      100, 1)
  
  # allNodes <- sapply(1:nrow(enrichInfo), function(nodeI) {
  #   newNode <- newXMLNode("node", attrs = c(id = enrichInfo$id[nodeI],
  #     label = enrichInfo$term[nodeI]), parent = graphNode,
  #     addFinalizer = T)
  #   newXMLNode("att", attrs = c(type = "string", name = "formattedLabel",
  #     value = enrichInfo$formattedLabel[nodeI]), parent = newNode)
  #   newXMLNode("att", attrs = c(type = "string", name = "overlapGenes",
  #     value = enrichInfo$overlapGenes[nodeI]), parent = newNode)
  #   newXMLNode("att", attrs = c(type = "string", name = "cluster",
  #     value = enrichInfo$cluster[nodeI]), parent = newNode)
  #   newXMLNode("att", attrs = c(type = "real", name = "FDR",
  #     value = enrichInfo$FDR[nodeI]), parent = newNode)
  #   newXMLNode("att", attrs = c(type = "real", name = "geneSetFraction",
  #     value = enrichInfo$geneSetFraction[nodeI]), parent = newNode)
  #   newXMLNode("att", attrs = c(type = "real", name = "querySetFraction",
  #     value = enrichInfo$querySetFraction[nodeI]), parent = newNode)
  #   newXMLNode("att", attrs = c(type = "integer", name = "nGenes",
  #     value = enrichInfo$nGenes[nodeI]), parent = newNode)
  #   newNode
  # })
  
  
  
  # if (plotForEachEnrichedTerm) {
  #   maxScore <- max(enrichInfo$maxOverlapGeneScore,na.rm=T, na.rm = T)
  #   if (!is.null(nonEnrichInfo)) {
  #     maxScore <- max(maxScore, max(nonEnrichInfo$maxOverlapGeneScore, na.rm = T))
  #   }
  #   tmp <- sapply(enrichInfo$overlapGenes, mygenOverlapGenePlot.gChart,
  #     c(0, maxScore), barModGenes, scoreLabel = scoreName)
  #   enrichInfo$image <- tmp[3, ]
  # }
  # else {
  #   nodeInfo.plot <- t(clusterPlots(enrichInfo, scoreName,
  #     barModGenes))
  #   colnames(nodeInfo.plot) <- c("w", "h", "image")
  #   nodeInfo.plot <- as.data.frame(nodeInfo.plot)
  #   nodeInfo.plot$w <- round(as.numeric(as.character(nodeInfo.plot$w)),
  #     2)
  #   nodeInfo.plot$h <- round(as.numeric(as.character(nodeInfo.plot$h)),
  #     2)
  #   nodeInfo.plot$image <- as.character(nodeInfo.plot$image)
  #   nodeInfo.plot$cluster <- rownames(nodeInfo.plot)
  #   nodeInfo.plot$id <- substring(nodeInfo.plot$cluster,
  #     2)
  #tmp <- strsplit(enrichInfo$overlapGenes, "\\|")
  #enrichInfo$overlapGenes <- sapply(tmp, paste, collapse = "| ")
  # }
  #toDoI <- which(enrichInfo$cluster != prunedCol)
  #toDoI <- toDoI[enrichInfo$FDR[toDoI] == min(enrichInfo$FDR[toDoI])]
  
  ### here come the edge matrrix
  if (is.null(edgeMat)) print("edgeMat is NULL")
  #below doesn't work...
  #if(is.null(edgeMat))  {stop("No fucking Go enrichment")}
  
  if (!is.null(edgeMat)) {
    nam = c("source","target","label","overlapCoeff","width")
    orig = c("nodeA","nodeB","label","coeff","size")
    m = match(edgeMat$nodeA,as.numeric(as.factor(names(geneSets))))
    src = names(geneSets)[m]
    m = match(edgeMat$nodeB,as.numeric(as.factor(names(geneSets))))
    trg = names(geneSets)[m]
    edgeMat$label = paste(src,"(overlap)",trg)
    m = match(names(edgeMat),orig)
    names(edgeMat) = nam[m]
  }
 
  #if (is.null(edgeMat)) edgeMat = ""
  #if(is.null(edgeMat))  {stop("No Go enrichment")}
  output = list(enrichInfo = enrichInfo,edgeMat = edgeMat)
  
  return(output)
}
# generates overlap gene (i.e. those that drive enrichment) barplots for each cluster
# enrichInfo - dataframe with info on enriched terms (one per row), with the following columns:
#              overlapGenes (genes that drive the enrichment of the term), cluster (cluster of the term)
# scoreLabel - score label to use in the barplots
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly in the
#              barplots, should they be in the top overlap genes
# scoreInfo - if provided (i.e. not NULL), a data.frame of scores with the following columns:
#              gene (gene IDs) and score (gene scores)
# maxGenes - maximum number of top overlap genes to show in a barplot
# posScores - if TRUE, enrichment is amongst genes with positive scores; if FALSE, negative scores
# plotCol - the colour of the bars, in hexidecimal format without the '#' character
# RETURNS a matrix of plot info where each column corresponds to a different cluster bar plot,
# and the rows correspond to plot width (1st), plot height (2nd) hand the google chart URL for the plot (3rd)
clusterPlots <- function(enrichInfo, scoreLabel, barModGenes=NULL, scoreInfo=NULL, maxGenes=10, posScores=T, plotCol="BEBEBE") {
  barWidth <- 10

  overlapGenes <- strsplit(enrichInfo$overlapGenes, "\\|")

  if (is.null(scoreInfo)) {
    oGenes <- unique(unlist(overlapGenes))
    genes <- strsplit(oGenes, "\\(")
    scores <- sapply(genes, function(vec) { vec[length(vec)] })
    genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
    scores <- unlist(strsplit(scores, ")"))
    scoreInfo <- data.frame(gene=genes, score=as.numeric(scores), geneStr=oGenes, stringsAsFactors=F)
  }
  else {
    scoreInfo$geneStr <- scoreInfo$gene
  }

  scoreInfo <- scoreInfo[order(scoreInfo$score, decreasing=posScores), ]
  uniGenes <- unique(scoreInfo$gene)
  i <- match(uniGenes, scoreInfo$gene)
  scoreInfo <- scoreInfo[i, ]

  clusters <- split(1:nrow(enrichInfo), enrichInfo$cluster)

  plotData <- lapply(clusters, function(clustI) {
    clusterGenes <- table(unlist(overlapGenes[clustI]))
    clusterGenes <- clusterGenes/length(clustI) * 100

    i <- match(names(clusterGenes), scoreInfo$geneStr)
    clusterGenes <- cbind(clusterGenes, scoreInfo$score[i])
    colnames(clusterGenes) <- c("% of gene sets", "score")
    rownames(clusterGenes) <- scoreInfo$gene[i]

    # sort the genes by score, then by % of gene sets they are in
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]
    orderI <- order(clusterGenes[, 1], decreasing=T)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (nrow(clusterGenes) > maxGenes) {
      clusterGenes <- clusterGenes[1:maxGenes, , drop=F]
    }

    # re-order by score
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (is.null(barModGenes)) {
      return(clusterGenes)
    }

    cbind(clusterGenes, rownames(clusterGenes) %in% barModGenes) })

  if (posScores) {
    scoreRange <- c(0, ceiling(max(scoreInfo$score, na.rm=T)))
  }
  else {
    scoreRange <- c(floor(min(scoreInfo$score, na.rm=T)), 0)
  }

  sapply(plotData, mygenLeadingEdgePlot.gChart, plotCol, scoreRange, barWidth, scoreLabel)
}

# generate a barplot of the common leading edge genes using google charts (can be visualized
# with the html img tag), bar length corresponds to score
# leadInfo - matrix/dataframe of leading edge genes; rownames = gene IDs, column 1 = % of gene sets, column 2 = score
#          - if a 3rd column is provided, it should provide TRUE/FALSE indicating whether or not the gene
#            should be marked
# plotCol - the colour of the bars; hexidecimal format without the '#' character
# scoreRange - a vector of the range of scores in the profile; 1st value = min,  2nd value = max
# barWidth - bar width in the plot
# scoreLabel - label for the score axis
# RETURNS a vector of plot info: plot width, plot height and the google chart URL for the plot
#
#####hiphop:::mclWrapper
# MCL clustering (http://micans.org/mcl/)
# requires MCL to be installed and runnable from any directory
# simMat - 3-column similarity matrix; each row specifies a pair of items followed by a similarity value (>= 0)
# runDir - the directory where the input and output files will be saved
# inflation - the inflation parameter
# RETURNS a list where each element in the list is a separate cluster, and the element is a vector of cluster items
# mclWrapper = function (simMat, runDir, inflation = 2)
# {
#   inFile <- paste(runDir, "mcl.in", sep = "/")
#   write.table(simMat, inFile, row.names = F, col.names = F,
#     sep = "\t", quote = F)
#   outFile <- paste(runDir, "mcl.out", sep = "/")
#   cmdStr <- paste("mcl", inFile, "--abc -I", inflation, "-o",
#     outFile, "&>/dev/null")
#   system(cmdStr)
#   clusters <- scan(outFile, what = character(), sep = "\n")
#   strsplit(clusters, "\t")
# }

mygenLeadingEdgePlot.gChart =
  function (leadInfo, plotCol, scoreRange, barWidth, scoreLabel = "Sensitivity")
  {
    dataRange <- scoreRange[2] - scoreRange[1]
    barLens <- round((leadInfo[, 2] - scoreRange[1])/dataRange *
        100)
    if (dataRange <= 1) {
      stepSize <- 0.5
    }
    else if (dataRange <= 5) {
      stepSize <- 2
    }
    else if (dataRange <= 20) {
      stepSize <- 5
    }
    else if (dataRange <= 50) {
      stepSize <- 10
    }
    else if (dataRange <= 100) {
      stepSize <- 20
    }
    else if (dataRange <= 500) {
      stepSize <- 100
    }
    else {
      stepSize <- 250
    }
    scoreLabel <- unlist(strsplit(scoreLabel, ""))
    scoreLabel[scoreLabel == " "] <- "+"
    scoreLabel <- paste(scoreLabel, collapse = "")
    tmpStep <- 100/length(barLens)
    labelPos <- round(seq(tmpStep/2, 100, by = tmpStep))
    w <- 150
    h <- barWidth * length(barLens) + 50
    if (scoreRange[1] < 0) {
      zeroLineStr <- paste("&chp=", round(abs(scoreRange[1])/dataRange,
        2), sep = "")
    }
    else {
      zeroLineStr <- ""
    }
    if (ncol(leadInfo) > 2 && any(leadInfo[, 3] == 1)) {
      barModStr <- paste("o,000000,", which(leadInfo[, 3] ==
          1) - 1, ",-1,5", sep = "")
      barModStr <- paste("&chm=", paste(barModStr, collapse = "|"),
        sep = "")
    }
    else {
      barModStr <- ""
    }
    c(w, h, paste("http://chart.apis.google.com/chart?chxt=x,x,y&chs=",
      w, "x", h, "&cht=bhg&chd=t:", paste(barLens, collapse = "|"),
      "&chco=", plotCol, "&chxl=1:|", scoreLabel, "|2:|", paste(rev(rownames(leadInfo)),
        collapse = "|"), "&chxp=1,50|2,", paste(labelPos,
          collapse = ","), "&chxr=0,", scoreRange[1], ",",
      scoreRange[2], ",", stepSize, "&chbh=", barWidth, ",1,0",
      zeroLineStr, barModStr, sep = ""))
  }
####
#this is equivalent to t(combn(maxVal,2))
##### required functions
# computes the number of unique pairs given the number of items to consider
# maxVal - the maximum number of items
# RETURNS a 2-column matrix where each row contains a different pair, specified with item indices
getUniquePairs = function (maxVal)
{
  firstI <- rep(1:(maxVal - 1), (maxVal - 1):1)
  secondI <- sapply(2:maxVal, function(x) {
    x:maxVal
  })
  cbind(firstI, unlist(secondI))
}

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
    max = 255,
    alpha = (100 - percent) * 255 / 100,
    names = name)
  
  ## Save the color
  invisible(t.col)
}

lighten <- function(color, factor=1.4){
col <- col2rgb(color)
col <- col*factor
col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
col
}


CLUST.COL <- c("#FF00CC","#33CCFF", "#33CC00", "#9900FF", "#FF9900", "#FFFF00", "#FFCCFF", "#FF0000", "#006600", "#009999", "#CCCC00", "#993300", "#CC99CC", "#6699CC","#CCCCFF", "#FFCC99", "#9966FF", "#CC6600", "#CCFFFF", "#99CC00", "#FF99FF", "#0066FF", "#66FFCC", "#99CCFF", "#9999CC", "#CC9900", "#CC33FF", "#006699", "#F5DF16", "#B5185E", "#99FF00", "#00FFFF", "#990000", "#CC0000", "#33CCCC", "#CC6666", "#996600", "#9999FF", "#3366FF")
CLUST.COL <- c(CLUST.COL,CLUST.COL)
coln=NULL
for(i in 1:39) coln[i] = lighten(CLUST.COL[i])
w = c(6,15,16,19,21,23)
coln[w]=CLUST.COL[w]
CLUST.COL2 = coln
CLUST.COL2[1]="lightpink"
CLUST.COL2[3]="lightgreen"
#CLUST.COL2[1]="pink1"
CLUST.COL2[4]="lavender"
CLUST.COL2[5]="darkorange"
CLUST.COL2[6]="dodgerblue"



# generate a barplot of the top-scoring overlap genes (driving enrichment) using google charts (can be
# visualized with the html img tag), bar length corresponds to score
# oGeneStr - |-separated string of overlap genes, and after each gene, its score is provided in parentheses,
#          genes are sorted by score in decreasing order
# scoreRange - a vector of the range of scores in the profile; 1st value = min,  2nd value = max
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly,
#          should they be in the top overlap genes
# barWidth - bar width in the plot
# maxGenes - maximum number of top overlap genes to return
# plotCol - the colour of the bars; hexidecimal format without the '#' character
# scoreLabel - label for the score axis
# RETURNS a vector of plot info: plot width, plot height and the google chart URL for the plot
mygenOverlapGenePlot.gChart <- function(oGeneStr, scoreRange, barModGenes=NULL, barWidth=10, maxGenes=10, plotCol="BEBEBE", scoreLabel="Sensitivity") {
  oGenes <- unlist(strsplit(oGeneStr, "\\|"))
  genes <- strsplit(oGenes, "\\(")
  scores <- sapply(genes, function(vec) { vec[length(vec)] })
  genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
  scores <- unlist(strsplit(scores, ")"))

  scoreMat <- data.frame(percent=rep(NA, length(genes)), score=as.numeric(scores), stringsAsFactors=F)
  rownames(scoreMat) <- genes
  if (nrow(scoreMat) > maxGenes) {
    scoreMat <- scoreMat[1:maxGenes, ]
  }

  if (!is.null(barModGenes)) {
    scoreMat$mod <- as.numeric(scoreMat$gene %in% barModGenes)
  }

  mygenLeadingEdgePlot.gChart(scoreMat, plotCol, scoreRange, barWidth, scoreLabel)
}


# generates overlap gene (i.e. those that drive enrichment) barplots for each cluster
# enrichInfo - dataframe with info on enriched terms (one per row), with the following columns:
#              overlapGenes (genes that drive the enrichment of the term), cluster (cluster of the term)
# scoreLabel - score label to use in the barplots
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly in the
#              barplots, should they be in the top overlap genes
# scoreInfo - if provided (i.e. not NULL), a data.frame of scores with the following columns:
#              gene (gene IDs) and score (gene scores)
# maxGenes - maximum number of top overlap genes to show in a barplot
# posScores - if TRUE, enrichment is amongst genes with positive scores; if FALSE, negative scores
# plotCol - the colour of the bars, in hexidecimal format without the '#' character
# RETURNS a matrix of plot info where each column corresponds to a different cluster bar plot,
# and the rows correspond to plot width (1st), plot height (2nd) hand the google chart URL for the plot (3rd)
clusterPlots <- function(enrichInfo, scoreLabel, barModGenes=NULL, scoreInfo=NULL, maxGenes=10, posScores=T, plotCol="BEBEBE") {
  barWidth <- 10

  overlapGenes <- strsplit(enrichInfo$overlapGenes, "\\|")

  if (is.null(scoreInfo)) {
    oGenes <- unique(unlist(overlapGenes))
    genes <- strsplit(oGenes, "\\(")
    scores <- sapply(genes, function(vec) { vec[length(vec)] })
    genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
    scores <- unlist(strsplit(scores, ")"))
    scoreInfo <- data.frame(gene=genes, score=as.numeric(scores), geneStr=oGenes, stringsAsFactors=F)
  }
  else {
    scoreInfo$geneStr <- scoreInfo$gene
  }

  scoreInfo <- scoreInfo[order(scoreInfo$score, decreasing=posScores), ]
  uniGenes <- unique(scoreInfo$gene)
  i <- match(uniGenes, scoreInfo$gene)
  scoreInfo <- scoreInfo[i, ]

  clusters <- split(1:nrow(enrichInfo), enrichInfo$cluster)

  plotData <- lapply(clusters, function(clustI) {
    clusterGenes <- table(unlist(overlapGenes[clustI]))
    clusterGenes <- clusterGenes/length(clustI) * 100

    i <- match(names(clusterGenes), scoreInfo$geneStr)
    clusterGenes <- cbind(clusterGenes, scoreInfo$score[i])
    colnames(clusterGenes) <- c("% of gene sets", "score")
    rownames(clusterGenes) <- scoreInfo$gene[i]

    # sort the genes by score, then by % of gene sets they are in
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]
    orderI <- order(clusterGenes[, 1], decreasing=T)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (nrow(clusterGenes) > maxGenes) {
      clusterGenes <- clusterGenes[1:maxGenes, , drop=F]
    }

    # re-order by score
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (is.null(barModGenes)) {
      return(clusterGenes)
    }

    cbind(clusterGenes, rownames(clusterGenes) %in% barModGenes) })

  if (posScores) {
    scoreRange <- c(0, ceiling(max(scoreInfo$score, na.rm=T)))
  }
  else {
    scoreRange <- c(floor(min(scoreInfo$score, na.rm=T)), 0)
  }

  sapply(plotData, mygenLeadingEdgePlot.gChart, plotCol, scoreRange, barWidth, scoreLabel)
}

# generate a barplot of the common leading edge genes using google charts (can be visualized
# with the html img tag), bar length corresponds to score
# leadInfo - matrix/dataframe of leading edge genes; rownames = gene IDs, column 1 = % of gene sets, column 2 = score
#          - if a 3rd column is provided, it should provide TRUE/FALSE indicating whether or not the gene
#            should be marked
# plotCol - the colour of the bars; hexidecimal format without the '#' character
# scoreRange - a vector of the range of scores in the profile; 1st value = min,  2nd value = max
# barWidth - bar width in the plot
# scoreLabel - label for the score axis
# RETURNS a vector of plot info: plot width, plot height and the google chart URL for the plot
mygenLeadingEdgePlot.gChart <- function(leadInfo, plotCol, scoreRange, barWidth, scoreLabel="Sensitivity") {
  # express bar lengths as values in [0, 100]
  #
  # ggadded aug3
  max = max(scoreRange,na.rm = T)
  min = min(scoreRange,na.rm = T)

  dataRange <- max - min
  #dataRange <- scoreRange[2] - scoreRange[1]
  barLens <- round((leadInfo[, 2] - scoreRange[1])/dataRange * 100)

  # determine the score step size
  if (dataRange <= 1) {
    stepSize <- 0.5
  }
  else if (dataRange <= 5) {
    stepSize <- 2
  }
  else if (dataRange <= 20) {
    stepSize <- 5
  }
  else if (dataRange <= 50) {
    stepSize <- 10
  }
  else if (dataRange <= 100) {
    stepSize <- 20
  }
  else if (dataRange <= 500) {
    stepSize <- 100
  }
  else {
    stepSize <- 250
  }

  # replace any spaces in scoreLabel with a +
  scoreLabel <- unlist(strsplit(scoreLabel, ""))
  scoreLabel[scoreLabel == " "] <- "+"
  scoreLabel <- paste(scoreLabel, collapse="")

  # determine the positions of the gene labels
  tmpStep <- 100/length(barLens)
  labelPos <- round(seq(tmpStep/2, 100, by=tmpStep))

  # compute plot size
  w <- 150
  h <- barWidth*length(barLens) + 50

  # specify the zero line if have negative values
  if (scoreRange[1] < 0) {
    zeroLineStr <- paste("&chp=", round(abs(scoreRange[1])/dataRange, 2), sep="")
  }
  else {
    zeroLineStr <- ""
  }

  # if a 3rd column is in leadInfo, use it to determine which gene bars to mark
  if (ncol(leadInfo) > 2 && any(leadInfo[, 3] == 1)) {
    barModStr <- paste("o,000000,", which(leadInfo[,3]==1) - 1, ",-1,5", sep="")
    barModStr <- paste("&chm=", paste(barModStr, collapse="|"), sep="")
  }
  else {
    barModStr <- ""
  }

  c(w, h, paste("http://chart.apis.google.com/chart?chxt=x,x,y&chs=", w, "x", h,
    "&cht=bhg&chd=t:", paste(barLens, collapse="|"),
    "&chco=", plotCol,
    "&chxl=1:|", scoreLabel, "|2:|", paste(rev(rownames(leadInfo)), collapse="|"),
    "&chxp=1,50|2,", paste(labelPos, collapse=","),
    "&chxr=0,", scoreRange[1], ",", scoreRange[2], ",", stepSize,
    "&chbh=", barWidth, ",1,0",
    zeroLineStr, barModStr, sep=""))
}
###
# computes enrichment using the hypergeometric test, and uses the resulting P values with
# the Benjamini Hochberg method to estimate FDR values
# querySet - character vector of genes in query set
# geneSets - named list of gene sets to test for significant overlap w/ the query set
# scoreMat - dataframe of gene scores
#          - first column = scores, gene column
#          - can be NULL
# uni - character vector of genes in the universe (i.e. background set)
#     - if NULL, must specify uniSize
# uniSize - the # of genes in the universe
# minSetSize, maxSetSize - min/max # of genes in geneSets (after restricting to the gene universe)
# RETURNS a dataframe of enrichment results, sorted by increasing FDR value. The columns are:
#         term = name of gene set
#         querySetFraction = the fraction of the query set that overlaps with the term set
#         geneSetFraction = the fraction of the term set that overlaps with the query set
#         foldEnrichment = the fold enrichment of the query set with the term genes
#         P = P value estimating the significance with which the query set is enriched with the term genes
#         FDR = FDR value estimating the significance of enrichment
#         overlapGenes = a |-separated list of genes in the overlap of the query set and the term set;
#                        if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#	  maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
#####hiphop:::hyperG
myhyperG = function (querySet, geneSets, uni, scoreMat, minSetSize = 1, 
  maxSetSize = 300, uniSize = NA) 
{
  if (!is.null(uni)) {
    geneSets <- lapply(geneSets, intersect, uni)
    lens <- sapply(geneSets, length)
    geneSets <- geneSets[lens >= minSetSize & lens <= maxSetSize]
    uniSize <- length(uni)
  }
  if (!is.null(scoreMat)) {
    scoreMat <- scoreMat[order(scoreMat$score, decreasing = T), 
      ]
    if (!is.null(uni)) {
      i <- match(uni, scoreMat$gene)
      scoreMat <- scoreMat[sort(i[!is.na(i)]), ]
    }
    scoreMat$score <- round(scoreMat$score, 2)
  }
  enrichInfo <- sapply(geneSets, function(geneSet) {
    overlapSet <- intersect(querySet, geneSet)
    pVal <- phyper(length(overlapSet) - 1, length(geneSet), 
      uniSize - length(geneSet), length(querySet), lower.tail = F)
    if (length(overlapSet) > 0) {
      overlapSet <- sort(overlapSet)
    }
    overlapSize <- length(overlapSet)
    if (is.null(scoreMat)) {
      maxScore <- NA
    }
    else {
      i <- sort(match(overlapSet, scoreMat$gene))
      maxScore <- scoreMat$score[i[1]]
      overlapSet <- paste(scoreMat$gene[i], "(", scoreMat$score[i], 
        ")", sep = "")
    }
    overlapSet <- paste(overlapSet, collapse = "|")
    bgRate <- length(geneSet)/uniSize
    foldEnrich <- overlapSize/length(querySet)/bgRate
    c(overlapSet, overlapSize/length(geneSet), foldEnrich, 
      pVal, maxScore, overlapSize/length(querySet))
  })
  enrichInfo <- t(enrichInfo)
  enrichCol <- data.frame(term = names(geneSets), querySetFraction = as.numeric(enrichInfo[, 
    6]), geneSetFraction = as.numeric(enrichInfo[, 2]), foldEnrichment = as.numeric(enrichInfo[, 
      3]), P = as.numeric(enrichInfo[, 4]), FDR = p.adjust(as.numeric(enrichInfo[, 
        4]), method = "BH"), overlapGenes = enrichInfo[, 1], 
    maxOverlapGeneScore = as.numeric(enrichInfo[, 5]), stringsAsFactors = F)
  rownames(enrichCol) <- NULL
  enrichCol = enrichCol[order(enrichCol$FDR), ]
}
#######hiphop:::overlapCoeff
#######the overlap of genesets for all combinations
#######If set X is a subset of Y or the converse then the overlap coefficient is equal to 1.
# compute the overlap coefficient given a pair of (gene) sets
# gsPairList - a list of two sets (each set is a vector of IDs)
# RETURNS the overlap coefficient
myoverlapCoeff = function (gsPairList) 
{
  length(intersect(gsPairList[[1]], gsPairList[[2]]))/min(length(gsPairList[[1]]), 
    length(gsPairList[[2]]))
}
######
mygenesNotInEnrichedTerm = function (queryGeneSets, enrichMat, scoreMat, termsToExclude, 
  fdrThresh = 0.1) 
{
  scoreMat <- as.matrix(scoreMat)
  enrichMat <- enrichMat[!(enrichMat$term %in% termsToExclude), 
    , drop = F]
  lens <- sapply(queryGeneSets, length)
  queryGeneSets <- queryGeneSets[lens > 0]
  oGenes <- strsplit(enrichMat$overlapGenes, "\\|")
  oGenes <- lapply(oGenes, function(genes) {
    genes <- strsplit(genes, "\\(")
    sapply(genes, function(vec) {
      vec[1]
    })
  })
  rowI <- split(1:nrow(enrichMat), enrichMat$filename)
  enrichI <- match(names(queryGeneSets), names(rowI))
  extraGenes <- queryGeneSets[is.na(enrichI)]
  queryGeneSets <- queryGeneSets[!is.na(enrichI)]
  rowI <- rowI[enrichI[!is.na(enrichI)]]
  tmp <- lapply(1:length(queryGeneSets), function(expI) {
    setdiff(queryGeneSets[[expI]], unlist(oGenes[rowI[[expI]]]))
  })
  names(tmp) <- names(queryGeneSets)
  extraGenes <- c(extraGenes, tmp)
  lens <- sapply(extraGenes, length)
  extraGenes <- extraGenes[lens > 0]
  if (length(extraGenes) > 0) {
    lens <- lens[lens > 0]
    extraGenes <- data.frame(filename = rep(names(extraGenes), 
      lens), gene = unlist(extraGenes), stringsAsFactors = F)
    i <- match(extraGenes$gene, rownames(scoreMat))
    i <- cbind(i, match(extraGenes$filename, colnames(scoreMat)))
    extraGenes$score <- round(scoreMat[i], 2)
    extraGenes <- extraGenes[order(extraGenes$score, decreasing = T), 
      ]
    i <- split(1:nrow(extraGenes), extraGenes$filename)
    extraGenes <- lapply(i, function(curRow) {
      tmp <- paste(extraGenes$gene[curRow], "(", extraGenes$score[curRow], 
        ")", sep = "")
      c(extraGenes$score[curRow[1]], paste(tmp, collapse = "|"))
    })
  }
  tmp <- lapply(1:length(queryGeneSets), function(expI) {
    curRow <- rowI[[expI]]
    sigI <- curRow[enrichMat$FDR[curRow] <= fdrThresh]
    unenrichedGenes <- setdiff(queryGeneSets[[expI]], unlist(oGenes[sigI]))
    curRow <- setdiff(curRow, sigI)
    if (length(curRow) == 0) {
      return(list(rowI = NULL, unenrichedGenes = NULL))
    }
    unenrichedGenes <- lapply(oGenes[curRow], function(genes) {
      intersect(unenrichedGenes, genes)
    })
    lens <- sapply(unenrichedGenes, length)
    unenrichedGenes <- unenrichedGenes[lens > 0]
    curRow <- curRow[lens > 0]
    if (length(curRow) == 0) {
      return(list(rowI = NULL, unenrichedGenes = NULL))
    }
    expI <- match(enrichMat$filename[curRow[1]], colnames(scoreMat))
    unenrichedGenes <- lapply(unenrichedGenes, function(curGenes) {
      geneI <- match(curGenes, rownames(scoreMat))
      geneStr <- scoreMat[geneI, expI]
      names(geneStr) <- curGenes
      geneStr <- round(sort(geneStr, decreasing = T), 2)
      geneStr <- paste(names(geneStr), "(", geneStr, ")", 
        sep = "")
      paste(geneStr, collapse = "|")
    })
    list(rowI = curRow, unenrichedGenes = unenrichedGenes)
  })
  unenrichedMat <- enrichMat[unlist(lapply(tmp, function(ob) {
    ob$rowI
  })), ]
  unenrichedMat$unenrichedGenes <- unlist(lapply(tmp, function(ob) {
    ob$unenrichedGenes
  }))
  if (length(extraGenes) > 0) {
    unenrichedMat <- unenrichedMat[c(rep(1, length(extraGenes)), 
      1:nrow(unenrichedMat)), ]
    toDoI <- 1:length(extraGenes)
    unenrichedMat$filename[toDoI] <- names(extraGenes)
    unenrichedMat$term[toDoI] <- "OTHER"
    unenrichedMat$overlapGenes[toDoI] <- sapply(extraGenes, 
      function(vec) {
        vec[2]
      })
    unenrichedMat$maxOverlapGeneScore[toDoI] <- as.numeric(sapply(extraGenes, 
      function(vec) {
        vec[1]
      }))
    unenrichedMat$unenrichedGenes[toDoI] <- unenrichedMat$overlapGenes[toDoI]
    if (!is.null(unenrichedMat$pruneOutcome)) {
      unenrichedMat$pruneOutcome[toDoI] <- "OTHER"
    }
    naCol <- setdiff(colnames(unenrichedMat), c("filename", 
      "term", "overlapGenes", "maxOverlapGeneScore", "unenrichedGenes"))
    colI <- match(naCol, colnames(unenrichedMat))
    unenrichedMat[toDoI, colI] <- NA
  }
  rownames(unenrichedMat) <- NULL
  unenrichedMat <- unenrichedMat[order(unenrichedMat$geneSetFraction, 
    decreasing = T), ]
  unenrichedMat[order(unenrichedMat$maxOverlapGeneScore, decreasing = T), 
    ]
}
####
# mygetNodes = function(xml){
#   library(dplyr)
#   library(xml2)
#   df <- bind_rows(lapply(xml_find_all(xml, "//node"), function(x) {
#     
#     # extract the attributes from the parent tag as a data.frame
#     parent <- data.frame(as.list(xml_attrs(x)), stringsAsFactors=FALSE)
#     
#     # make a data.frame out of the attributes of the kids
#     kids <- bind_rows(lapply(xml_children(x), function(x) as.list(xml_attrs(x))))
#     #kids2 = bind_rows(lapply(xml_children(xml_children(x)),function(x) as.list(xml_attrs(x))))
#     # combine them
#     cbind.data.frame(parent, kids, stringsAsFactors=FALSE)
#   }))
#   w = which(is.na(df$value))
#   if(length(w) > 0) df2 = df[-w,] else (df2 = df)
#   lens = length(table(df2$id))
#   spl = split(df2$value,df2$name)
#   w = which(sapply(spl,length) == lens)
#   if(length(w) > 0) spl = spl[w]
#   df3 = do.call(data.frame,spl)
#   w = which(df3$formattedLabel%in%df2$value)
#   m = match(df3$formattedLabel[w],df2$value)
#   df3$label[w]=df2$label[m]
#   df3$id[w]=df2$id[m]
#   ###assumes you got overlap genes
#   nrows = ncol(df3)
#   df3 = df3[,c(nrows,nrows-1,1:(nrows-2))]
#   s = lapply(df3,as.character)
#   w = which(names(s) %in% c("FDR","querySetFraction","geneSetFraction","nGenes","size"))
#   s[w]= lapply(s[w],as.numeric)
#   ds= data.frame(s,stringsAsFactors = F)
#   ds
# }


# mygetEdges = function(xml){
#   library(tidyr)
#   library(xml2)
#   df <- bind_rows(lapply(xml_find_all(xml, "//edge"), function(x) {
#     
#     # extract the attributes from the parent tag as a data.frame
#     parent <- data.frame(as.list(xml_attrs(x)), stringsAsFactors=FALSE)
#     
#     # make a data.frame out of the attributes of the kids
#     kids <- bind_rows(lapply(xml_children(x), function(x) as.list(xml_attrs(x))))
#     #kids2 = bind_rows(lapply(xml_children(xml_children(x)),function(x) as.list(xml_attrs(x))))
#     # combine them
#     cbind.data.frame(parent, kids, stringsAsFactors=FALSE)
#   }))
#   df2 = spread(df,name,value = value)
#   df2$source = as.numeric(df2$source)
#   df2$target = as.numeric(df2$target)
#   df2$width = as.numeric(df2$width)
#   df2
# }
####
mygenebarplot = function(overlapGenes){
  s = strsplit(overlapGenes,"\\|")
  s2 = lapply(s,strsplit,"\\(")
  s3 = lapply(s2,sapply,strsplit,"\\)")
  s4 = lapply(s3,"t")
  s5 = lapply(s4,as.data.frame,stringsAsFactors=F)
  s5 = lapply(s5,function(x) {
    names(x) = c("gene","score")
    x
  }
  )
  s6 = lapply(s5,function(x) {
    x$gene = unlist(x$gene)
    x$score = unlist(x$score)
    x$score = as.numeric(x$score)
    x = x %>% arrange(score)
    x
  }
  )
  nrows = sapply(s6,nrow)
  w = which(nrows > 10)
  if(length(w) > 0) {
    s6[w] = lapply(s6[w],function(x) x = x[1:10,])
  }
  s6
}
