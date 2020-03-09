#library(shiny)
#library(shinydashboard)
#library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(DT)


p10 <- function(matI,x,sig = 1,ylab = "fitness defect",xlab = "gene",las = 2,font = 3,cex = 1.2,cex.main = 1.8,cex.axis = 1.2,cex.lab=1.5,...)  {
  col = ifelse (matI[, x] > sig, 1,      ifelse (matI[, x] < -sig, 3, 2))
  w <- which(matI[, x]  > sig | matI[, x]  < -sig, arr.ind = T)
  palette(mycolors)
  posw = ifelse (w > nrow(matI) - 0.1 * nrow(matI) ,
    2,
    ifelse (w <  nrow(matI) + 0.1 * nrow(matI), 4 ,  2))
  plot(
    matI[, x],
    col = col ,cex.lab = 1.5,
    main = paste(colnames(matI)[x]),
    ylab = ylab, xaxt = "n",
    xlab = xlab,
    las = las,
    cex.axis = cex.axis,
    cex.main = cex.main,
    ...
  )
  if (length(w != 0))
    text(
      w,
      matI[w, x],
      names(w),
      pos = posw,
      cex = cex,
      font = font,
      ...
    )
  abline(
    h = sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
  abline(
    h = -sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
}

###########
# myxcor = function(gene,mat,method = "spearman") {
#   ord = order(mat[gene,])
#   xcor<- cor(mat[gene,ord],t(mat[,ord]), use="pairwise.complete.obs", method=method)
#   xcor = t(xcor)
#   xcor = xcor[order(xcor[,1],decreasing = T),]
#   xcor
# }

###########
myplotpts = function(mat,rows,group,ylab = "fitness defect",xlab = "", main = rows,cex = 1.2,cex.lab=1.5,cex.main=1.8,font=3,...){
  color = function(x) {
    col = as.numeric(as.factor(x))
    col
  }
  palette(mycolors)
  group = factor(group)
  lev = levels(group)
  plot(mat[rows,]~group,las=2,ylab = ylab,xlab=xlab, main = main,cex.lab=cex.lab,cex=cex,font=font,cex.main=cex.main,...)
  points(mat[rows,]~group,col=color(group),cex.lab = cex.lab,cex=cex,...)
  mat[rows,]
}
mycolors = c(
  "darkorange",
  "dodgerblue",
  "limegreen",
  "navy",
  "mediumpurple"  ,
  "royalblue3",
  "darkolivegreen4",
  "firebrick",
  "cyan",
  "hotpink",
  "plum4",
  "blue",
  "magenta2",
  "skyblue",
  "green",
  "red",
  "steelblue",
  "tomato",
  "purple",
  "gold",
  "steelblue",
  "red",
  "lightpink",
  "purple",
  "darkred",
  "lightblue"
)


# setnicepar = function(mar = c(3, 3, 2, 1), mgp = c(2, 0.4, 0),
#   tck = -.01, cex.axis = 0.9,
#   las = 1, mfrow = c(1, 1), ...) {
#   par(mar = mar, mgp = mgp, tck = tck,
#     cex.axis = cex.axis, las = las,
#     mfrow = mfrow, ...)
# }

#setnicepar(pch = 19,mfrow = c(1,1),pch = 19,las = 2)
#par(mar = c(3,3,2,1))
###########
mymeltdf = function(mat, row, df = phiphop){
  require(reshape2)
  
  mx = melt(mat[row,], value.name = "fitness_defect")
  mx$gene = row
  mx$screen = rownames(mx)
  
  mx = mx[,c("gene","screen","fitness_defect")]
  m = match(mx$screen,df$name)
  table(is.na(m))
  mx$id = df$id[m]
  mx$drug = df$drug[m]
  mx$dose = df$conc[m]
  mx$name = df$name[m]
  mx$drug_conc = df$drug_conc[m]
  
  #mx$type = df$drug[m]
  #mx$signature = df$name[m]
  
  #mx$shape = ifelse(mx$hit == 1,17,19)
  mx
  
}
###########
myplot1cnt = function(mat,vect,row,refer = "ctrl", cex = 0.4,main = row,...){
  g = grep(refer,vect)
  if(is.factor(vect)) vect = as.character(vect)
  lens = length(g)
  names(vect) = colnames(mat)
  fac = c(vect[g],sort(vect[-g]))
  palette(mycolors)
  r = rle(as.character(fac))
  
  med = tapply(mat[row,names(fac)],fac,median)
  m = match(r$values,names(med))
  med = med[m]
  names(med) = r$values
  lens = r$lengths
  cum = cumsum(r$lengths)
  d = c(1,cum+1)
  colg = as.numeric(as.factor(fac))
  ucolg = unique(colg)
  plot(mat[row,names(fac)],col = colg,cex = cex,pch = 17,xaxt = "n",ylab = "counts",xlab = "",main = main,...)
  for(i in 1:length(d)-1) segments(d[i],med[i],cum[i],med[i],col = ucolg[i],lwd=3)
  lens = rle(fac)$lengths
  lens2 = cumsum(lens)
  factor = relevel(factor(fac),ref = "ctrl")
  
  axis(
    1,
    at = lens2 - lens/2,
    labels =  levels(factor),
    cex = 0.7, las = 2
  )
  med
}
