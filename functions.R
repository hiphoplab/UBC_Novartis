


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

mybarheight = function(leadInfo){
  
  scoreRange = range(leadInfo$score)
  
  dataRange <- scoreRange[2] - scoreRange[1]
  
  barLens <- round((leadInfo[, 2] - scoreRange[1])/dataRange *
      100)
  w <- 150
  barWidth = 15
  h <- barWidth * length(barLens) + 50
  h
}

mybar = function(df){
  
  scoreRange = range(df$score)
  
  dataRange <- scoreRange[2] - scoreRange[1]
  
  barLens <- round((df[, 4] - scoreRange[1])/dataRange *
                     100)
  w <- 150
  barWidth = 15
  h <- barWidth * length(barLens) + 50
  h
}
