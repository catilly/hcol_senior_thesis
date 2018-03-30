arg = commandArgs(trailingOnly=T)
PAR.FILE = toString(arg[1])
PAR.METHOD = toString(arg[2])

if (PAR.METHOD == "XPBLUP"){
  library(plyr)
  dat = read.table(PAR.FILE, header=TRUE, row.names=NULL)
  agg <- ddply(dat, .(CORR), function(df) cor(df$PHENO, df$SCORE)) # Outputs correlation!
  out <- data.frame(corrs.to.test = agg$CORR, variable = "crosspop.corr", value = agg$V1, method = PAR.METHOD)
} else {
  library(reshape2)
  dat = read.table(PAR.FILE, row.names=NULL)[,-1]
  dat.filt <- dat[, c(1, grep(".corr", colnames(dat)))]
  out <- melt(dat.filt, id.vars = "corrs.to.test")
  out$method = PAR.METHOD
}

write.table( out , row.names = FALSE, quote=F )

q()
