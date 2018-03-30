library(snpStats)
library(dplyr)
source("simulate_expression_prediction_v2.R")

arg = commandArgs(trailingOnly=T)
QUANT.DATA = toString(arg[1])
PAR.METHOD = toString(arg[2])

set.seed( 1 )

#################################################################
# Load data

# Load in quantification data (assumed to be all from the same chromosome!)
quant <- read.table(QUANT.DATA, header = TRUE)

# Load in data from relevant chromosome to speed things up
plink <- read.plink(paste0("./Data/Genotype/1000G.all.deduped.", quant[1,]$Chr, ".bed"),
                    paste0("./Data/Genotype/1000G.all.deduped.", quant[1,]$Chr, ".bim"),
                    paste0("./Data/Genotype/1000G.all.deduped.", quant[1,]$Chr, ".fam"))
pop = read.table("./Data/Genotype/1000G.POP") # Read in population data

#################################################################
# Predict!
# This script trains models on only the EUR data and predicts in all of the YRI data.
# It also performs cross-validation to obtain within-EUR accuracy.
# It is suitable only for the top SNP, PRS, BLUP, lasso, and EN methods.

# Initialize parameters and vectors to store corresponding scores
method = PAR.METHOD
window.size = 1000000 # Specify window size
k.flds = 5 # Folds to use for k-fold CV
rarity = "all" # If "common," SNPs with MAF <1% are filtered out before training

withinpop.cv.corr = rep(NA, nrow(quant))
crosspop.corr = rep(NA, nrow(quant))

withinpop.cv.MSE = rep(NA, nrow(quant))
crosspop.MSE = rep(NA, nrow(quant))

# Loop through genes
for (i in 1:nrow(quant)){
  
  # Grab the relevant expression values for
  E <- t(quant[i, -c(1:4)]); colnames(E) <- "E"
  
  # Extract SNPs that are near the gene of interest
  X = extract.X.from.plink(plink, quant[i,]$Coord - window.size/2, quant[i,]$Coord + window.size/2)
  
  # If there's only one SNP, append a column of zeros so that glmnet won't freak out
  while ((ncol(X) < 2) || (is.null(ncol(X)))){
    X = cbind(X, 0)
  }
  
  # Check if we want only the common variants in population 1
  if (rarity == "common"){
    MAF = apply( X[which(pop[,2] != "YRI"),], 2, function(x) min( mean(x-1)/2, 1-mean(x-1)/2 ) )
    SNP.subset <- which(MAF > 0.01) # Note that we define common as >1% MAF here
    
    X <- X[, SNP.subset, drop = FALSE]
    
    # Check again if there are enough SNPs
    while ((ncol(X) < 2) || (is.null(ncol(X)))){
      X = cbind(X, 0)
    }
  }
  
  # Merge the expression and genotype data; format
  df = transform(merge(X, E, by = "row.names", all.x = TRUE), row.names=Row.names, Row.names=NULL)
  df <- df %>% select(E, everything())
  
  # Subset by population
  df1 = df[which(pop[,2] != "YRI"),]
  df2 = df[which(pop[,2] == "YRI"),]
  
  # Use k-fold CV to evaluate accuracy
  # https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
  folds <- cut(seq(1, nrow(df1)), breaks = k.flds, labels = FALSE)
  
  # Shuffle data
  df1 = df1[sample(nrow(df1)),]
  
  corr = rep(NA, k.flds)
  MSE = rep(NA, k.flds)
  
  for (k in 1:k.flds){
    
    test.indexes = which(folds == k, arr.ind = TRUE)
    
    test = df1[test.indexes,]; train = df1[-test.indexes,]
    model = fit.gene.expression.model(train, method = method)
    corr[k] = score.gene.expression(model, test, metric = "corr")
    MSE[k] = score.gene.expression(model, test, metric = "MSE")
  }
  
  # Average across the k folds
  withinpop.cv.corr[i] = mean(na.omit(corr))
  withinpop.cv.MSE[i] = mean(na.omit(MSE))
  
  # Last set of models to get cross-population accuracy
  # Cross-population accuracy, causal observed
  test = df2; train = df1
  model = fit.gene.expression.model(train, method = method)
  crosspop.corr[i] = score.gene.expression(model, test, metric = "corr")
  crosspop.MSE[i] = score.gene.expression(model, test, metric = "MSE")
}

output <- data.frame(corrs.to.test = quant$TargetID,
                     withinpop.cv.corr = withinpop.cv.corr,
                     crosspop.corr = crosspop.corr,
                     withinpop.cv.MSE = withinpop.cv.MSE,
                     crosspop.MSE = crosspop.MSE)
write.table( output , quote=F )

q()
