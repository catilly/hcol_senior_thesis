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
# This script trains XP-EN in all of the EUR data + 87.5% of the YRI data;
# it predicts in the remaining 12.5%.
# It can be used with top SNP, PRS, lasso, EN, or BLUP, but we focus mainly on combined EN in our work.

# Initialize parameters and vectors to store corresponding scores
method = PAR.METHOD
window.size = 1000000 # Specify window size
rarity = "all" # If "common," SNPs with MAF <1% are filtered out before training

crosspop.corr = rep(NA, nrow(quant))
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
  
  # Cross-population accuracy, causal observed
  # Sample some of the pop2 individuals, then combine with pop1 to use as a training set
  train2.ix = sample.int(n = nrow(df2), size = floor(.875*nrow(df2)), replace = F)
  train = rbind(df1, df2[train2.ix,])
  test = df2[-train2.ix,]
  
  model = fit.gene.expression.model(train, method = method, nfolds = 5)
  crosspop.corr[i] = score.gene.expression(model, test, metric = "corr", lambda = "lambda.min")
  crosspop.MSE[i] = score.gene.expression(model, test, metric = "MSE", lambda = "lambda.min")
  }

output <- data.frame(corrs.to.test = quant$TargetID,
                     crosspop.corr = crosspop.corr,
                     crosspop.MSE = crosspop.MSE)
write.table( output , quote=F )

q()
