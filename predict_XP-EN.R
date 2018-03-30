library(snpStats)
library(dplyr)
source("simulate_expression_prediction_v2.R")

arg = commandArgs(trailingOnly=T)
QUANT.DATA = toString(arg[1])
C1.LOC = toString(arg[2])

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
# This script trains XP-EN in 87.5% of the YRI data and predicts in the remaining 12.5%.

# Initialize parameters and vectors to store corresponding scores
window.size = 1000000 # Specify window size
k.flds = 5 # Folds to use for k-fold CV
rarity = "all" # If "common," SNPs with MAF <1% are filtered out before training

crosspop.corr = rep(NA, nrow(quant))
crosspop.MSE = rep(NA, nrow(quant))

best.scaling.factor = rep(NA, nrow(quant))

# Loop through genes
for (i in 1:nrow(quant)){
  
  # Grab the relevant expression values
  E <- t(quant[i, -c(1:4)]); colnames(E) <- "E"
  
  # Get C1
  C1 <- scan(paste0(C1.LOC, "/C1.chr", quant[i,]$Chr, ".TargetID_", quant[i,]$TargetID, ".txt"), what = "character")
  
  # Extract SNPs that are near the gene of interest
  pos.start = quant[i,]$Coord - window.size/2
  pos.end = quant[i,]$Coord + window.size/2
  X = extract.X.from.plink(plink, pos.start, pos.end)
  
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
  
  # Cross-population accuracy
  # Split pop2 into test and train
  train.ix = sample.int(n = nrow(df2), size = floor(.875*nrow(df2)), replace = F)
  train = df2[train.ix,]
  test = df2[-train.ix,]
  
  model = fit.gene.expression.model(train, method = "XPEN", C1 = C1)
  
  crosspop.corr[i] = score.gene.expression(model, test, metric = "corr", lambda = "lambda.min")
  crosspop.MSE[i] = score.gene.expression(model, test, metric = "MSE", lambda = "lambda.min")
}

output <- data.frame(corrs.to.test = quant$TargetID,
                     crosspop.corr = crosspop.corr,
                     crosspop.MSE = crosspop.MSE)
write.table( output , quote=F )

q()