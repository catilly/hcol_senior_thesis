library(snpStats)
library(dplyr)
source("simulate_expression_prediction_v2.R")

arg = commandArgs(trailingOnly=T)
QUANT.DATA = toString(arg[1])
PAR.OUTDIR = toString(arg[2])

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
# Ready files for use in XP-BLUP!
# Assumes we will train XP-BLUP in 87.5% of the YRI data and predict in the remaining 12.5%.
# Also outputs files needed to conduct association testing in EUR.

# Initialize parameters and vectors to store corresponding scores
window.size = 1000000 # Specify window size
rarity = "all" # If "common," SNPs with MAF <1% are filtered out before training

# Loop through the probe targets in our specified chromosome,
# generating a SNP list, a pop1 phenotype file,
# a pop2 training phenotype file, and a pop2 testing phenotype file

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
  
  # Write locus SNP list
  locus.SNP.list = colnames(X)
  write.table( locus.SNP.list ,
               file = file.path(PAR.OUTDIR,
                                paste("locus.SNP.list.chr", quant[1,]$Chr,
                                      ".TargetID_", quant[i,]$TargetID,
                                      ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
  
  # Merge the expression and genotype data; format
  df = transform(merge(X, E, by = "row.names", all.x = TRUE), row.names=Row.names, Row.names=NULL)
  df <- df %>% select(E, everything())
  
  # Subset by population
  df1 = df[which(pop[,2] != "YRI"),]
  df2 = df[which(pop[,2] == "YRI"),]
  
  # Pop 1 phenotypes (used for defining C1 SNPs)
  pop1.pheno = data.frame(row.names(df1), row.names(df1), df1[,1])
  write.table( pop1.pheno ,
               file = file.path(PAR.OUTDIR, paste("pop1.pheno.chr", quant[1,]$Chr,
                                                  ".TargetID_", quant[i,]$TargetID,
                                                  ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
  
  # Pop 2 phenotypes (used for training/testing)
  sample <- sample.int(n = nrow(df2), size = floor(.875*nrow(df2)), replace = F)
  pop2.train.pheno = data.frame(row.names(df2[sample, ]), row.names(df2[sample, ]), df2[sample, ][,1])
  write.table( pop2.train.pheno ,
               file = file.path(PAR.OUTDIR, paste("pop2.train.pheno.chr", quant[1,]$Chr,
                                                  ".TargetID_", quant[i,]$TargetID,
                                                  ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
  
  pop2.test.pheno = data.frame(row.names(df2[-sample, ]), row.names(df2[-sample, ]), df2[-sample, ][,1])
  write.table( pop2.test.pheno ,
               file = file.path(PAR.OUTDIR, paste("pop2.test.pheno.chr", quant[1,]$Chr,
                                                  ".TargetID_", quant[i,]$TargetID,
                                                  ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
}

