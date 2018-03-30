library(snpStats)

extract.X.from.plink = function(plink, pos.start, pos.end){
  # Converts genotype output from snpStats::read.plink() to a data frame.
  # Args:
  #   plink: Object loaded using snpStats::read.plink().
  #   pos.start, pos.end: Start and end positions of the block from which SNPs should be extracted.
  # Returns:
  #   X: Genotype matrix formatted as a data frame.
  #   Rows are individuals and columns are SNPs.
  #   1, 2, 3 indicate 00, 01/10, and 11 genotypes respectively; 0 indicates a missing value.
  
  # Filter for SNPs that are between >= pos.start and < pos.end
  genotypes = plink$genotypes[, which((plink$map[,'position'] >= pos.start) & (plink$map[,'position'] < pos.end))]
  
  # Convert to data frame
  X = matrix(as.numeric(genotypes@.Data), nrow = nrow(genotypes), ncol = ncol(genotypes))
  X = data.frame(X)
  rownames(X) <- rownames(genotypes@.Data); colnames(X) <- colnames(genotypes@.Data)
  
  return(X)
}


simulate.gene.expression = function(X1, X2, hsq1, hsq2, corr, n.causal, hide.causal = FALSE, diff.causal = FALSE, sample.causal.from.subset = NA){
  # Appends simulated gene expression levels to genotype data frames from two populations.
  # Args:
  #   X1, X2: Genotype matrices from populations 1 and 2.
  #   hsq1, hsq2: Heritability of gene expression for populations 1 and 2.
  #   corr: Correlation between expression effect sizes between populations 1 and 2.
  #   n.causal: Number of SNPs that will contribute to expression.
  #   hide.causal: If TRUE, removes the SNPs used to generate expression from the output data frame.
  #   diff.causal: If TRUE, different casual variants are selected for populations 1 and 2.
  #   sample.causal.from.subset: If provided, the index of the causal variant
  #   will be sampled from only this set of indices.
  # Returns:
  #   list(df1, df1.hidden, df2, df2.hidden, causal.ix1, causal.ix2): List of data frames for each population.
  #   Also returns the column index of the (non-hidden) data frames in which the causal variant is located
  #   If hide.causal is TRUE, include data frames where the causal variant is hidden.
  #   The first column of each data frame contains the simulated expression levels.
  #
  # NOTE: Setting hsq1 = hsq2 assumes that within-pop genotypic variance has been scaled to be the same.
  
  library(MASS)
  
  N1 = nrow(X1); N2 = nrow(X2) # Sample size
  M = ncol(X2) # Number of markers (SNPs)
  
  # Generate expression effect sizes for n.causal SNPs
  b = mvrnorm(n = n.causal, mu = c(0, 0), Sigma = matrix(c(1, corr, corr, 1), 2, 2))
  dim(b) <- c(n.causal, 2) # Ensure output is a matrix, not a vector
  
  # While loop: we resample if we end up picking a SNP with no variance in one of the populations
  Eg1 = rep(NA, N1); Eg2 = rep(NA, N2)
  while (anyNA(Eg1) || anyNA(Eg2)){
    
    # Generate genetic component of expression
    if (anyNA(sample.causal.from.subset)){
      # Randomly choose n.causal SNPs
      causal.ix1 = sample(ncol(X1), size = n.causal, replace = FALSE)
    } else {
      # Sample n.causal SNPs from a provided subset of indices
      causal.ix1 = sample(sample.causal.from.subset, size = n.causal, replace = FALSE)
    }
    
    if (diff.causal == TRUE){
      causal.ix2 = sample(ncol(X1), size = n.causal, replace = FALSE)
    } else {
      causal.ix2 = causal.ix1
    }
    
    Eg1 = scale(as.matrix(X1[, causal.ix1]) %*% b[, 1])
    Eg2 = scale(as.matrix(X2[, causal.ix2]) %*% b[, 2])
  }
  
  # Generate overall expression by adding random noise
  E1 = scale(Eg1 * sqrt(hsq1) + rnorm(N1, 0, sqrt(1 - hsq1)))
  E2 = scale(Eg2 * sqrt(hsq2) + rnorm(N2, 0, sqrt(1 - hsq2)))
  
  # If hide.causal, delete the causal SNPs after having used them to generate expression
  if (hide.causal == TRUE){
    X1.hidden = X1[, -causal.ix1]
    X2.hidden = X2[, -causal.ix2]
    
    # Concatenate simulated expression levels to input dataframes
    df1 = as.data.frame(cbind(E1, X1)); colnames(df1)[1] <- "E"
    df1.hidden = as.data.frame(cbind(E1, X1.hidden)); colnames(df1.hidden)[1] <- "E"
    df2 = as.data.frame(cbind(E2, X2)); colnames(df2)[1] <- "E"
    df2.hidden = as.data.frame(cbind(E2, X2.hidden)); colnames(df2.hidden)[1] <- "E"
    
    return(list(df1, df1.hidden, df2, df2.hidden, causal.ix1 + 1, causal.ix2 + 1))
    
  } else {
    # Concatenate simulated expression levels to input dataframes
    df1 = as.data.frame(cbind(E1, X1)); colnames(df1)[1] <- "E"
    df2 = as.data.frame(cbind(E2, X2)); colnames(df2)[1] <- "E"
    
    return(list(df1, df2, causal.ix1 + 1, causal.ix2 + 1)) 
  }
}

fit.gene.expression.model = function(df.train, method = NA, nfolds = 5, C1 = numeric(0)){
  # Fits model that predicts gene expression using genotype, given training data.
  # Args:
  #   df.train: Training data frame, where first column records expression levels
  #   and the rest of the data frame is a genotype matrix.
  #   method: One of "multiple," "marginal", "PRS", "BLUP", "lasso", "EN," or "XPEN."
  #   nfolds: Default 10. Specifies number of folds to be used for CV when fitting cv.glmnet models.
  #   C1: If using the XPEN method, this specifies the set of SNPs that belong in C1
  #       (i.e., they were significant predictors in the reference population).
  # Returns:
  #   model: Fitted model.
  
  if (method == "multiple") {
    
    # Fit multiple regression
    f <- paste(colnames(df.train)[1], "~ .")
    model = do.call("lm", list(as.formula(f), data = as.name("df.train")))
    
  } else if (method == "marginal") {
    
    # Take dot product of gene expression with all other columns to get correlation (since they are all scaled)
    predictors = scale(as.matrix(df.train[, 2:ncol(df.train)]))
    corrs = t(df.train[, 1]) %*% predictors
    
    # Fit single linear regression with SNP that has greatest (absolute) correlation with gene expresion
    f <- paste(colnames(df.train)[1], "~", colnames(df.train)[which.max(abs(corrs)) + 1])
    model = do.call("lm", list(as.formula(f), data=as.name("df.train")))
    
    ## Initialize high p-value
    #best.p = 1
    
    ## Loop through all columns in X.train to select most significant SNP
    #for (snp in 2:ncol(df.train)){
    
    #f <- paste(colnames(df.train)[1], "~", colnames(df.train)[snp])
    #cur.model = do.call("lm", list(as.formula(f), data=as.name("df.train")))
    
    #if (summary(cur.model)$coefficients[,4][2] < best.p){
    #model = cur.model
    #best.p = summary(cur.model)$coefficients[,4][2]
    #}
    #}
    
  } else if (method == "random"){
    
    # Fit single linear regression with a random SNP
    rand.ix <- sample(2:ncol(df.train), size = 1)
    f <- paste(colnames(df.train)[1], "~", colnames(df.train)[rand.ix])
    model = do.call("lm", list(as.formula(f), data=as.name("df.train")))
    
  } else if (method == "PRS"){
    
    # Get OLS coefficients quickly (no intercept)
    # https://brenocon.com/blog/2012/03/cosine-similarity-pearson-correlation-and-ols-coefficients/
    predictors = as.matrix(df.train[, 2:ncol(df.train)])
    OLS.coef = (t(df.train[, 1]) %*% predictors) / apply(predictors, 2, function(x) sum(x^2))
    
    # Return OLS coefficients as the "model"
    model = OLS.coef
    
  } else if (method == "BLUP") {
    
    library(glmnet)
    model = cv.glmnet(as.matrix(df.train[, 2:ncol(df.train)]), df.train[, 1], nfolds = nfolds, family = "gaussian", alpha = 0)
    
  } else if (method == "lasso") {
    
    library(glmnet)
    model = cv.glmnet(as.matrix(df.train[, 2:ncol(df.train)]), df.train[, 1], nfolds = nfolds, family = "gaussian", alpha = 1)
    
  } else if (method == "EN") {
    
    library(glmnet)
    model = cv.glmnet(as.matrix(df.train[, 2:ncol(df.train)]), df.train[, 1], nfolds = nfolds, family = "gaussian", alpha = 0.5)
    
  } else if (method == "XPEN") {
    
    library(glmnet)
    # Set a very high initial cross-validated MSE
    cur.best.cvm = 1000000
    model = NA
    
    # First check if C1 is length 0 to save some computational time
    if (length(C1) == 0){
      model = cv.glmnet(as.matrix(df.train[, 2:ncol(df.train)]), df.train[, 1], nfolds = k.flds, family = "gaussian", alpha = 0.5,
                        lambda = exp(seq(log(0.001), log(5), length.out=100)))
    } else {
      for (scaling.factor in c(0.01, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000)){
        
        # Check which SNPs are in C1
        # We shrink the SNPs that aren't in C1 more by a factor of scaling.factor
        penalty = colnames(df.train) %in% C1
        penalty[penalty == 0] <- scaling.factor
        
        # Fit EN with alpha = 0.5
        cur.model = cv.glmnet(as.matrix(df.train[, 2:ncol(df.train)]), df.train[, 1], nfolds = k.flds, penalty.factor = penalty, family = "gaussian", alpha = 0.5,
                          lambda = exp(seq(log(0.001), log(5), length.out=100)))
        
        # Get model with the best scaling.factor
        if (min(cur.model$cvm) < cur.best.cvm){
          model = cur.model
          cur.best.cvm = min(cur.model$cvm)
        }
      }
    }
    
    # Sometimes no good model is selected, so in that case just fit the vanilla EN
    if (is.na(model)){
      model = cv.glmnet(as.matrix(train[, 2:ncol(train)]), train[, 1], nfolds = k.flds, family = "gaussian", alpha = 0.5,
                        lambda = exp(seq(log(0.001), log(5), length.out=100)))
    }
    
  } else {
    
    print("Specify prediction method.")
    
  }
  
  return(model)
}

score.gene.expression = function(model, df.test, metric = "MSE", lambda = "lambda.min"){
  # Evaluates performance of model on a testing dataset.
  # Args:
  #   model: Fitted model.
  #   df.test: Testing data frame, where first column records expression levels
  #   and the rest of the data frame is a genotype matrix.
  #   metric: One of "MSE", "corr", or "Rsq"
  #   lambda: One of "lambda.1se" or "lambda.min," to be passed into the predict function for cv.glmnet methods
  # Returns:
  #   score: Evaluation of model on test set using the specified metric.
  # Note that if predictions result in NAs, they're removed before scoring.
  # To compute R2, we use the square of the Pearson correlation between observed and predicted
  
  # Check class of model that was just passed in; predict
  if (class(model) == "cv.glmnet"){
    E.pred = predict(model, newx = as.matrix(df.test[, 2:ncol(df.test)]), s = lambda)
  } else if (class(model) == "matrix"){
    E.pred = as.matrix(df.test[, 2:ncol(df.test)]) %*% t(model) 
  } else
    E.pred = predict(model, newdata = df.test)
  
  if (metric == "MSE"){
    score = mean((E.pred - df.test[, 1])^2, na.rm = TRUE)
    
  } else if (metric == "corr"){
    score = cor(E.pred, df.test[, 1], use = "na.or.complete")
    if (is.na(score)){
      score <- 0
    }
  } else if (metric == "Rsq"){
    score = cor(E.pred, df.test[, 1], use = "na.or.complete")^2
    if (is.na(score)){
      score <- 0
    }
  } else
    print ("Specify scoring metric.")
  
  return(score)
}