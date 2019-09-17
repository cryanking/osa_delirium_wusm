

bcf_combined <- function(data, y_var="Y", treatment_var="T", multicore=FALSE, ...) {

  ## no missing data allowed
  if(any(is.na(data))) { stop("missing data not supported - consider imputing first")}
  
  ## detect the y variable type
  if(length(unique(data[,y_var]))) {var_type <- "binom"} else {var_type <- "cont"}
  
  
  ## replace this with turning treatment into an integer first (maybe factor -> least common)
  sort_order <- c( sort(which(data[,treatment_var]==1)) , sort(which(data[,treatment_var]==0)) )
  
  num_treated <- sum(treatment_var==1)
  
  data <- data[sort_order , ]
  
  y.train <- data[,y_var]
  t.train <- data[,treatment_var]
  cov.train <- data[,-c(treatment_var,y_var) ]
  qx <-  colMeans(pbart(y.train=t.train, x.train=cov.train 
                         , nskip = 300, ndpost=1000
                         , sparse=TRUE)$yhat.train	)
  antilogit <- function(x) {exp(x)/(1+exp(x))}
  qx_main <- antilogit(qx) 
  cov.train <- cbind(cov.train,qx_main )
  
  ## branching
  
  if(multicore) {
    if(var_type == 'binom') {
      res <- mc.pbart_bcf(x.train=cov.train,   y.train=y.train,  numtreated = num_treated, ...)
    } else {
      res <- mc.wbart_bcf(x.train=cov.train,   y.train=y.train,  numtreated = num_treated, ...)
    }
  
  } else {
    if(var_type == 'binom') {
      res <- pbart_bcf(x.train=cov.train,   y.train=y.train,  numtreated = num_treated, ...)
    } else {
      res <- wbart_bcf(x.train=cov.train,   y.train=y.train,  numtreated = num_treated, ...)
    }  
  
  }

## add a pscore matching estimate and plain bart estimate as well?
  return(res)
}
