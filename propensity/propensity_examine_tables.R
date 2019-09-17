
library('tidyverse')
library('magrittr')
library('BART')
library('ggplot2')
library('gridExtra')
library('spatstat')


#' Compute empirical cumulative distribution
#'
#' The empirical cumulative distribution function (ECDF) provides an alternative
#' visualisation of distribution. Compared to other visualisations that rely on
#' density (like [geom_histogram()]), the ECDF doesn't require any
#' tuning parameters and handles both continuous and categorical variables.
#' The downside is that it requires more training to accurately interpret,
#' and the underlying visual tasks are somewhat more challenging.
#'
#' @inheritParams layer
#' @inheritParams geom_point
#' @param na.rm If `FALSE` (the default), removes missing values with
#'    a warning.  If `TRUE` silently removes missing values.
#' @param n if NULL, do not interpolate. If not NULL, this is the number
#'   of points to interpolate with.
#' @param pad If `TRUE`, pad the ecdf with additional points (-Inf, 0)
#'   and (Inf, 1)
#' @section Computed variables:
#' \describe{
#'   \item{x}{x in data}
#'   \item{y}{cumulative density corresponding x}
#' }
#' @export
#' @examples
#' df <- data.frame(
#'   x = c(rnorm(100, 0, 3), rnorm(100, 0, 10)),
#'   g = gl(2, 100)
#' )
#' ggplot(df, aes(x)) + stat_ecdf(geom = "step")
#'
#' # Don't go to positive/negative infinity
#' ggplot(df, aes(x)) + stat_ecdf(geom = "step", pad = FALSE)
#'
#' # Multiple ECDFs
#' ggplot(df, aes(x, colour = g)) + stat_ecdf()
stat_ecdf <- function(mapping = NULL, data = NULL,
                      geom = "step", position = "identity",
                      weight =  NULL, 
                      ...,
                      n = NULL,
                      pad = TRUE,
                      na.rm = FALSE,
                      show.legend = NA,
                      inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatEcdf,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      n = n,
      pad = pad,
      na.rm = na.rm,
      weight = weight,
      ...
    )
  )
}

StatEcdf <- ggproto("StatEcdf", Stat,
                    compute_group = function(data, scales, weight, n = NULL, pad = TRUE) {
                      # If n is NULL, use raw values; otherwise interpolate
                      if (is.null(n)) {
                        x <- unique(data$x)
                      } else {
                        x <- seq(min(data$x), max(data$x), length.out = n)
                      }
                      
                      if (pad) {
                        x <- c(-Inf, x, Inf)
                      }
                      y <- ewcdf(data$x, weights=data$weight/sum(data$weight))(x)
                      
                      data.frame(x = x, y = y)
                    },
                    
                    default_aes = aes(y = stat(y)),
                    
                    required_aes = c("x")
)


standardized_difference_w <- function(data, treatmentname, variable, weightname=NULL) {
  treatment <- data[,treatmentname]
  w1<- sum(data[treatment,variable] * data[treatment,weightname], na.rm=TRUE)/sum(data[treatment,weightname], na.rm=TRUE)
  w2 <- sum(data[!treatment,variable] * data[!treatment,weightname], na.rm=TRUE)/sum(data[!treatment,weightname], na.rm=TRUE)
  s1 <- sum(data[treatment,weightname], na.rm=TRUE) / var(data[treatment,weightname], na.rm=TRUE) * sum( data[treatment,weightname] * (data[treatment,variable] - w1)^2, na.rm=TRUE)
  s2 <- sum(data[!treatment,weightname], na.rm=TRUE) / var(data[!treatment,weightname], na.rm=TRUE) * sum( data[!treatment,weightname] * (data[!treatment,variable] - w2)^2, na.rm=TRUE)

  return( ( w1 - w2 ) * sqrt(2) / sqrt(s1+s2) )
}

standardized_difference <- function(data, treatmentname, variable) {
    treatment <- data[,treatmentname]
return( ( mean(data[treatment,variable], na.rm=TRUE) - mean(data[!treatment,variable], na.rm=TRUE) ) * sqrt(2) / sqrt(var(data[treatment,variable], na.rm=TRUE)+var(data[!treatment,variable], na.rm=TRUE) ) )
}



## code for the population variance and inclusion rate (equiv to variance ==0) for each predictor -> table for propensity and predictive models
load(file="imputation_folders/icu_pop/1/imputed_with_propensity.Rdata")
 load(file="osa_data/pre_imputation_preop.Rdata" )


options(tibble.width = Inf)
options(dplyr.print_max = 100)
sink('osa_results/propensity_desc.txt')

## hisplay the propensity scores by treatment group


plot.data <- local.imputed
plot.data %<>% rename(propensity = predicted_osa)

plot.data %<>% mutate(prop_weight = new_osa/pmax( pmin( pnorm(propensity), .95) , .05) + (1-new_osa)/(1-pmax( pmin( pnorm(propensity), .95) , .05) ) )

## bring back some variables excluded from the imputation setup just to see if they balanced
# plot.data <- left_join(plot.data, filtered_cpap %>% select(one_of(c("PatientID", "BMI") ) ), by="PatientID")



print("what fraction of data by osa are within stable bounds (say .95 /.05) and would require clipping in IPTW?")

plot.data %>% group_by(new_osa) %>% summarize(z5=mean( abs(propensity) < qnorm(.95)), n5=sum( abs(propensity) < qnorm(.95)) , z10=mean( abs(propensity) < qnorm(.9)), n10=sum( abs(propensity) < qnorm(.9))) %>% print(.)

print(" mean weight in each group")
plot.data%>% group_by(new_osa) %>% summarize( w=mean(new_osa/pnorm(propensity) + (1-new_osa)/(1-pnorm(propensity))) ) %>% print(.)

print(' how much treated data is outside decent overlap? ')
# N(x>a) / N(y>a) > 3

treated_ecdf <- plot.data %>% filter(new_osa) %>% select(propensity) %>% unlist(.) %>% ecdf(.)
untreated_ecdf <- plot.data %>% filter(!new_osa) %>% select(propensity) %>% unlist(.) %>% ecdf(.)

n_treated <- sum(plot.data$new_osa)
n_untreated <- sum(!plot.data$new_osa)

interesting_x <- seq(from = -qnorm(.95), to= qnorm(.95), length.out=100)

count_ratio <- (1.-treated_ecdf(interesting_x)) * n_treated / (1.-untreated_ecdf(interesting_x)) / n_untreated

# plot(interesting_x, count_ratio)

# par(mfrow=c(1,2))
# plot(interesting_x, count_ratio)
# p

print(paste("fraction of data below where treated start to outnumber untreated by 5 or 3", treated_ecdf(interesting_x[which.min(count_ratio < 5)]) %>% round(.,3) ,
treated_ecdf(interesting_x[which.min(count_ratio < 3)]) %>% round(.,3) ))

# predict(osa_model, irmi_imputed)

## picture of overlap
  ## there is pretty good overlap
  
 p<-ggplot(plot.data, aes(x=propensity, fill=new_osa, color=new_osa)) +
  geom_histogram(position="identity", alpha=0.5)

  ## the units aregument doesn't work - the width and height are actually in pixels
  ggsave(file="osa_results/propensity_overlap.png", plot=p, units="in", device=png, dpi=300, width=600, height=600, limitsize = FALSE)


  
## loop through standardized differences on binary and continuous variables

binary_variables <- c("HTN", "SEX" , "CAD" , "CAD_PRIORMI" , "CHF" , "PPM_ICD" , "CV_TIA_STROKE" , "PAD" , "DVT" , "PE" , "CKD" , "DM" , "PHTN" , "COPD" , "ASTHMA"  , "CIRRHOSIS" , "CANCER_HX" , "GERD" , "ANEMIA" , "COOMBS_POS" , 'CHF_Diastolic_Function', "DEMENTIA" , 'AFIB' , 'Outpatient_Insulin')
continuous_variables <- c("Age" , "WEIGHT" ,  "BMI" , "CCI"  , "LVEF"  , "Neck" ,  "PreOp_Diastolic" , "PreOp_Systolic", "rsi_1", "urban_lodds", "Total.", "employed_lodds",  "poverty_fraction", "emergency") #, "StopBang_Total" , "StopBang_Observed" , "StopBang_Pressure" , "StopBang_Snore" , "StopBang_Tired"
ordered_set <- c('FUNCTIONAL_CAPACITY', 'ASA' , 'VALVULAR_DISEASE' , "Dialysis_History", "case_year")
multinomial_variables <- c("Surg_Type"  , "RACE"  , "SMOKING_EVER" , 'PAP_Type')
multinomial_variables <- setdiff(multinomial_variables, binary_variables)
multinomial_variables <- setdiff(multinomial_variables, ordered_set)
multinomial_variables <- setdiff(multinomial_variables, multinomial_variables)
plot.data %<>% mutate_at( ordered_set, as.numeric)
plot.data %<>% mutate_at( binary_variables, as.numeric)

std_outs <- matrix(NA, nrow=length(c(binary_variables, continuous_variables)), ncol=2)
for( this.var in seq_along(c(binary_variables, continuous_variables) ) ) {

  local.name <- c(binary_variables, continuous_variables)[this.var]
  std_outs[this.var,1] <- standardized_difference(data=plot.data, treatmentname="new_osa", variable=local.name)
  std_outs[this.var,2] <- standardized_difference_w(data=plot.data, treatmentname="new_osa", variable=local.name, weightname='prop_weight')

}
rownames(std_outs) <- c(binary_variables, continuous_variables)
std_outs[,1] <- round(std_outs[,1], 2)
std_outsC <- matrix(NA, nrow=length(c(binary_variables, continuous_variables)), ncol=2)
std_outsC[,1] <- as.character(std_outs[,1])
std_outsC[,2] <- formatC(std_outs[,2], format = "e", digits = 1)
rownames(std_outsC) <- c(binary_variables, continuous_variables)
print("standardized differences in covariates before and after propensity weighting")
print(std_outsC)


## pictures for continuous, ordered, multinomial  

name_vector <- c(ordered_set, continuous_variables, multinomial_variables)
picture_holder <- vector('list', length(name_vector)*2)
for( this.var in seq_along(name_vector ) ) {
picture_holder[[(this.var*2 -1) ]] <- ggplot(plot.data, aes_string(name_vector[this.var],  weight = NULL, color='new_osa'))  +   stat_ecdf()
picture_holder[[(this.var*2) ]] <- ggplot(plot.data, aes_string(name_vector[this.var],  weight = 'prop_weight', color='new_osa'))  +   stat_ecdf()
}
ml <- marrangeGrob(picture_holder, nrow=4, ncol=2, layout_matrix = matrix(seq_len(4 *2), nrow = 4, ncol = 2, byrow=TRUE))
ggsave("osa_results/ecdf_propensity.pdf", ml)


if(FALSE) {

## MCMC convergence diagnostics: they are grossly bad
## on inspection of longer runs, the acf is low (~7%) at 50*75 -> 100k samples = 27 independent obs
## the effectiveSize ratio based on ar0 is about .07
## however, the z's from independent chains are also way overdispersed (IQR about 4.6), suggesting that the ar0 tends to underestimate the variance of the sample estimates(overestimate the neff or underestimate the marginal variance)
## ar0 based estimates are known to behave badly in many cases
## it is the case that single-sample probs will mix much more slowly than composites from the entire sample; they often will have only a few variables that can change and f(x_i) doesn't have the clt over sample effect with at most 1 tree changing per sample
## it's a little strange that the 
pdf()
    i <- floor(seq(1, 10000, length.out=50))
    auto.corr <- acf(osa_model$yhat.train[seq(nrow(osa_model$yhat.train)/3) , i], plot=FALSE, lag.max=100)
    max.lag <- max(auto.corr$lag[ , 1, 1])
    

    these_acf <- matrix(NA, nrow=length(i), ncol=dim(auto.corr$acf)[1] )
    for( h in seq_along(i)) {
    these_acf[h,] <- auto.corr$acf[, h, h]
    }
    
    ## the acf function is noticing a legit slow start in eg case 1225 - is it related to how low the baseline prob is? if i start it after iter 2000 it's normal, so maybe it's just a burn-in thing. need to take a look if it's ALWAYS a beginning thing and if it's similar in the different chains.
    ## other low values don't necessarily do that
    ## mayble also it's a function of the sparse options and lowish N trees -> no relevant variables changing (it does have a long slow curve) 
    tacf <- acf(osa_model$yhat.train[seq(from=2000, to=nrow(osa_model$yhat.train)/3) , i[7] , drop=FALSE])
    
    j <- seq(-0.5, 0.4, length.out=10)
    for(h in 1:10) {
        if(h==1)
            plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
                 type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
                 ylab='acf', xlab='lag')
        else
            lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
                 type='h', col=h)
    }
dev.off()

# tacf <- acf(osa_model$yhat.train[seq(from=1, to=nrow(osa_model$yhat.train)/3) , i[10] , drop=FALSE], lag.max=100)
# tacf <- acf(osa_model$yhat.train[seq(to=nrow(osa_model$yhat.train), from=nrow(osa_model$yhat.train)*2/3) , i[10] , drop=FALSE], lag.max=100, plot=FALSE)
# tacf <- acf(osa_model$yhat.train[seq(from=nrow(osa_model$yhat.train)/3, to=nrow(osa_model$yhat.train)*2/3) , i[10] , drop=FALSE], lag.max=100, plot=FALSE)
# tacf <- acf(osa_model$yhat.train[seq(to=nrow(osa_model$yhat.train)/3, from=1) , i[10] , drop=FALSE], lag.max=100, plot=FALSE)


pdf()
plot(osa_model$yhat.train[,i[10]])
abline(v=nrow(osa_model$yhat.train)/3 * seq(2), col='red')
dev.off()

pdf()
plot(osa_model$yhat.train[,i[15]])
abline(v=nrow(osa_model$yhat.train)/3 * seq(2), col='red')
dev.off()


    geweke <- gewekediag(osa_model$yhat.train)
    temp_data <- gewekediag(osa_model$yhat.train[seq(to=nrow(osa_model$yhat.train)/3, from=1) ,] , frac1=.2, frac2=.2)
    temp_ar<-spectrum0ar(osa_model$yhat.train[seq(nrow(osa_model$yhat.train)/3), ])$spec
    deflation <-  apply(osa_model$yhat.train[seq(nrow(osa_model$yhat.train)/3), ] ,2, var) / temp_ar
    summary(deflation)
    summary(temp_data$z)
    k <- 50
    
    mean(abs(temp_data$z) > qnorm(.95)) 

plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',)

mean(abs(geweke$z) > qnorm(.95)) 
mean(abs(geweke$z) > qnorm(.99)) 
mean(abs(geweke$z) > qnorm(.999)) 

dev.off()


## something similar for cpap


}
