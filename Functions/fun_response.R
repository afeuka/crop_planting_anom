### Title: Functional response function for crop planting analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes:

# dat_df <- dat_clean
# scale_cov_name <- "take.hog.intens.sc"
# beta_long <- beta_long_all
# covsx <- covsx_all

fun_response <- function(dat_df, #dat_clean or dat_clean_op
                         scale_cov_name, #scaled covariate column name from dat_df
                         covsx, #covariate matrix, in order fit to modek
                         beta_long,#pivoted mcmc samples for beta coefficients
                         # tau,#vector of mcmc samples for tau regression error
                         x_incr=0.5,#'by' argument for creating sequence of x's
                         intercept=FALSE){ #logical to include intercept in calculation (FALSE for spatial RE's)
  require(tidyverse)
  
  dat_df <- as.data.frame(dat_df)
  anom_cov <- data.frame(sc=seq(range(dat_df[,scale_cov_name])[1],
                                 range(dat_df[,scale_cov_name])[2],by=x_incr))
  anom_cov$bt <- anom_cov$sc *attr(dat_df[,scale_cov_name], 'scaled:scale') + 
    attr(dat_df[,scale_cov_name], 'scaled:center')
  
  beta_cov <- beta_long %>% filter(cov==covsx$name[covsx$cov==scale_cov_name]) %>% select(value)
  beta_cov <- unlist(c(beta_cov),use.names = F)
  
  if(intercept){
    
    beta_int <- beta_long %>% filter(cov=="Intercept") %>% select(value)
    beta_int <- unlist(c(beta_int),use.names = F)
    
    est_cov_sc <- sapply(1:length(beta_cov),function(i){
      beta_int[i] + beta_cov[i]*anom_cov$sc
      # rnorm(nrow(anom_cov),beta_int[i] + beta_cov[i]*anom_cov$sc,sqrt(1/tau[i]))
    })
  } else {
    
    est_cov_sc <- sapply(1:length(beta_cov),function(i){
      beta_cov[i]*anom_cov$sc
      # rnorm(nrow(anom_cov),beta_cov[i]*anom_cov$sc,sqrt(1/tau[i]))
    })
  }

  est_cov_sum <- data.frame(mn=rowMeans(est_cov_sc),
                             lci=sapply(1:nrow(est_cov_sc),function(i)quantile(est_cov_sc[i,],probs=0.025)),
                             uci=sapply(1:nrow(est_cov_sc),function(i)quantile(est_cov_sc[i,],probs=0.975)))
  est_cov_sum <- cbind.data.frame(anom_cov,est_cov_sum)
  return(est_cov_sum)
}
