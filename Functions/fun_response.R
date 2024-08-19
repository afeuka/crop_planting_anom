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
                         beta_long){ #pivoted mcmc samples for beta coefficients
  require(tidyverse)
  
  dat_df <- as.data.frame(dat_df)
  anom_cov <- data.frame(sc=seq(range(dat_df[,scale_cov_name])[1],
                                 range(dat_df[,scale_cov_name])[2],by=0.5))
  anom_cov$bt <- anom_cov$sc *attr(dat_df[,scale_cov_name], 'scaled:scale') + 
    attr(dat_clean$take.hog.intens.sc, 'scaled:center')

  beta_int <- beta_long %>% filter(cov=="Intercept") %>% select(value)
  beta_int <- unlist(c(beta_int),use.names = F)
  
  beta_cov <- beta_long %>% filter(cov==covsx$name[covsx$cov==scale_cov_name]) %>% select(value)
  beta_cov <- unlist(c(beta_cov),use.names = F)
  est_cov_sc <- sapply(1:length(beta_cov),function(i)beta_int[i] + beta_cov[i]*anom_cov$sc)
  est_cov_sum <- data.frame(mn=rowMeans(est_cov_sc),
                             lci=sapply(1:nrow(est_cov_sc),function(i)quantile(est_cov_sc[i,],probs=0.025)),
                             uci=sapply(1:nrow(est_cov_sc),function(i)quantile(est_cov_sc[i,],probs=0.975)))
  est_cov_sum <- cbind.data.frame(anom_cov,est_cov_sum)
  return(est_cov_sum)
}
