### Title: Functional response function for crop planting analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes:

# dat_df <- dat_op$dat_clean
# scale_cov_name <- "crp.nfsp.sc"
# beta_long <- beta_long_all
# s_long <- s_long_all
# covsx <- dat_op$covsx
# x_incr<-0.5
# spatial=TRUE

fun_response <- function(dat_df, #dat_clean or dat_clean_op
                         scale_cov_name, #scaled covariate column name from dat_df
                         covsx, #covariate matrix, in order fit to modek
                         beta_long,#pivoted mcmc samples for beta coefficients
                         s, #samples[,"s"] for spatial random effect
                         # tau,#vector of mcmc samples for tau regression error
                         x_incr=0.5,#'by' argument for creating sequence of scaled x's
                         spatial=FALSE){ #logical to include intercept in calculation (FALSE for spatial RE's)
  require(tidyverse)
  
  dat_df <- as.data.frame(dat_df)
  anom_cov <- data.frame(sc=seq(range(dat_df[,scale_cov_name])[1],
                                 range(dat_df[,scale_cov_name])[2],
                                by=x_incr))
  anom_cov$bt <- anom_cov$sc *attr(dat_df[,scale_cov_name], 'scaled:scale') + 
    attr(dat_df[,scale_cov_name], 'scaled:center')
  
  crp_cov <- data.frame(sc=seq(range(dat_df[,'crp.prop.sc'])[1],
                                range(dat_df[,'crp.prop.sc'])[2],
                                length=length(anom_cov$sc)))
  crp_cov$bt <- crp_cov$sc *attr(dat_df[,'crp.prop.sc'], 'scaled:scale') + 
    attr(dat_df[,'crp.prop.sc'], 'scaled:center')

  nfsp_cov <- data.frame(sc=seq(range(dat_df[,'prop.nfsp.sc'])[1],
                               range(dat_df[,'prop.nfsp.sc'])[2],
                               length=length(anom_cov$sc)))
  nfsp_cov$bt <- nfsp_cov$sc *attr(dat_df[,'prop.nfsp.sc'], 'scaled:scale') + 
    attr(dat_df[,'prop.nfsp.sc'], 'scaled:center')
  
  beta_cov <- beta_long %>% filter(cov==covsx$name[covsx$cov==scale_cov_name]) %>% select(value)
  beta_cov <- unlist(c(beta_cov),use.names = F)
  
  beta_crp <- beta_long %>% filter(cov==covsx$name[covsx$cov=='crp.prop.sc']) %>% select(value)
  beta_crp <- unlist(c(beta_crp),use.names = F)
  
  beta_nfsp <- beta_long %>% filter(cov==covsx$name[covsx$cov=='prop.nfsp.sc']) %>% select(value)
  beta_nfsp <- unlist(c(beta_nfsp),use.names = F)
  
  ncounties <- length(unique(dat_df$GEOID))
  
  s_long <- s %>%
    pivot_longer(cols=all_of(1:ncounties),values_to="value",names_to="s")
  
  s_long$county_idx <- rep(1:ncounties,nrow(s))
  s_long <- s_long %>% select(-s) %>% arrange(county_idx)

  s_long$tot_idx <- rep(1:(max(s_long$samp_idx)*max(s_long$chain_idx)),max(s_long$county_idx))

  if(spatial){
    # est_cov_sc <- sapply(1:length(s_long),function(s){
    #   sapply(1:length(beta_cov),function(i){
    #     beta_cov[i]*anom_cov$sc + s_long$value
    #   })
    # })
    
    est_cov_sc <- array(NA,dim=c(nrow(anom_cov),length(beta_cov),length(unique(s_long$county_idx))))
 
    if(scale_cov_name!="crp.nfsp.sc"){
      for(i in 1:nrow(s_long)){
        est_cov_sc[,s_long$tot_idx[i],s_long$county_idx[i]] <- 
          beta_cov[s_long$tot_idx[i]]*anom_cov$sc + s_long$value[i]
      }
    } else {
      for(i in 1:nrow(s_long)){
        est_cov_sc[,s_long$tot_idx[i],s_long$county_idx[i]] <- 
          beta_crp[s_long$tot_idx[i]]*crp_cov$sc +
          beta_nfsp[s_long$tot_idx[i]]*nfsp_cov$sc +
          beta_cov[s_long$tot_idx[i]]*anom_cov$sc + s_long$value[i] 
      }
    }

  } else {
    beta_int <- beta_long %>% filter(cov=="Intercept") %>% select(value)
    beta_int <- unlist(c(beta_int),use.names = F)
    
    est_cov_sc <- sapply(1:length(beta_cov),function(i){
      beta_int[i] + beta_cov[i]*anom_cov$sc
      # rnorm(nrow(anom_cov),beta_int[i] + beta_cov[i]*anom_cov$sc,sqrt(1/tau[i]))

    })
  }
  
  quant_lci <- function(i){quantile(i,probs=0.025)}
  quant_uci <- function(i){quantile(i,probs=0.975)}
  
  est_cov_sc[,1:10,1]
  
  mn = apply(est_cov_sc,c(1,3),mean)
  md = apply(est_cov_sc,c(1,3),median)
  lci=apply(est_cov_sc,c(1,3),quant_lci)
  uci=apply(est_cov_sc,c(1,3),quant_uci)
  
  mn_long <- mn %>% as_data_frame() %>% 
    mutate(cov_sc=anom_cov$sc,
           cov_bt=anom_cov$bt) %>% 
    pivot_longer(cols=all_of(1:length(unique(s_long$county_idx))),
                 values_to="mn",names_to="county_idx")
  
  md_long <- md %>% as_data_frame() %>% 
    mutate(cov_sc=anom_cov$sc,
           cov_bt=anom_cov$bt) %>% 
    pivot_longer(cols=all_of(1:length(unique(s_long$county_idx))),
                 values_to="md",names_to="county_idx")
  
  lci_long <- lci %>% as_data_frame() %>% 
    mutate(cov_sc=anom_cov$sc,
           cov_bt=anom_cov$bt) %>% 
    pivot_longer(cols=all_of(1:length(unique(s_long$county_idx))),
                 values_to="lci",names_to="county_idx")
  
  uci_long <- uci %>% as_data_frame() %>% 
    mutate(cov_sc=anom_cov$sc,
           cov_bt=anom_cov$bt) %>% 
    pivot_longer(cols=all_of(1:length(unique(s_long$county_idx))),
                 values_to="uci",names_to="county_idx")
  
  est_cov_sum <- mn_long%>% left_join(md_long) %>% left_join(lci_long) %>% left_join(uci_long)
  # est_cov_sum <- cbind.data.frame(anom_cov,est_cov_sum)
  return(est_cov_sum)
}
