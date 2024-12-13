### Title: Model validation for multi-chain crop analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes: 

library(nimble)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

source("./crop_planting_anom/Functions/clean_crop_dat.R")

nb=FALSE
temporal=FALSE

if(nb==FALSE & temporal==FALSE){
  subfolder <- "Farm Bill Sites No NB No Temp"
} else if (nb==FALSE & temporal==TRUE){
  subfolder <- "Farm Bill Sites No NB"
} else if (nb==TRUE & temporal==FALSE){
  subfolder <- "Farm Bill Sites No Temp"
} else {
  subfolder <- "Farm Bill Sites"
}


load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

count_all <- read.csv("./Data/county_state_names.csv")
count_fb <- read.csv("./Data/FB_Projects_Counties_long.csv")
count_fb <- count_fb %>% 
  rename(ST_ABBR=State,
         CNTY_NAME=County) %>% 
  select(-X)
count_fb$FB <-TRUE
count_all <- count_all %>% 
  mutate(CNTY_NAME=toupper(CNTY_NAME))
count_all <- count_all %>% left_join(count_fb)
count_all$FB[is.na(count_all$FB)] <- FALSE

dat_orig_fb <- dat_orig %>% left_join(count_all %>% 
                                        select(-CNTY_NAME) %>% #use geo id to join, not names (. and spaces cause differences)
                                        mutate(state_name=toupper(ST_NAME),
                                               GEOID=factor(GEOID))) %>% 
  filter(FB==TRUE)

# commod_idx <- 1
pval<- numeric()
r2_v <- numeric()
for(commod_idx in 1:length(commod_names_c)){
  if(commod_idx!=4){
    # only counties with pigs -------------
    ## load samples -----------------
    load(paste0("./Model outputs/",subfolder,"/",
                commod_names[commod_idx],"_op_region_multi.RData"))
    
    ## load data------------
    dat_op  <- clean_crop_dat(dat_orig=dat_orig_fb,
                              nb=nb,
                              temporal=temporal,
                              commod_name = commod_names[commod_idx], 
                              only_pigs = T)
    
    dat_clean_op <- dat_op$dat_clean
    
    ##model checks ---------------------
    ypred <- samples_all[,grep("ypred",colnames(samples_all))]
    ypred_sum <- data.frame(mn=colMeans(ypred),
                            md=apply(ypred,2,median),
                            lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
                            uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
                            obs = dat_clean_op$plant.anom)
    
    ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
    ypred_sum$std_resid <- scale(ypred_sum$resid)
    
    ypred_sum = cbind(ypred_sum,dat_clean_op)
    ypred_sum <- ypred_sum %>%
      left_join(dat_clean_op %>% dplyr::select(GEOID) %>% distinct())
    
    ###posterior predictive distribution ----------------
    ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
                                      values_to="value",
                                      names_to="typ"))+
      geom_histogram(aes(x=value,fill=typ),alpha=0.5,col="transparent")+
      scale_fill_discrete(name="",labels=c("Post pred",
                                           "Data"))+
      ggtitle(paste(commod_names_t[commod_idx]," - Posterior predictions - Counties with pigs"))
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/data_post_dist_op.jpeg"),
           device="jpeg",height=5,width=7,units="in")
    
    ##simulated model stats --------------------
    ypred_stat <- data.frame(mn=rowMeans(ypred),
                             sd=apply(ypred,1,sd),
                             idx=1:nrow(ypred))
    ggplot()+
      geom_histogram(data=ypred_stat,aes(x=mn))+
      geom_vline(xintercept=mean(dat_clean_op$plant.anom),col="red",lty=2)+
      ggtitle(paste(commod_names_t[commod_idx]," - Posterior data mean - Counties with pigs"))
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/pval_mn_op.jpeg"),
           device="jpeg",width=7,height=5,units="in")
    ggplot()+
      geom_histogram(data=ypred_stat,aes(x=sd))+
      geom_vline(xintercept=sqrt(var(dat_clean_op$plant.anom)),col="red",lty=2)+
      ggtitle(paste(commod_names_t[commod_idx]," - Posterior data SD - Counties with pigs"))
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/pval_sd_op.jpeg"),
           device="jpeg",width=7,height=5,units="in")
    
    ### bayesian p-value ---------------
    tau <- cbind.data.frame(value=samples_all[,"tau"],
                            chain_idx=samples_all$chain_idx,
                            samp_idx=samples_all$samp_idx)
    tau$sd <- sqrt(1/tau$value)
    
    mu <- samples_all[,grep("mu",colnames(samples_all))]
    mu <- as.matrix(mu)
    ypred <- as.matrix(ypred)
    
    
    dat_ll <- sapply(1:nrow(mu),function(i){
      sum(dnorm(dat_clean_op$plant.anom,mu[i,],sd=tau$sd[i],log=T))
    })
    pred_ll <- sapply(1:nrow(mu),function(i){
      sum(dnorm(ypred[i,],mu[i,],sd=tau$sd[i],log=T))
    })
    
    pval[commod_idx] <- sum(dat_ll>pred_ll)/nrow(mu)
    
    #r2 ------------
    ssr <- sum(ypred_sum$resid^2)
    sst <- sum((ypred_sum$obs-mean(ypred_sum$obs))^2)
    r2_v[commod_idx] <- 1- ssr/sst
    
    ggplot(ypred_sum)+
      geom_point(aes(x=obs,y=md))+
      geom_abline(intercept=0,slope = 1)+
      ggtitle(paste("R2=",round(r2_v[commod_idx],2)))
    
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/r2_plot_op.jpeg"),
           device="jpeg",width=7,height=5,units="in")
    
    # ### standardized residuals--------------------
    # ggplot(ypred_sum)+
    #   geom_point(aes(x=mn,y=std_resid,size=long),alpha=0.7)+
    #   geom_hline(yintercept=0,col="blue",lty=2)
    # 
    # covs_outliers <- ypred_sum %>% filter(std_resid<=-5) %>% 
    #   pivot_longer(cols=c(grep("sc",colnames(ypred_sum)),
    #                       grep("lat",colnames(ypred_sum)),
    #                       grep("long",colnames(ypred_sum))),
    #                names_to="cov",values_to="value") %>% 
    #   group_by(cov) %>% 
    #   summarise(min=min(value),
    #             max=max(value),
    #             typ="outlier")
    # 
    # covs_all <- ypred_sum %>% 
    #   pivot_longer(cols=c(grep("sc",colnames(ypred_sum)),
    #                       grep("lat",colnames(ypred_sum)),
    #                       grep("long",colnames(ypred_sum))),
    #                names_to="cov",values_to="value") %>% 
    #   group_by(cov) %>% 
    #   summarise(min=min(value),
    #             max=max(value),
    #             typ="all")
    # 
    # covs_resid <- full_join(covs_outliers,covs_all) 
    # 
    # ggplot(covs_resid)+
    #   geom_errorbar(aes(x=cov,ymin=min,ymax=max,col=typ),
    #                 position = position_dodge(width=0.5))+
    #   theme(axis.text.x=element_text(angle=90,hjust=1))
  }
}
pval_df <- data.frame(Crop=commod_names_t,pval,r2_v)
write.csv(pval_df,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Model outputs/",subfolder,"/Plots/Combined Figures/pval_table.csv"))

