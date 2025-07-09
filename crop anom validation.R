### Title: Model validation for multi-chain crop analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes: 

library(nimble)
library(tidyverse)

# setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

source("./crop_planting_anom/Functions/clean_crop_dat.R")

mod_typ <- "spatial" #"region_multi"

nb=FALSE
temporal=TRUE

load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

commod_idx <- 1
pval <- r2_v <- samp_mn<-samp_sd<- numeric()
ypred_sum_list <- ypred_stat_list<-list()

for(commod_idx in 1:length(commod_names_c)){
  # all counties-------------
  dat <- clean_crop_dat(dat_orig=dat_orig,
                        nb=nb,
                        temporal=temporal,
                        commod_name = commod_names[commod_idx], only_pigs = F)
  
  dat_clean <- dat$dat_clean
  
  if("take.hog.intens.5yeartrend.sc"%in%colnames(dat_clean) & mod_typ!="spatial"){
    subfolder <- "Take Trend"
  } else if("take.hog.intens.prev.sc"%in%colnames(dat_clean) & mod_typ!="spatial") {
    subfolder <- "Take Previous" 
  } else if("crp.nfsp.sc"%in%colnames(dat_clean) & mod_typ=="spatial") {
    subfolder <- "CRPxPigs"
  } else if(mod_typ=="spatial") {
    subfolder <- "Spatial" 
  } else {
    stop("Must include take covariate. Check clean_crop_dat.R")
  }
  
  # 
  # ## load samples -----------------
  # load(paste0("./Model outputs/",subfolder,"/",
  #             commod_names[commod_idx],"_region_multi.RData"))
  # 
  # ##model checks ---------------------
  # ypred <- samples_all[,grep("ypred",colnames(samples_all))]
  # ypred_sum <- data.frame(mn=colMeans(ypred),
  #                         lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
  #                         uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
  #                         obs = dat_clean$plant.anom)
  # 
  # ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
  # ypred_sum$std_resid <- scale(ypred_sum$resid)
  # 
  # ypred_sum = cbind(ypred_sum,dat_clean)
  # ypred_sum <- ypred_sum %>%
  #   left_join(dat_clean %>% dplyr::select(GEOID) %>% distinct())
  # 
  # ###posterior predictive distribution ----------------
  # ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
  #                                   values_to="value",
  #                                   names_to="typ"))+
  #   geom_histogram(aes(x=value,fill=typ),alpha=0.5)+
  #   scale_fill_discrete(name="",labels=c("Post pred",
  #                                        "Data"))+
  #   ggtitle(paste(commod_names_t[commod_idx]," - Posterior predictions"))
  # ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/data_post_dist.jpeg"),
  #        device="jpeg",height=5,width=7,units="in")
  # 
  # ##simulated model stats --------------------
  # sd <- function(x){sqrt(var(x))}
  # ypred_stat <- data.frame(mn=rowMeans(ypred),
  #                          sd=apply(ypred,1,sd),
  #                          idx=1:nrow(ypred))
  # ggplot()+
  #   geom_histogram(data=ypred_stat,aes(x=mn))+
  #   geom_vline(xintercept=mean(dat_clean$plant.anom),col="red",lty=2)+
  #   ggtitle(paste(commod_names_t[commod_idx]," - Posterior data mean"))
  # ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/pval_mn.jpeg"),
  #        device="jpeg",width=7,height=5,units="in")
  # 
  # ggplot()+
  #   geom_histogram(data=ypred_stat,aes(x=sd))+
  #   geom_vline(xintercept=sqrt(var(dat_clean$plant.anom)),col="red",lty=2)+
  #   ggtitle(paste(commod_names_t[commod_idx]," - Posterior data SD"))
  # ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/pval_sd.jpeg"),
  #        device="jpeg",width=7,height=5,units="in")
  # 
  ### standardized residuals--------------------
  # ggplot(ypred_sum)+
  #   geom_point(aes(x=mn,y=std_resid,size=long),alpha=0.7)+
  #   geom_hline(yintercept=0,col="blue",lty=2)
  # 
  # covs_outliers <- ypred_sum %>% filter(abs(std_resid)<=5) %>% 
  #   pivot_longer(cols=c(grep("ever.pigs",colnames(ypred_sum)),
  #                       grep("sc",colnames(ypred_sum)),
  #                       grep("lat",colnames(ypred_sum)),
  #                       grep("long",colnames(ypred_sum))),
  #                names_to="cov",values_to="value") %>% 
  #   group_by(cov) %>% 
  #   summarise(min=min(value),
  #             max=max(value),
  #             typ="outlier")
  # 
  # covs_all <- ypred_sum %>% 
  #   pivot_longer(cols=c(grep("ever.pigs",colnames(ypred_sum)),
  #                       grep("sc",colnames(ypred_sum)),
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
  
  # rm(samples_all)
  
  # only counties with pigs -------------
  ## load samples -----------------
  load(paste0("./Model outputs/",subfolder,"/",
              commod_names[commod_idx],"_op_",mod_typ,".RData"))
  
  ## load data------------
  dat_op  <- clean_crop_dat(dat_orig=dat_orig,
                            nb=nb,
                            temporal=temporal,
                            commod_name = commod_names[commod_idx], only_pigs = T)
  
  dat_clean_op <- dat_op$dat_clean
  samp_mn[commod_idx] <- mean(dat_clean_op$plant.anom)
  samp_sd[commod_idx] <- sqrt(var(dat_clean_op$plant.anom))
  
  ##model checks ---------------------
  ypred <- samples_all[,grep("ypred",colnames(samples_all))]
  ypred_sum <- data.frame(mn=colMeans(ypred),
                          lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
                          uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
                          obs = dat_clean_op$plant.anom)

  ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
  ypred_sum$std_resid <- scale(ypred_sum$resid)
  
  ypred_sum = cbind(ypred_sum,dat_clean_op)
  ypred_sum <- ypred_sum %>%
    left_join(dat_clean_op %>% dplyr::select(GEOID) %>% distinct())
  
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
  
  ###posterior predictive distribution ----------------
  ypred_sum <- ypred_sum %>% pivot_longer(cols=c("mn","obs"),
                             values_to="value",
                             names_to="typ")
  ypred_sum$crop <- commod_names[commod_idx]
  ypred_sum_list[[commod_idx]] <- ypred_sum
  
  ##simulated model stats --------------------
  ypred_stat <- data.frame(mn=rowMeans(ypred),
                           sd=apply(ypred,1,sd),
                           idx=1:nrow(ypred))
  
  ypred_stat$crop <- commod_names[commod_idx]
  ypred_stat_list[[commod_idx]] <- ypred_stat
  
  #r2 ------------
  ssr <- sum(ypred_sum$resid^2)
  sst <- sum((ypred_sum$plant.anom-mean(ypred_sum$plant.anom))^2)
  r2_v[commod_idx] <- 1- ssr/sst
  
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
pval_df <- data.frame(Crop=commod_names_t,pval)
r2_df <- data.frame(Crop=commod_names_t,r2_v)
fit_df <- pval_df %>% left_join(r2_df)

# if(!dir.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Model outputs/",
#                       subfolder,"/Plots/Combined Figures"))){
#   dir.create(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Model outputs/",
#                          subfolder,"/Plots/Combined Figures"))
# }
if(!dir.exists(paste0("./Model outputs/",
                      subfolder,"/Plots/Combined Figures"))){
  dir.create(paste0("./Model outputs/",
                    subfolder,"/Plots/Combined Figures"))
}
write.csv(fit_df,paste0("./Model outputs/",
                         subfolder,"/Plots/Combined Figures/pval_r2_table.csv"))

#posterior predictive vs value 
ypred_stat <- do.call("rbind.data.frame",ypred_stat_list)
samp_mn <- data.frame(mn=samp_mn,crop=commod_names)
samp_sd <- data.frame(sd=samp_sd,crop=commod_names)
ggplot()+
  geom_histogram(data=ypred_stat,aes(x=mn))+
  geom_vline(data=samp_mn,aes(xintercept=mn),col="red",lty=2)+
  facet_wrap(.~crop)+
  ggtitle("Posterior data mean - Counties with pigs")
ggsave(filename=paste0("./Model outputs/",
                       subfolder,"/Plots/Combined Figures/pval_mn.jpeg"),
       device="jpeg",width=7,height=5,units="in")
ggplot()+
  geom_histogram(data=ypred_stat,aes(x=sd))+
  geom_vline(data=samp_sd,aes(xintercept=sd),col="red",lty=2)+
  facet_wrap(.~crop)+
  ggtitle("Posterior data SD - Counties with pigs")
ggsave(filename=paste0("./Model outputs/",
                       subfolder,"/Plots/Combined Figures/pval_sd.jpeg"),
       device="jpeg",width=7,height=5,units="in")


#posterior predictive vs data
ypred_sum <- do.call("rbind.data.frame",ypred_sum_list)
ggplot(ypred_sum)+
  geom_histogram(aes(x=value,fill=typ),alpha=0.5,col="transparent")+
  facet_wrap(.~crop)+
  scale_fill_discrete(name="",labels=c("Post pred",
                                       "Data"))+
  ggtitle("Posterior predictions - Counties with pigs")

ggsave(filename=paste0("./Model outputs/",
                       subfolder,"/Plots/Combined Figures/data_post_dist_op.jpeg"),
       device="jpeg",height=5,width=7,units="in")
