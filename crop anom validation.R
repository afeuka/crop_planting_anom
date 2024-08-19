### Title: Model validation for multi-chain crop analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes: 

library(nimble)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)
commod_idx <- 1

for(commod_idx in 1:length(commod_names_c)){
  # all counties-------------
  ## load samples -----------------
  load(paste0("./Model outputs/",commod_names[commod_idx],"_region_rs_multi.RData"))
  
  dat <- dat_orig %>% 
    filter(commodity_desc==commod_names_c[commod_idx]) %>% 
    filter(!is.nan(plant.anom) & 
             !is.na(plant.anom) & 
             !is.na(GEOID) & 
             year>=2009 
    )
  pig_cov <- "take.hog.intens" 
  
  covs_all <- data.frame(
    cov=c("plant.anom",
          "GEOID",
          "division_grp",
          # "pig.last.year.sc",
          "temp5trend.sc",
          "precip5trend.sc",
          "reg.roi5trend.sc",
          "plant.anom.nb.sc",
          "plant.anom.prev.sc",
          # "take5trend.sc",
          "take.hog.intens.sc",
          "prop.nfsp.sc",
          "crp.prop.sc"
    ),
    name=c("Planting anomaly",
           "County",
           "Ecoregion",
           # "Pig presence previous year",
           "Temperature 5 yr trend",
           "Precipitation 5 yr trend",
           "ROI 5 yr trend",
           "Neighboring planting anomaly",
           "Previous year's planting anomaly",
           # "Take 5 yr trend",
           "Take per hog intensity",
           "Prop. of county with pigs",
           "Prop. CRP land"
    ))
  pig_cov_name <- covs_all$name[which(covs_all$cov==paste0(pig_cov,".sc"))]
  
  covsx_all <- covs_all[-(1:3),]
  
  dat_clean <- na.omit(dat[,covs_all$cov])
  dat_clean$county_idx <- as.numeric(factor(dat_clean$GEOID))
  dat_clean$region_idx <- as.numeric(factor(dat_clean$division_grp))
  
  reg_count_idx <- dat_clean %>% group_by(county_idx) %>% 
    summarise(reg_count_idx=unique(region_idx)) %>% 
    arrange(county_idx)
  
  xmat <- dat_clean[,covsx_all$cov]
  xmat <- sapply(xmat,as.numeric)
  
  if(pig_cov=="ever.pigs" | pig_cov=="nfsp.level" | pig_cov=="hog.intensity"){
    tidx <- which(covsx_all$cov=="take5trend.sc")
    xmat <- xmat[,-tidx]
    covsx_all <- covsx_all[-tidx,]
  }
  
  ##model checks ---------------------
  ypred <- samples_all[,grep("ypred",colnames(samples_all))]
  ypred_sum <- data.frame(mn=colMeans(ypred),
                          lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
                          uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
                          obs = dat_clean$plant.anom)
  
  ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
  ypred_sum$std_resid <- scale(ypred_sum$resid)
  
  ypred_sum = cbind(ypred_sum,dat_clean)
  ypred_sum <- ypred_sum %>%
    left_join(dat %>% dplyr::select(GEOID,lat,long) %>% distinct())
  
  ###posterior predictive distribution ----------------
  ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
                                    values_to="value",
                                    names_to="typ"))+
    geom_histogram(aes(x=value,fill=typ),alpha=0.5)+
    scale_fill_discrete(name="",labels=c("Post pred",
                                         "Data"))+
    ggtitle(paste(commod_names_t[commod_idx]," - Posterior predictions"))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"_data_post_dist.jpeg"),
         device="jpeg",height=5,width=7,units="in")

  ##simulated model stats --------------------
  sd <- function(x){sqrt(var(x))}
  ypred_stat <- data.frame(mn=rowMeans(ypred),
                           sd=apply(ypred,1,sd),
                           idx=1:nrow(ypred))
  ggplot()+
    geom_histogram(data=ypred_stat,aes(x=mn))+
    geom_vline(xintercept=mean(dat_clean$plant.anom),col="red",lty=2)+
    ggtitle(paste(commod_names_t[commod_idx]," - Posterior data mean"))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/pval_mn.jpeg"),
         device="jpeg",width=7,height=5,units="in")

  ggplot()+
    geom_histogram(data=ypred_stat,aes(x=sd))+
    geom_vline(xintercept=sqrt(var(dat_clean$plant.anom)),col="red",lty=2)+
    ggtitle(paste(commod_names_t[commod_idx]," - Posterior data SD"))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/pval_sd.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
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
  
  rm(samples_all)
  
  # only counties with pigs -------------
  ## load samples -----------------
  load(paste0("./Model outputs/",commod_names[commod_idx],"_op_region_rs_multi.RData"))
  
  ## load data------------
  dat_op <- dat %>% filter(year<2023) %>% 
    filter(!is.nan(plant.anom) & !is.na(plant.anom) & 
             !is.na(GEOID) & year>=2009 &
             (pig.last.year==1 | pig.in.year==1 )#& long>=-112
    )%>% 
    distinct() 
  
  pig_cov <- "take.hog.intens" 
  
  covs_op <- data.frame(
    cov=c("plant.anom",
          "GEOID",
          "division_grp",
          # "pig.last.year.sc",
          "temp5trend.sc",
          "precip5trend.sc",
          "reg.roi5trend.sc",
          "plant.anom.nb.sc",
          "plant.anom.prev.sc",
          # "take5trend.sc",
          "take.hog.intens.sc",
          "prop.nfsp.sc",
          "crp.prop.sc"
    ),
    name=c("Planting anomaly",
           "County",
           "Ecoregion",
           # "Pig presence previous year",
           "Temperature 5 yr trend",
           "Precipitation 5 yr trend",
           "ROI 5 yr trend",
           "Neighboring planting anomaly",
           "Previous year's planting anomaly",
           # "Take 5 yr trend",
           "Take per hog intensity",
           "Prop. of county with pigs",
           "Prop. CRP land"
    ))
  pig_cov_name <- covs_op$name[which(covs_op$cov==paste0(pig_cov,".sc"))]
  
  covsx_op <- covs_op[-(1:3),]
  
  dat_clean_op <- na.omit(dat_op[,covs_op$cov])
  dat_clean_op$county_idx <- as.numeric(factor(dat_clean_op$GEOID))
  dat_clean_op$region_idx <- as.numeric(factor(dat_clean_op$division_grp))
  
  reg_count_idx <- dat_clean_op %>% group_by(county_idx) %>% 
    summarise(reg_count_idx=unique(region_idx)) %>% 
    arrange(county_idx)
  
  xmat <- dat_clean_op[,covsx_op$cov]
  xmat <- sapply(xmat,as.numeric)
  
  if(pig_cov=="ever.pigs" | pig_cov=="nfsp.level" | pig_cov=="hog.intensity"){
    tidx <- which(covsx_op$cov=="take5trend.sc")
    xmat <- xmat[,-tidx]
    covsx_op <- covsx_op[-tidx,]
  }
  
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
    left_join(dat %>% dplyr::select(GEOID,lat,long) %>% distinct())
  
  ###posterior predictive distribution ----------------
  ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
                                    values_to="value",
                                    names_to="typ"))+
    geom_histogram(aes(x=value,fill=typ),alpha=0.5,col="transparent")+
    scale_fill_discrete(name="",labels=c("Post pred",
                                         "Data"))+
    ggtitle(paste(commod_names_t[commod_idx]," - Posterior predictions - Counties with pigs"))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"_data_post_dist_op.jpeg"),
         device="jpeg",height=5,width=7,units="in")
  
  ##simulated model stats --------------------
  ypred_stat <- data.frame(mn=rowMeans(ypred),
                           sd=apply(ypred,1,sd),
                           idx=1:nrow(ypred))
  ggplot()+
    geom_histogram(data=ypred_stat,aes(x=mn))+
    geom_vline(xintercept=mean(dat_clean_op$plant.anom),col="red",lty=2)+
    ggtitle(paste(commod_names_t[commod_idx]," - Posterior data mean - Counties with pigs"))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/pval_mn_op.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  ggplot()+
    geom_histogram(data=ypred_stat,aes(x=sd))+
    geom_vline(xintercept=sqrt(var(dat_clean_op$plant.anom)),col="red",lty=2)+
    ggtitle(paste(commod_names_t[commod_idx]," - Posterior data SD - Counties with pigs"))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/pval_sd_op.jpeg"),
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
