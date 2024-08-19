### Title: Plotting multi-chain crop analysis
### Author: Abbey Feuka
### Date: 15AUG24
### Notes: 

library(nimble)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)
# commod_idx <- 1

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
  # xmat <- cbind(xmat,xmat[,which(covsx_all$cov=="crp.prop.sc")]*
  #                 xmat[,which(covsx_all$cov==paste0(pig_cov,".sc"))])
  # xmat <- cbind(xmat,xmat[,which(covsx_all$cov=="take5trend.sc")]*
  #                 xmat[,which(covsx_all$cov==paste0(pig_cov,".sc"))])
  
  if(pig_cov=="ever.pigs" | pig_cov=="nfsp.level" | pig_cov=="hog.intensity"){
    tidx <- which(covsx_all$cov=="take5trend.sc")
    xmat <- xmat[,-tidx]
    covsx_all <- covsx_all[-tidx,]
  }
  
  ##trace plots -------------------------
  tau <- cbind.data.frame(samples_all[,"tau"],
                chain_idx=samples_all$chain_idx,
                samp_idx=samples_all$samp_idx)
  lscale <- cbind.data.frame(samples_all[,"lscale"],
               chain_idx=samples_all$chain_idx,
               samp_idx=samples_all$samp_idx)
  
  tau_long <- tau %>%  
    pivot_longer(grep("tau",colnames(tau)),values_to="value",names_to="beta")
  lscale_long <- lscale %>%  
    pivot_longer(grep("lscale",colnames(lscale)),values_to="value",names_to="beta")
  
  ggplot(tau_long)+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    ylab("tau")+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/tau_trace.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  ggplot(lscale_long)+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    ylab("logistic scale")+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/lscale_trace.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  beta <- samples_all[,grepl("beta",colnames(samples_all)) & 
                        !grepl("beta_county",colnames(samples_all)) &
                        !grepl("beta_region",colnames(samples_all))]
  beta <- cbind(beta,
                chain_idx=samples_all$chain_idx,
                samp_idx=samples_all$samp_idx)
  # beta$idx <- samples_all[,"idx"]
  
  
  # beta_cor<- cor(beta)
  # # beta_cor <- cor(beta %>% select(-idx))
  # beta_cor_idx <- as.data.frame(which(abs(beta_cor)>0.7 &beta_cor!=1,arr.ind = T))
  # beta_cor_idx[,1] <- c("int",covsx_all$cov)[beta_cor_idx[,1]]
  # beta_cor_idx[,2] <- c("int",covsx_all$cov)[beta_cor_idx[,2]]
  # beta_cor_idx
  
  # beta$idx_samp <- 1:nrow(beta)
  beta_long_all <- beta %>%  
    pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
  beta_long_all$cov <- rep(c("Intercept",covsx_all$name),nrow(beta))
  
  ggplot(beta_long_all %>% filter(beta!="beta[1]"))+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    geom_hline(yintercept=0,col="blue",lty=2)+
    facet_wrap(.~cov)+
    ggtitle(paste(commod_names_t[commod_idx],"- Fixed Effects"))+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/beta_trace.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  ##model checks ---------------------
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
  #   left_join(dat %>% dplyr::select(GEOID,lat,long) %>% distinct())
  
  ###posterior predictive distribution ----------------
  # ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
  #                                   values_to="value",
  #                                   names_to="typ"))+
  #   geom_density(aes(x=value,fill=typ),alpha=0.5)+
  #   scale_fill_discrete(name="",labels=c("Post pred",
  #                                        "Data"))
  
  ##simulated model stats --------------------
  # ypred_stat <- data.frame(mn=rowMeans(ypred),
  #                          sd=sapply(1:nrow(ypred),function(i)sqrt(var(ypred[i,]))),
  #                          idx=1:nrow(ypred))
  # ggplot()+
  #   geom_density(data=ypred_stat,aes(x=mn))+
  #   geom_vline(xintercept=mean(dat_clean$plant.anom),col="red",lty=2)
  # ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/pval_mn.jpeg"),
  #        device="jpeg",width=7,height=5,units="in")
  # 
  # ggplot()+
  #   geom_density(data=ypred_stat,aes(x=sd))+
  #   geom_vline(xintercept=sqrt(var(dat_clean$plant.anom)),col="red",lty=2)
  # ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/pval_sd.jpeg"),
  #        device="jpeg",width=7,height=5,units="in")
  
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
  
  ## fixed effects --------------------------
  beta_sum_all <- beta_long_all %>% 
    group_by(cov) %>% 
    summarise(mn=mean(value),
              lci=quantile(value,probs=0.025),
              uci=quantile(value,probs=0.975))
  
  beta_sum_all$signif <- ifelse(sign(beta_sum_all$lci)==sign(beta_sum_all$uci),TRUE,FALSE)
  
  ggplot(beta_sum_all %>% filter(cov!="Intercept"))+
    geom_point(aes(y=cov,x=mn,alpha=signif),size=2.5)+
    geom_errorbar(aes(y=cov,xmin=lci,xmax=uci,alpha=signif),width=0,lwd=2)+
    geom_vline(xintercept = 0,col="blue",lty=2)+
    scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
    ylab("Covariate")+xlab("Coefficient estimate")+
    ggtitle(paste(commod_names_t[commod_idx],"- All counties"))+
    theme(text=element_text(size=15))
  
  ggsave(filename = paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/all_counties_betas.jpeg"),
         device = "jpeg",
         width=8,height=8,units="in")
  
  ## random effects ---------------------------------
  ### regions--------------------
  reg_col<- c("Hot Continental"="#F8766D", 
              "Warm Continental"="#DD9600",
              "Subtropical"="springgreen4",
              "Tropical/Subtropical Steppe"="#39B600",
              "Prairie"="darksalmon",
              "Marine"="#00B0F6",
              "Temperate Desert"="goldenrod4",
              "Temperate Steppe"="hotpink4",
              "Tropical/Subtropical Desert"="chocolate1",
              "Tropical/Subtropical"="#A3A500",
              "Mediterranean"="paleturquoise4")
  
  beta_region <- samples_all[,grep("beta_region",colnames(samples_all))]
  beta_region_sum <-data.frame(mn=colMeans(beta_region),
                               lci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.025)),
                               uci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.975)),
                               division_grp=unique(dat_clean$division_grp))
  
  beta_region_sum$signif <- ifelse(sign(beta_region_sum$lci)==sign(beta_region_sum$uci),TRUE,FALSE)
  
  ggplot(beta_region_sum)+
    geom_errorbar(aes(y=division_grp,col=division_grp,
                      xmin=lci,xmax=uci,alpha=signif),width=0,lwd=2)+
    geom_point(aes(y=division_grp,x=mn,alpha=signif,col=division_grp),size=3)+
    geom_vline(xintercept = 0,col="blue",lty=2)+
    scale_color_manual(values=reg_col)+
    scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
    guides(col="none")+
    ylab("Ecoregion")+xlab("Random intercept estimate")+
    theme(text=element_text(size=15),
          axis.line = element_blank(),
          axis.ticks = element_blank())+
    ggtitle(paste(commod_names_t[commod_idx],"- All counties"))
  
  ggsave(filename = paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/all_regions_re.jpeg"),
         device = "jpeg",
         width=8,height=8,units="in")
  
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
  # xmat <- cbind(xmat,xmat[,which(covsx_op$cov=="crp.prop.sc")]*
  #                 xmat[,which(covsx_op$cov==paste0(pig_cov,".sc"))])
  # xmat <- cbind(xmat,xmat[,which(covsx_op$cov=="take5trend.sc")]*
  #                 xmat[,which(covsx_op$cov==paste0(pig_cov,".sc"))])
  
  if(pig_cov=="ever.pigs" | pig_cov=="nfsp.level" | pig_cov=="hog.intensity"){
    tidx <- which(covsx_op$cov=="take5trend.sc")
    xmat <- xmat[,-tidx]
    covsx_op <- covsx_op[-tidx,]
  }
  
  ## trace plots -------------------------
  tau <- cbind.data.frame(samples_all[,"tau"],
                          chain_idx=samples_all$chain_idx,
                          samp_idx=samples_all$samp_idx)
  lscale <- cbind.data.frame(samples_all[,"lscale"],
                             chain_idx=samples_all$chain_idx,
                             samp_idx=samples_all$samp_idx)
  
  tau_long <- tau %>%  
    pivot_longer(grep("tau",colnames(tau)),values_to="value",names_to="beta")
  lscale_long <- lscale %>%  
    pivot_longer(grep("lscale",colnames(lscale)),values_to="value",names_to="beta")
  
  ggplot(tau_long)+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    ylab("tau")+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/tau_trace_op.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  ggplot(lscale_long)+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    ylab("logistic scale")+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/lscale_trace_op.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  beta <- samples_all[,grepl("beta",colnames(samples_all)) & 
                        !grepl("beta_county",colnames(samples_all)) &
                        !grepl("beta_region",colnames(samples_all))]
  beta <- cbind(beta,
                chain_idx=samples_all$chain_idx,
                samp_idx=samples_all$samp_idx)
  # beta$idx <- samples_all[,"idx"]
  
  
  # beta_cor<- cor(beta)
  # # beta_cor <- cor(beta %>% select(-idx))
  # beta_cor_idx <- as.data.frame(which(abs(beta_cor)>0.7 &beta_cor!=1,arr.ind = T))
  # beta_cor_idx[,1] <- c("int",covsx_all$cov)[beta_cor_idx[,1]]
  # beta_cor_idx[,2] <- c("int",covsx_all$cov)[beta_cor_idx[,2]]
  # beta_cor_idx
  
  # beta$idx_samp <- 1:nrow(beta)
  beta_long_all <- beta %>%  
    pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
  beta_long_all$cov <- rep(c("Intercept",covsx_all$name),nrow(beta))
  
  ggplot(beta_long_all %>% filter(beta!="beta[1]"))+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    geom_hline(yintercept=0,col="blue",lty=2)+
    facet_wrap(.~cov)+
    ggtitle(paste(commod_names_t[commod_idx],"- Fixed Effects"))+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/beta_trace_op.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  ##model checks ---------------------
  # ypred <- samples_op[,grep("ypred",colnames(samples_op))]
  # ypred_sum <- data.frame(mn=colMeans(ypred),
  #                         lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
  #                         uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
  #                         obs = dat_clean_op$plant.anom)
  # 
  # ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
  # ypred_sum$std_resid <- scale(ypred_sum$resid)
  # 
  # ypred_sum = cbind(ypred_sum,dat_clean_op)
  # ypred_sum <- ypred_sum %>% 
  #   left_join(dat %>% dplyr::select(GEOID,lat,long) %>% distinct())
  # 
  # ###posterior predictive distribution ----------------
  # ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
  #                                   values_to="value",
  #                                   names_to="typ"))+
  #   geom_density(aes(x=value,fill=typ),alpha=0.5,col="transparent")+
  #   scale_fill_discrete(name="",labels=c("Post pred",
  #                                        "Data"))
  # 
  # ##simulated model stats --------------------
  # ypred_stat <- data.frame(mn=rowMeans(ypred),
  #                          sd=sapply(1:nrow(ypred),function(i)sqrt(var(ypred[i,]))),
  #                          idx=1:nrow(ypred))
  # ggplot()+
  #   geom_density(data=ypred_stat,aes(x=mn))+
  #   geom_vline(xintercept=mean(dat_clean_op$plant.anom),col="red",lty=2)
  # 
  # ggplot()+
  #   geom_density(data=ypred_stat,aes(x=sd))+
  #   geom_vline(xintercept=sqrt(var(dat_clean_op$plant.anom)),col="red",lty=2)
  # 
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
  
  ##fixed effects --------------------------
  beta_sum <- beta_long_all %>% 
    group_by(cov) %>% 
    summarise(mn=mean(value),
              lci=quantile(value,probs=0.025),
              uci=quantile(value,probs=0.975))
  
  beta_sum$signif <- ifelse(sign(beta_sum$lci)==sign(beta_sum$uci),TRUE,FALSE)
  
  ggplot(beta_sum %>% filter(cov!="Intercept"))+
    geom_point(aes(y=cov,x=mn,alpha=signif),size=2.5)+
    geom_errorbar(aes(y=cov,xmin=lci,xmax=uci,alpha=signif),width=0,lwd=2)+
    geom_vline(xintercept = 0,col="blue",lty=2)+
    scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
    ylab("Covariate")+xlab("Coefficient estimate")+
    ggtitle("Counties with pigs")+
    theme(text=element_text(size=15))
  ggsave(filename = paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/op_betas.jpeg"),
         device = "jpeg",
         width=8,height=8,units="in")
  
  ## random effects ---------------------------
  beta_region <- samples_all[,grep("beta_region",colnames(samples_all))]
  beta_region_sum <-data.frame(mn=colMeans(beta_region),
                               lci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.025)),
                               uci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.975)),
                               division_grp=unique(dat_clean$division_grp))
  
  beta_region_sum$signif <- ifelse(sign(beta_region_sum$lci)==sign(beta_region_sum$uci),TRUE,FALSE)
  
  ggplot(beta_region_sum)+
    geom_errorbar(aes(y=division_grp,col=division_grp,
                      xmin=lci,xmax=uci,alpha=signif),width=0,lwd=2)+
    geom_point(aes(y=division_grp,x=mn,alpha=signif,col=division_grp),size=3)+
    geom_vline(xintercept = 0,col="blue",lty=2)+
    scale_color_manual(values=reg_col)+
    scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
    guides(col="none")+
    ylab("Ecoregion")+xlab("Random intercept estimate")+
    theme(text=element_text(size=15),
          axis.line = element_blank(),
          axis.ticks = element_blank())+
    ggtitle(paste(commod_names_t[commod_idx],"- All counties"))
  
  ggsave(filename = paste0("./Model outputs/Plots/",commod_names_t[commod_idx],"/op_regions_re.jpeg"),
         device = "jpeg",
         width=8,height=8,units="in")
}
