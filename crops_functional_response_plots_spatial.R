### Title: Functional response plots for crop planting analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes:

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/crop_planting_anom")

library(tidyverse)

source("./Functions/fun_response.R")
source("./Functions/clean_crop_dat.R")

mod_typ="spatial"
subfolder='Spatial'
model_dir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/"

#original data --------------------
load(paste0(model_dir,"/Data/all_crops_anom_scaled_2023_2009_anom.RData"))
dat_orig <- dat 

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)


#load beta table --------------
beta_tab<- read.csv(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Model outputs/",
                           subfolder,"/Plots/Combined Figures/op_betas_all_table.csv"))
beta_tab <- beta_tab %>% select(-X) 
beta_tab <- beta_tab %>% 
  filter((grepl("pig",cov) & signif==TRUE))

g<-list()
for(i in 1:nrow(beta_tab)){
  dat_op <- clean_crop_dat(dat_orig=dat,commod_name = tolower(beta_tab$crop[i]),
                           nb=FALSE,
                           temporal=TRUE,
                           only_pigs = T)
  
  load(paste0(model_dir,"/Model outputs/",subfolder,"/",
              tolower(beta_tab$crop[i]),"_op_",mod_typ,".RData"))
  
  beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                        !grepl("beta_county",colnames(samples_all)) &
                        !grepl("beta_region",colnames(samples_all))]
  beta <- cbind(beta,
                chain_idx=samples_all$chain_idx,
                samp_idx=samples_all$samp_idx)
  
  beta_long_all <- beta %>%
    pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
  beta_long_all$cov <- rep(c("Intercept",dat_op$covsx$name),nrow(beta))
  
  s <- samples_all[,(grepl("s",colnames(samples_all)) | grepl("chain_idx",colnames(samples_all))) &
                     !grepl("lscale",colnames(samples_all)) &
                     !grepl("tau_s",colnames(samples_all))]
  # tau <- samples_all[,"tau"]
  
  fun_take_op <- fun_response(dat_df=dat_op$dat_clean,
                              scale_cov_name=dat_op$covsx$cov[dat_op$covsx$name==beta_tab$cov[i]],
                              covsx = dat_op$covsx,
                              beta_long=beta_long_all,
                              s=s,
                              # tau=tau,
                              x_incr=0.5,
                              spatial=T)

  
  g[[i]] <- ggplot(fun_take_op)+
    geom_line(aes(x=cov_bt,y=mn,col=county_idx))+
    guides(col="none")+
    xlab(dat_op$covsx$name[dat_op$covsx$name==beta_tab$cov[i]]) +
    geom_hline(yintercept=0,lty=2)+
    ylab("County-level  planting anomaly")+
    ggtitle(paste("Counties with pigs - ",beta_tab$crop[i]))+
    theme(text=element_text(size=15))
  
}

cowplot::plot_grid(g[[1]],g[[2]])


ggsave(filename=paste0(model_dir,"/Model outputs/",subfolder,"/Plots/Combined Figures/func_response.jpeg"),
       width=5,height=5,units="in")
