### Title: Functional response plots for crop planting analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes:
library(tidyverse)

source("./crop_planting_anom/Functions/fun_response.R")
source("./crop_planting_anom/Functions/clean_crop_dat.R")

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

#original data --------------------
load("./Data/all_crops_anom_scaled_2023_2009_anom.RData") #dat

# corn ------------------------
## all counties ------------------
dat_corn <- clean_crop_dat(dat_orig=dat,commod_name = "corn",only_pigs = F)

load("./Model outputs/corn_region_rs_multi.RData")

beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                      !grepl("beta_county",colnames(samples_all)) &
                      !grepl("beta_region",colnames(samples_all))]
beta <- cbind(beta,
              chain_idx=samples_all$chain_idx,
              samp_idx=samples_all$samp_idx)

beta_long_all <- beta %>%
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long_all$cov <- rep(c("Intercept",dat_corn$covsx$name),nrow(beta))


fun_take <- fun_response(dat_df=dat_corn$dat_clean,
             scale_cov_name="take.hog.intens.sc",
             covsx = dat_corn$covsx,
             beta_long=beta_long_all)

ggplot(fun_take)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Take by hog intensity") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level corn planting anomaly")+
  ggtitle(paste("All counties - Corn"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Corn/all_counties_take_fr.jpeg"),
       width=5,height=5,units="in")

## only pigs-----------------------
dat_corn_op <- clean_crop_dat(dat_orig=dat,commod_name = "corn",only_pigs = T)

load("./Model outputs/corn_op_region_rs_multi.RData")

beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                      !grepl("beta_county",colnames(samples_all)) &
                      !grepl("beta_region",colnames(samples_all))]
beta <- cbind(beta,
              chain_idx=samples_all$chain_idx,
              samp_idx=samples_all$samp_idx)

beta_long_all <- beta %>%
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long_all$cov <- rep(c("Intercept",dat_corn_op$covsx$name),nrow(beta))

fun_take_op <- fun_response(dat_df=dat_corn_op$dat_clean,
                         scale_cov_name="take.hog.intens.sc",
                         covsx = dat_corn_op$covsx,
                         beta_long=beta_long_all)

ggplot(fun_take_op)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Take by hog intensity") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level corn planting anomaly")+
  ggtitle(paste("Counties with pigs - Corn"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Corn/op_take_fr.jpeg"),
       width=5,height=5,units="in")

# cotton -------------------
## all counties---------------------
### nfsp ------------------
dat_cotton <- clean_crop_dat(dat_orig=dat,commod_name = "cotton",only_pigs = F)

load("./Model outputs/cotton_region_rs_multi.RData")

beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                      !grepl("beta_county",colnames(samples_all)) &
                      !grepl("beta_region",colnames(samples_all))]
beta <- cbind(beta,
              chain_idx=samples_all$chain_idx,
              samp_idx=samples_all$samp_idx)

beta_long_all <- beta %>%
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long_all$cov <- rep(c("Intercept",dat_cotton$covsx$name),nrow(beta))

fun_nfsp <- fun_response(dat_df=dat_cotton$dat_clean,
                         scale_cov_name="prop.nfsp.sc",
                         covsx = dat_cotton$covsx,
                         beta_long=beta_long_all)

ggplot(fun_nfsp)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Proportion of NFSP-estimated hog range") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level cotton planting anomaly")+
  ggtitle(paste("All counties - Cotton"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Cotton/all_counties_nfsp_fr.jpeg"),
       width=5,height=5,units="in")

### crp -------------------------
fun_crp <- fun_response(dat_df=dat_cotton$dat_clean,
                         scale_cov_name="crp.prop.sc",
                         covsx = dat_cotton$covsx,
                         beta_long=beta_long_all)

ggplot(fun_crp)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Proportion of CRP land") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level cotton planting anomaly")+
  ggtitle(paste("All counties - Cotton"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Cotton/all_counties_crp_fr.jpeg"),
       width=5,height=5,units="in")

## only pigs ----------------------
dat_cotton_op <- clean_crop_dat(dat_orig=dat,commod_name = "cotton",only_pigs = T)

load("./Model outputs/cotton_op_region_rs_multi.RData")

beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                      !grepl("beta_county",colnames(samples_all)) &
                      !grepl("beta_region",colnames(samples_all))]
beta <- cbind(beta,
              chain_idx=samples_all$chain_idx,
              samp_idx=samples_all$samp_idx)

beta_long_all <- beta %>%
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long_all$cov <- rep(c("Intercept",dat_cotton_op$covsx$name),nrow(beta))


fun_nfsp_op <- fun_response(dat_df=dat_cotton_op$dat_clean,
                         scale_cov_name="prop.nfsp.sc",
                         covsx = dat_cotton_op$covsx,
                         beta_long=beta_long_all)

ggplot(fun_nfsp_op)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Proportion of NFSP-estimated hog range") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level cotton planting anomaly")+
  ggtitle(paste("Counties with pigs - Cotton"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Cotton/op_nfsp_fr.jpeg"),
       width=5,height=5,units="in")

### crp -------------------------
fun_crp_op <- fun_response(dat_df=dat_cotton_op$dat_clean,
                        scale_cov_name="crp.prop.sc",
                        covsx = dat_cotton_op$covsx,
                        beta_long=beta_long_all)

ggplot(fun_crp_op)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Proportion of CRP land") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level cotton planting anomaly")+
  ggtitle(paste("Counties with pigs - Cotton"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Cotton/op_crp_fr.jpeg"),
       width=5,height=5,units="in")

# soybeans ----------------------------
## all counties -----------------------
dat_soy <- clean_crop_dat(dat_orig=dat,commod_name = "soybeans",only_pigs = F)

load("./Model outputs/soybeans_region_rs_multi.RData")

beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                      !grepl("beta_county",colnames(samples_all)) &
                      !grepl("beta_region",colnames(samples_all))]
beta <- cbind(beta,
              chain_idx=samples_all$chain_idx,
              samp_idx=samples_all$samp_idx)

beta_long_all <- beta %>%
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long_all$cov <- rep(c("Intercept",dat_soy$covsx$name),nrow(beta))


fun_roi <- fun_response(dat_df=dat_soy$dat_clean,
                         scale_cov_name="reg.roi5trend.sc",
                         covsx = dat_soy$covsx,
                         beta_long=beta_long_all)

ggplot(fun_roi)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("5-yr ROI trend") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level soybean planting anomaly")+
  ggtitle(paste("All counties - Soybeans"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Soybeans/all_counties_roi_fr.jpeg"),
       width=5,height=5,units="in")

## only pigs -------------------
dat_soy_op <- clean_crop_dat(dat_orig=dat,commod_name = "soybeans",only_pigs = T)

load("./Model outputs/soybeans_op_region_rs_multi.RData")

beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                      !grepl("beta_county",colnames(samples_all)) &
                      !grepl("beta_region",colnames(samples_all))]
beta <- cbind(beta,
              chain_idx=samples_all$chain_idx,
              samp_idx=samples_all$samp_idx)

beta_long_all <- beta %>%
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long_all$cov <- rep(c("Intercept",dat_soy_op$covsx$name),nrow(beta))


fun_roi_op <- fun_response(dat_df=dat_soy_op$dat_clean,
                        scale_cov_name="reg.roi5trend.sc",
                        covsx = dat_soy_op$covsx,
                        beta_long=beta_long_all)

ggplot(fun_roi_op)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("5-yr ROI trend") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level soybean planting anomaly")+
  ggtitle(paste("Counties with pigs - Soybeans"))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/Soybeans/op_roi_fr.jpeg"),
       width=5,height=5,units="in")
