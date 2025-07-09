### Title: Data cleaning function for crop planting analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes:


clean_crop_dat <- function(dat_orig, #original data, uncleaned
                           commod_name, #lower case commodity/crop name ("corn","wheat",etc)
                           nb=T,
                           temporal=T,
                           only_pigs=F #logical, only include counties with pigs
                           # pig.cov #column name of covariate used for pig inference
                           ){
  require(tidyverse)
  
  dat <- dat_orig %>% 
    filter(year<2023 & year >=2009) %>% 
    filter(commodity_desc==toupper(commod_name)) %>% 
    filter(!is.nan(plant.anom) & 
             !is.na(plant.anom) & 
             !is.na(GEOID) & year>=2009) %>% 
    distinct() 
  
  if(only_pigs){
    dat <- dat %>% filter(!is.nan(plant.anom) & 
                            !is.na(plant.anom) & 
                            !is.na(GEOID)  &
                            pig.last.year==1) 
  }
  
  #crp interaction
  dat$crp.nfsp <- dat$prop_nfsp * dat$crp_prop
  dat$crp.nfsp.sc <- scale(dat$crp.nfsp)

  covs_all <- data.frame(
    cov=c("plant.anom",
          "GEOID",
          "division_grp",
          # "pig.last.year.sc",
          "temp5trend.sc",
          "precip5trend.sc",
          "reg.roi5trend.sc",
          if(nb){"plant.anom.nb.sc"},
          if(temporal){"plant.anom.prev.sc"},
          # "take5trend.sc",
          # "take.hog.intens.prev.sc",
          "take.hog.intens.5yeartrend.sc",
          "prop.nfsp.sc",
          "crp.prop.sc",
          "crp.nfsp.sc"
    ),
    name=c("Planting anomaly",
           "County",
           "Ecoregion",
           # "Pig presence previous year",
           "Temperature 5 yr trend",
           "Precipitation 5 yr trend",
           "ROI 5 yr trend",
           if(nb){"Neighboring planting anomaly"},
           if(temporal){"Previous year's planting anomaly"},
           # "Take 5 yr trend",
           # "Take per hog intensity previous year",
           "Take per wild pig intensity 5 yr trend",
           "Prop. of county with wild pigs",
           "Prop. CRP land",
           "CRP x pigs"
    ))
  # pig_cov_name <- covs_all$name[which(covs_all$cov==paste0(pig_cov,".sc"))]
  
  covsx_all <- covs_all[-(1:3),]
  
  # if(pig_cov=="ever.pigs" | pig_cov=="nfsp.level" | pig_cov=="hog.intensity"){
  #   tidx <- which(covsx_all$cov=="take5trend.sc")
  #   xmat <- xmat[,-tidx]
  #   covsx_all <- covsx_all[-tidx,]
  # }
  
  dat_clean <- na.omit(dat[,covs_all$cov])
  dat_clean$county_idx <- as.numeric(factor(dat_clean$GEOID))
  dat_clean$region_idx <- as.numeric(factor(dat_clean$division_grp))
  
  xmat <- dat_clean[,covsx_all$cov]
  xmat <- sapply(xmat,as.numeric)
  
  covsx_all <- covs_all[-(1:3),]
  
  dat_clean <- na.omit(dat_clean[,covs_all$cov])
  dat_clean$county_idx <- as.numeric(factor(dat_clean$GEOID))
  dat_clean$region_idx <- as.numeric(factor(dat_clean$division_grp))
  
  reg_count_idx <- dat_clean %>% group_by(county_idx) %>% 
    summarise(reg_count_idx=unique(region_idx)) %>% 
    arrange(county_idx)
  
  xmat <- dat_clean[,covsx_all$cov]
  xmat <- sapply(xmat,as.numeric)
  
  return(list(dat_clean=dat_clean,
         covsx=covsx_all,
         xmat=xmat,
         reg_count_idx=reg_count_idx))
}
