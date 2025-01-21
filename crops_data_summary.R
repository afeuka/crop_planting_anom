### Title: Summary statistics for crop planting
### Author: Abbey Feuka
### Date: 07NOV24

library(tidyverse)

source("./Functions/clean_crop_dat.R")

#all counties-------------
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

commod_idx <- 1
counties_df <- data.frame(Crop=commod_names_t,
                         nCounties=NA,
                         nEcoRegions=NA,
                         n=NA)

nb=FALSE
temporal=TRUE

for(commod_idx in 1:length(commod_names_c)){
  dat_op  <- clean_crop_dat(dat_orig=dat_orig,
                            nb=nb,
                            temporal=temporal,
                            commod_name = commod_names[commod_idx], only_pigs = T)
  
  dat_clean_op <- dat_op$dat_clean
  # covsx_op  <- dat_op$covsx
  # xmat <- dat_op$xmat
  # reg_count_idx <- dat_op$reg_count_idx
  
  if("take.hog.intens.5yeartrend.sc"%in%colnames(dat_clean) & mod_typ!="spatial"){
    subfolder <- "Take Trend"
  } else if("take.hog.intens.prev.sc"%in%colnames(dat_clean) & mod_typ!="spatial") {
    subfolder <- "Take Previous" 
  } else if(mod_typ=="spatial") {
    subfolder <- "Spatial" 
  } else {
    stop("Must include take covariate. Check clean_crop_dat.R")
  }
  
  counties_df$nCounties[commod_idx] <- length(unique(dat_clean_op$GEOID))
  counties_df$nEcoRegions[commod_idx] <- length(unique(dat_clean_op$division_grp))
  counties_df$n[commod_idx] <-  nrow(dat_clean_op)
}

write.csv(counties_df,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Model outputs/",subfolder,"/Plots/Combined Figures/data_summary_table.csv"))
