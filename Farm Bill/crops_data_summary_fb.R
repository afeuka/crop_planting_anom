### Title: Summary statistics for crop planting
### Author: Abbey Feuka
### Date: 07NOV24

library(tidyverse)

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


#all counties-------------
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


counties_df <- data.frame(Crop=commod_names_t,
                          nCounties=NA,
                          nEcoRegions=NA,
                          n=NA)
for(commod_idx in 1:length(commod_names_c)){
  # if(commod_idx!=4){ #not enough peanut FB counties
    # only counties with pigs ----------------------------
    dat_op  <- clean_crop_dat(dat_orig=dat_orig_fb,
                              nb=nb,
                              temporal=temporal,
                              commod_name = commod_names[commod_idx], 
                              only_pigs = T)
    
    dat_clean_op <- dat_op$dat_clean
    covsx_op  <- dat_op$covsx
    xmat <- dat_op$xmat
    reg_count_idx <- dat_op$reg_count_idx
    

  
  
  counties_df$nCounties[commod_idx] <- length(unique(dat_clean_op$GEOID))
  counties_df$nEcoRegions[commod_idx] <- length(unique(dat_clean_op$division_grp))
  counties_df$n[commod_idx] <-  nrow(dat_clean_op)
  # }
}

write.csv(counties_df,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/Model outputs/",subfolder,"/Plots/Combined Figures/data_summary_table.csv"))
