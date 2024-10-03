library(tidyverse)
library(rgdal)
library(lubridate)
library(foreign)
library(sf)
library(dplyr)

setwd("~/Desktop/git_RiskModel/RiskModel/All_Code_2024")

input.folder <- "./model_outputs"
classes <- c("M" ,"L", "XL") 

for (c in 1:length(classes)) {
  class <- classes[c]
  years <- c("2017", "2018", "2019", "2021", "2022") 
  
  output.folder <- paste0(input.folder, "/summary_outputs_", class)
  if (!dir.exists(output.folder)) {dir.create(output.folder)}
  
  n.boot <- 1000 #from original model 
  n.whales <- 350 #from original model
  
  risk_model_summary_table_withAvoidance <- NULL
  risk_model_summary_table_noAvoidance <- NULL
  
  
  for (k in 1:length(years)) {
    year <- years[k]
    files <- list.files(paste0(input.folder, "/Class", class, "/bootstrap_model_output_", year), full.names = T)
    length(files)
    
    tot.strike.avoid.slow.none <- mat.or.vec(n.boot, 12)
    tot.strike.mort.slow.none <- mat.or.vec(n.boot, 12)
    
    tot.strike.avoid.slow.all <- mat.or.vec(n.boot, 12)
    tot.strike.mort.slow.all <- mat.or.vec(n.boot, 12)
    
    for (i in 1:length(files)) {
      whales <- read.csv(files[i], header=T)
      boot <- str_sub(files[i], -13,-10)
      if (substr(boot,4,4) == "_") {
        boot <- as.numeric(substr(boot,1,3))
      } else {
        boot <- as.numeric(boot)
      }
      month <- as.numeric(str_sub(files[i], -6,-5))
      
      print(boot)
      print(month)
      
      ## AVOIDANCE VS. NO AVOIDANCE FOR SENSITIVITY (.mort)
      #With Avoidance
      if (class == "M") {
        tot.strike.avoid.slow.none[boot, month] <-  length(whales$mort.avoid.slow.none[whales$mort.avoid.slow.none > 0]) #this is unique number of deaths - including avoidance (status quo scenario)
      } else {
        tot.strike.avoid.slow.none[boot, month] <-  length(whales$mort.avoid.slow.none[whales$mort.avoid.slow.none > 0]) #this is unique number of deaths - including avoidance (status quo scenario)
      }
      
      tot.strike.avoid.slow.all[boot, month] <-  length(whales$mort.avoid.slow.all[whales$mort.avoid.slow.all > 0]) #this is unique number of deaths - including avoidance(slow all scenario)
      
      #### WITHOUT Avoidance ####
      if (class == "M") {
        tot.strike.mort.slow.none[boot, month] <-  length(whales$mort.n.slow.none[whales$mort.n.slow.none > 0]) #this is unique number of deaths - including avoidance (ssz scenario)
      } else {
        tot.strike.mort.slow.none[boot, month] <-  length(whales$mort.n.slow.none[whales$mort.n.slow.none > 0]) #this is unique number of deaths - without avoidance (status quo scenario)
      }
      
      tot.strike.mort.slow.all[boot, month] <-  length(whales$mort.n.slow.all[whales$mort.n.slow.all > 0]) #this is unique number of deaths - without avoidance (slow all scenario)
      
    }
    
    ## WITH avoidance ##
    ann.tot.avoid.slow.none <- rowSums(tot.strike.avoid.slow.none)
    ann.tot.avoid.slow.all <- rowSums(tot.strike.avoid.slow.all)
    
    tot.strike.df.slow.none <- data.frame(tot.strike.avoid.slow.none/n.whales)
    names(tot.strike.df.slow.none) <- month.abb
    tot.strike.df.slow.none$boot <- c(1:n.boot)
    
    tot.strike.df.slow.all <- data.frame(tot.strike.avoid.slow.all/n.whales)
    names(tot.strike.df.slow.all) <- month.abb
    tot.strike.df.slow.all$boot <- c(1:n.boot)
    
    write.csv(tot.strike.df.slow.none, paste0(output.folder, "/Avoid_tot_strike_Slow_None_", class, "_", year,  ".csv"), row.names = FALSE)
    
    write.csv(tot.strike.df.slow.all, paste0(output.folder, "/Avoid_tot_strike_Sim_Slow_All", class, "_", year,  ".csv"), row.names = FALSE)
    
    mean.mortality.slow.none <- data.frame(mean(ann.tot.avoid.slow.none))
    colnames(mean.mortality.slow.none) <- "mean_mortality"
    mean.mortality.slow.none$sd_mortality <- sd(ann.tot.avoid.slow.none)
    mean.mortality.slow.none$scenario <- "real world"
    
    mean.mortality.slow.all <- data.frame(mean(ann.tot.avoid.slow.all))
    colnames(mean.mortality.slow.all) <- "mean_mortality"
    mean.mortality.slow.all$sd_mortality <- sd(ann.tot.avoid.slow.all)
    mean.mortality.slow.all$scenario <- "slow all simulation"
    
    
    summary.df <- rbind(mean.mortality.slow.none, mean.mortality.slow.all)
    summary.df$year <- year
    summary.df$avoidance <- "avoidance"
    
    risk_model_summary_table_withAvoidance <- rbind(risk_model_summary_table_withAvoidance, summary.df)
    
    ## WITHOUT avoidance ## 
    ann.tot.mort.slow.none <- rowSums(tot.strike.mort.slow.none)
    ann.tot.mort.slow.all <- rowSums(tot.strike.mort.slow.all)
    
    tot.strike.df.mort.slow.none <- data.frame(tot.strike.mort.slow.none/n.whales)
    names(tot.strike.df.mort.slow.none) <- month.abb
    tot.strike.df.mort.slow.none$boot <- c(1:n.boot)
    
    tot.strike.df.mort.slow.all <- data.frame(tot.strike.mort.slow.all/n.whales)
    names(tot.strike.df.mort.slow.all) <- month.abb
    tot.strike.df.mort.slow.all$boot <- c(1:n.boot)
    
    write.csv(tot.strike.df.mort.slow.none, paste0(output.folder, "/No_Avoid_tot_strike_Slow_None_", class, "_", year,  ".csv"), row.names = FALSE)
    
    write.csv(tot.strike.df.mort.slow.all, paste0(output.folder, "/No_Avoid_tot_strike_Sim_Slow_All", class, "_", year,  ".csv"), row.names = FALSE)
    
    mean.mortality.mort.slow.none <- data.frame(mean(ann.tot.mort.slow.none))
    colnames(mean.mortality.mort.slow.none) <- "mean_mortality"
    mean.mortality.mort.slow.none$sd_mortality <- sd(ann.tot.mort.slow.none)
    mean.mortality.mort.slow.none$scenario <- "real world"
    
    mean.mortality.mort.slow.all <- data.frame(mean(ann.tot.mort.slow.all))
    colnames(mean.mortality.mort.slow.all) <- "mean_mortality"
    mean.mortality.mort.slow.all$sd_mortality <- sd(ann.tot.mort.slow.all)
    mean.mortality.mort.slow.all$scenario <- "slow all simulation"
    
    summary.df.noAvoid <- rbind(mean.mortality.mort.slow.none, mean.mortality.mort.slow.all)
    summary.df.noAvoid$year <- year
    summary.df.noAvoid$avoidance <- "no avoidance"
    
    risk_model_summary_table_noAvoidance <- rbind(risk_model_summary_table_noAvoidance, summary.df.noAvoid)
    
    
  }
  all_data_summaries <- rbind(risk_model_summary_table_withAvoidance,risk_model_summary_table_noAvoidance)
  
  write.csv(all_data_summaries, paste0(output.folder, "risk_model_summary_table_", class,  ".csv"), row.names = F)
  
}

