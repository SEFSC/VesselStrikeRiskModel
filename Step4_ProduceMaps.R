## Step 5 - Spatial Model for plotting

library(tidyverse)
library(rgdal)
library(RStoolbox)
library(lubridate)
library(foreign)
library(sf)


############## Get clean vessel data which includes all possible imputed values ##################

setwd("~/Desktop/git_RiskModel/RiskModel/All_Code_2024")

#read the cleaned vessel data
vessel.data <- readRDS("./Data/ais_rds/all_vessels_clean_distinct_2017_2022.rds")

#there are duplicate vessel ids - so get rid of dupes
vessel.data <- vessel.data %>% distinct(vessels_id, .keep_all = TRUE)

#read in the list of track files - summarized tracks by transit/vessel
track.files <- dir(path = "./Data/ais_rds/", pattern = glob2rx("narw_tracks_*.rds"))
n.track.files <- length(track.files)

#AIS grid shapefile
grid <- readOGR(dsn="./Data/shapefiles", layer="ec22_ease2_n_10km_pl")

#NARW Density grid converted to polygon shapefile - with regional designations for dive behavior and mid-Atlantic
#cells
grid.density <- readRDS("./Data/shapefiles/new_NARW_density_shp_v3.RDS")
grid.density <- data.frame(grid.density %>% select(-"geometry"))

#Get vessel categories for summaries and groupings of some rare vessel types
vessel.cats <- read.csv("./Data/vessel_categories.csv")

LethalityCurveGLM <- readRDS("./Data/LethalityCurveGLM_26Apr24.RDS")


###################################################################################################
################## ENCOUNTER RISK SCRIPTS FROM MARTIN ET AL 2016 ########################################
#####  Martin, Julien et al. (2016), Data from: A quantitative framework for investigating risk of deadly 
#####   collisions between marine wildlife and boats, Dryad, Dataset, https://doi.org/10.5061/dryad.vv150
#########################################################################################################

#function to estimate lambda_e [eqn 1 in the text], animal speed is assumed fixed
getLambda <- function(vb, vm, r, S)
  #vb: boat speed; vm: animal speed; r= critical distance of encounter (rc in eqn 1); S = Area
{
  getJ <- Vectorize(function(vb, vm) {
    alpha <- (2*vm*vb)/(vm^2+vb^2)
    f.theta <- function(theta, alpha) sqrt(1 - alpha * cos(theta))/(2 * pi)
    integrate(f.theta, lower = 0, upper = 2 * pi, alpha = alpha)$value})
  
  (2*r/S)*sqrt(vm^2 + vb^2)*getJ(vb, vm)
  
}

f.vm <- function(v, ...) dweibull(v, ...)

getLambda2 <- function(vb, f.vm=f.vm, r=r, S=S, ...)
{
  getJ <- function(vb, vm) {
    alpha <- ((2*vm*vb)/((vm^2)+(vb^2)))
    f.theta <- function(theta, alpha) sqrt(1 - alpha * cos(theta))/(2 * pi)
    integrate(f.theta, lower = 0, upper = 2 * pi, alpha = alpha)$value}
  
  FUN <- Vectorize(function(x) sqrt(x^2 + vb^2)*getJ(vb, x)*f.vm(x, ...))
  (2*r)/S *integrate(FUN, lower=0, upper=Inf)$value
}

##########################################
#Probability of death given strike speed #
##########################################
#"Fate.LPG ~ Vess.Speed + vess.cat_f + spe.HB + Vess.Speed * spe.HB"

p.strike.mort <- function(v.boat.knots, reclass_Vessel_Size, spe.HB) {
  
  Vess.Speed =  v.boat.knots #vessel speed (knots)
  vess.cat_f = reclass_Vessel_Size
  spe.HB = "Not Humpback" #humpback or not 
  
  newdata <- data.frame(Vess.Speed,vess.cat_f, spe.HB)
  
  p.strike.mort  <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
}

#############################################
#Probability of avoidance given vessel speed#
#############################################
### NEW FUNCTION AFTER MCKENNA ##############

p.avoid <- function(v.boat.ms, PropSuctDepth, depth) {
  ShipSpeed <- v.boat.ms #speed in m/s-1
  
  #Descent Rate - norm distribution
  max.rx.descent <- 2.0 #feeding= based on an average descent speed of 1.40 m s–1, range to 2.0
  min.rx.descent <- 0.81 #0.81 = descent speed range from Baumgartner et al. 2003
  rx.descent <- (runif(1) * (max.rx.descent - min.rx.descent)) + min.rx.descent
  DescentRate <- rx.descent
  
  # Reaction Distance - norm distribution
  max.rx.dist <- 1200
  min.rx.dist <- 10
  rx.dist <- (runif(1) * (max.rx.dist - min.rx.dist)) + min.rx.dist
  ReactionDistance <- rx.dist
  
  WhaleDepth <- DescentRate * (ReactionDistance/ShipSpeed)
  
  ## Depth limit
  PropSuctDepth <- PropSuctDepth
  
  PosDepth <- depth*-1
  
  WhaleHeight <- 3.0 #Moore et al. (2005) Jour. of Cet. Res. & Man.
  
  DepthLimit <- PropSuctDepth + WhaleHeight
  
  p.avoid <- ifelse(WhaleDepth > PropSuctDepth & PosDepth > DepthLimit, 1, 0)
  return(p.avoid)
}

p.avoid.boot <- function(v.boat.ms, PropSuctDepth, depth) {
  n.boot <- 100
  avoid.boot <- rep(NA, n.boot)
  ShipSpeed <- v.boat.ms #speed in m/s-1
  
  #Descent Rate - norm distribution
  max.rx.descent <- 2.0 #feeding= based on an average descent speed of 1.40 m s–1, range to 2.0
  min.rx.descent <- 0.81 #0.81 = descent speed range from Baumgartner et al. 2003
  
  for (i in 1:n.boot) {  
    rx.descent <- (runif(1) * (max.rx.descent - min.rx.descent)) + min.rx.descent
    DescentRate <- rx.descent
    
    # Reaction Distance - norm distribution
    max.rx.dist <- 1200
    min.rx.dist <- 10
    rx.dist <- (runif(1) * (max.rx.dist - min.rx.dist)) + min.rx.dist
    ReactionDistance <- rx.dist
    
    WhaleDepth <- DescentRate * (ReactionDistance/ShipSpeed)
    
    ## Depth limit
    PropSuctDepth <- PropSuctDepth
    
    PosDepth <- depth*-1
    
    WhaleHeight <- 3.0 #Moore et al. (2005) Jour. of Cet. Res. & Man.
    
    DepthLimit <- PropSuctDepth + WhaleHeight
    
    avoid.boot[i] <- ifelse(WhaleDepth > PropSuctDepth & PosDepth > DepthLimit, 1, 0)
  }
  p.avoid <- mean(avoid.boot)
  return(p.avoid)
}

############ Probability at surface given DiveArea region #############
depth_props_by_region <- readRDS("./Data/depth_props_by_region_weighted_final_Apr2024.RDS")


prob.surface <- function(DiveArea, draft) {
  
  if (DiveArea == 0 | DiveArea == 2) {
    loc <- "Mid-Atlantic U.S."
  } 
  
  if (DiveArea == 1) {
    loc = "Southeast U.S."
  }
  
  if (DiveArea == 3 | DiveArea == 4 ) {
    loc = "Northeast U.S."
  }
  
  if (DiveArea == 5) {
    loc =  "Cape Cod Bay"
  }
  
  mean.draft <- draft
  
  mean_prop <- depth_props_by_region[which(depth_props_by_region$location == loc & depth_props_by_region$layer == mean.draft),]
  mean_prop <- mean_prop$weighted_mean
  
  # 1 = above X m, 0 = below X m
  
  # Simulate the number of successes based on the mean proportion
  #n_successes <- rbinom(n = 100, size = 1, prob = mean_prop)
  #psurface <- sample(n_successes, size=1)
  psurface <- mean_prop
  return(psurface)
  
}  

###### FIXED PARAMS AND CONVERSIONS ######
ms2knot <- 1.94384  # convert m/s to knots
area <- 1e+08       # grid cell area 10 x 10 km in m^2

whale.len <- 13.5   # same generate distribution of whale lengths - or select from a distribution of known lengths
whale.len.sd <- 0.61 # Fortune et al. 2020 mean and sd for adult females
whale.width <- 3.5  # scale to length

radius.whale <- sqrt((whale.len * whale.width)/pi)  

#fixed mean whale speed
v.whale <- 0.389    # 0.389 m/s  #fixed speed for whales
#reflects average speed in se US...should look at literature for better values

#Random distribution for whale speed from weibull per 
shape.whale <- 1.48 # shape parameter (k) for the weibull distribution for whale speed
scale.whale <- 0.43 # scale parameter (L) for the weibull distribution for whale speed
#summary(rweibull(n=100000, shape=shape.whale, scale =scale.whale));
#hist(rweibull(n=1000, shape=shape.whale, scale =scale.whale));


###################### LOAD VESSEL TRACK DATA AND APPLY RISK EQUATIONS ###############
############ This first section of code reads in track and vessel data and calculates mortality risk
############ on a cell by cell basis and saves out monthly files for a given scenario

### customize suffix for different scenario outputs
output.folder <- "./model_outputs"
if (!dir.exists(output.folder)) {dir.create(output.folder)}

output.suffix <- "_riskMaps_baseline"

months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

for (i in 1:n.track.files) {
  # for (i in 1:12) {
  
  start.time <- Sys.time()
  
  track.file <- track.files[i]
  tracks <- readRDS(paste("./Data/ais_rds/", track.file, sep = "/"))
  yr <- substr(track.file, 13, 16)
  mon.num <- i %% 12 
  if (mon.num == 0) {mon.num <- 12}
  
  #merge in the vessels
  tracks.vessel <- merge(tracks, vessel.data, by = "vessels_id", all.x = TRUE)
  
  #filter missing LOA
  tracks.vessel.final <- tracks.vessel %>% filter(!is.na(report_loa))
  
  tracks.vessel.final <- left_join(tracks.vessel.final, vessel.cats)
  
  ###calcualte EE params
  tracks.vessel.final$v.boat.ms <- tracks.vessel.final$w.mu.speed / ms2knot #vessel speed in m/s
  
  ### Get whale density grid
  grid.density.mon <- data.frame(grid.density[, names(grid.density) %in% c("Id",months[mon.num], "DiveArea","ModelArea","depth15")])
  grid.density.mon[is.na(grid.density.mon)] <- 0
  grid.density.mon$gid <- as.numeric(as.character(grid.density.mon$Id))
  grid.density.mon$n.whale <- grid.density.mon[, names(grid.density.mon) == months[mon.num]]
  grid.region <- grid.density.mon[,names(grid.density.mon) %in% c("gid", "region", "n.whale", "DiveArea", "ModelArea","depth15")]
  tracks.vessel.final <- left_join(tracks.vessel.final, grid.region)
  # tracks.vessel.final <- tracks.vessel.final %>% rowwise() %>% mutate(p.mort = p.strike.mort(w.mu.speed))
  
  ## New mortality
  ## vessel classes: S/M 26-65 ft, L 65-350 ft, XL > 350 ft
  tracks.vessel.final <-mutate(tracks.vessel.final, reclass_Vessel_Size = case_when(
    report_loa > 106.7 ~ "XL",
    report_loa >19.8 &  report_loa <=106.7 ~ "L", 
    report_loa >7.9 &  report_loa <=19.8 ~ "M",
    report_loa <= 7.9 ~ "S",
  ))
  
  ###drop S vessels
  tracks.vessel.final <- tracks.vessel.final[tracks.vessel.final$reclass_Vessel_Size != "S",]
  
  ## draft/prop suction depth by class
  tracks.vessel.final<-mutate(tracks.vessel.final, mean.draft = case_when(
    reclass_Vessel_Size == "M" ~ 5,
    reclass_Vessel_Size == "L" ~ 5,
    reclass_Vessel_Size == "XL" ~ 15,
  ))
  
  ####steps to account for missing AIS data for medium size class
  #convert length from meters to feet
  tracks.vessel.final$loa_ft <- tracks.vessel.final$report_loa*3.28084
  
  ###LINK TO PERMIT:AIS vessel counts by econ_vessel_type and loa_ft
  tracks.vessel.final<-mutate(tracks.vessel.final, size_5ft = case_when(
    loa_ft >= 26 & loa_ft < 30 ~ "26-30",       
    loa_ft >= 30 & loa_ft < 35 ~ "30-35",        
    loa_ft >= 35 & loa_ft < 40 ~ "35-40",
    loa_ft >= 40 & loa_ft < 45 ~ "40-45",
    loa_ft >= 45 & loa_ft < 50 ~ "45-50",
    loa_ft >= 50 & loa_ft < 55 ~ "50-55",
    loa_ft >= 55 & loa_ft < 60 ~ "55-60",
    loa_ft >= 60 & loa_ft < 65 ~ "60-65",
    loa_ft >= 65 ~ ">65"
  ))
  
  #vessel count corrections by size/type 
  vess.corrections.init <- read.csv("./data/permit_vesscount_corrections_14Aug2023.csv", header = TRUE)
  
  #adds collapse.type to the tracks.vessel.final data set
  tracks.vessel.final <- left_join(tracks.vessel.final, vess.corrections.init[,c(1,2,5)])
  
  #### COLLAPSE ACROSS VESS CORRECTIONS COLLAPSE TPYE AND CALCULATE CORRECTION RATIOS
  vess.corrections.collapse <- vess.corrections.init %>% group_by(collapse.type, size_5ft) %>% 
    summarize(permit.count = sum(permit.vess.count), ais.count = sum(ais.vess.count)) %>%
    mutate(permit.to.ais = permit.count/ais.count)
  
  vess.corrections.collapse$permit.to.ais[vess.corrections.collapse$collapse.type == "Pilot"] <- 1
  vess.corrections.collapse$permit.to.ais[vess.corrections.collapse$collapse.type == "Exempt"] <- 1
  
  ### Join permit.to.ais into tracks.vessel.final - joins by 5ft class and collapse type
  tracks.vessel.final <- left_join(tracks.vessel.final, vess.corrections.collapse)
  
  ### permit.to.ais is NA for all large and xl vessels - set it to 1
  tracks.vessel.final$permit.to.ais[is.na(tracks.vessel.final$permit.to.ais)] <- 1
  
  tracks.vessel.final$PropSuctDepth <- tracks.vessel.final$mean.draft
  tracks.vessel.final$PropSuctDepth[tracks.vessel.final$reclass_Vessel_Size == "XL"] <- tracks.vessel.final$mean.draft[tracks.vessel.final$reclass_Vessel_Size == "XL"] * 2
  
  #faster p.mort
  newdata <- tracks.vessel.final[, c("w.mu.speed", "reclass_Vessel_Size")]
  colnames(newdata) <- c("Vess.Speed", "vess.cat_f")
  newdata$spe.HB <- "Not Humpback"
  tracks.vessel.final$p.mort <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
  
  #for consistency with pavoid and psurf - take a random draw from p.mort
  #tracks.vessel.final <- tracks.vessel.final %>% rowwise() %>% mutate(mort.yes = rbinom(n = 1, size = 1, p.mort))
  tracks.vessel.final <- tracks.vessel.final %>% rowwise() %>% mutate(mort.yes = p.mort)
  
  ##random whale params
  tracks.vessel.final$whale.len <-  rnorm(nrow(tracks.vessel.final), mean = whale.len, sd = whale.len.sd)
  tracks.vessel.final$v.whale <- rweibull(nrow(tracks.vessel.final), shape=shape.whale, scale =scale.whale)
  tracks.vessel.final$whale.width <- tracks.vessel.final$whale.len * 3.5/14.5 #based on original width/length ratio
  tracks.vessel.final$radius.whale <- sqrt((tracks.vessel.final$whale.len * tracks.vessel.final$whale.width)/pi)
  
  #### tracks.vessel.final.sim is the output data frame for modifying vessel track speed, etc.  
  
  ######## BASELINE - no changes to vessel traffic characteristics
  tracks.vessel.final.sim <- tracks.vessel.final 
  #######
  
  ################### Simulation: SPEED_ALL: Reduce all traffic >10 knots to around 10 knots
  #Speed restriction - coastwide - set all tracks > 10 to 10 [with error +/- 0.5 knot]
  #tracks.vessel.final.sim <- tracks.vessel.final
  #tracks.vessel.final.sim$w.mu.speed[tracks.vessel.final.sim$w.mu.speed > 10] <-  runif(length(tracks.vessel.final.sim$w.mu.speed[tracks.vessel.final.sim$w.mu.speed > 10]),9.5,10.5)
  #tracks.vessel.final.sim$v.boat.ms <- tracks.vessel.final.sim$w.mu.speed / ms2knot
  ###########################################################################
  
  #calculate encounter rate based on these factors
  tracks.vessel.final.sim <- tracks.vessel.final.sim %>% rowwise() %>% mutate(ER.fixed =
                                                                                getLambda(
                                                                                  vb = v.boat.ms,
                                                                                  vm = v.whale,
                                                                                  S = area,
                                                                                  r = radius.whale
                                                                                ))
  
  #apply prob mortality given vessel speed
  tracks.vessel.final.sim <- tracks.vessel.final.sim %>% mutate(EE.fixed = ER.fixed * tot.time.s) %>% mutate(EE.Mi = EE.fixed * mort.yes)
  
  #get surface presence
  tracks.vessel.final.sim <- tracks.vessel.final.sim %>% rowwise() %>% mutate(surface.yes = prob.surface(DiveArea, mean.draft))

  tracks.vessel.final.sim <- tracks.vessel.final.sim %>% rowwise() %>% mutate(avoid.no = 1 - p.avoid.boot(v.boat.ms, PropSuctDepth, depth15))
  
 
  #multiply to get likelihood of mortality - for one whale in the cell - for each vessel track
  tracks.vessel.final.sim <- tracks.vessel.final.sim %>% mutate(EE.M = EE.Mi * surface.yes * avoid.no,
                                                                EE.M.noavoid = EE.Mi * surface.yes)
  
  #for the corrected - treating the correction factor as multiplying the number of tracks of each type
  #that is - more than one trial...so this is a binomial with n.trials = permit.to.ais and want to 
  #know likelihood of one or more "successes" or moralities...so
  #pbinom(number of successes, number of trials, probability, lower.tail = FALSE so >x)
  #for classes with permit.to.ais = 1, this is equal to EE.M
  tracks.vessel.final.sim$EE.M.c <- pbinom(0, as.integer(tracks.vessel.final.sim$permit.to.ais), tracks.vessel.final.sim$EE.M, lower.tail = FALSE)
  tracks.vessel.final.sim$EE.M.c.na <- pbinom(0, as.integer(tracks.vessel.final.sim$permit.to.ais), tracks.vessel.final.sim$EE.M.noavoid, lower.tail = FALSE)
  
  
  #code speed gt10
  tracks.vessel.final.sim$speed.gt10 <- 0
  tracks.vessel.final.sim$speed.gt10[tracks.vessel.final.sim$w.mu.speed > 10] <- 1
  
  #get total risk of one or more strikes occurring - assume each track is an independent probability
  #so get the product of all the 1-EE.M.c to get the probability of no mortalities
  #then 1-that product to get prob of one or more.
  
  risk.cell <- tracks.vessel.final.sim %>% group_by(gid, reclass_Vessel_Size, speed.gt10) %>% summarize(n.segments = n(), mu.speed = mean(w.mu.speed), 
                                                                                                        prod.EE.M = prod(1-EE.M.c),
                                                                                                        prod.EE.M.na = prod(1-EE.M.c.na))
  ##this should be the probability of one or more strikes given the presence of 1 whale in the cell
  risk.cell$p.Strike <- 1-risk.cell$prod.EE.M
  risk.cell$p.Strike.na <- 1-risk.cell$prod.EE.M.na
  
  #join in the density
  risk.cell <- left_join(risk.cell, grid.density.mon)
  
  #density is probability a whale will be there - so this gives the overall risk of a whale strike occurring
  #could consider scaling this to 350 total whales...
  risk.cell$n.mort <- risk.cell$p.Strike * risk.cell$n.whale
  risk.cell$n.mort.na <- risk.cell$p.Strike.na * risk.cell$n.whale
  
  ##remove inshore cells and shallow cells
  ##drop grid cells with depth of > 3m
  risk.cell <- risk.cell[risk.cell$depth15 <= -3,]
  #sum(risk.cell$n.mort)
  
  #drop inshore cells see file ./data/excluded_inshore_cells.csv
  risk.cell <-  risk.cell[!( risk.cell$gid == 1886034 |  risk.cell$gid == 1911259 |  risk.cell$gid == 1913061 | 
                               risk.cell$gid == 1913062 |  risk.cell$gid == 1914863 |  risk.cell$gid == 1916664 | 
                               risk.cell$gid == 1923873),]
  
  #save out the monthly file
  out.file.name <- paste0(output.folder, month.abb[mon.num], "_", yr, output.suffix, ".csv", sep = "")
  write.csv(risk.cell, out.file.name, row.names = FALSE) 
  print(paste(yr, month.abb[mon.num], sep = "_"))
  print(Sys.time() - start.time)
  
} #month/track file



############ This section reads output files created with the code section above to make shapefiles for risk maps
output.folder <- "./model_outputs"
output.suffix <- "_riskMaps_baseline"

### First, get average (across years) baseline mortality by vessel size for each month
yr.list <- c(2017, 2018, 2019, 2021, 2022)

for (i in 1:12) {
  
  ## load files and average across years
  file.pattern <- paste0("riskMaps", month.abb[i],"_*", output.suffix, ".csv")
  risk.files <- dir(path = output.folder, pattern = glob2rx(file.pattern))
  
  for (yr.num in 1:(length(risk.files))) {
    
    file.name <- paste0(output.folder, "riskMaps",month.abb[i], "_", yr.list[yr.num], output.suffix, ".csv")
    yr.data <- read.csv(file.name, header = TRUE)
    
    yr.data <- yr.data[,c("gid", "reclass_Vessel_Size", "speed.gt10","n.mort", "n.mort.na")]
    names(yr.data) [4] <- paste0("mort", yr.list[yr.num])
    names(yr.data) [5] <- paste0("mort.na.", yr.list[yr.num])
    #yr.data$mort <- yr.list[yr.num]
    
    if (yr.num == 1) {
      yr.data.all <- yr.data
    }
    
    if (yr.num > 1) {
      yr.data.all <- full_join(yr.data.all, yr.data)
    }
  }
  
  yr.data.all[is.na(yr.data.all)] <- 0 
  
  
  #mean.mort <- rowMeans(yr.data.all[,4:8])
  #yr.data.all$mean.mort <- mean.mort
  #yr.data.all$mean.mort <- yr.data.all$mort2017
  
  mean.mort <- rowMeans(yr.data.all[,c("mort2017", "mort2018", "mort2019", "mort2021", "mort2022")])
  mean.mort.na <- rowMeans(yr.data.all[,c("mort.na.2017", "mort.na.2018", "mort.na.2019", "mort.na.2021", "mort.na.2022")])
  
  yr.data.all$mean.mort <- mean.mort
  yr.data.all$mean.mort.na <- mean.mort.na
  
  ##collapse across gt10
  yr.data.summary <- yr.data.all %>% group_by(gid, reclass_Vessel_Size) %>% summarize(mort = sum(mean.mort), mort.na = sum(mean.mort.na))
  
  
  if (i == 1) {
    base.M.out <- yr.data.summary[yr.data.summary$reclass_Vessel_Size == "M", c("gid", "mort", "mort.na")]
    names(base.M.out) <- c("gid", month.abb[i], paste0(month.abb[i], "_na"))
    
    base.L.out <- yr.data.summary[yr.data.summary$reclass_Vessel_Size == "L", c("gid", "mort", "mort.na")]
    names(base.L.out) <- c("gid", month.abb[i], paste0(month.abb[i], "_na"))
    
    base.XL.out <- yr.data.summary[yr.data.summary$reclass_Vessel_Size == "XL", c("gid", "mort", "mort.na")]
    names(base.XL.out) <- c("gid", month.abb[i], paste0(month.abb[i], "_na"))
    
  }
  
  if (i > 1) {
    base.M <- yr.data.summary[yr.data.summary$reclass_Vessel_Size == "M", c("gid", "mort", "mort.na")]
    names(base.M) <- c("gid", month.abb[i], paste0(month.abb[i], "_na"))
    
    base.M.out <- full_join(base.M.out, base.M)
    
    base.L <- yr.data.summary[yr.data.summary$reclass_Vessel_Size == "L", c("gid", "mort", "mort.na")]
    names(base.L) <-c("gid", month.abb[i], paste0(month.abb[i], "_na"))
    
    base.L.out <- full_join(base.L.out, base.L)
    
    base.XL <- yr.data.summary[yr.data.summary$reclass_Vessel_Size == "XL", c("gid", "mort", "mort.na")]
    names(base.XL) <- c("gid", month.abb[i], paste0(month.abb[i], "_na"))
    
    base.XL.out <- full_join(base.XL.out, base.XL)
    
  }
  
  base.M.out[is.na(base.M.out)] <- 0
  base.L.out[is.na(base.L.out)] <- 0
  base.XL.out[is.na(base.XL.out)] <- 0
  
  
} #end i - Month  

####join to grid input shapefiles for mapping
grid@data$gid <- as.numeric(grid@data$Id)

grid.M.out <- grid
grid.M.out@data <- left_join(grid.M.out@data, base.M.out)

grid.L.out <- grid
grid.L.out@data <- left_join(grid.L.out@data, base.L.out)

grid.XL.out <- grid
grid.XL.out@data <- left_join(grid.XL.out@data, base.XL.out)

writeOGR(grid.M.out, dsn = output.folder, layer = "baseline_risk_m_v2", driver = "ESRI Shapefile",
         overwrite_layer=TRUE)

writeOGR(grid.L.out, dsn = output.folder, layer = "baseline_risk_lg_v2", driver = "ESRI Shapefile",
         overwrite_layer=TRUE)

writeOGR(grid.XL.out, dsn = output.folder, layer = "baseline_risk_xl_v2", driver = "ESRI Shapefile",
         overwrite_layer=TRUE)



##############################################################################################
###################### Maps/plotting for manuscript ##########################################
##############################################################################################
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(dplyr)

output.folder <- "~/Desktop/git_RiskModel/RiskModel/All_Code_2024/model_outputs"

world.sfs <- ne_countries(scale="medium", returnclass = "sf")

baseline.risk.OGV <- st_read(paste0(output.folder, "/baseline_risk_xl_v2.shp"))
head(baseline.risk.OGV)

baseline.risk.OGV.sf <- st_transform(baseline.risk.OGV, crs=4329)

baseline.risk.OGV.sf.january <- baseline.risk.OGV.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jan)
names(baseline.risk.OGV.sf.january)[names(baseline.risk.OGV.sf.january) == 'Jan'] <- 'risk'
baseline.risk.OGV.sf.january$month <- "Jan"

baseline.risk.OGV.sf.july <- baseline.risk.OGV.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jul)
names(baseline.risk.OGV.sf.july)[names(baseline.risk.OGV.sf.july) == 'Jul'] <- 'risk'
baseline.risk.OGV.sf.july$month <- "Jul"

baseline.risk.OGV.sf.months <- rbind(baseline.risk.OGV.sf.january, baseline.risk.OGV.sf.july)
baseline.risk.OGV.sf.months$size_class <- "OGV"

########### LARGE ############################

baseline.risk.L <- st_read(paste0(output.folder, "/baseline_risk_lg_v2.shp"))
head(baseline.risk.L)

baseline.risk.L.sf <- st_transform(baseline.risk.L, crs=4329)


baseline.risk.L.sf.january <- baseline.risk.L.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jan)
names(baseline.risk.L.sf.january)[names(baseline.risk.L.sf.january) == 'Jan'] <- 'risk'
baseline.risk.L.sf.january$month <- "Jan"

baseline.risk.L.sf.july <- baseline.risk.L.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jul)
names(baseline.risk.L.sf.july)[names(baseline.risk.L.sf.july) == 'Jul'] <- 'risk'
baseline.risk.L.sf.july$month <- "Jul"


baseline.risk.L.sf.months <- rbind(baseline.risk.L.sf.january, baseline.risk.L.sf.july)
baseline.risk.L.sf.months$size_class <- "Large"

########### MEDIUM ############################

baseline.risk.M <- st_read(paste0(output.folder, "/baseline_risk_m_v2.shp"))
head(baseline.risk.M)

baseline.risk.M.sf <- st_transform(baseline.risk.M, crs=4329)

baseline.risk.M.sf.january <- baseline.risk.M.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jan)
names(baseline.risk.M.sf.january)[names(baseline.risk.M.sf.january) == 'Jan'] <- 'risk'
baseline.risk.M.sf.january$month <- "Jan"

baseline.risk.M.sf.july <- baseline.risk.M.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jul)
names(baseline.risk.M.sf.july)[names(baseline.risk.M.sf.july) == 'Jul'] <- 'risk'
baseline.risk.M.sf.july$month <- "Jul"


baseline.risk.M.sf.months <- rbind(baseline.risk.M.sf.january, baseline.risk.M.sf.july)

baseline.risk.M.sf.months$size_class <- "Small/Medium"

baseline.risk.allClasses <- rbind(baseline.risk.OGV.sf.months, baseline.risk.L.sf.months, baseline.risk.M.sf.months)

baseline.risk.allClasses$size_class <- factor(baseline.risk.allClasses$size_class, levels = c(
  "OGV", "Large", "Small/Medium"))

baseline.risk.allClasses$risk_times1000 <- baseline.risk.allClasses$risk*1000

mybreaks.times1000 = c(0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 7.0) #120

library(scales); q_colors = 20; v_colors = turbo(q_colors)

v_colors <- c("#30123BFF", "#455ED2FF",   "#18DEC1FF", "#B5F836FF", "#EFCD3AFF" ,"#FD8A26FF" , "#C82803FF",  "#7A0403FF")

#[1] "#30123BFF" "#3F3994FF" "#455ED2FF" "#4681F7FF" "#3AA2FCFF" "#23C3E4FF" "#18DEC1FF" "#2CF09EFF" "#5BFB72FF" "#8EFF49FF"
#[11] "#B5F836FF" "#D6E635FF" "#EFCD3AFF" "#FCB036FF" "#FD8A26FF" "#F36215FF" "#E14209FF" "#C82803FF" "#A51301FF" "#7A0403FF"

names(v_colors) <- levels(baseline.risk.allClasses$risk_factor)

baseline.risk.allClasses <- baseline.risk.allClasses %>%
  mutate(risk_factor = case_when(
    risk_times1000 >= 0 & risk_times1000 < 0.000001 ~ '0 - 0.000001',
    risk_times1000 >= 0.000001 & risk_times1000 < 0.00001 ~ '0.000001 - 0.00001',
    risk_times1000 >= 0.00001 & risk_times1000 < 0.0001 ~ '0.00001 - 0.0001',
    risk_times1000 >= 0.0001 & risk_times1000 < 0.001 ~ '0.0001 - 0.001',
    risk_times1000 >= 0.001 & risk_times1000 < 0.01 ~ '0.001 - 0.01',
    risk_times1000 >= 0.01 & risk_times1000 < 0.1 ~ '0.01 - 0.1',
    risk_times1000 >= 0.1 & risk_times1000 < 1.0 ~ '0.1 - 1.0',
    risk_times1000 >= 1.0 ~ '1.0 - 76.0'))

baseline.risk.allClasses$risk_factor <- factor(baseline.risk.allClasses$risk_factor, levels=
                                                 c('0 - 0.000001', '0.000001 - 0.00001', '0.00001 - 0.0001', '0.0001 - 0.001',
                                                   '0.001 - 0.01', '0.01 - 0.1','0.1 - 1.0', '1.0 - 76.0'))



baseline.risk.allClasses$Month <- baseline.risk.allClasses$month
baseline.risk.allClasses$Month[baseline.risk.allClasses$month == "Jan"] <- "January"
baseline.risk.allClasses$Month[baseline.risk.allClasses$month == "Jul"] <- "July"

risk_all_facetGrid <- ggplot()+
  geom_sf(data=baseline.risk.allClasses, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  facet_grid(size_class~Month) +
  theme(strip.placement = "outside") +
  #,  strip.position = c("top", "left")) +
  theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  
#risk_all.1000.3

risk_all_facetGrid

ggsave(paste0(output.folder, "riskMap_all_1000_facetGrid.png"), risk_all_facetGrid, dpi=300, height=12, width=10)

risk_spatial_OGV_jan <- baseline.risk.allClasses %>%
  filter(size_class == "OGV") %>%
  filter(Month == "January")

risk_OGV_jan <- ggplot()+
  geom_sf(data=risk_spatial_OGV_jan, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_OGV_jul <- baseline.risk.allClasses %>%
  filter(size_class == "OGV") %>%
  filter(Month == "July")

risk_OGV_jul <- ggplot()+
  geom_sf(data=risk_spatial_OGV_jul, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_L_jan <- baseline.risk.allClasses %>%
  filter(size_class == "Large") %>%
  filter(Month == "January")

risk_L_jan <- ggplot()+
  geom_sf(data=risk_spatial_L_jan, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_L_jul <- baseline.risk.allClasses %>%
  filter(size_class == "Large") %>%
  filter(Month == "July")

risk_L_jul <- ggplot()+
  geom_sf(data=risk_spatial_L_jul, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
# theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_M_jan <- baseline.risk.allClasses %>%
  filter(size_class == "Small/Medium") %>%
  filter(Month == "January")

risk_M_jan <- ggplot()+
  geom_sf(data=risk_spatial_M_jan, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month)# +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_M_jul <- baseline.risk.allClasses %>%
  filter(size_class == "Small/Medium") %>%
  filter(Month == "July")

risk_M_jul <- ggplot()+
  geom_sf(data=risk_spatial_M_jul, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk.OGV <- ggarrange(risk_OGV_jan, risk_OGV_jul, ncol=2, common.legend = T, legend="none")
ggsave(paste0(output.folder, "risk.OGV.png"), risk.OGV, dpi=300, height=3, width=6)

risk.L <- ggarrange(risk_L_jan, risk_L_jul, ncol=2, common.legend = T, legend="none")
ggsave(paste0(output.folder, "risk.L.png"), risk.L, dpi=300, height=3, width=6)

risk.M <- ggarrange(risk_M_jan, risk_M_jul, ncol=2, common.legend = T, legend="none")
ggsave(paste0(output.folder, "risk.M.png"), risk.M, dpi=300, height=3, width=6)

ggsave(paste0(output.folder, "risk.legend.png"), risk_OGV_jan, dpi=300, height=3, width=6)



######### NO AVOIDANCE #######################

#OGV
baseline.risk.OGV.sf.january.na <- baseline.risk.OGV.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jan_na)
names(baseline.risk.OGV.sf.january.na)[names(baseline.risk.OGV.sf.january.na) == 'Jan_na'] <- 'risk'
baseline.risk.OGV.sf.january.na$month <- "Jan_na"

baseline.risk.OGV.sf.july.na <- baseline.risk.OGV.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jul_na)
names(baseline.risk.OGV.sf.july.na)[names(baseline.risk.OGV.sf.july.na) == 'Jul_na'] <- 'risk'
baseline.risk.OGV.sf.july.na$month <- "Jul_na"


baseline.risk.OGV.sf.months.na <- rbind(baseline.risk.OGV.sf.january.na, baseline.risk.OGV.sf.july.na)
baseline.risk.OGV.sf.months.na$size_class <- "OGV"

# Large 
baseline.risk.L.sf.january.na <- baseline.risk.L.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jan_na)
names(baseline.risk.L.sf.january.na)[names(baseline.risk.L.sf.january.na) == 'Jan_na'] <- 'risk'
baseline.risk.L.sf.january.na$month <- "Jan_na"

baseline.risk.L.sf.july.na <- baseline.risk.L.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jul_na)
names(baseline.risk.L.sf.july.na)[names(baseline.risk.L.sf.july.na) == 'Jul_na'] <- 'risk'
baseline.risk.L.sf.july.na$month <- "Jul_na"


baseline.risk.L.sf.months.na <- rbind(baseline.risk.L.sf.january.na, baseline.risk.L.sf.july.na)
baseline.risk.L.sf.months.na$size_class <- "Large"

# Medium 
baseline.risk.M.sf.january.na <- baseline.risk.M.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jan_na)
names(baseline.risk.M.sf.january.na)[names(baseline.risk.M.sf.january.na) == 'Jan_na'] <- 'risk'
baseline.risk.M.sf.january.na$month <- "Jan_na"

baseline.risk.M.sf.july.na <- baseline.risk.M.sf %>%
  select(Id, gridcode, Shape_Leng, Shape_Area, ec22_study, gid, Jul_na)
names(baseline.risk.M.sf.july.na)[names(baseline.risk.M.sf.july.na) == 'Jul_na'] <- 'risk'
baseline.risk.M.sf.july.na$month <- "Jul_na"


baseline.risk.M.sf.months.na <- rbind(baseline.risk.M.sf.january.na, baseline.risk.M.sf.july.na)

baseline.risk.M.sf.months.na$size_class <- "Small/Medium"

baseline.risk.allClasses.na <- rbind(baseline.risk.OGV.sf.months.na, baseline.risk.L.sf.months.na, baseline.risk.M.sf.months.na)

baseline.risk.allClasses.na$size_class <- factor(baseline.risk.allClasses.na$size_class, levels = c(
  "OGV", "Large", "Small/Medium"))

baseline.risk.allClasses.na$risk_times1000 <- baseline.risk.allClasses.na$risk*1000

mybreaks.times1000 = c(0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 7.0) #120

library(scales); q_colors = 20; v_colors = turbo(q_colors)

v_colors <- c("#30123BFF", "#455ED2FF",   "#18DEC1FF", "#B5F836FF", "#EFCD3AFF" ,"#FD8A26FF" , "#C82803FF",  "#7A0403FF")

#[1] "#30123BFF" "#3F3994FF" "#455ED2FF" "#4681F7FF" "#3AA2FCFF" "#23C3E4FF" "#18DEC1FF" "#2CF09EFF" "#5BFB72FF" "#8EFF49FF"
#[11] "#B5F836FF" "#D6E635FF" "#EFCD3AFF" "#FCB036FF" "#FD8A26FF" "#F36215FF" "#E14209FF" "#C82803FF" "#A51301FF" "#7A0403FF"

names(v_colors) <- levels(baseline.risk.allClasses.na$risk_factor)

baseline.risk.allClasses.na <- baseline.risk.allClasses.na %>%
  mutate(risk_factor = case_when(
    risk_times1000 >= 0 & risk_times1000 < 0.000001 ~ '0 - 0.000001',
    risk_times1000 >= 0.000001 & risk_times1000 < 0.00001 ~ '0.000001 - 0.00001',
    risk_times1000 >= 0.00001 & risk_times1000 < 0.0001 ~ '0.00001 - 0.0001',
    risk_times1000 >= 0.0001 & risk_times1000 < 0.001 ~ '0.0001 - 0.001',
    risk_times1000 >= 0.001 & risk_times1000 < 0.01 ~ '0.001 - 0.01',
    risk_times1000 >= 0.01 & risk_times1000 < 0.1 ~ '0.01 - 0.1',
    risk_times1000 >= 0.1 & risk_times1000 < 1.0 ~ '0.1 - 1.0',
    risk_times1000 >= 1.0 ~ '1.0 - 128.0'))

baseline.risk.allClasses.na$risk_factor <- factor(baseline.risk.allClasses.na$risk_factor, levels=
                                                    c('0 - 0.000001', '0.000001 - 0.00001', '0.00001 - 0.0001', '0.0001 - 0.001',
                                                      '0.001 - 0.01', '0.01 - 0.1','0.1 - 1.0', '1.0 - 128.0'))



baseline.risk.allClasses.na$Month <- baseline.risk.allClasses.na$month
baseline.risk.allClasses.na$Month[baseline.risk.allClasses.na$month == "Jan"] <- "January"
baseline.risk.allClasses.na$Month[baseline.risk.allClasses.na$month == "Jul"] <- "July"

risk_all_facetGrid.na <- ggplot()+
  geom_sf(data=baseline.risk.allClasses.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  facet_grid(size_class~Month) +
  theme(strip.placement = "outside") +
  #,  strip.position = c("top", "left")) +
  theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  
#risk_all.1000.3

risk_all_facetGrid.na

ggsave(paste0(output.folder, "riskMap_noAvoid_all_1000_facetGrid.png"), risk_all_facetGrid.na, dpi=300, height=12, width=10)

risk_spatial_OGV_jan.na <- baseline.risk.allClasses.na %>%
  filter(size_class == "OGV") %>%
  filter(Month == "Jan_na")

risk_OGV_jan.na <- ggplot()+
  geom_sf(data=risk_spatial_OGV_jan.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_OGV_jul.na <- baseline.risk.allClasses.na %>%
  filter(size_class == "OGV") %>%
  filter(Month == "Jul_na")

risk_OGV_jul.na <- ggplot()+
  geom_sf(data=risk_spatial_OGV_jul.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_L_jan.na <- baseline.risk.allClasses.na %>%
  filter(size_class == "Large") %>%
  filter(Month == "Jan_na")

risk_L_jan.na <- ggplot()+
  geom_sf(data=risk_spatial_L_jan.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_L_jul.na <- baseline.risk.allClasses.na %>%
  filter(size_class == "Large") %>%
  filter(Month == "Jul_na")

risk_L_jul.na <- ggplot()+
  geom_sf(data=risk_spatial_L_jul.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
# theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_M_jan.na <- baseline.risk.allClasses.na %>%
  filter(size_class == "Small/Medium") %>%
  filter(Month == "Jan_na")

risk_M_jan.na <- ggplot()+
  geom_sf(data=risk_spatial_M_jan.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month)# +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk_spatial_M_jul.na <- baseline.risk.allClasses.na %>%
  filter(size_class == "Small/Medium") %>%
  filter(Month == "Jul_na")

risk_M_jul.na <- ggplot()+
  geom_sf(data=risk_spatial_M_jul.na, aes(fill=risk_factor), col=NA) +
  #scale_fill_viridis_d(option="turbo",na.value = "gray90",  name ="Risk", breaks=mybreaks.times1000) + 
  scale_fill_manual(name="Mortality Risk", values=v_colors) +
  geom_sf(data=world.sfs, col="black", fill="gray75") +
  coord_sf(ylim=c(25,46), xlim=c(-82, -57))+
  theme(axis.text=element_text(size=24))  +
  theme_minimal() +
  #facet_grid(size_class~Month) +
  theme(strip.placement = "outside") #+
#,  strip.position = c("top", "left")) +
#theme(axis.text=element_text(size=12), strip.text = element_text(size=12))  

risk.OGV.na <- ggarrange(risk_OGV_jan.na, risk_OGV_jul.na, ncol=2, common.legend = T, legend="none")
ggsave(paste0(output.folder, "risk.na.OGV.png"), risk.OGV.na, dpi=300, height=3, width=6)

risk.L.na <- ggarrange(risk_L_jan.na, risk_L_jul.na, ncol=2, common.legend = T, legend="none")
ggsave(paste0(output.folder, "risk.na.L.png"), risk.L.na, dpi=300, height=3, width=6)

risk.M.na <- ggarrange(risk_M_jan.na, risk_M_jul.na, ncol=2, common.legend = T, legend="none")
ggsave(paste0(output.folder, "risk.na.M.png"), risk.M.na, dpi=300, height=3, width=6)

ggsave(paste0(output.folder, "risk.na.legend.png"), risk_OGV_jan.na, dpi=300, height=3, width=6)
