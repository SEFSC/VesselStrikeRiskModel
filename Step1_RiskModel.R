library(tidyverse)
library(rgdal)
#library(RStoolbox)
library(lubridate)
library(foreign)
library(sf)
library(dplyr)
library(parallel)
#install.packages("foreach")
#install.packages("doParallel")
library("foreach")
library("doParallel")

setwd("~/Desktop/git_RiskModel/RiskModel/All_Code_2024")

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
summary(rweibull(n=100000, shape=shape.whale, scale =scale.whale));hist(rweibull(n=1000, shape=shape.whale, scale =scale.whale));

LethalityCurveGLM <- readRDS("./Data/LethalityCurveGLM_26Apr24.RDS") # Garrison et al. 2024

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

p.strike.mort <- function(v.boat.knots, reclass_Vessel_Size, spe.HB) {
  
  Vess.Speed =  v.boat.knots #vessel speed (knots)
  vess.cat_f = reclass_Vessel_Size
  spe.HB = "Not Humpback" #humpback or not 
  
  newdata <- data.frame(Vess.Speed,vess.cat_f, spe.HB)
  
  p.strike.mort  <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
}


####### Probability of Avoidance ###########

p.avoid <- function(v.boat.ms, PropSuctDepth, depth, class) {
  ShipSpeed <- v.boat.ms #speed in m/s-1
  
  class <- class
  
  if (class == "M" | class == "L") {
    depth.lim <- 5
  } else if (class == "XL") {
    depth.lim <- 15
    #PropSuctDepth <- 1*mean.draft
  }
  
  WhaleDepth.start <- (runif(1) * (depth.lim - 0))  #select random whale start depth within limits of draft
  
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
  #print(ReactionDistance)
  
  WhaleDepth <- WhaleDepth.start + (DescentRate * (ReactionDistance/ShipSpeed))
  
  ## Depth limit
  PropSuctDepth <- PropSuctDepth
  
  PosDepth <- depth*-1
  
  WhaleHeight <- 3.0 #Moore et al. (2005) Jour. of Cet. Res. & Man.
  
  DepthLimit <- PropSuctDepth + WhaleHeight
  p.avoid <- ifelse(PosDepth > DepthLimit & WhaleDepth > PropSuctDepth, 1, 0)
}


##### Probability in the Upper 5 or 15 m 
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
  psurface <- rbinom(n = 1, size = 1, prob = mean_prop)
  return(psurface)
  
}

############## Get clean vessel data which includes all possible imputed values ##################

vessel.data <- readRDS("./Data/ais_rds/all_vessels_clean_distinct_2017_2022.rds")

#AIS grid shapefile
grid <- readOGR(dsn="./Data/shapefiles", layer="ec22_ease2_n_10km_pl")


grid.density <- readRDS("./Data/shapefiles/new_NARW_density_shp_v3.RDS")

#Get vessel categories for summaries and groupings of some rare vessel types
vessel.cats <- read.csv("./Data/vessel_categories.csv")



####################################################################################

###################### LOAD VESSEL TRACK DATA AND APPLY RISK EQUATIONS ###############
o.folder <- "./model_outputs"
### create directories if needed
if (!dir.exists(o.folder)) {dir.create(o.folder)}

library(magrittr)
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

years <- c("2017", "2018", "2019", "2021", "2022") 

for (k in 1:length(years)) {
  
  year <- years[k]
  
  track.files <- dir(path = "./Data/ais_rds/", pattern = glob2rx(paste0("narw_tracks_", year, "*.rds")))
  n.track.files <- length(track.files)
  
  boat.classes <- c("L", "XL") 
  
  for(y in 1:length(unique(boat.classes))) {
    
    months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
    
    #stuff for the bootstrap
    n.whales <- 350
    n.boot <- 1000
    
    tot.strike.avoid.slow.none <- mat.or.vec(n.boot, 12)
    tot.strike.avoid.slow.all <- mat.or.vec(n.boot, 12)
    
    tot.strike.EE.M.slow.none <-  mat.or.vec(n.boot, 12)
    tot.strike.EE.M.slow.all <-  mat.or.vec(n.boot, 12)
    
    tot.risk.out <- data.frame(id = c(1:n.track.files), yr = NA, mon = NA, tot.risk = NA, tot.whales = NA)
    
    class <- boat.classes[y]  
    output.folder <- paste0(o.folder, "/Class", class)
    if (!dir.exists(output.folder)) {dir.create(output.folder)}
    output.bootstrap.folder <- paste0(output.folder, "/bootstrap_model_output_", year, "/")
    if (!dir.exists(output.bootstrap.folder)) {dir.create(output.bootstrap.folder)}
    
    start_month <- 1
    end_month <- 12
    
    ###This loop cycles through AIS data for months of one year
    #now parallel instead of for loop
    foreach(mon.num=seq(start_month, end_month), .combine="c", .packages = c("dplyr", "tidyr")) %dopar% {
      #for (i in c(1:12)) { #12
      library(dplyr)
      library(tidyr)
      
      start.time <- Sys.time()
      
      mon.num=mon.num
      track.file <- track.files[mon.num]
      tracks <- readRDS(paste("./Data/ais_rds", track.file, sep = "/"))
      
      yr <- substr(track.file, 13, 16)
      #mon.num <- month %% 12 
      if (mon.num == 0) {mon.num <- 12}
      
      #merge in the vessels
      tracks.vessel <- merge(tracks, vessel.data, by = "vessels_id", all.x = TRUE)
      
      #filter missing LOA and <10m reported length
      tracks.vessel <- tracks.vessel %>% dplyr::filter(!is.na(report_loa)) %>% dplyr::filter(report_loa > 10)
      
      tracks.vessel.final <- tracks.vessel
      
      #merge vessel categories to standardize name and group rare vessel types
      tracks.vessel.final <- left_join(tracks.vessel.final, vessel.cats)
      
      tracks.vessel.final$radius.boat <- sqrt((tracks.vessel.final$report_loa * tracks.vessel.final$beam)/pi)
      
      tracks.vessel.final<-mutate(tracks.vessel.final, reclass_Vessel_Size = case_when(
        report_loa >= 106.68 ~ "XL",
        report_loa >= 19.812 &  report_loa < 106.68 ~ "L",
        report_loa >= 10.668 &  report_loa < 19.812 ~ "M",
        report_loa < 10.668 ~ "S"
      )) 
      
      
      tracks.vessel.final <- tracks.vessel.final[which(tracks.vessel.final$reclass_Vessel_Size == class),]
      
      ## draft by class
      if (class == "M" | class == "L") {
        mean.draft <- 5
        PropSuctDepth <- 1*mean.draft
      } else if (class == "XL") {
        mean.draft <- 15
        PropSuctDepth <- 2*mean.draft 
      }
      print(mean.draft)
      print(PropSuctDepth)
      
      ###calculate EE params
      
      ### Note: next we would like to put random parameters in here too..so that each track has a different whale speed, size, etc. That way you still end up dealing with uncertainty in those parameters.
      tracks.vessel.final$whale.len <-  rnorm(nrow(tracks.vessel.final), mean = whale.len, sd = whale.len.sd)
      tracks.vessel.final$v.whale <- rweibull(nrow(tracks.vessel.final), shape=shape.whale, scale =scale.whale)
      tracks.vessel.final$whale.width <- tracks.vessel.final$whale.len * 3.5/14.5 #based on original width/length ratio
      tracks.vessel.final$radius.whale <- sqrt((tracks.vessel.final$whale.len * tracks.vessel.final$whale.width)/pi)
      
      ### Get whale density grid
      grid.density.mon <- data.frame(grid.density[, names(grid.density) %in% c("Id",months[mon.num], "region", "DiveArea","MidAtlArea", "depth", "depth15", "minZ15")])
      grid.density.mon[is.na(grid.density.mon)] <- 0
      grid.density.mon$gid <- as.numeric(as.character(grid.density.mon$Id))
      grid.density.mon$n.whale <- grid.density.mon[, names(grid.density.mon) == months[mon.num]]
      #remove cells with depths less than 3 m (height of adult whale)
      grid.density.mon <- grid.density.mon[!(grid.density.mon$depth15 >= 3.0),]
      grid.region <- grid.density.mon[,names(grid.density.mon) %in% c("gid", "region", "n.whale", "DiveArea", "MidAtlArea", "depth", "depth15", "minZ15")]
      
      tracks.vessel.final <- left_join(tracks.vessel.final, grid.region)
      
      #drop inshore cells see file ./data/excluded_inshore_cells.csv
      tracks.vessel.final <- tracks.vessel.final[!(tracks.vessel.final$gid == 1886034 | tracks.vessel.final$gid == 1911259 | tracks.vessel.final$gid == 1913061 | tracks.vessel.final$gid == 1913062 | tracks.vessel.final$gid == 1914863 | tracks.vessel.final$gid == 1916664 | tracks.vessel.final$gid == 1923873),]
      
      tracks.vessel.final$day <- format(as.Date(tracks.vessel.final$s.date, format="%Y-%m-%d"), "%d")
      
      tracks.vessel.final$uniqueID <- seq(1:nrow(tracks.vessel.final))
      
      ## speed simulation
      ## first save original speed (status quo/slow none)
      tracks.vessel.final$w.mu.speed.slow.none <- tracks.vessel.final$w.mu.speed
      
      #nrow(tracks.vessel.final) 
       
      ### Select vessels above 10 knots  (exclude Exempt vessels)
      fast.vessels <- tracks.vessel.final[tracks.vessel.final$w.mu.speed >10 & tracks.vessel.final$econ_vessel_type != "Exempt",]
      
      n.fast.tracks <- dim(fast.vessels) [1]
      ## remove these from the data - because we will put them back in later
      fast.vessels$w.mu.speed.slow.all <- runif(n.fast.tracks, 9.5, 10.0)
      
      tracks.vessel.final.sim <- tracks.vessel.final[!(tracks.vessel.final$w.mu.speed >10 & tracks.vessel.final$econ_vessel_type != "Exempt"),] 
      tracks.vessel.final.sim$w.mu.speed.slow.all <- tracks.vessel.final.sim$w.mu.speed.slow.none

      tracks.vessel.final <- rbind(tracks.vessel.final.sim,fast.vessels)
      
      tracks.vessel.final$v.boat.ms.slow.none <- tracks.vessel.final$w.mu.speed.slow.none/ms2knot
      tracks.vessel.final$v.boat.ms.slow.all <- tracks.vessel.final$w.mu.speed.slow.all/ms2knot
      
      ### end speed sims ###
      
      #test speed sim 
      #hist(tracks.vessel.final$w.mu.speed)
      #hist(tracks.vessel.final$w.mu.speed.slow.none) # should be the same as above
      #hist(tracks.vessel.final$w.mu.speed.slow.all)
      
      #mean(tracks.vessel.final$w.mu.speed)
      #mean(tracks.vessel.final$w.mu.speed.slow.none) # should be the same as above
      #mean(tracks.vessel.final$w.mu.speed.slow.all) # should be lower than above
      
      
      
      #calculate probability of mortality given vessel strike for each track
      #tracks.vessel.final.func <- tracks.vessel.final %>% rowwise() %>% mutate(p.mort.old = p.strike.mort.old(w.mu.speed))
      newdata <- tracks.vessel.final[, c("w.mu.speed.slow.none", "reclass_Vessel_Size")]
      colnames(newdata) <- c("Vess.Speed", "vess.cat_f")
      newdata$spe.HB <- "Not Humpback"
      tracks.vessel.final$p.mort.slow.none <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
   
      newdata <- tracks.vessel.final[, c("w.mu.speed.slow.all", "reclass_Vessel_Size")]
      colnames(newdata) <- c("Vess.Speed", "vess.cat_f")
      newdata$spe.HB <- "Not Humpback"
      tracks.vessel.final$p.mort.slow.all <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
      
      #hist(tracks.vessel.final$p.mort.slow.none)
      #hist(tracks.vessel.final$p.mort.slow.all)
      
      #mean(tracks.vessel.final$p.mort.slow.none)
      #mean(tracks.vessel.final$p.mort.slow.all)
      
      
      #tracks.vessel.final.func <- tracks.vessel.final %>% dplyr::rowwise() %>% 
      # dplyr::mutate(p.mort = p.strike.mort(v.boat.knots=w.mu.speed, 
      #                                        reclass_Vessel_Size=reclass_Vessel_Size, 
      #                                       spe.HB="Not Humpback"))
      
      
      #calculate encounter rate for each track given the presence of a random whale - this step takes some time
      tracks.vessel.final <- tracks.vessel.final[!(is.na(tracks.vessel.final$vessels_id)),]
      tracks.vessel.final <- tracks.vessel.final %>% dplyr::rowwise() %>% dplyr::mutate(ER.fixed.slow.none =
                                                                                          getLambda(
                                                                                            vb = v.boat.ms.slow.none,
                                                                                            vm = v.whale,
                                                                                            S = area,
                                                                                            r = radius.whale
                                                                                          ))
      
      
      tracks.vessel.final <- tracks.vessel.final %>% dplyr::rowwise() %>% dplyr::mutate(ER.fixed.slow.all =
                                                                                          getLambda(
                                                                                            vb = v.boat.ms.slow.all,
                                                                                            vm = v.whale,
                                                                                            S = area,
                                                                                            r = radius.whale
                                                                                          ))
      
      
      # hist(tracks.vessel.final$ER.fixed.slow.none)
      #hist(tracks.vessel.final$ER.fixed.slow.all)
      
      #mean(tracks.vessel.final$ER.fixed.slow.none)
      #mean(tracks.vessel.final$ER.fixed.slow.all)
      
      
      #calculate mortality probability - does not include avoidance or surface behavior at this point
      #Time is dependent on speed --> we should look to make this distance
      
      tracks.vessel.final <- tracks.vessel.final %>% dplyr::mutate(EE.fixed.slow.none = ER.fixed.slow.none * tot.time.s) %>% dplyr::mutate(EE.Mi.slow.none = EE.fixed.slow.none * p.mort.slow.none)
      
      tracks.vessel.final <- tracks.vessel.final %>% dplyr::mutate(EE.fixed.slow.all = ER.fixed.slow.all * tot.time.s) %>% dplyr::mutate(EE.Mi.slow.all = EE.fixed.slow.all * p.mort.slow.all)

      # Test encounter rate 
      
      #hist(tracks.vessel.final$ER.fixed.slow.none)
      #hist(tracks.vessel.final$ER.fixed.slow.all)
      
      #mean(tracks.vessel.final$ER.fixed.slow.none)
      #mean(tracks.vessel.final$ER.fixed.slow.all)
      
      #  hist(tracks.vessel.final$EE.Mi.slow.none)
      #hist(tracks.vessel.final$EE.Mi.slow.all)
      
      #mean(tracks.vessel.final$EE.Mi.slow.none)
      #mean(tracks.vessel.final$EE.Mi.slow.all)
      
      print(Sys.time() - start.time)
      ### ~1.5 minutes for this step
      
      start <- Sys.time()
      ### START BOOTSTRAP LOOP
      for (boot in 1:n.boot) {
        
        ### 1) Make a data.frame to hold individual whales
        whales <- data.frame(whale = c(1:n.whales))
        
        ### 2) Put the whales in the cells - could add a time step here to move the whales on a random walk
        ### selection of cells weighted by predicted density
        
        ran.cells <- grid.density.mon[sample(nrow(grid.density.mon), n.whales, prob = grid.density.mon[,2], replace = TRUE), ]
        whales$gid <- ran.cells$gid
        whales$DiveArea <- ran.cells$DiveArea
        
        ### 3) Calculate risk for each whale - individuals stays in the same cell for a whole month
        whales$tot.EE.slow.none <- NA
        whales$tot.EE.M.slow.none <- NA
        whales$avoid.tracks.slow.none <- NA #number of positive avoidances
        whales$mort.n.slow.none <- NA #number of positive moralities
        whales$mort.avoid.slow.none <- NA #number of morts after avoided
        
        whales$tot.EE.slow.all<- NA
        whales$tot.EE.M.slow.all <- NA
        #whales$tot.tracks.e <- NA
        whales$avoid.tracks.slow.all <- NA #number of positive avoidances
        whales$mort.n.slow.all <- NA #number of positive moralities
        whales$mort.avoid.slow.all <- NA #number of morts after avoided
        
        for(w in 1:n.whales) {
          
          #get the cell for whale w
          #get the tracks - currently this is all tracks for a whole month
          tracks.cell <- data.frame(tracks.vessel.final[tracks.vessel.final$gid == whales$gid[w],])
          n.tracks <- dim(tracks.cell) [1]
          
          if (n.tracks> 0) {
            
            ##apply the p.surface function - get a random probability at the surface depending on region
            
            p.surf <- mapply(prob.surface, tracks.cell$DiveArea, mean.draft)  
            
            tracks.cell$EE.M.slow.none <- tracks.cell$EE.Mi.slow.none * p.surf
            tracks.cell$EE.M.slow.all <- tracks.cell$EE.Mi.slow.all * p.surf
            
            #generate hits on a per track basis - this accounts for potential variation in whale size, speed, etc.
            hityes.slow.none <- rbinom(size = 1, n = n.tracks, prob = tracks.cell$EE.M.slow.none)
            hityes.slow.all <- rbinom(size = 1, n = n.tracks, prob = tracks.cell$EE.M.slow.all)
            
            
            avoidyes.slow.none <- mapply(p.avoid, tracks.cell$v.boat.ms.slow.none, PropSuctDepth, tracks.cell$depth15, class=class)
            avoidyes.slow.all <- mapply(p.avoid, tracks.cell$v.boat.ms.slow.all, PropSuctDepth, tracks.cell$depth15, class=class)
            
          }
          
          if (n.tracks == 0) {
            hityes.slow.none <- 0
            avoidyes.slow.none <- 0
            
            hityes.slow.all <- 0
            avoidyes.slow.all <- 0 
          }
          
          
          whales$tot.EE.slow.none[w] <- sum(tracks.cell$EE.fixed.slow.none)
          whales$tot.EE.M.slow.none[w] <- sum(tracks.cell$EE.M.slow.none)
          whales$tot.tracks.slow.none[w] <- n.tracks
          whales$avoid.tracks.slow.none[w] <- sum(avoidyes.slow.none) #number of tracks avoided
          whales$mort.n.slow.none[w] <- sum(hityes.slow.none) #number of positive mortalities
          whales$mort.avoid.slow.none[w] <- sum(hityes.slow.none * (1-avoidyes.slow.none)) #1 means the track was avoided
          
          whales$tot.EE.slow.all[w] <- sum(tracks.cell$EE.fixed.slow.all)
          whales$tot.EE.M.slow.all[w] <- sum(tracks.cell$EE.M.slow.all)
          whales$tot.tracks.slow.all[w] <- n.tracks
          whales$avoid.tracks.slow.all[w] <- sum(avoidyes.slow.all) #number of tracks avoided
          whales$mort.n.slow.all[w] <- sum(hityes.slow.all) #number of positive mortalities
          whales$mort.avoid.slow.all[w] <- sum(hityes.slow.all * (1-avoidyes.slow.all)) #1 means the track was avoided
          
        } #whale
        
        ### Save output from the bootstrap for all whales
        boot_no <- boot
        write.csv(whales, paste0(output.bootstrap.folder,"whales_boot_", class, "_", formatC(boot_no,width = 3, flag = "0"),"_mon", formatC(mon.num, width = 2, flag = "0"), ".csv"), row.names = FALSE)
        
        #print(boot)
        
        
        
        
      } #boot
      print(Sys.time() - start)
    } #month/track file
  }
}

stopCluster(cl)


## Small/Medium size class - account for correction factor
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

years <- c("2017", "2018", "2019", "2021", "2022")  

set.seed(1234)

for (k in 1:length(years)) {
  
  year <- years[k]
  
  track.files <- dir(path = "./Data/ais_rds/", pattern = glob2rx(paste0("narw_tracks_", year, "*.rds")))
  n.track.files <- length(track.files)
  
  class <- "M"
  
  #stuff for the bootstrap
  n.whales <- 350
  n.boot <- 1000
  months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  
  ###This loop cycles through AIS data for months of one year [k]
  
  output.folder <- paste0(o.folder, "/Class", class)
  if (!dir.exists(output.folder)) {dir.create(output.folder)}
  output.bootstrap.folder <- paste0(output.folder, "/bootstrap_model_output_", year, "/")
  if (!dir.exists(output.bootstrap.folder)) {dir.create(output.bootstrap.folder)}
  
  ##Parallel by month
  library(parallel) 
  num_cores <- 6 #CHECK COMPUTER CORES - recommended to do at least 1 or 2 fewer cores than your computer has, depending on what you have you can drop this down
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  start_month <- 1
  end_month <- 12
  
  start <- Sys.time()
  #(parallel)
  foreach(mon.num=seq(start_month, end_month), .combine="c", .packages = c("dplyr", "tidyr")) %dopar% {
    #for (mon.num in 1:12) {
    
    track.file <- track.files[mon.num] 
    tracks <- readRDS(paste("./Data/ais_rds", track.file, sep = "/"))
    
    
    yr <- substr(track.file, 13, 16)
    #mon.num <- mon.num %% 12 
    #if (mon.num == 0) {mon.num <- 12}
    
    #merge in the vessels
    tracks.vessel <- merge(tracks, vessel.data, by = "vessels_id", all.x = TRUE)
    
    #filter missing LOA and <7.925m reported length 
    tracks.vessel <- tracks.vessel %>% dplyr::filter(!is.na(report_loa)) %>% dplyr::filter(report_loa >= 7.925)
    
    tracks.vessel.final <- tracks.vessel
    
    #fix vessels with missing report_vessel_type
    tracks.vessel.final$report_vessel_type[tracks.vessel.final$report_vessel_type == ""] <- "Other"
    
    #merge vessel categories to standardize name and group rare vessel types
    suppressMessages(
      tracks.vessel.final <- left_join(tracks.vessel.final, vessel.cats)
    )
    
    #Classify size to match Garrison et al. 2024 categories:
    #S: < 40 feet (12.192 m) - contains vessels from 26-40 feet in length
    #M: >= 40 ft (12.192 m) - <65 ft (19.812 m)
    #L: >= 65 ft (19.812 m) - <350 ft (106.68 m)
    #XL >= 350 ft (106.68 m)
    
    tracks.vessel.final<-mutate(tracks.vessel.final, reclass_Vessel_Size = case_when(
      report_loa >= 106.68 ~ "XL",
      report_loa >= 19.812 &  report_loa < 106.68 ~ "L",
      report_loa >= 12.192 &  report_loa < 19.812 ~ "M",
      report_loa < 12.192 ~ "S"
    )) 
    
    
    ###select small and medium vessels
    tracks.vessel.final <- tracks.vessel.final[tracks.vessel.final$reclass_Vessel_Size %in% c("S", "M"),]
    
    tracks.vessel.final$loa_ft <- tracks.vessel.final$report_loa*3.28084
    summary(tracks.vessel.final$loa_ft)
    
    #note that about 1/3 of the tracks are in vessels <40m in length
    tracks.vessel.final<-mutate(tracks.vessel.final, sim_slow_class = case_when(
      loa_ft >= 40 ~ "40_65",
      loa_ft < 40 ~ "<40"
    ))
    
    ## set vessel mean draft and propsuctdepth for avoidance model
    mean.draft <- 5
    PropSuctDepth <- 1*mean.draft
    
    ###LINK TO PERMIT:AIS vessel counts by econ_vessel_type and loa_ft
    tracks.vessel.final<-mutate(tracks.vessel.final, size_5ft = case_when(
      loa_ft >= 26 & loa_ft < 30 ~ "26-30",       
      loa_ft >= 30 & loa_ft < 35 ~ "30-35",        
      loa_ft >= 35 & loa_ft < 40 ~ "35-40",
      loa_ft >= 40 & loa_ft < 45 ~ "40-45",
      loa_ft >= 45 & loa_ft < 50 ~ "45-50",
      loa_ft >= 50 & loa_ft < 55 ~ "50-55",
      loa_ft >= 55 & loa_ft < 60 ~ "55-60",
      loa_ft >= 60 & loa_ft < 65 ~ "60-65"
    ))

    
    #need to update to include smaller classes 
    vess.corrections.init <- read.csv("./Data/permit_vesscount_corrections_14Aug2023.csv", header = TRUE)
    
    suppressMessages(
      tracks.vessel.final <- left_join(tracks.vessel.final, vess.corrections.init[,c(1,2,5)])
    )
    
    #### COLLAPSE ACROSS VESS CORRECTIONS COLLAPSE TPYE AND CALCULATE CORRECTION RATIOS
    suppressMessages(
      vess.corrections.collapse <- vess.corrections.init %>% group_by(collapse.type, size_5ft) %>% 
        summarize(permit.count = sum(permit.vess.count), ais.count = sum(ais.vess.count)) %>%
        mutate(permit.to.ais = permit.count/ais.count)
    )
    
    vess.corrections.collapse$permit.to.ais[vess.corrections.collapse$collapse.type == "Pilot"] <- 1
    vess.corrections.collapse$permit.to.ais[vess.corrections.collapse$collapse.type == "Exempt"] <- 1
    
    ###calculate EE params
    tracks.vessel.final$whale.len <-  rnorm(nrow(tracks.vessel.final), mean = whale.len, sd = whale.len.sd)
    tracks.vessel.final$v.whale <- rweibull(nrow(tracks.vessel.final), shape=shape.whale, scale =scale.whale)
    tracks.vessel.final$whale.width <- tracks.vessel.final$whale.len * 3.5/14.5 #based on original width/length ratio
    tracks.vessel.final$radius.whale <- sqrt((tracks.vessel.final$whale.len * tracks.vessel.final$whale.width)/pi)
    
    ### Get whale density grid
    grid.density.mon <- data.frame(grid.density[, names(grid.density) %in% c("Id",months[mon.num], "region", "DiveArea","MidAtlArea", "depth", "depth15", "minZ15")])
    grid.density.mon[is.na(grid.density.mon)] <- 0
    grid.density.mon$gid <- as.numeric(as.character(grid.density.mon$Id))
    grid.density.mon$n.whale <- grid.density.mon[, names(grid.density.mon) == months[mon.num]]
    grid.region <- grid.density.mon[,names(grid.density.mon) %in% c("gid", "region", "n.whale", "DiveArea", "MidAtlArea", "depth", "depth15", "minZ15")]
    
    suppressMessages(
      tracks.vessel.final <- left_join(tracks.vessel.final, grid.region)
    )
    ##drop grid cells with depth of > 3m
    tracks.vessel.final <- tracks.vessel.final[tracks.vessel.final$depth15 <= -3,]
    
    #drop inshore cells see file ./data/excluded_inshore_cells.csv
    tracks.vessel.final <- tracks.vessel.final[!(tracks.vessel.final$gid == 1886034 | tracks.vessel.final$gid == 1911259 | tracks.vessel.final$gid == 1913061 | tracks.vessel.final$gid == 1913062 | tracks.vessel.final$gid == 1914863 | tracks.vessel.final$gid == 1916664 | tracks.vessel.final$gid == 1923873),]
    tracks.vessel.final$day <- format(as.Date(tracks.vessel.final$s.date, format="%Y-%m-%d"), "%d")
    
    ## speed simulation
    ## first save original speed (status quo)
    tracks.vessel.final$w.mu.speed.o <- tracks.vessel.final$w.mu.speed
    
    ### Select vessels above 10 knots
    fast.vessels <- tracks.vessel.final[tracks.vessel.final$w.mu.speed >10 & tracks.vessel.final$econ_vessel_type != "Exempt",]
    n.fast.tracks <- dim(fast.vessels) [1]
    
    ## remove these from the data - because we will put them back in later
    tracks.vessel.final.sim <- tracks.vessel.final[!(tracks.vessel.final$w.mu.speed >10 & tracks.vessel.final$econ_vessel_type != "Exempt"),] 
    
    
    ##add the speeds you will simulate back into the data so they are in the original data - unchanged for those remaining in the dataset
    ##because they are not in above 10 knots
    tracks.vessel.final.sim$v.boat.ms.o <- tracks.vessel.final.sim$w.mu.speed.o/ms2knot
    tracks.vessel.final.sim$w.mu.speed.sim <- tracks.vessel.final.sim$w.mu.speed.o
    tracks.vessel.final.sim$v.boat.ms.sim <- tracks.vessel.final.sim$w.mu.speed.o/ms2knot
    
    ### Select a random new speed for all tracks in fast.vessels
    fast.vessels$w.mu.speed.sim <- runif(n.fast.tracks, 9.5, 10.0)
    fast.vessels$v.boat.ms.sim <- fast.vessels$w.mu.speed.sim/ms2knot
    fast.vessels$v.boat.ms.o <- fast.vessels$w.mu.speed.o/ms2knot
    
    #### done with speeds - bind fast.vessels back into tracks
    tracks.vessel.final.sim <- rbind(tracks.vessel.final.sim, fast.vessels)
    
    #### calculate strike prob stuff for original speed
    ### original vessel speed
    newdata <- tracks.vessel.final.sim[, c("w.mu.speed.o", "reclass_Vessel_Size")]
    colnames(newdata) <- c("Vess.Speed", "vess.cat_f")
    newdata$spe.HB <- "Not Humpback"
    tracks.vessel.final.sim$p.mort.o <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
    
    #calculate encounter rate for each track given the presence of a random whale - for original vessel speed
    tracks.vessel.final.sim <- tracks.vessel.final.sim %>% dplyr::rowwise() %>% 
      dplyr::mutate(ER.fixed.o = getLambda(vb = v.boat.ms.o, vm = v.whale, S = area, r = radius.whale))
    
    tracks.vessel.final.sim <- tracks.vessel.final.sim %>% dplyr::mutate(EE.fixed.o = ER.fixed.o * tot.time.s) %>% 
      dplyr::mutate(EE.Mi.o = EE.fixed.o * p.mort.o)
    
    ### get p.mort for simulated speeds
    newdata <- tracks.vessel.final.sim[, c("v.boat.ms.sim", "reclass_Vessel_Size")]
    colnames(newdata) <- c("Vess.Speed", "vess.cat_f")
    newdata$Vess.Speed <- newdata$Vess.Speed * ms2knot
    newdata$spe.HB <- "Not Humpback"
    tracks.vessel.final.sim$p.mort.sim <- predict(LethalityCurveGLM, newdata = newdata, se.fit = F, type = "response")
    
    #calculate encounter rate for each track given the presence of a random whale - for simulated vessel speed
    tracks.vessel.final.sim <- tracks.vessel.final.sim %>% dplyr::rowwise() %>% 
      dplyr::mutate(ER.fixed.sim = getLambda(vb = v.boat.ms.sim, vm = v.whale, S = area, r = radius.whale))
    
    tracks.vessel.final.sim <- tracks.vessel.final.sim %>% dplyr::mutate(EE.fixed.sim = ER.fixed.sim * tot.time.s) %>% 
      dplyr::mutate(EE.Mi.sim = EE.fixed.sim * p.mort.sim)
    
    
    #### BOOTSTRAP FOR WHALE SIMULATION FOLLOWS    
    ### START BOOTSTRAP LOOP
    
    for (boot in 1:n.boot) {
      
      ### 1) Make a data.frame to hold individual whales
      whales <- data.frame(whale = c(1:n.whales))
      
      ### 2) Put the whales in the cells - could add a time step here to move the whales on a random walk
      
      ran.cells <- grid.density.mon[sample(nrow(grid.density.mon), n.whales, prob = grid.density.mon[,2], replace = TRUE), ]
      whales$gid <- ran.cells$gid
      whales$DiveArea <- ran.cells$DiveArea
      whales$tot.tracks <- NA
      whales$tot.tracks.gt10 <- NA
      ### 3) Calculate risk for each whale - individuals stays in the same cell for a whole month
      
      #status quo - uses o speed for all vessels
      whales$tot.EE.M.slow.none <- NA
      whales$avoid.tracks.slow.none <- NA #number of positive avoidances
      whales$mort.n.slow.none <- NA #number of positive moralities
      whales$mort.avoid.slow.none <- NA #number of morts after avoided
      
      #slow all scenario - uses sim speed for all vessels
      whales$tot.EE.M.slow.all <- NA
      whales$avoid.tracks.slow.all <- NA #number of positive avoidances
      whales$mort.n.slow.all <- NA #number of positive moralities
      whales$mort.avoid.slow.all <- NA #number of morts after avoided
      
      
      for(w in 1:n.whales) {
        
        #get the cell for whale w
        #get the tracks - currently this is all tracks for a whole month...
        tracks.cell <- data.frame(tracks.vessel.final.sim[tracks.vessel.final.sim$gid == whales$gid[w],])
        n.tracks <- dim(tracks.cell) [1]
        whales$tot.tracks[w] <- n.tracks
        whales$tot.tracks.gt10[w] <- length(tracks.cell$w.mu.speed.o[tracks.cell$w.mu.speed.o > 10])

        if (n.tracks> 0) {
          
          ####expand tracks.cell to account for missing AIS data - replicate the vessels (by size and collapse.type)
          ####sampled from tracks.vessel.final.sim - to represent all possible vessels
          
          ##1) Merge in the ratio - should automagically do by collapse-type and 
          tracks.cell <- left_join(tracks.cell, vess.corrections.collapse, by = c("collapse.type", "size_5ft"))
          
          cell.vess.types <- tracks.cell %>% distinct(collapse.type, size_5ft, permit.to.ais) %>% mutate(permit.to.ais.int = round(permit.to.ais))
          
          tracks.cell$p.surf <- mapply(prob.surface, tracks.cell$DiveArea, mean.draft)  
          
          
          tracks.cell$EE.M.o <-  tracks.cell$EE.Mi.o * tracks.cell$p.surf
          tracks.cell$EE.M.sim <- tracks.cell$EE.Mi.sim * tracks.cell$p.surf
          
          hityes <- mat.or.vec(n.tracks, 2) #status quo, slow all
          avoidyes <- mat.or.vec(n.tracks, 2) #status quo, slow all
          hit.avoid <- mat.or.vec(n.tracks, 2) #status quo, slow all
          
          tracks.cell$EE.M.slow.none   <- NA
          tracks.cell$EE.M.slow.all <- NA
          
          
          for (track.cell.i in 1:n.tracks) {
            
            permit.ratio <- cell.vess.types$permit.to.ais.int[cell.vess.types$collapse.type ==  tracks.cell$collapse.type[track.cell.i] & cell.vess.types$size_5ft == tracks.cell$size_5ft[track.cell.i] ]
            
            
            #speed <= 10 all speeds same (o) for all scenarios, one draw
            if (tracks.cell$w.mu.speed.o[track.cell.i] <= 10) {
              
              tracks.cell$EE.M.slow.none[track.cell.i]   <- tracks.cell$EE.M.o[track.cell.i]
              tracks.cell$EE.M.slow.all[track.cell.i] <- tracks.cell$EE.M.o[track.cell.i]
              
              ###scale to the permit ratio by taking permit.ratio draws and summing
              hityes.scale <- rbinom(size = 1, n = permit.ratio, prob = tracks.cell$EE.M.o[track.cell.i])
              
              avoidyes.scale <- rep(NA, permit.ratio)
              for (avoidyes.i in 1:permit.ratio) {
                avoidyes.scale[avoidyes.i] <- mapply(p.avoid, tracks.cell$v.boat.ms.o[track.cell.i], PropSuctDepth, tracks.cell$depth15[track.cell.i], class=class)
              }
              
              hit.avoid.scale <- hityes.scale * (1-avoidyes.scale)
              
              hityes[track.cell.i,] <- sum(hityes.scale)
              avoidyes[track.cell.i,] <- sum(avoidyes.scale)/permit.ratio 
              hit.avoid[track.cell.i,] <- sum(hit.avoid.scale) #This is the total number of actual hits
              #### END SCALE TO PERMIT
            }
            
            ##speed >= 10 all speeds same (o), different for slow everything scenario
            if (tracks.cell$w.mu.speed.o[track.cell.i] > 10 ) {
              
              tracks.cell$EE.M.slow.none[track.cell.i]   <- tracks.cell$EE.M.o[track.cell.i]
              tracks.cell$EE.M.slow.all[track.cell.i] <- tracks.cell$EE.M.sim[track.cell.i]
              
              ###scale to the permit ratio by taking permit.ratio draws and summing
              hityes.scale <- rbinom(size = 1, n = permit.ratio, prob = tracks.cell$EE.M.o[track.cell.i])
              
              avoidyes.scale <- rep(NA, permit.ratio)
              for (avoidyes.i in 1:permit.ratio) {
                avoidyes.scale[avoidyes.i] <- mapply(p.avoid, tracks.cell$v.boat.ms.o[track.cell.i], PropSuctDepth, tracks.cell$depth15[track.cell.i], class=class)
              }
              
              hit.avoid.scale <- hityes.scale * (1-avoidyes.scale)
              
              hityes[track.cell.i,c(1, 2)] <- sum(hityes.scale)
              avoidyes[track.cell.i,c(1, 2)] <- sum(avoidyes.scale)/permit.ratio 
              hit.avoid[track.cell.i,c(1, 2)] <- sum(hit.avoid.scale) #This is the total number of actual hits
              
              #slow all scenario uses sim speed
              hityes.scale <- rbinom(size = 1, n = permit.ratio, prob = tracks.cell$EE.M.sim[track.cell.i])
              
              avoidyes.scale <- rep(NA, permit.ratio)
              for (avoidyes.i in 1:permit.ratio) {
                avoidyes.scale[avoidyes.i] <- mapply(p.avoid, tracks.cell$v.boat.ms.sim[track.cell.i], PropSuctDepth, tracks.cell$depth15[track.cell.i], class=class)
              }
              
              hit.avoid.scale <- hityes.scale * (1-avoidyes.scale)
              
            } #vessel type if
          } #track.cell.i
          
          whales$tot.EE.M.slow.none[w] <- sum(tracks.cell$EE.M.slow.none)
          whales$avoid.tracks.slow.none[w] <- sum(avoidyes[,1]) #number of tracks avoided
          whales$mort.n.slow.none[w] <- sum(hityes[,1]) #number of positive mortalities
          whales$mort.avoid.slow.none[w] <- sum(hit.avoid[,1])
          
          whales$tot.EE.M.slow.all[w] <- sum(tracks.cell$EE.M.slow.all)
          whales$avoid.tracks.slow.all[w] <- sum(avoidyes[,2]) #number of tracks avoided
          whales$mort.n.slow.all[w] <- sum(hityes[,2]) #number of positive mortalities
          whales$mort.avoid.slow.all[w] <- sum(hit.avoid[,2])
          
        } #tracks > 0 if
        
      } #whale
      
      ##NAs turn to zeros
      whales[is.na(whales)] <- 0
      
      boot_no <- boot
      write.csv(whales, paste0(output.bootstrap.folder,"whales_boot_", class, "_", formatC(boot_no,width = 3, flag = "0"),"_mon", formatC(mon.num, width = 2, flag = "0"), ".csv"), row.names = FALSE)
      
    } #boot
    
  } #end foreach month
  
  print(paste0("Year:", year, " Runtime: ", Sys.time() - start))
  stopCluster(cl)
  
} #years        

