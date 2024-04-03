
library(tidyverse)
library(ncdf4)
library(raster)
library(lhs)

setwd('your_wd')
load('WHAM_INFERNO_EO_data.RData')

.resultpath <- 'your_results_dir'

MOD.crop    <- crop(MOD.JULES, extent(c(-180, 180, -60, 90)))

###########################################################################################

### Latin hypercube of parameters

###########################################################################################

#############################################

### Hypercube sampling

#############################################

for(seed_val in 1987:2006) {
  
  set.seed(seed_val)

A     <- maximinLHS(1000, 20)
B     <- matrix(nrow = nrow(A), ncol = ncol(A))

### Burned area per PFT
B[,1] <- qunif(A[,1], 0.85, 2.55) #tree
B[,2] <- qunif(A[,2], 0.85, 2.55) #tree
B[,3] <- qunif(A[,3], 1.6, 4.8)  #grass
B[,4] <- qunif(A[,4], 1.6, 4.8)  #pasture
B[,5] <- qunif(A[,5], 1.6, 4.8)  #grass
B[,6] <- qunif(A[,6], 1.6, 4.8)  #pasture
B[,7] <- qunif(A[,7], 1.35, 4.05) #shrub

### Ignition-fire scaling
B[, 8]<- qunif(A[, 8], 0.01, 0.05) # ig scaling
B[, 9]<- qunif(A[, 9], 10, 50) # arson scaling
B[, 10]<-qunif(A[, 10], 3.85, 11.55) # lightning scaling
B[, 11]<-qunif(A[, 11], 400, 900) # ontology scaling

### Suppression parameters
B[,12] <- qunif(A[,12], min = 0, max = 0.05)
B[,13] <- qunif(A[,13], min = 0, max = 0.1)
B[,14] <- qunif(A[,14], min = 0.8, max = 1)


### Fragmentation parameters
B[,15] <- qunif(A[,15], 4.55, 13.4) # Road
B[,16] <- qunif(A[,16], 1, 2.5) # Logging
B[,17] <- qunif(A[,17], min = 0.1, max = 0.4) #previous fire
B[,18] <- qunif(A[,18], min = 0.1, max = 0.4) #previous fire

B[,19] <- qunif(A[,19], min = 0.5, max = 1.5) #Scaling factor for managed pasture
B[,20] <- qunif(A[,20], min = 0.5, max = 1.5) #Scaling factor for managed vegetation

B     <- B[B[, 1] < B[, 3], ]
B     <- B[B[, 2] < B[, 3], ]
B     <- B[B[, 1] < B[, 5], ]
B     <- B[B[, 2] < B[, 5], ]
B     <- B[B[, 12] < B[, 13], ]
B     <- B[B[, 13] < B[, 14], ]

B     <- B[apply(B, 1, function(x) {!any(x < 0)}), ]

PFT_pars <- data.frame(B)
colnames(PFT_pars) <- c('Tree1', 'Tree2', 'Grass1', 'Pasture1', 'Grass2','Pasture2', 'Shrub', 
                        'igs', 'arson', 'lightning', 'scaling',
                        'Pre', 'Trans', 'Intense', 
                        'Road', 'Logging', 'fire_thresh', 'fire_slope', 
                        'Managed_pasture', 'Managed_veg')

###########################################################################################

### Calculate WHAM unmanaged fire

###########################################################################################

fires        <- list()
igs          <- list()
WHAM_ba      <- list()
lightning_ba <- list()
Flam_total   <- list()


calibration_res <- vector(mode = 'list', length = length(nrow(PFT_pars)))



for(z in 1:nrow(PFT_pars)) {
  
  
  
  ###########################################################################################
  
  ### Update suppression
  
  ###########################################################################################
  
  ######################################################
  ### Suppression model
  ######################################################
  
  Sup.frame$lin_extinct <- ifelse(Sup.frame$Sup_index == 1, PFT_pars$Pre[z], 
                                  ifelse(Sup.frame$Sup_index == 2, PFT_pars$Trans[z], 
                                         ifelse(Sup.frame$Sup_index == 3, PFT_pars$Intense[z], 0)))
  
  E.lin <- lm(lin_extinct ~ ., 
              data = Sup.frame[, c(5, 7:8, 10)])
  
  
  ######################################################
  ### Compile maps
  ######################################################
  
  Suppression <- list()
  
  for(year in 1:nlayers(Arson)) {
    
    Map.dat           <- data.frame(Pre[[year]][], Trans[[year]][], Intense[[year]][], Post[[year]][])
    colnames(Map.dat) <- c('Pre', 'Trans', 'Intense', 'Post')
    
    lin.pred              <- JULES.mask
    values(lin.pred)      <- predict(E.lin, Map.dat)
    lin.pred[lin.pred >1] <- 1
    lin.pred[lin.pred <0] <- 0
    
    Suppression[[year]] <- lin.pred
    
  }
  
  ###########################################################################################
  
  ### Update BA per pft
  
  ###########################################################################################
  
  PFT.key <- list('broadleaf deciduous'          = c('tree', PFT_pars$Tree1[z]), 
                  'broadleaf evergreen tropical' = c('tree', PFT_pars$Tree1[z]), 
                  'broadleaf evergreen temperate'= c('tree', PFT_pars$Tree1[z]), 
                  'needleleaf deciduous'         = c('tree', PFT_pars$Tree2[z]), 
                  'needleleaf evergreen'         = c('tree', PFT_pars$Tree2[z]), 
                  'C3 grass'                     = c('grass', PFT_pars$Grass1[z]), 
                  'C3 crop'                      = c('grass', 0), 
                  'C3 pasture'                   = c('grass', PFT_pars$Pasture1[z]), 
                  'C4 grass'                     = c('grass', PFT_pars$Grass2[z]), 
                  'C4 crop'                      = c('grass', 0),
                  'C4 pasture'                   = c('grass', PFT_pars$Pasture2[z]), 
                  'deciduous shrub'              = c('shrub', PFT_pars$Shrub[z]), 
                  'evergreen shrub'              = c('shrub', PFT_pars$Shrub[z]))
  
  PFT.ba <- as.numeric(unlist(lapply(PFT.key, function(x) x[2])))
  PFT.ba <- brick(lapply(PFT.ba, function(x) x/area(JULES.mask))) #ba frac per fire per pft
  
  
  ### adjust for roads
  thresh<- PFT_pars$Road[z]
  Roads            <- 1 - brick(lapply(1:nlayers(raw_Roads), function(x) log(raw_Roads[[x]]) /thresh))
  Roads[Roads < 0] <- 0
  
  for(i in 1:nlayers(Arson)) {
    
    
    ### Annual ignitions
    fires[[i]]   <- list()
    igs[[i]]     <- list()
    WHAM_ba[[i]] <- list()
    lightning_ba[[i]] <- list()
    Flam_total[[i]]   <- list()
    
    ### Adjust BA per pft for roads
    PFT.roads <- (PFT.ba) * (Roads[[i]])
    
    for(j in c(1:12)) {
      
      ### Monthly anthropogenic fires   
      igs[[i]][[j]]     <- ((PFT_pars$igs[z]/12) * (1 - Suppression[[i]])) + (Escaped[[i]]/12) + Arson[[i]]/PFT_pars$arson[z] # background rate + escaped + arson
      igs[[i]][[j]]     <- igs[[i]][[j]] * area(JULES.mask) # ig per cell
      
      #Monthly Flam per cell
      Flam.tot               <- brick(lapply(PFT[[i]][1:13], function(x){
                                 x[[j]]})) * brick(lapply(Flam[[i]], function(x) x[[j]]))
      
      Flam.tot               <- sum(Flam.tot)
      Flam.tot[Flam.tot > 1] <- 1
      Flam_total[[i]][[j]]   <- Flam.tot
      
      
      ### Adjust BA per pft for Logging
      logged_for <- Logging[[i]] / sum(brick(lapply(PFT[[i]][1:3], function(x) {
        x[[j]]
      }))) 
      
      logged_for[logged_for > 1] <- 1
      PFT.for                    <- PFT.roads
      
      for(k in 1:3) {
        PFT.for[[k]] <- PFT.for[[k]] + (PFT_pars$Logging[z]*logged_for) / area(JULES.mask)
      }
      
      ### BA per pft per cell - adjusted for roads
      BA_PFT.tot             <- sum(brick(lapply(PFT[[i]][1:13], 
                                                 function(x) x[[j]])) * PFT.for)

      ### Monthly BA
      ### constant rate of unmanaged fire
      fires[[i]][[j]]      <-  (igs[[i]][[j]]) * Flam.tot ### fires per cell
      WHAM_ba[[i]][[j]]    <-  BA_PFT.tot * (fires[[i]][[j]]) ### ba frac per cell
      
      ### scaling factor -> ABM is giving fires not ignitions!
      WHAM_ba[[i]][[j]]    <- WHAM_ba[[i]][[j]] * PFT_pars$scaling[z] ### highly sensitive parameter
      
      ### lightning
      strikes                <- (30 * lightning[[j]] * area(JULES.mask)) * (1- Suppression[[i]]) #strikes per box
      lightning_ba[[i]][[j]] <- PFT_pars$lightning[z] * strikes * Flam.tot * BA_PFT.tot 
      
      
      gc()
      
      
      
    }
    
    print(i)
    
  }
  
  
  #################################################################################
  
  ### Add managed and lightning
  
  #################################################################################
  
  tot_WHAM <- list()
  comb_BA  <- list()
  
  
  for(i in 1:length(WHAM_ba)) {
    
    tot_WHAM[[i]] <- list()
    comb_BA[[i]]  <- list()
    
    Managed.adj <- (Arable[[i]] + (Vegetation[[i]] * PFT_pars$Managed_veg[z]) + 
                          (Pasture[[i]] * PFT_pars$Managed_pasture[z]))
    
    for(j in c(1:12)) {
      
      tot_WHAM[[i]][[j]]  <- WHAM_ba[[i]][[j]] + Managed.adj/12
      comb_BA[[i]][[j]]   <- lightning_ba[[i]][[j]] + tot_WHAM[[i]][[j]]
      
      
    }
    
    print(i)
    
  }
  
  
  ### previous fires
  comb.annual <- brick(lapply(comb_BA, function(x) sum(brick(x))))
  UM.annual   <- brick(lapply(WHAM_ba, function(x) sum(brick(x)))) + brick(lapply(lightning_ba, function(x) sum(brick(x))))
  M.annual    <- comb.annual - UM.annual
  
  k2 <- PFT_pars$fire_thresh[z]
  s  <- PFT_pars$fire_slope[z]
  
  ca <- brick(lapply(1:25, function(x) comb.annual[[x]] - Arable[[x]]))
  
  for(i in 1:nlayers(UM.annual)) {
    
    UM.annual[[i]][ca[[i]] > k2] <- UM.annual[[i]][ca[[i]] > k2] - ((UM.annual[[i]][ca[[i]] > k2] - k2) * (1-s))
    comb.annual[[i]]             <- M.annual[[i]] + UM.annual[[i]]
    
  }
  
  comb.annual[is.na(comb.annual[[1]]) & !is.na(MOD.JULES)] <- 0
  comb.annual                                              <- crop(comb.annual, extent(c(-180, 180, -60, 90)))
  
  calibration_res[[z]]              <- list()
  
  calibration_res[[z]]$BA_ssa       <- sum((area(Sentinel.JULES) * comb.annual[[25]])[filt], na.rm = T)
  
  calibration_res[[z]]$BA_global    <- summary(sapply(1:14, 
                                                      function(i) sum((area(JULES.mask) * comb.annual[[i+11]])[], na.rm = T)))
  
  calibration_res[[z]]$cor_ssa      <- cor.test(sqrt(Sentinel.JULES)[filt], sqrt(comb.annual[[25]])[filt])
  
  calibration_res[[z]]$cor_global   <- summary(unlist(lapply(1:14, function(i) {cor.test(sqrt(MOD.crop[[i]])[], 
                                                                                         sqrt(comb.annual[[i+11]])[])$estimate})))

  
  print(calibration_res[[z]])
  print(PFT_pars[z, ])
  
  if(z %% 100 == 0 | z == nrow(B)) {
    
    save.image(paste0(.resultpath, 'Seed_', seed_val, '_GFED5_Thesis.RData'))
    
  }
  
  print('done a run')
  
  }
  
}

