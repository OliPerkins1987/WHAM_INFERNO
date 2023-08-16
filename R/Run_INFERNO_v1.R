
library(tidyverse)
library(ncdf4)
library(raster)
library(lhs)

load("C:/Users/Oli/Documents/PhD/wham_coupled/RCode/Tidy/INFERNO_Data.RData")

.resultpath <- 'D:/WHAM_coupled/calibration_results/Annual/Seed_'

for(seed_val in 1987:2004) {

###########################################################################################

### Latin hypercube of parameters

###########################################################################################

set.seed(seed_val)

A     <- maximinLHS(1000, 15)
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
B[, 8]<- qunif(A[, 8], 0.5, 2)  # ig scaling
B[, 9]<- qunif(A[, 9], 3.85, 11.55) # lightning scaling

### Suppression parameters
B[,10] <- qunif(A[,10], min = 0.5, max = 2) #suppression

### Fragmentation parameters
B[,11] <- qunif(A[,11], 4.55, 13.4) # Road
B[,12] <- qunif(A[,12], 1.6, 4.8)   # Crop
B[,13] <- qunif(A[,13], 1.6, 4.8)   # Crop
B[,14] <- qunif(A[,14], min = 0.1, max = 0.4) #previous fire
B[,15] <- qunif(A[,15], min = 0.1, max = 0.4) #previous fire

B     <- B[apply(B, 1, function(x) {!any(x < 0)}), ]

PFT_pars           <- data.frame(B)
colnames(PFT_pars) <- c('Tree1', 'Tree2', 'Grass1', 'Pasture1', 'Grass2','Pasture2', 'Shrub', 
                        'igs', 'lightning', 'suppression', 'Road', 'Crop1', 'Crop2', 'fire_thresh', 'fire_slope')


###########################################################################################

### Calculate WHAM unmanaged fire

###########################################################################################

fires        <- list()
igs          <- list()
ba           <- list()
lightning_ba <- list()
Flam_total   <- list()


calibration_res <- vector(mode = 'list', length = length(nrow(PFT_pars)))


for(z in 1:nrow(PFT_pars)) {
  
  
  ###########################################################################################
  
  ### Update BA per pft
  
  ###########################################################################################
  
  PFT.key <- list('broadleaf deciduous'          = c('tree', PFT_pars$Tree1[z]), 
                  'broadleaf evergreen tropical' = c('tree', PFT_pars$Tree1[z]), 
                  'broadleaf evergreen temperate'= c('tree', PFT_pars$Tree1[z]), 
                  'needleleaf deciduous'         = c('tree', PFT_pars$Tree2[z]), 
                  'needleleaf evergreen'         = c('tree', PFT_pars$Tree2[z]), 
                  'C3 grass'                     = c('grass', PFT_pars$Grass1[z]), 
                  'C3 crop'                      = c('grass', PFT_pars$Crop1[z]), 
                  'C3 pasture'                   = c('grass', PFT_pars$Pasture1[z]), 
                  'C4 grass'                     = c('grass', PFT_pars$Grass2[z]), 
                  'C4 crop'                      = c('grass', PFT_pars$Crop2[z]),
                  'C4 pasture'                   = c('grass', PFT_pars$Pasture2[z]), 
                  'deciduous shrub'              = c('shrub', PFT_pars$Shrub[z]), 
                  'evergreen shrub'              = c('shrub', PFT_pars$Shrub[z]))
  
  PFT.ba <- as.numeric(unlist(lapply(PFT.key, function(x) x[2])))
  PFT.ba <- brick(lapply(PFT.ba, function(x) x/area(JULES.mask))) #ba frac per fire per pft
  
  
  ### adjust for roads
  #thresh           <- PFT_pars$Road[z]
  #Roads            <- 1 - brick(lapply(1:nlayers(raw_Roads), function(x) log(raw_Roads[[x]]) /thresh))
  #Roads[Roads < 0] <- 0
  
  for(i in 1:25) {
    
    
    ### Annual ignitions
    fires[[i]]        <- list()
    igs[[i]]          <- list()
    ba[[i]]           <- list()
    lightning_ba[[i]] <- list()
    Flam_total[[i]]   <- list()
    
    ### Adjust BA per pft for roads
    #PFT.roads <- (PFT.ba) * (Roads[[i]])
    
    ### Suppression
    NS.adjust <- (PFT_pars$suppression[z] * frac.NS[[i]])

    
    for(j in c(1:12)) {
      
      
      ### Monthly anthropogenic fires   
      igs[[i]][[j]]   <- (PFT_pars$igs[z] * Pop.igs[[i]]/12) * NS.adjust  # ig per km-2
      igs[[i]][[j]]   <- igs[[i]][[j]] * area(JULES.mask)                 # ig per cell
      
      #Monthly Flam per cell
      Flam.tot               <- brick(lapply(PFT[[i]][1:13], function(x){
                                 x[[j]]})) * brick(lapply(Flam[[i]], function(x) x[[j]]))
      
      Flam.tot               <- sum(Flam.tot)
      Flam.tot[Flam.tot > 1] <- 1
      Flam_total[[i]][[j]]   <- Flam.tot
      
      
      ### BA per pft per cell - adjusted for roads
      BA_PFT.tot             <- sum(brick(lapply(PFT[[i]][1:13], 
                                                 function(x) x[[j]])) * PFT.ba) * JULES.mask

      ### Monthly BA
      ### constant rate of unmanaged fire
      fires[[i]][[j]] <-  (igs[[i]][[j]]) * Flam.tot     ### fires per cell
      ba[[i]][[j]]    <-  BA_PFT.tot * (fires[[i]][[j]]) ### ba frac per cell

      ### lightning
      strikes                <- (30 * lightning[[j]] * area(JULES.mask)) * (NS.adjust) #strikes per box
      lightning_ba[[i]][[j]] <- PFT_pars$lightning[z] * strikes * Flam.tot * BA_PFT.tot 
      
      
      gc()
      
      
      
    }
    
    print(i)
    
  }
  
  ### previous fires
  BA.annual <- brick(lapply(ba, function(x) {sum(brick(x))})) + brick(lapply(lightning_ba, function(x) {sum(brick(x))}))
  
  k2 <- PFT_pars$fire_thresh[z]
  s  <- PFT_pars$fire_slope[z]
  
  for(i in 1:nlayers(BA.annual)) {
    
    BA.annual[[i]][BA.annual[[i]] > k2] <- BA.annual[[i]][BA.annual[[i]] > k2] - ((BA.annual[[i]][BA.annual[[i]] > k2] - k2) * (1-s))
    
  }
  
  calibration_res[[z]]              <- list()
  
  calibration_res[[z]]$BA_ssa       <- sum((area(Sentinel.JULES) * BA.annual[[25]])[filt], na.rm = T)
  
  calibration_res[[z]]$BA_global    <- summary(sapply(1:14, 
                                        function(i) sum((area(JULES.mask) * BA.annual[[i+11]])[], na.rm = T)))
  
  calibration_res[[z]]$cor_ssa      <- cor.test(sqrt(Sentinel.JULES)[filt], sqrt(BA.annual[[25]])[filt])
  
  calibration_res[[z]]$cor_global   <- summary(unlist(lapply(1:14, function(i) {cor.test(sqrt(MOD.JULES[[i]])[], 
                                          sqrt(BA.annual[[i+11]])[])$estimate})))

  print(calibration_res[[z]])
  print(PFT_pars[z, ])
  
  if(z %% 100 == 0 | z == nrow(B)) {
    
    save.image(paste0(.resultpath, seed_val, '_GFED5_Thesis.RData'))
    
  }
  
  print('done a run')
  
  
  }
  

}


