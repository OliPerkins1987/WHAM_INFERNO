
library(rPref)
library(raster)
library(tidyverse)
library(ncdf4)
library(lhs)

###########################################

### Process GFED

###########################################

setwd('C:/Users/Oli/Documents/PhD/wham_coupled/results/thesis/Calibration/Spatial')
comb.pareto <- read.csv('Pareto_pars.csv')

setwd('C:/Users/Oli/Documents/PhD/wham_coupled/RCode/tidy')
load('WHAM_INFERNO_Data.RData')
Igs      <- 12*(Arson/34.44) + Escaped + 0.04025 ### annual unmanaged fires
Road.res <- lapply(lapply(1:nrow(comb.pareto), function(z) {
               1 - brick(lapply(1:nlayers(raw_Roads), 
                  function(x) log(raw_Roads[[x]]) /comb.pareto$Road[z]))
}), function(y) {y[y<0] <- 0
                y})

setwd('F:/PhD/Data/Unprocessed secondary/Fires/GFED5/Combined')
MOD.GFED5  <- brick('Total_GFED5_025deg.tif')
MOD.JULES  <- resample(MOD.GFED5/area(MOD.GFED5), JULES.mask)
MOD.JULES[JULES.mask <= 0] <- NA


##################################################################################################

### Process WHAM

##################################################################################################

setwd('C:/Users/Oli/Documents/PhD/wham_coupled/results/thesis/Calibration/Spatial/Annual')

Tot.GFED     <- lapply(list.files(pattern = 'Total.*.tif'), brick)
UM.GFED      <- lapply(list.files(pattern = 'Unmanaged.*.tif'), brick)
Managed.GFED <- lapply(list.files(pattern = 'Managed.*.tif'), brick)
UM_anthro    <- lapply(list.files(pattern = 'Anthro.*.tif'), brick)
Lightning    <- lapply(list.files(pattern = 'Light.*.tif'), brick)
Flam.tot     <- brick('Flam_tot.tif')
Sup.pareto   <- brick('Suppression_pareto.tif')

result <- list()
UM     <- list()
Man    <- list()
Road   <- list()

for(i in 1:25) {
  
  result[[i]]   <- list()
  UM[[i]]       <- list()
  Man[[i]]      <- list()
  #Road[[i]]     <- list()
  
  for(j in 1:length(Tot.GFED)) {
    
    result[[i]][[j]]  <- Tot.GFED[[j]][[i]]
    UM[[i]][[j]]      <- UM.GFED[[j]][[i]]
    Man[[i]][[j]]     <- Managed.GFED[[j]][[i]]
    #Road[[i]][[j]]    <- Road.res[[j]][[i]]
  }
  
}

result.rast <- brick(lapply(result, function(j) {mean(brick(j))}))
UM.rast     <- brick(lapply(UM, function(j) {mean(brick(j))}))
M.rast      <- brick(lapply(Man, function(j) {mean(brick(j))}))
Road.rast   <- brick(lapply(Road, function(j) {mean(brick(j))}))

#########################
### Continents 
########################

setwd('C:/Users/Oli/Documents/PhD/Model development/Data/Secondary data/World_Continents')
continents <- rgdal::readOGR('4a7d27e1-84a3-4d6a-b4c2-6b6919f3cf4b202034-1-2zg7ul.ht5ut.shp')

WHAM.co<- raster::extract(result.rast*area(result.rast), continents, fun = sum, na.rm = T)
GFED.co<- raster::extract(MOD.JULES*area(result.rast), continents, fun = sum, na.rm = T)
Man.co <- raster::extract(M.rast*area(result.rast), continents, fun = sum, na.rm = T)
UM.co  <- raster::extract(UM.rast*area(result.rast), continents, fun = sum, na.rm = T)


###########################################################################################

### 1) Basic plots & analysis

###########################################################################################

### basic map plot
map.df <- as.data.frame(as(brick(list(Total = result.rast[[c(1)]], 
                Managed = M.rast[[c(1)]], 
                Unmanaged = UM.rast[[c(1)]])), "SpatialPixelsDataFrame"))

map.df2<- as.data.frame(as(brick(list(Total = result.rast[[c(25)]], 
          Managed = M.rast[[c(25)]], Unmanaged = UM.rast[[c(25)]])), 
          "SpatialPixelsDataFrame"))

### Figure 7

rbind(mutate(map.df, Year = 1990), mutate(map.df2, Year= 2014)) %>%
  pivot_longer(1:3, values_to = 'Burned fraction') %>%
  ggplot(aes(x = x, y = y, fill = `Burned fraction`)) +
  theme_classic() +
  geom_raster() + 
  scale_fill_viridis_c(na.value = 'white', 
        trans = scales::pseudo_log_trans(sigma = 0.05), 
        labels = c(0, 0.1, 0.25, 0.5, 1), 
        breaks = c(0, 0.1, 0.25, 0.5, 1)) + 
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(Year~name)+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="bottom", 
        text = element_text(size = 14), 
        legend.text = element_text(angle = 65, hjust = 0.25, vjust = 0.5))

setwd('C:/Users/Oli/Documents/PhD/Thesis/Draft content/ABM Coupled/Charts')
#ggsave('WHAM_map.tiff', width = 12, height = 7)

### basic TS plot
ts.df <- data.frame(WHAM_INFERNO = sapply(12:25, function(z) {sum((result.rast[[z]]*area(JULES.mask))[], na.rm = T)}), 
                    GFED5 =   sapply(1:14, function(i) {sum((MOD.JULES[[i]] * area(MOD.JULES))[], na.rm = T)}),
                    UM    = sapply(12:25, function(z) {sum((UM.rast[[z]]*area(JULES.mask))[], na.rm = T)}),
                    Man   = sapply(12:25, function(z) {sum((M.rast[[z]]*area(JULES.mask))[], na.rm = T)}),
                    Year = 2001:2014)

### Figure 9
ts.df %>%
  pivot_longer(c(WHAM_INFERNO, GFED5), names_to = 'Method') %>%
   ggplot(aes(x = Year, y = value/10000, colour = Method)) + geom_line(size = 1) +
    scale_colour_viridis_d() + theme_classic() + 
     theme(text = element_text(size = 14)) + ylab('Burned area (Mha)')

#ggsave('WHAM_GFED_ts.tiff', width = 9, height = 6)

##############################
### Continent-scale plot
##############################

Man.co   <- t(Man.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
  setNames(c(continents$CONTINENT, 'Year'))

UM.co    <-  t(UM.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
  setNames(c(continents$CONTINENT, 'Year'))

cont.ba          <- rbind(Man.co, UM.co)
cont.ba$Type     <- c(rep('Managed', times = nrow(cont.ba)/2), 
                     rep('Unmanaged', times = nrow(cont.ba)/2))

### Figure 8
cont.ba %>% pivot_longer(continents$CONTINENT) %>%
  filter(!name %in% c('Antarctica', 'Oceania')) %>%
  ggplot(aes(x = Year, y = value/10000, colour = Type)) + geom_line(size = 1) +
  facet_grid(.~name) +
  ylab('Burned area (Mha)') + scale_colour_viridis_d() +
  theme_classic() + theme(text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

#ggsave('WHAM_ts.tiff', width = 13, height = 9)

#################################
## summary stats
#################################

### Global managed and unmanaged
summary(sapply(1:25, function(z) {sum((UM.rast[[z]]*area(JULES.mask))[], na.rm = T)}))
summary(sapply(1:25, function(z) {sum((M.rast[[z]]*area(JULES.mask))[], na.rm = T)}))

### cont breakdown

cont.ba %>%
  pivot_longer(c(1:4, 6, 8)) %>%
  group_by(name, Year) %>%
   mutate(Total = sum(value)) %>% ungroup() %>%
   group_by(name, Type) %>%
    summarise(.tot = mean(value), 
              .frac= mean(value) / mean(Total))

### Global trend
ts.df %>%
  summarise_all(function(z) {z[14] - z[1]})

### Continent temporal trend
cont.ba %>%
  group_by(Type) %>%
  summarise_if(is.numeric, function(z) {z[25] - z[12]})

#############################################################################################

### 2) Drivers of unmanaged fire by continent

##############################################################################################

Sup.co <- raster::extract(Sup.pareto*area(result.rast), continents, fun = sum, na.rm = T)
Flam.co<- raster::extract(Flam.tot*area(result.rast), continents, fun = sum, na.rm = T)
Igs.co <- raster::extract(Igs*area(result.rast), continents, fun = sum, na.rm = T)
Road.co <- raster::extract((1-Road.rast)*area(result.rast), continents, fun = sum, na.rm = T)

### continents
WHAM.co  <-  t(WHAM.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
  setNames(c(continents$CONTINENT, 'Year'))

GFED.co  <-  t(GFED.co) %>% data.frame() %>% mutate(Year = 2001:2020) %>%
  setNames(c(continents$CONTINENT, 'Year'))

Road.co  <-  t(Road.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
              setNames(c(continents$CONTINENT, 'Year'))

Sup.co  <-  t(Sup.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
  setNames(c(continents$CONTINENT, 'Year'))

Flam.co  <-  t(Flam.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
  setNames(c(continents$CONTINENT, 'Year'))

Igs.co   <- t(Igs.co) %>% data.frame() %>% mutate(Year = 1990:2014) %>%
  setNames(c(continents$CONTINENT, 'Year'))


###############################################

### plots

################################################

### Figure 10

cont.ba          <- rbind(WHAM.co, GFED.co, Road.co, Sup.co, Flam.co, Igs.co)
cont.ba$Method <- c(rep('WHAM-INFERNO', times = nrow(WHAM.co)), rep('GFED5', times = nrow(GFED.co)), 
                    rep('Road_density', times = nrow(WHAM.co)), rep('Suppression', times = nrow(WHAM.co)), 
                    rep('Flammability', times = nrow(WHAM.co)), rep('NFires', times = nrow(WHAM.co)))

cont.ba %>% pivot_longer(continents$CONTINENT) %>%
  filter(!name %in% c('Antarctica', 'Oceania')) %>%
  filter(name %in% c('Africa', 'Asia', 'South America') & Year %in% 2001:2014) %>%
  filter(Method %in% c('WHAM-INFERNO', 'GFED5')) %>%
  ggplot(aes(x = Year, y = value/10000, colour = Method)) + geom_line(size = 1) +
  facet_grid(.~name) + ylab('Burned area (Mha)') + scale_colour_viridis_d() +
  theme_classic() + theme(text = element_text(size = 14))

#ggsave('WHAM_GFED_by_conts.tiff', width = 9, height = 5)

cont.ba          <- rbind(Road.co, Sup.co, Flam.co, Igs.co)
cont.ba$Variable <- c(rep('Road density', times = nrow(cont.df)), rep('Suppression', times = nrow(cont.df)), 
                      rep('Flammability', times = nrow(cont.df)), rep('Nfires', times = nrow(cont.df)))

cont.ba %>% filter(Year %in% 2001:2014) %>%
  pivot_longer(continents$CONTINENT) %>%
  filter(name %in% c('Africa', 'Asia', 'South America')) %>%
  filter(!Variable %in% c('WHAM-INFERNO', 'GFED5')) %>%
  group_by(name, Variable) %>% mutate(value = value /value[1]) %>%
  ggplot(aes(x = Year, y = value, colour = Variable)) + geom_line(size = 1) +
  facet_grid(.~name) + ylab('Dimensionless scale (2001 = 1)') + scale_colour_viridis_d() +
  theme_classic() + theme(text = element_text(size = 14))

#ggsave('Cont_drivers.tiff', width = 9, height = 5)


#######################################

### drivers of fire by cont: correlations

#######################################

cont.ba        <- rbind(UM.co, GFED.co, Road.co, Sup.co, Flam.co, Igs.co)
cont.ba$Method <- c(rep('WHAM-INFERNO', times = nrow(cont.df)), rep('GFED5', times = nrow(GFED.co)), 
                    rep('Road_density', times = nrow(cont.df)), rep('Suppression', times = nrow(cont.df)), 
                    rep('Flammability', times = nrow(cont.df)), rep('Nfires', times=  nrow(cont.df)))

for(cont in c('Africa', 'Asia', 'South America')) {
  
  teg <- cont.ba[c(cont, 'Year', 'Method')]
  
  print(cor.test(teg[[cont]][teg$Method == 'WHAM-INFERNO'], 
        teg[[cont]][teg$Method == 'Nfires']))
  
}




#################################################################################################

## 3) Error analysis

#################################################################################################

I <- list()
for(i in 1:14) {
  
  I[[i]] <- list()
  
  for(j in 1:length(Tot.GFED)) {
    
    I[[i]][[j]] <- Tot.GFED[[j]][[i+11]] - MOD.JULES[[i]] 
    
  }
  
}

error.rast <- brick(lapply(I, function(z) {mean((brick(z)))}))

###################

### plot

##################

error.error<- as.data.frame(as(brick(list(error_2001 = error.rast[[1]], 
                                          error_2014 = error.rast[[14]])), "SpatialPixelsDataFrame"))

error.error %>%
  pivot_longer(c(error_2001, error_2014), values_to = 'Burned fraction') %>%
  ggplot(aes(x = x, y = y, fill = `Burned fraction`)) +
  geom_raster() + scale_fill_viridis_c(na.value = 'white', option = 'E') + 
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(.~name)+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        text = element_text(size = 12)) +
  theme(legend.position="bottom")

ggsave('error_maps.tiff', width = 9, height = 5)

### weighted error plots
error.UM <- (error.rast * ((UM.rast[[12:25]]) /result.rast[[12:25]]))
error.M  <- (error.rast * ((M.rast[[12:25]]) /result.rast[[12:25]]))

error.error <- as.data.frame(as(brick(list(error_Unmanaged = mean(error.UM), 
                                           error_Managed = mean(error.M))), "SpatialPixelsDataFrame"))

error.error %>%
  pivot_longer(c(error_Unmanaged, error_Managed), values_to = 'Burned fraction') %>%
  ggplot(aes(x = x, y = y, fill = `Burned fraction`)) +
  geom_raster() + scale_fill_viridis_c(na.value = 'white', option = 'E') + 
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(.~name)+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        text = element_text(size = 12), 
        legend.text = element_text(angle = 90, hjust = 0.8)) +
  theme(legend.position="bottom")

ggsave('weighted_error_maps.tiff', width = 9, height = 5)

### WHAM GFED plot
error.2014 <- as.data.frame(as(brick(list(WHAM_INFERNO = result.rast[[14]], GFED5 = MOD.JULES[[14]])), "SpatialPixelsDataFrame"))
error.2001 <- as.data.frame(as(brick(list(WHAM_INFERNO = result.rast[[1]], GFED5 = MOD.JULES[[1]])), "SpatialPixelsDataFrame"))
error.df.plot      <- rbind(error.2014, error.2001)
error.df.plot$Year <- c(rep('2014', times = nrow(error.2014)), 
                        rep('2001', times = nrow(error.2001)))

error.df.plot %>%
  pivot_longer(c(WHAM_INFERNO, GFED5), values_to = 'Burned fraction') %>%
  ggplot(aes(x = x, y = y, fill = `Burned fraction`)) +
  geom_raster() + 
  scale_fill_viridis_c(na.value = 'white', 
                       trans = scales::pseudo_log_trans(sigma = 0.05), 
                       labels = c(0, 0.1, 0.25, 0.5, 1), 
                       breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(Year~name)+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        text = element_text(size = 14), 
        legend.text = element_text(angle = 90, hjust = 0.8)) +
  theme(legend.position="bottom")

ggsave('GFED_WHAM_maps.tiff', width = 13, height = 8)

###########################################################################################

### What is error associated with?

###########################################################################################

setwd('C:/Users/Oli/Documents/PhD/Model development/Data/Processed secondary/Merged with JULES')
base_dir           <- getwd()

setwd(paste0(base_dir, '/Socio economic/GDP and HDI'))
HDI            <- brick('HDI_JULES.tif')
GDP            <- brick('GDP_JULES.tif')

setwd('C:/Users/Oli/Documents/PhD/Model development/Data/wham_dynamic')
Pop_dense      <- brick('Pop_dense.nc')
Past.lc        <- brick('Pasture.nc')
Range.lc       <- brick('Rangeland.nc')

setwd(paste0(base_dir, '/Climate/JULES biophysical/NPP'))
NPP            <- brick(unlist(lapply(list.files(), raster)))

setwd(paste0(base_dir, '/Climate/JULES biophysical/ET'))
ET             <- brick(unlist(lapply(list.files(), raster)))


load(url("https://github.com/OliPerkins1987/Fire_GBM/blob/main/Data/Plot_data.RData?raw=true"))

error.df            <- data.frame(error = abs(as.data.frame(error.rast, long = T)[, 2]))

error.df$Trans      <- as.data.frame(trans[[12:25]], long = T)[, 2]
error.df$Pre        <- as.data.frame(pre[[12:25]], long = T)[, 2]
error.df$Industrial <- as.data.frame(intense[[12:25]], long = T)[, 2]
error.df$Post       <- as.data.frame(Post[[12:25]], long = T)[, 2]
error.df$ET         <- as.data.frame(ET[[12:25]], long = T)[, 2]
error.df$NPP        <- as.data.frame(NPP[[12:25]], long = T)[, 2]
error.df$HDI        <- as.data.frame(HDI[[12:25]], long = T)[, 2]
error.df$GDP        <- log(as.data.frame(GDP[[12:25]], long = T)[, 2])
error.df$Pop_dense  <- log(as.data.frame(Pop_dense[[12:25]], long = T)[, 2])
error.df$Cropland   <- as.data.frame(Cropland[[12:25]], long = T)[, 2]
error.df$Pasture    <- as.data.frame(Past.lc[[12:25]], long = T)[, 2]
error.df$Rangeland  <- as.data.frame(Range.lc[[12:25]], long = T)[, 2]

error.df$UM.error  <- as.data.frame(abs(error.rast * ((UM.rast[[12:25]]) /result.rast[[12:25]])), long = T)[, 2]
error.df$M.error   <- as.data.frame(abs(error.rast * ((M.rast[[12:25]]) /result.rast[[12:25]])), long = T)[, 2]

cor(error.df[complete.cases(error.df), ], method = 'spearman')


fit1     <- lm(log1p(error) ~. - GDP - UM.error - M.error, data = error.df[complete.cases(error.df), ])
fit.UM   <- lm(log1p(UM.error) ~. - GDP - M.error - error, data = error.df[complete.cases(error.df), ])
fit.M    <- lm(M.error ~. - GDP - UM.error - error, data = error.df[complete.cases(error.df), ])


#error.df$Pop        <- as.data.frame(Population)
library(tree)
fit.tree  <- tree(log1p(error) ~ ., data= error.df[complete.cases(error.df), ])
fit.tree2 <- tree(error ~., data= error.df[complete.cases(error.df), ])

error.pred <- predict(fit.tree2, error.df)
1 - (sqrt(mean((error.pred - error.df$error)^2, na.rm = T)) / sd(error.df$error -mean(error.df$error, na.rm = T), na.rm = T))

library(caret)
model <- caret::train(
  error ~ ., 
  error.df[complete.cases(error.df), -c(14, 15)],
  method = "rpart", 
  trControl = trainControl(
    method = "cv", number = 10,
    verboseIter = TRUE
  )
)



