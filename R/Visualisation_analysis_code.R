
library(raster)
library(tidyverse)
library(rPref)
library(DescTools)

load("visualisation_data.RData")

##################################################################################

### 1) Evaluation of models

##################################################################################

############################
### Calculate pareto spaces
############################

col.key <- c('GFED_cor', 'Global_BA', 'SSA_cor', 'SSA_BA', 'I.thresh_glob', 'I.thresh_SSA', 
             'I.glob', 'I.SSA', 'Tranche', 'I_global')
cor.res <- list()

for(i in 1:length(cor.comb)) {
  
  cor.error           <- cor.comb[[i]][cor.comb[[i]]$I.thresh_glob == 1, colnames(cor.comb[[i]]) %in% col.key]
  cor.error$ID        <- (1:nrow(cor.comb[[i]]))[cor.comb[[i]]$I.thresh_glob == 1]
  cor.error$Global_BA <- abs(cor.error$Global_BA - 
                               (802.5*10000)) / (802.5*10000) 
  cor.error$GFED_cor  <- 1 - cor.error$GFED_cor / max(cor.error$GFED_cor)
  
  front <- psel(cor.error, low(Global_BA) * low(GFED_cor), show_level = TRUE, top_level = 1)
  
  ### Var importance
  
  comb.pareto                <- cor.comb[[i]]
  comb.pareto$Tranche        <- ifelse(comb.pareto$I.thresh_glob == 1,'NROY', 'Ruled out')
  comb.pareto$Tranche[front$ID] <- 'Pareto'
  
  cor.comb[[i]]<- comb.pareto
  
  cor.res[[i]] <- comb.pareto
  
}


############################
### Figure 3
############################

for(i in 1:length(cor.res)) {
  cor.res[[i]] <- mutate(cor.res[[i]], 'Model' = 
                           c('WI-EO', 'INFERNO_V1', 'WI-JULES')[i])
}

cor.res %>% plyr::rbind.fill() %>%
  mutate(Global_BA= Global_BA/10000) %>%
  mutate(Tranche = factor(Tranche, levels = c('Ruled out', 'NROY', 'Pareto'), ordered = T)) %>%
  #pivot_longer(c(GFED_cor, Global_BA)) %>%
  ggplot(aes(x = Model, y = Global_BA, fill = Tranche)) + geom_boxplot() +
  ylab("Burned area (Mha)") + theme_classic() + theme(text = element_text(size= 12)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  scale_fill_viridis_d()

cor.res %>% plyr::rbind.fill() %>%
  mutate(Global_BA= Global_BA/10000) %>%
  mutate(Tranche = factor(Tranche, levels = c('Ruled out', 'NROY', 'Pareto'), ordered = T)) %>%
  #pivot_longer(c(GFED_cor, Global_BA)) %>%
  ggplot(aes(x = Model, y = GFED_cor, fill = Tranche)) + geom_boxplot() +
  ylab("Correlation (r)") + theme_classic() + theme(text = element_text(size= 12)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  scale_fill_viridis_d()


##################################################################################

### 2) Global maps

##################################################################################


##########################################

## Figure 4

##########################################

map.df <- as.data.frame(as(brick(list(WHAM_Empirical = mean(BA_Empirical[[c(1:25)]]),
                                      WHAM_INFERNO = mean(BA_INFERNO[[c(1:25)]]),
                                      Unmanaged_EO = mean(UM_Empirical[[1:25]]),
                                      Unmanaged_INFERNO   = mean(UM_INFERNO[[1:25]]),
                                      Managed_EO   = mean(M_Empirical[[1:25]]),
                                      Managed_INFERNO   = mean(M_INFERNO[[1:25]]))), "SpatialPixelsDataFrame"))

map.df %>%
  pivot_longer(3:6) %>%
  #mutate(value = ifelse(value < 0.00001, NA, value)) %>%
  mutate(Model = ifelse(grepl('INFERNO', name), 'WI-JULES', 'WI-EO')) %>%
  mutate(name = ifelse(grepl('Unmanaged', name), 'Unmanaged', 'Managed')) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  theme_classic() +
  geom_raster() + borders(colour = 'black') + 
  scale_fill_gradientn(colours = c('white', 'blue', 'red'), na.value = 'White', 
                       trans = scales::pseudo_log_trans(sigma = 0.05), 
                       labels = c(0, 0.1, 0.25, 0.5, 1), 
                       breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(Model ~ name)+
  labs(fill = 'Burned fraction') + 
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="bottom", 
        text = element_text(size = 14), 
        legend.text = element_text(angle = 65, hjust = 0.25, vjust = 0.5))


############################

### Figure 6

############################

map.df <- as.data.frame(as(brick(list(GFED5 = mean(MOD.JULES[[c(1:14)]]),
  WI_EO = mean(BA_Empirical[[c(12:25)]]),
  WI_JULES = mean(BA_INFERNO[[c(12:25)]]))), "SpatialPixelsDataFrame"))


map.df %>%
  pivot_longer(1:3, names_to = 'Model') %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  theme_classic() +
  geom_raster() + 
  borders(colour = 'black') + 
  scale_fill_gradientn(colours = c('white', 'blue', 'red'), na.value = 'White', 
                       trans = scales::pseudo_log_trans(sigma = 0.05), 
                       labels = c(0, 0.1, 0.25, 0.5, 1), 
                       breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(Model ~ .)+
  labs(fill = 'Burned fraction') + 
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="bottom", 
        text = element_text(size = 14), 
        legend.text = element_text(angle = 65, hjust = 0.25, vjust = 0.5))


#ggsave('Burned_area_GFEDv2.tiff', height = 7.5, width = 5)

#######################################################################

## 3) Time series

#######################################################################

ts.df <- data.frame(WHAM_Empirical = sapply(1:25, function(z) {sum((BA_Empirical[[z]]*area(JULES.mask))[], na.rm = T)}),
                    WHAM_INFERNO   = sapply(1:25, function(z) {sum((BA_INFERNO[[z]]*area(JULES.mask))[], na.rm = T)}),
                    INFERNO_v1     =  sapply(1:25, function(z) {sum((BA_OG[[z]]*area(JULES.mask))[], na.rm = T)}),
                    UM_Empirical   =  sapply(1:25, function(z) {sum((UM_Empirical[[z]]*area(JULES.mask))[], na.rm = T)}),
                    UM_INFERNO     =  sapply(1:25, function(z) {sum((UM_INFERNO[[z]]*area(JULES.mask))[], na.rm = T)}),
                    M_Empirical   =  sapply(1:25, function(z) {sum((M_Empirical[[z]]*area(JULES.mask))[], na.rm = T)}),
                    M_INFERNO     =  sapply(1:25, function(z) {sum((M_INFERNO[[z]]*area(JULES.mask))[], na.rm = T)}),
                    Year = 1990:2014)


##################################################################################

### Regional plots

##################################################################################

cont.ba        <- rbind(WHAM_Empirical_WB,WHAM_INFERNO_WB)
cont.ba$Model  <- c(rep('WI-EO', times = 95), rep('WI-JULES', times = 95))
cont.ba$Method <- rep(c(rep('Burned_area', times= 25), rep('GFED5', times = 20), 
                        rep('Managed fire', times = 25), rep('Unmanaged fire', times = 25)), times = 2)

cont.ba$Method <- ifelse(cont.ba$Method == 'Burned_area' & cont.ba$Model == 'WI-EO', 
                         'WI-EO', cont.ba$Method)
cont.ba$Method <- ifelse(cont.ba$Method == 'Burned_area' & cont.ba$Model == 'WI-JULES', 
                         'WI-JULES', cont.ba$Method)

####################

### Figure 5

####################


cont.ba %>% pivot_longer(1:7) %>%
  mutate(name = gsub('\\.', '_', name)) %>%
  mutate(name = gsub('___', ' & ', name)) %>%
  mutate(name = gsub('_', ' ', name)) %>%
  mutate(name = recode(name, 'Latin America & Caribbean' = 'Lat. Am & Car', 
                       'Middle East & North Africa' = 'MENA', 
                       "Europe & Central Asia" = "Eu. & Central Asia")) %>%
  filter(grepl('Unmanaged', Method) | grepl('Managed', Method)) %>%
  ggplot(aes(x = Year, y = value/10000, colour = Model, linetype = Method)) + geom_line(size = 1) +
  facet_grid(.~name) + ylab('Burned area (Mha)') + scale_colour_viridis_d() +
  labs('linetype' = 'Fire type') + 
  theme_classic() + theme(text = element_text(size = 13)) + 
  scale_x_continuous(breaks = c(1990, 2000, 2010))

####################################################

### Figure 7

####################################################

cont.ba %>% pivot_longer(1:7) %>%
  mutate(name = gsub('\\.', '_', name)) %>%
  mutate(name = gsub('___', ' & ', name)) %>%
  mutate(name = gsub('_', ' ', name)) %>%
  filter(Year %in% 2001:2014) %>%
  mutate(name = recode(name, 'Latin America & Caribbean' = 'Lat. Am & Car', 
                       'Middle East & North Africa' = 'MENA')) %>%
  filter(Method %in% c('WI-JULES', 'WI-EO', 'GFED5')) %>%
  ggplot(aes(x = Year, y = value/10000, colour = Method)) + geom_line(size = 1) +
  facet_grid(.~name) + ylab('Burned area (Mha)') + scale_colour_viridis_d() +
  theme_classic() + theme(text = element_text(size = 14))

#########################################################################

### Autocorrelation

#########################################################################

detrend <- function(x) {
  df  <- data.frame(x = x, t = 1:length(x))
  tmp <- lm(x ~ t, data = df)
  print(summary(tmp))
  x   <- x - tmp$fitted.values
  x
}


### global

ts.df <- data.frame(GFED5 =  sapply(1:14, function(i) {sum((MOD.JULES[[i]] * area(MOD.JULES))[], na.rm = T)}),
                    WI_EO = sapply(12:25, function(z) {sum((BA_Empirical[[z]]*area(JULES.mask))[], na.rm = T)}),
                    WI_JULES   = sapply(12:25, function(z) {sum((BA_INFERNO[[z]]*area(JULES.mask))[], na.rm = T)}),
                    INFERNO_baseline     =  sapply(12:25, function(z) {sum((BA_OG[[z]]*area(JULES.mask))[], na.rm = T)}),
                    UM_Empirical   =  sapply(12:25, function(z) {sum((UM_Empirical[[z]]*area(JULES.mask))[], na.rm = T)}),
                    UM_INFERNO     =  sapply(12:25, function(z) {sum((UM_INFERNO[[z]]*area(JULES.mask))[], na.rm = T)}),
                    M_Empirical   =  sapply(12:25, function(z) {sum((M_Empirical[[z]]*area(JULES.mask))[], na.rm = T)}),
                    M_INFERNO     =  sapply(12:25, function(z) {sum((M_INFERNO[[z]]*area(JULES.mask))[], na.rm = T)}),
                    Year = 2001:2014)

glob.acf <- ts.df %>% summarise_at(1:8, detrend)
glob.acf %>% summarise_at(2:8, function(x) {cor(x, (.)$GFED5)})

glob.acf %>% mutate(x = 2001:2014) %>%
  pivot_longer(1:4) %>%
ggplot(aes(x = x, y= value/10000, colour = name)) + geom_line(size = 1.5)



###################################################################

### 4) Drivers of burned area

###################################################################


pix <- brick(list(raster::which.max(abs(Emp.pix)), 
                  raster::which.max(abs(Inf.pix)), 
                  MOD.JULES[[14]]))

pix         <- as.data.frame(pix, xy = TRUE)
pix$layer.1 <- ifelse(pix$layer.14 < 0.001, NA, pix$layer.1)
pix$layer.2 <- ifelse(pix$layer.14 < 0.001, NA, pix$layer.2)

colnames(pix)[3:4] <- c('WI-EO', 'WI-JULES')

#######################################

### Figure 8a)

#######################################

pix %>%
  pivot_longer(3:4) %>%
  filter(!is.na(value)) %>%
  mutate(value = recode(value, '1' = 'Ignitions', 
                        '2' = 'Road density', '3' = 'Flammability', 
                        '4' = 'Fire suppression', '5' = 'Cropland conversion')) %>%
  mutate(Model = ifelse(name == 'layer.1', 'WHAM-EO', 'WHAM-INFERNO')) %>%
  ggplot(aes(x = x, y = y, fill = factor(value))) +
  theme_classic() +
  geom_raster() + 
  scale_fill_manual(values = c('Red', 'Yellow', 'Grey', 'Green', 'Blue'), 
                    na.value = 'White') + 
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(. ~ name)+
  labs(fill = 'Most correlated independent variable:') + 
  borders() + 
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="bottom", 
        text = element_text(size = 11))


#######################################################

### Figure 8b)

#######################################################

ba_delta <- brick(list('WI-EO' = (UM_Empirical[[25]] - UM_Empirical[[12]]), 
                       'WI-JULES' = (UM_INFERNO[[25]] - UM_INFERNO[[12]])))

ba_delta <- as.data.frame(ba_delta, xy = TRUE)


ba_delta %>%
  pivot_longer(3:4) %>%
  mutate(name = gsub('\\.', '-', name)) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  theme_classic() +
  geom_raster() + 
  scale_fill_gradientn(colours = c('blue', 'gray96', 'red'), na.value = 'White', 
                       values = scales::rescale(c(-0.4, -0.1, 0, 0.1, 0.4)), 
                       breaks = c(-0.4, -0.1, 0.1, 0.4),
                       limits = c(-0.4, 0.4)) + 
  scale_y_continuous(limits = c(-60, 85)) +
  ylab('') + xlab('') + facet_grid(. ~ name)+
  labs(fill = "Change in burned area fraction (2001-2014)") + 
  borders() + 
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="bottom", 
        text = element_text(size = 12))

####################################################################################

### Table 2 - regression coefficients

####################################################################################

### WI-EO

EO.regression.dat %>% group_by(x, y) %>%
  mutate_at(c(4, 6:10), function(k) {c(NA, diff(k))}) -> EO.IAV.dat

EO.regression.dat %>% filter(Year %in% c(2001, 2014)) %>%
  group_by(x, y) %>% mutate_at(c(4, 6:10), function(k) {c(NA, diff(k))}) -> EO.Trend.dat

EO.IAV.dat$Road <- (0-EO.IAV.dat$Road)
EO.Trend.dat$Road<- (0-EO.Trend.dat$Road)

EO.IAV.lm   <- lm(BA_Empirical ~ Igs + Flammability + Cropland + Road + Suppression + 0, 
                data = EO.IAV.dat)

EO.Trend.lm  <- lm(BA_Empirical ~ Igs + Flammability + Cropland + Road + Suppression + 0, 
                data = EO.Trend.dat)

### WI-JULES

JULES.regression.dat %>% group_by(x, y) %>%
  mutate_at(c(4, 6:10), function(k) {c(NA, diff(k))}) -> JULES.IAV.dat

JULES.regression.dat %>% filter(Year %in% c(2001, 2014)) %>%
  group_by(x, y) %>% mutate_at(c(4, 6:10), function(k) {c(NA, diff(k))}) -> JULES.Trend.dat

JULES.IAV.dat$Road <- (0-JULES.IAV.dat$Road)
JULES.Trend.dat$Road<- (0-JULES.Trend.dat$Road)

JULES.IAV.lm   <- lm(BA_Empirical ~ Igs + Flammability + Cropland + Road + Suppression + 0, 
                  data = JULES.IAV.dat)

JULES.Trend.lm <- lm(BA_Empirical ~ Igs + Flammability + Cropland + Road + Suppression + 0, 
                   data = JULES.Trend.dat)
