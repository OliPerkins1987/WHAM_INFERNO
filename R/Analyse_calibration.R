
library(tidyverse)
library(raster)
library(rPref)
library(DescTools)


###########################################

### GFED eval

###########################################

setwd('C:/Users/Oli/Documents/PhD/Thesis/Draft content/ABM Coupled/Results')

cor.comb <- lapply(list.files(pattern = '*.csv'), read.csv)
cor.comb <- lapply(cor.comb, function(x) {mutate(x, I.thresh_glob = abs(((Global_BA/10000) - 802.5) / sqrt(106.7325^2 + 68.33^2)) <=3)})
cor.comb <- lapply(cor.comb, function(x) {if(nrow(x) > 10000) {x[1:10000, ]} else{x}})


#######################################################

### Error trade offs

#######################################################

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

#########################

### Table 1

#########################

lapply(cor.res[c(1, 5:6)], function(x) {
  
  x %>% group_by(Tranche) %>% 
    summarise(cor.mean = mean(GFED_cor), BA.mean = mean(Global_BA), 
              cor.min  = min(GFED_cor), BA.min = min(Global_BA),
              cor.max  = max(GFED_cor), BA.max = max(Global_BA), 
              count    = n())
  
})

### Are correlations significantly different?
std <- function(x) sd(x)/sqrt(length(x))
ZTest(FisherZ(as.numeric(cor.res[[2]]$GFED_cor[cor.res[[2]]$Tranche == 'Pareto'])), 
  FisherZ(as.numeric(cor.res[[3]]$GFED_cor[cor.res[[3]]$Tranche == 'Pareto'])), 
    alternative = 'greater', 
  sd_pop = std(FisherZ(as.numeric(cor.res[[2]]$GFED_cor[cor.res[[2]]$Tranche == 'Pareto']) - 
                         std(FisherZ(as.numeric(cor.res[[3]]$GFED_cor[cor.res[[3]]$Tranche == 'Pareto']))))))


######################################################

### Filter variables

######################################################

### prepare results

cor.dist    <- cor.res
res.kruskal <- list()

for(j in 1:length(cor.res)) {

### NROY
comb.filt  <- cor.res[[j]]
res.t_test <- list()
res.kruskal[[j]] <- list()

for(i in 1:length(colnames(cor.res[[j]])[!colnames(cor.res[[j]]) %in% col.key])) {

  res.t_test[[i]] <- t.test(comb.filt[comb.filt$Tranche == 'Pareto', i], 
                            comb.filt[comb.filt$Tranche == 'NROY', i])
  
  res.kruskal[[j]][[i]] <- kruskal.test(cor.res[[j]][, i] ~ cor.res[[j]]$Tranche)$p.value
  
  }

cor.dist[[j]] <- cor.res[[j]] %>% 
  dplyr::select(c(which(sapply(res.t_test, function(z) {z$p.value})<(0.05/20)), ncol(cor.res[[j]]))) %>%
  pivot_longer(colnames(comb.filt)[which(sapply(res.t_test, function(z) {z$p.value[1]})<(0.05/20))])

}

############################

### Analyse NROY vs ruled out

############################

### get t-values
setNames(sapply(res.t_test[which(sapply(res.t_test, 
      function(z) {z$p.value})<(0.05/20))], function(z) {z$statistic}), 
         colnames(cor.comb[[1]])[which(sapply(res.t_test, function(z) {z$p.value})<(0.05/20))])


### Figure 1
cor.res[[1]] %>%
  dplyr::select(c(which(sapply(res.kruskal[[1]], function(z) {z[1]})<(0.05/20)), 26)) %>%
  pivot_longer(colnames(cor.res[[1]])[which(sapply(res.kruskal[[1]], function(z) {z[1]})<(0.05/20))]) %>%
  filter(!(grepl('Managed', name) & value > 1.5)) %>%
  mutate(name = recode(name, 'arson' = 'Arson', 'Grass1' = 'Grass C3', 'Grass2' = 'Grass C4', 
   'igs' = 'Background fires', 'lightning' = 'Lightning ignitions', 'Managed_pasture' = 'Managed pasture', 
   'Managed_veg' = 'Managed vegetation', 'Pasture2' = 'Pasture C4', 'Road' = 'Road density', 
   'scaling' = 'Fire-ignitions-scaling', 'Tree1' = 'Broadleaf Tree', 'Tree2' = 'Needleaf Tree')) %>%
  group_by(name) %>% mutate(value = value/mean(value)) %>% ungroup() %>%
  ggplot(aes(x = name, y = value, fill = Tranche)) + geom_boxplot() +
  theme_classic() + theme(text = element_text(size= 12)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  ylab('Parameter values (Overall mean = 1)') + scale_fill_viridis_d()



setwd('C:/Users/Oli/Documents/PhD/Thesis/Draft content/ABM Coupled/Charts')
#ggsave('Par_drivers.tiff', width = 12, height = 7)


############################################

### NROY vs pareto (rerun results above!)

############################################

cor.dist    <- lapply(1:length(cor.res), function(x) {cor.dist[[x]] %>% mutate(Model = as.character(x))}) %>% dplyr::bind_rows()


### Figure 2
cor.dist %>% 
  filter(Model %in% c('1', '5', '6')) %>%
  mutate(Model = recode(Model, '1' = 'WHAM-INFERNO',
   '5' = 'INFERNO_v1', '6'= 'INFERNO_Road')) %>%
  filter(Tranche != 'Ruled out') %>% filter(name != 'Managed_pasture') %>%
  group_by(name, Model) %>% mutate(value = value/mean(value)) %>% ungroup() %>% 
  mutate(name = recode(name, 'suppression' = 'Suppression', 
     'Grass1' = 'Grass C3', 'igs' = "B'round fires",
      'Pasture2' = 'Pasture C4', 'Road' = 'Road density', 'Shrub' = 'Shrub')) %>%
  ggplot(aes(x = value, fill = Tranche)) + 
  geom_histogram(aes(y = ..density..), colour = 'black', position = position_dodge(), bins = 15) +
  theme_classic() + theme(text = element_text(size = 14)) +
  xlab('Scaled parameter values (Overall mean = 1)') + 
  ylab('Probability density') + facet_grid(name~Model)+ scale_fill_viridis_d()#Model

ggsave('Par_hists.tiff', width = 9, height = 6)



################################################################

### Weighted correlation

################################################################

### can we weight logging by its impact on burned area?
impact.BA <- list()
impact.cor<- list()
impact.rng<- list()

cor.w   <- cor.comb[[1]][!(cor.comb[[1]]$Managed_pasture > 1.5 |cor.comb[[1]]$Managed_veg > 1.5) & cor.comb[[1]]$Tranche != 'Ruled out', ]
cor.teg <- cor.w[cor.w$Tranche != 'Ruled out', ]

for(i in 1:20) {
  
  impact.BA[[i]] <- cor(cor.teg[, i], cor.teg$Global_BA)
  impact.cor[[i]]<- cor(cor.teg[, i], cor.teg$GFED_cor)
  impact.rng[[i]]<-  as.numeric(cor.w %>% #group_by(Tranche == 'Pareto') %>% 
    dplyr::select(i) %>%
     summarise_if(is.numeric, function(x) {max(x) - min(x)}) %>%
       summarise_if(is.numeric, function(x) {x[2]/x[1]}) %>% ungroup())
  
}

dat.w         <- data.frame(cor.BA = setNames(unlist(impact.BA), 
                                            colnames(cor.comb[[1]])[1:20]))
dat.w$cor.cor <- unlist(impact.cor)
dat.w$cor.rng <- unlist(impact.rng)

dat.w %>%
  mutate(Parameter = row.names(dat.w)) %>%
  mutate(cor.weight = cor.cor/cor.BA) %>%
  mutate(Category = recode(Parameter, 'arson' = 'Unmanaged fires', 'fire_slope' = 'Previous fire', 
                           'fire_thresh' = 'Previous fire', 'Grass1' = 'Grass PFT', 'Grass2' = 'Grass PFT', 
                           'igs' = 'Unmanaged fires', 'Intense' = 'Suppression','lightning' = 'Unmanaged fires', 
                           'Managed_pasture' = 'Managed fire', 'Managed_veg' = 'Managed fire', 'Shrub' = 'Shrub PFT',
                           'Pasture1' = 'Pasture PFT', 'Pasture2' = 'Pasture PFT', 'Pre' = 'Suppression', 
                           'Trans' = 'Suppression', 'scaling' = 'Unmanaged fires', 'Tree1' = 'Tree PFT',
                           'Tree2' = 'Tree PFT')) %>%
  group_by(Category) %>% summarise_if(is.numeric, function(z) {mean(abs(z))}) %>%
  mutate(cor.weight = cor.weight/max(cor.weight)) %>%
  pivot_longer(c(cor.BA, cor.cor, cor.weight), names_to = 'Measure') %>%
  ggplot(aes(x = Category, y = value, fill = Measure)) + 
  geom_col(position = position_dodge(), colour = 'black') +
  theme_classic() + xlab('Variable Category') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        text = element_text(size = 14)) + scale_fill_viridis_d() +
         ylab('Abolute correlation coefficient (r)')

ggsave('Weighted_correlation.tiff', width = 9, height = 5)

