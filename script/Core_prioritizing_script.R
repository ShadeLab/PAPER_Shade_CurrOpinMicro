library(tidyverse)
library(reshape2)
library(vegan)
theme_set(theme_light())

setwd('~/Documents/git/PAPER_Shade_CurrOpinMicro/')

#-----------------------------------------------------------------------------------------
##Prioritizing core microbiome based on TIME

#Example  - switchgrass dataset (Grady et al. (2019))

nReads=1000                                                            # input dataset needs to be rarified and the rarifaction depth included 
otu <- readRDS('data/switchgrassOTUtable.rds')
map <- readRDS('data/switchgrassMAPtable.rds')

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame

# Ranking OTUs based on their occupancy
# For caluclating raking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sequence_name, abun, -otu) %>%
  left_join(map, by = 'sequence_name') %>%
  group_by(otu, sampling_date) %>%
  summarise(time_freq=sum(abun>0)/length(abun),            # frequency of detection between time points
            coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 1 only if occupancy 1 with specific time, 0 if not
  group_by(otu) %>%
  summarise(sumF=sum(time_freq),
            sumG=sum(coreTime),
            nS=length(sampling_date)*2,           
            Index=(sumF+sumG)/nS)                 # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

for(i in otu_ranked$otu){
  # calculating BC dissimilarity based on the 1st ranked OTU
  otu_start=otu_ranked$otu[1]                   
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  # calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. Can be set to the entire length of OTUs in the dataset, however it might take some time if more than 5000 OTUs are included.
  for(i in 2:500){                              
    otu_add=otu_ranked$otu[i]                       
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  # calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition,df_full, by='x_names')
} 

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))

#Creating plot of Bray-Curtis similarity
ggplot(BC_ranked[1:100,], aes(x=factor(BC_ranked$rank[1:100], levels=BC_ranked$rank[1:100]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])), lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+14, y=.1, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))+3, y=.5, label=paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])),")",sep=''), color="blue")

#Creating occupancy abundance plot
occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'

#Fitting neutral model (Burns et al., 2016 (ISME J) - functions are in the sncm.fit.R)
spp=t(otu)
taxon=as.vector(rownames(otu))

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

ggplot() +
  geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
  geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='blue', size=1.8) +
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy")

#Creating a plot of core taxa occupancy by time point
core <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method="total", MARGIN=2)

plotDF <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
  gather(sequence_name, relabun, -otu) %>%
  left_join(map, by = 'sequence_name') %>%
  left_join(otu_ranked, bu='otu') %>%
  filter(otu %in% core) %>% 
  group_by(otu, sampling_date) %>%
  summarise(time_freq=sum(relabun>0)/length(relabun),        
            coreTime=ifelse(time_freq == 1, 1, 0),      
            detect=ifelse(time_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels=otu_ranked$otu[1:34])

ggplot(plotDF,aes(x=otu, time_freq,fill=factor(sampling_date))) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu))) +
  theme(axis.text = element_text(size=6)) +
  labs(x='Ranked OTUs', y='Occupancy by site')

#-----------------------------------------------------------------------------------------
##Prioritizing core microbiome based on the GENOTYPE

#Example data - Bowsher et al. (2019)

nReads=10000
otu <- readRDS('data/mimulusOTUtable.rds')
map <- readRDS('data/mimulusMAPtable.rds')

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance 

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(SampleID, abun, -otu) %>%
  left_join(map, by = 'SampleID') %>%
  group_by(otu, Origin) %>%
  summarise(genotype_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
            coreGenotype=ifelse(genotype_freq == 1, 1, 0)) %>% # 1 only if occupancy 1 with specific genotype, 0 if not
  group_by(otu) %>%
  summarise(sumF=sum(genotype_freq),
            sumG=sum(coreGenotype),
            nS=length(Origin)*2,
            Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

BCaddition <- NULL
for(i in otu_ranked$otu){
  otu_start=otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  
  for(i in 2:500){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition,df_full, by='x_names')
} 

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)
ggplot(BC_ranked[1:200,], aes(x=factor(BC_ranked$rank[1:200], levels=BC_ranked$rank[1:200]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=6, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])), lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+10, y=.08, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))-4, y=.08, label=paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])),")", sep=''), color="blue")

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'

spp=t(otu)
taxon=as.vector(rownames(otu))

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

ggplot() +
  geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
  geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='blue', size=1.8) +
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy")

core <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method="total", MARGIN=2)

plotDF <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
  gather(SampleID, relabun, -otu) %>%
  left_join(map, by = 'SampleID') %>%
  left_join(otu_ranked, bu='otu') %>%
  filter(otu %in% core) %>% 
  arrange(desc(rank))%>%
  group_by(Origin,Genotype, otu) %>%
    summarise(genotype_freq=sum(relabun>0)/length(relabun),        # frequency of detection between time points
              coreGenotype=ifelse(genotype_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
              detect=ifelse(genotype_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels=otu_ranked$otu[1:167])
ggplot(plotDF,aes(x=otu, genotype_freq, group=Genotype,fill=Genotype)) +    geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(~Origin) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu))) +
  theme(axis.text = element_text(size=6)) +
  labs(x='Ranked ZOTUs', y='Occupancy by genotype and origin')

#----------------------------------------------------------------------------------------------------
##Prioritizing core microbiome based on the SITE

#Example data - Stopnisek and Ashley (bioXiv, 2019)

nReads=10000
otu <- readRDS('data/beanOTUtable.rds')
map <- readRDS('data/beanMAPtable.rds')

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance 

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sample_ID, abun, -otu) %>%
  left_join(map, by = 'sample_ID') %>%
  group_by(otu, site) %>%
  summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
            coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
            detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
  group_by(otu) %>%
  summarise(sumF=sum(plot_freq),
            sumG=sum(coreSite),
            nS=length(site)*2,
            Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

BCaddition <- NULL
for(i in otu_ranked$otu){
  otu_start=otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  
  for(i in 2:500){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition,df_full, by='x_names')
} 

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

ggplot(BC_ranked[1:200,], aes(x=factor(BC_ranked$rank[1:200], levels=BC_ranked$rank[1:200]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])), lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))-4, y=.08, label=paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])),')',sep=''), color="blue")

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'

spp=t(otu)
taxon=as.vector(rownames(otu))

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

ggplot() +
  geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
  geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='blue', size=1.8) +
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy")

core <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method="total", MARGIN=2)

plotDF <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
  gather(sample_ID, relabun, -otu) %>%
  left_join(map, by = 'sample_ID') %>%
  left_join(otu_ranked, bu='otu') %>%
  filter(otu %in% core) %>% 
  group_by(site, otu) %>%
  summarise(plot_freq=sum(relabun>0)/length(relabun),        # frequency of detection between time points
            coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
            detect=ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels=otu_ranked$otu[1:174])
plotDF$group <- 1
plotDF$group[plotDF$otu %in% otu_ranked$otu[87:174]] <- 2

ggplot(plotDF,aes(x=otu, plot_freq, group=site,fill=site)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu))) +
  theme(axis.text = element_text(size=6)) +
  labs(x='Ranked OTUs', y='Occupancy by site')

