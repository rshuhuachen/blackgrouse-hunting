###### Here we redo all genetic analyses but taking only one chick per brood ######

library(readxl); library(tidyverse)

### First of all: combine dataframes to select genotypes of only one chick per brood

chicks <- read.csv("data/cleandata/Unsplit.microsat.chicks.noLOCUS1+13+14.csv")
#chicks.stru <- read.table("data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru")
pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")

### combine with
chicks.df <- left_join(chicks, pops[,c(2,3)], by = c("pop" = "pop_num"))

# combine with file that has brood numbers
load("data/rawdata/Fulldata_chicks.RData")

msats.chick <- hunted.chick
length(unique(msats.chick$ID)) #1370 id's

# select unrelated chicks
msats.chick.unrelated <- msats.chick
msats.chick.unrelated <- msats.chick.unrelated %>% group_by(Brood) %>%
  slice_sample(n=1) #randomly select one chick per brood

## select columns for structure

msats.chick.unrelated <- msats.chick.unrelated[,c(13,5,14:41)] #select id, site id and genotypes
names(msats.chick.unrelated)[2] <- "pop"
names(msats.chick.unrelated)[1] <- "id"

msats.chick.unrelated <- msats.chick.unrelated[,-c(3,4,27,28,29,30)] #filter out BG10 and BG20
write.csv(msats.chick.unrelated, file = "data/cleandata/Unsplit.microsat.unrelated.chicks.noLOCUS1+13+14.csv",row.names=F, quote =F)

### Create file for STRUCTURE 

# replace NA with -9 for structure
msats.chick.unrelated.stru <- as.data.frame(msats.chick.unrelated)
msats.chick.unrelated.stru[is.na(msats.chick.unrelated.stru)] <- -9
msats.chick.unrelated.stru[(msats.chick.unrelated.stru == 0)] <- -9
msats.chick.unrelated.stru$pop[which(msats.chick.unrelated.stru$pop == -9)] <- NA

# split data to have two rows per id
#unrelated chick
msats.chick.unrelated.stru.a <- msats.chick.unrelated.stru[c(F,T)]
msats.chick.unrelated.stru.a <- cbind(msats.chick.unrelated.stru[c(1)], msats.chick.unrelated.stru.a)

msats.chick.unrelated.stru.b <- msats.chick.unrelated.stru[c(T,F)]
msats.chick.unrelated.stru.b <- cbind(msats.chick.unrelated.stru[c(2)], msats.chick.unrelated.stru.b)
msats.chick.unrelated.stru.b <- msats.chick.unrelated.stru.b[,c(2,1,3:ncol(msats.chick.unrelated.stru.b))]

#rename col names to bind

names(msats.chick.unrelated.stru.b) <- sub("a", "", names(msats.chick.unrelated.stru.b), fixed = T)
names(msats.chick.unrelated.stru.a) <- sub("b", "", names(msats.chick.unrelated.stru.a), fixed = T)

#bind col together
msats.chick.unrelated.stru.tot <- rbind(msats.chick.unrelated.stru.a, msats.chick.unrelated.stru.b)

#reorder
msats.chick.unrelated.stru.tot <- msats.chick.unrelated.stru.tot %>% arrange(id) %>% arrange(pop)


write.table(msats.chick.unrelated.stru.tot, file = "data/cleandata/Microsat.unrelated.chicks.noLOCUS1+13+14.forstructure.stru",col.names=F, quote=F, row.names=F, sep = " ")


