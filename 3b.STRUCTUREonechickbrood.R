###### STRUCTURE b - taking only one chick per brood ######

library(readxl); library(tidyverse)
# first make a table for GenAlEx
chicks <- read.csv("data/cleandata/Unsplit.microsat.chicks.noLOCUS1+13+14.csv")

chicks.stru <- read.table("data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru")

pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")

### combine with
chicks.df <- left_join(chicks, pops[,c(2,3)], by = c("pop" = "pop_num"))
chicks.hunted <- subset(chicks.df, hunt == "hunted")
chicks.unhunted <- subset(chicks.df, hunt == "unhunted")
