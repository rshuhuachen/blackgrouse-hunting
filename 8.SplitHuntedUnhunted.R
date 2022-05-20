##### Script for creating dataframes splitting into hunted/unhunted ####

### Load in unsplit files which are necessary for GenAlEx (spatial autocor) as well as structure files for STRUCTURE

#males
males <- read.csv("data/cleandata/Unsplit.microsat.males.noLOCUS1+13.csv")
males.stru <- read.table("data/cleandata/Microsat.males.noLOCUS1+13.forstructure.stru")

#females
females <- read.csv("data/cleandata/Unsplit.microsat.females.noLOCUS1+13.csv")
females.stru <- read.table("data/cleandata/Microsat.females.noLOCUS1+13.forstructure.stru")

#chicks
chicks <- read.csv("data/cleandata/Unsplit.microsat.chicks.noLOCUS1+13+14.csv")
chicks.stru <- read.table("data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru")

#unrelated chicks
unr.chicks <- read.csv("data/cleandata/Unsplit.microsat.unrelated.chicks.noLOCUS1+13+14.csv")
unr.chicks.stru <- read.table("data/cleandata/Microsat.unrelated.chicks.noLOCUS1+13+14.forstructure.stru")

# pop/hunting info
pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")
hunted.pops <- subset(pops, hunt == "hunted")
unhunted.pops <- subset(pops, hunt == "unhunted")

#combine dataframes to include hunted status

males.hunted <- subset(males, pop %in% hunted.pops$pop_num)
males.unhunted <- subset(males, pop %in% unhunted.pops$pop_num)

males.hunted.stru <- subset(males.stru, V2 %in% hunted.pops$pop_num)
males.unhunted.stru <- subset(males.stru, V2 %in% unhunted.pops$pop_num)

females.hunted <- subset(females, pop %in% hunted.pops$pop_num)
females.unhunted <- subset(females, pop %in% unhunted.pops$pop_num)

females.hunted.stru <- subset(females.stru, V2 %in% hunted.pops$pop_num)
females.unhunted.stru <- subset(females.stru, V2 %in% unhunted.pops$pop_num)

chicks.hunted <- subset(chicks, pop %in% hunted.pops$pop_num)
chicks.unhunted <- subset(chicks, pop %in% unhunted.pops$pop_num)

chicks.hunted.stru <- subset(chicks.stru, V2 %in% hunted.pops$pop_num)
chicks.unhunted.stru <- subset(chicks.stru, V2 %in% unhunted.pops$pop_num)

unr.chicks.hunted <- subset(unr.chicks, pop %in% hunted.pops$pop_num)
unr.chicks.unhunted <- subset(unr.chicks, pop %in% unhunted.pops$pop_num)

unr.chicks.hunted.stru <- subset(unr.chicks.stru, V2 %in% hunted.pops$pop_num)
unr.chicks.unhunted.stru <- subset(unr.chicks.stru, V2 %in% unhunted.pops$pop_num)


# write out files
write.csv(males.hunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.males.hunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.csv(males.unhunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.males.unhunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.table(males.hunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.males.hunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")
write.table(males.unhunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.males.unhunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")

write.csv(females.hunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.females.hunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.csv(females.unhunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.females.unhunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.table(females.hunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.females.hunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")
write.table(females.unhunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.females.unhunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")

write.csv(chicks.hunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.chicks.hunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.csv(chicks.unhunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.chicks.unhunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.table(chicks.hunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.chicks.hunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")
write.table(chicks.unhunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.chicks.unhunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")

write.csv(unr.chicks.hunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.unrelated.chicks.hunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.csv(unr.chicks.unhunted, "data/cleandata/SplitHuntedUnhunted/Unsplit.microsat.unrelated.chicks.unhunted.noLOCUS1+13.csv",row.names=F, quote =F)
write.table(unr.chicks.hunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.unrelated.chicks.hunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")
write.table(unr.chicks.unhunted.stru, "data/cleandata/SplitHuntedUnhunted/Microsat.unrelated.chicks.unhunted.noLOCUS1+13.forstructure.stru",row.names=F, quote =F,sep = " ")

### Format for GenAlEx (add site names rather than num)

# data on coordinates and adult + chick data

#adults.coor2 <- read.csv("P:\\Black Grouse PhD\\Projects\\Hunting_Microsats\\Data\\Microsat_data_hunted_adults_CLEANED.csv")
#chicks.coor <- read.csv("P:\\Black Grouse PhD\\Projects\\Hunting_Microsats\\Data\\Microsat_data_hunted_chicks_CLEANED.csv")

#males hunted
# add site names rather than numbers and coordinates
males.hunted.genalex <- left_join(males.hunted, pops[,-3], by = c("pop" = "pop_num"))
males.hunted.genalex <- males.hunted.genalex[,c(1,27,3:26,28,29)]
names(males.hunted.genalex)[2] <- "SITE"
names(males.hunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
males.hunted.genalex[is.na(males.hunted.genalex)] <- 0
males.hunted.genalex[(males.hunted.genalex == -9)] <- 0

summary(as.factor(males.hunted.genalex$SITE))

#males unhunted
# add site names rather than numbers and coordinates
males.unhunted.genalex <- left_join(males.unhunted, pops[,-3], by = c("pop" = "pop_num"))
males.unhunted.genalex <- males.unhunted.genalex[,c(1,27,3:26,28,29)]
names(males.unhunted.genalex)[2] <- "SITE"
names(males.unhunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
males.unhunted.genalex[is.na(males.unhunted.genalex)] <- 0
males.unhunted.genalex[(males.unhunted.genalex == -9)] <- 0

summary(as.factor(males.unhunted.genalex$SITE))

#females hunted
# add site names rather than numbers and coordinates
females.hunted.genalex <- left_join(females.hunted, pops[,-3], by = c("pop" = "pop_num"))
females.hunted.genalex <- females.hunted.genalex[,c(1,27,3:26,28,29)]
names(females.hunted.genalex)[2] <- "SITE"
names(females.hunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
females.hunted.genalex[is.na(females.hunted.genalex)] <- 0
females.hunted.genalex[(females.hunted.genalex == -9)] <- 0

summary(as.factor(females.hunted.genalex$SITE))

#females unhunted
# add site names rather than numbers and coordinates
females.unhunted.genalex <- left_join(females.unhunted, pops[,-3], by = c("pop" = "pop_num"))
females.unhunted.genalex <- females.unhunted.genalex[,c(1,27,3:26,28,29)]
names(females.unhunted.genalex)[2] <- "SITE"
names(females.unhunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
females.unhunted.genalex[is.na(females.unhunted.genalex)] <- 0
females.unhunted.genalex[(females.unhunted.genalex == -9)] <- 0

summary(as.factor(females.unhunted.genalex$SITE))

#chicks hunted
# add site names rather than numbers and coordinates
chicks.hunted.genalex <- left_join(chicks.hunted, pops[,-3], by = c("pop" = "pop_num"))
chicks.hunted.genalex <- chicks.hunted.genalex[,c(1,25,3:24,26,27)]
names(chicks.hunted.genalex)[2] <- "SITE"
names(chicks.hunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
chicks.hunted.genalex[is.na(chicks.hunted.genalex)] <- 0
chicks.hunted.genalex[(chicks.hunted.genalex == -9)] <- 0

summary(as.factor(chicks.hunted.genalex$SITE))

#chicks unhunted
# add site names rather than numbers and coordinates
chicks.unhunted.genalex <- left_join(chicks.unhunted, pops[,-3], by = c("pop" = "pop_num"))
chicks.unhunted.genalex <- chicks.unhunted.genalex[,c(1,25,3:24,26,27)]
names(chicks.unhunted.genalex)[2] <- "SITE"
names(chicks.unhunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
chicks.unhunted.genalex[is.na(chicks.unhunted.genalex)] <- 0
chicks.unhunted.genalex[(chicks.unhunted.genalex == -9)] <- 0

summary(as.factor(chicks.unhunted.genalex$SITE))

#unr.chicks hunted
# add site names rather than numbers and coordinates
unr.chicks.hunted.genalex <- left_join(unr.chicks.hunted, pops[,-3], by = c("pop" = "pop_num"))
unr.chicks.hunted.genalex <- unr.chicks.hunted.genalex[,c(1,25,3:24,26,27)]
names(unr.chicks.hunted.genalex)[2] <- "SITE"
names(unr.chicks.hunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
unr.chicks.hunted.genalex[is.na(unr.chicks.hunted.genalex)] <- 0
unr.chicks.hunted.genalex[(unr.chicks.hunted.genalex == -9)] <- 0

summary(as.factor(unr.chicks.hunted.genalex$SITE))

#unr.chicks unhunted
# add site names rather than numbers and coordinates
unr.chicks.unhunted.genalex <- left_join(unr.chicks.unhunted, pops[,-3], by = c("pop" = "pop_num"))
unr.chicks.unhunted.genalex <- unr.chicks.unhunted.genalex[,c(1,25,3:24,26,27)]
names(unr.chicks.unhunted.genalex)[2] <- "SITE"
names(unr.chicks.unhunted.genalex)[1]<- "CODE"

### Change missing data to 0 for GenAIEx 
unr.chicks.unhunted.genalex[is.na(unr.chicks.unhunted.genalex)] <- 0
unr.chicks.unhunted.genalex[(unr.chicks.unhunted.genalex == -9)] <- 0

summary(as.factor(unr.chicks.unhunted.genalex$SITE))
## Then format exactly as GenAIEx requires in excel manually

write.csv(males.hunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_males_hunted.csv",quote=F, row.names=F)
write.csv(males.unhunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_males_unhunted.csv",quote=F, row.names=F)
write.csv(females.hunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_females_hunted.csv",quote=F, row.names=F)
write.csv(females.unhunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_females_unhunted.csv",quote=F, row.names=F)
write.csv(chicks.hunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_chicks_hunted.csv",quote=F, row.names=F)
write.csv(chicks.unhunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_chicks_unhunted.csv",quote=F, row.names=F)
write.csv(unr.chicks.hunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_unrelated_chicks_hunted.csv",quote=F, row.names=F)
write.csv(unr.chicks.unhunted.genalex, file = "data/cleandata/SplitHuntedUnhunted/GenAlEx_unrelated_chicks_unhunted.csv",quote=F, row.names=F)

