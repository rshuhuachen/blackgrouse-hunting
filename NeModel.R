################### Modelling Ne ################

## load packages
library(tidyverse)

## load in data and only pick years 2005 and 2007
## only for adult males and females, seperately, as well as hunting seperately

males <- fread("data/Unsplit.microsat.males.noLOCUS1+13.csv")
females <- fread("data/Unsplit.microsat.females.noLOCUS1+13.csv") 

# include year and hunted/unhunted in dataframes
load("/data/Microsats_hunted_adults_combined.25.3.22_CLEAN.RData")

males.df <- right_join(hunted.ad[,c(1,3,8)], males)
females.df <- right_join(hunted.ad[,c(1,3,8)], females)

# select only 2005 and 2007
males.hunted.df <- subset(males.df,  hunt == "hunted" & (year == 2005 | year == 2007 ))
males.unhunted.df <- subset(males.df,  hunt == "unhunted" & (year == 2005 | year == 2007 ))
nrow(males.hunted.df) + nrow(males.unhunted.df)

females.hunted.df <- subset(females.df,  hunt == "hunted" & (year == 2005 | year == 2007 ))
females.unhunted.df <- subset(females.df,  hunt == "unhunted" & (year == 2005 | year == 2007 ))
nrow(females.hunted.df) + nrow(females.unhunted.df)

## write out as tables
write.table(males.hunted.df, "P:/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/MLNe/Males.hunted/MLNe.males.hunted.txt", sep = " ", quote = F, row.names=F)
write.table(males.unhunted.df, "P:/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/MLNe/Males.unhunted/MLNe.males.unhunted.txt", sep = " ", quote = F, row.names=F)
write.table(females.hunted.df, "P:/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/MLNe/Females.hunted/MLNe.females.hunted.txt", sep = " ", quote = F, row.names=F)
write.table(females.unhunted.df, "P:/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/MLNe/Females.unhunted/MLNe.females.unhunted.txt", sep = " ", quote = F, row.names=F)
