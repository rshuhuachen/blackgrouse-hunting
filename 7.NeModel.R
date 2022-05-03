################### Modelling Ne ################

## load packages
library(tidyverse); library(RLDNe); library(data.table)

## load in data and only pick years 2005 and 2007
## only for adult males and females, seperately, as well as hunting seperately

males <- fread("data/cleandata/Unsplit.microsat.males.noLOCUS1+13.csv")
females <- fread("data/cleandata/Unsplit.microsat.females.noLOCUS1+13.csv") 

# include year and hunted/unhunted in dataframes
load("data/rawdata/Fulldata_adults.RData")

males.df <- right_join(hunted.ad[,c(1,3,8)], males)
females.df <- right_join(hunted.ad[,c(1,3,8)], females)

# just for contemporary method, choose all years
males.hunted.df <- subset(males.df,  hunt == "hunted")
males.unhunted.df <- subset(males.df,  hunt == "unhunted")
nrow(males.hunted.df) + nrow(males.unhunted.df)

females.hunted.df <- subset(females.df,  hunt == "hunted")
females.unhunted.df <- subset(females.df,  hunt == "unhunted")
nrow(females.hunted.df) + nrow(females.unhunted.df)


# # select only 2005 and 2007
# males.hunted.df <- subset(males.df,  hunt == "hunted" & (year == 2005 | year == 2007 ))
# males.unhunted.df <- subset(males.df,  hunt == "unhunted" & (year == 2005 | year == 2007 ))
# nrow(males.hunted.df) + nrow(males.unhunted.df)
# 
# females.hunted.df <- subset(females.df,  hunt == "hunted" & (year == 2005 | year == 2007 ))
# females.unhunted.df <- subset(females.df,  hunt == "unhunted" & (year == 2005 | year == 2007 ))
# nrow(females.hunted.df) + nrow(females.unhunted.df)

## write out as tables
write.csv(males.hunted.df, "data/neestimator/NeEstimator.males.hunted.txt", quote = F, row.names=F)
write.csv(males.unhunted.df, "data/neestimator//NeEstimator.males.unhunted.txt", quote = F, row.names=F)
write.csv(females.hunted.df, "data/neestimator//NeEstimator.females.hunted.txt", quote = F, row.names=F)
write.csv(females.unhunted.df, "data/neestimator//NeEstimator.females.unhunted.txt", quote = F, row.names=F)

#### NeEstimator #####

#males all
gpfile_maleshunted <- write_genepop_zlr(loci = males.df[,5:ncol(males.df)], pops=males.df$pop,
                                        ind.ids = males.df$id, folder = "", filename="genepop_output_malesall.txt",
                                        missingVal = NA, diploid = T, ncode =2)
gpfile_femaleshunted <- write_genepop_zlr(loci = females.df[,5:ncol(females.df)], pops=females.df$pop,
                                        ind.ids = females.df$id, folder = "", filename="genepop_output_femalesall.txt",
                                        missingVal = NA, diploid = T, ncode =2)


#males hunted
setwd("data/neestimator/")

gpfile_maleshunted <- write_genepop_zlr(loci = males.hunted.df[,5:ncol(males.hunted.df)], pops=males.hunted.df$pop,
                                        ind.ids = males.hunted.df$id, folder = "", filename="genepop_output_maleshunted.txt",
                                        missingVal = NA, diploid = T, ncode =2)

param_files_maleshunted<- NeV2_LDNe_create(input_file = gpfile_maleshunted$Output_File ,param_file = "Ne_params_maleshunted.txt", NE_out_file = "Ne_out_maleshunted.txt")

run_LDNe(LDNe_params = param_files_maleshunted$param_file)

Ne_estimates_maleshunted<-readLDNe_tab(path = param_files_maleshunted$Ne_out_tab)

#males unhunted
gpfile_malesunhunted <- write_genepop_zlr(loci = males.unhunted.df[,5:ncol(males.unhunted.df)], pops=males.unhunted.df$pop,
                                        ind.ids = males.unhunted.df$id, folder = "", filename="genepop_output_malesunhunted.txt",
                                        missingVal = NA, diploid = T, ncode =2)

param_files_malesunhunted<- NeV2_LDNe_create(input_file = gpfile_malesunhunted$Output_File ,param_file = "Ne_params_malesunhunted.txt", NE_out_file = "Ne_out_malesunhunted.txt")

run_LDNe(LDNe_params = param_files_malesunhunted$param_file)

Ne_estimates_malesunhunted<-readLDNe_tab(path = param_files_malesunhunted$Ne_out_tab)

#females hunted
gpfile_femaleshunted <- write_genepop_zlr(loci = females.hunted.df[,5:ncol(females.hunted.df)], pops=females.hunted.df$pop,
                                        ind.ids = females.hunted.df$id, folder = "", filename="genepop_output_femaleshunted.txt",
                                        missingVal = NA, diploid = T, ncode =2)

param_files_femaleshunted<- NeV2_LDNe_create(input_file = gpfile_femaleshunted$Output_File ,param_file = "Ne_params_femaleshunted.txt", NE_out_file = "Ne_out_femaleshunted.txt")

run_LDNe(LDNe_params = param_files_femaleshunted$param_file)

Ne_estimates_femaleshunted<-readLDNe_tab(path = param_files_femaleshunted$Ne_out_tab)

#females unhunted
gpfile_femalesunhunted <- write_genepop_zlr(loci = females.unhunted.df[,5:ncol(females.unhunted.df)], pops=females.unhunted.df$pop,
                                          ind.ids = females.unhunted.df$id, folder = "", filename="genepop_output_femalesunhunted.txt",
                                          missingVal = NA, diploid = T, ncode =2)

param_files_femalesunhunted<- NeV2_LDNe_create(input_file = gpfile_femalesunhunted$Output_File ,param_file = "Ne_params_femalesunhunted.txt", NE_out_file = "Ne_out_femalesunhunted.txt")

run_LDNe(LDNe_params = param_files_femalesunhunted$param_file)

Ne_estimates_femalesunhunted<-readLDNe_tab(path = param_files_femalesunhunted$Ne_out_tab)

### Analyse results
Ne_estimates_maleshunted
Ne_estimates_malesunhunted
Ne_estimates_femaleshunted
Ne_estimates_femalesunhunted
