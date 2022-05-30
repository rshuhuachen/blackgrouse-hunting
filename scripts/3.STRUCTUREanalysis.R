###### STRUCTURE ANALYSIS HUNTING MICROSATS ######

library(data.table); library(ParallelStructure)
library(dplyr); library(tidyr); library(stringr)
library(pophelper)

#### NB: the raw results from structre are not uploaded, but can be fully recreated
#### using the workflow below

##### Running structure for males #####

males <- fread("data/cleandata/Microsat.adults.noLOCUS1+13.forstructure.stru") #located in /data in github directory

infile <- "data/cleandata/Microsat.males.noLOCUS1+13.forstructure.stru"
# system("mkdir data/Results_stru_males")
outpath <- "data/structure/Results_stru_males/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- "1,2,3,4,5,6,7,8,9,10,11,12" #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "/data/structure/hunt_jobs_adults.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_adults.txt', 
                                      n_cpu=45, infile=infile,
                                      outpath=outpath,numinds = nrow(males)/2,
                                      numloci=ncol(males)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)

# use admixture model (noadmix = 0), not using prior information regarding location (locprior =0) 
# not using prior population info (usepopinfo = 0)
# Having multiple family members in the sample also violates the model assumptions which could lead to overestimation
# of K, but currently we do have multiple family members. Despite that, K=1 is still optimum

##### Running structure for females #####
females <- fread("data/cleandata/Microsat.females.noLOCUS1+13.forstructure.stru")
infile <- "data/cleandata/Microsat.females.noLOCUS1+13.forstructure.stru"
# system("mkdir data/structure/Results_stru_females")
outpath <- "data/structure/Results_stru_females/"

ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_adults.txt', 
                                      n_cpu=45, infile=infile,
                                      outpath=outpath,numinds = nrow(females)/2,
                                      numloci=ncol(females)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)

##### Running structure for chicks #####

chicks <- fread("data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru")

infile <- "data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru"
# system("mkdir data/structure/Results_stru_chicks")
outpath <- "data/structure/Results_stru_chicks/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 199 #number of broods

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- "1,2,4,5,10,11,12" #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "/data/structure/hunt_jobs_chicks.txt")

# file path to structure

STR_path='/usr/local/bin/'


# Run Parallel Structure

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='/data/structure/hunt_jobs_chicks.txt', 
                                      n_cpu=45, infile=infile,
                                      outpath=outpath,numinds = nrow(chicks)/2,
                                      numloci=ncol(chicks)-2, noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)

# use admixture model (noadmix = 0), not using prior information regarding location (locprior =0) 
# not using prior population info (usepopinfo = 0)
# Having multiple family members in the sample also violates the model assumptions which could lead to overestimation
# of K, but currently we do have multiple family members. Despite that, K=1 is still optimum

#### Collecting output and inferring K with highest likelihood
# collect output

#males

source("scripts/STRUCTURE_functions.R")

## males 
make_directories("data/structure/Results_stru_males/")
load_and_clumpp("data/structure/Results_stru_males/Run_files/")
load_and_K("data/structure/Results_stru_males/Run_files/")
ks_ad_male <- load_and_K("data/structure/Results_stru_males/Run_files/")
save(ks_ad_male, file = "data/structure/ks_ad_male.R")
optimal_k(ks_ad_male) 
both_k(ks_ad_male) 

## females 
make_directories("data/structure/Results_stru_females/")
load_and_clumpp("data/structure/Results_stru_females/Run_files/")
load_and_K("data/structure/Results_stru_females/Run_files/")
ks_ad_fem <- load_and_K("data/structure/Results_stru_females/Run_files/")
save(ks_ad_fem, file = "data/structure/ks_ad_fem.R")
optimal_k(ks_ad_fem) 
both_k(ks_ad_fem) 


## chicks 
make_directories("data/structure/Results_stru_chicks/")
load_and_clumpp("data/structure/Results_stru_chicks/Run_files/")
load_and_K("data/structure/Results_stru_chicks/Run_files/")
ks_ad_chicks <- load_and_K("data/structure/Results_stru_chicks/Run_files/")
save(ks_ad_chicks, file = "data/structure/ks_ad_chick_2.R")
optimal_k(ks_ad_chicks) #8
both_k(ks_ad_chicks) #5, 8

#### Plotting barcharts ####

## Males: best log likelihood K = 9, highest delta K = 11
path_to_structure_out <- "data/structure/Results_stru_males/Run_files/"
path_to_struc_file <- "data/cleandata/Microsat.males.noLOCUS1+13.forstructure.stru"
all_files <- list.files(path_to_structure_out, pattern = "^results")
struc_out_paths <- paste0(path_to_structure_out, all_files)
slist <- readQ(files=struc_out_paths, filetype = "structure")
fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(slist,fn1))
ggthemr(palette = "solarized", layout = "clean",
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
struc_file <- path_to_struc_file
pops <- fread(struc_file) %>%
  select(V2) %>%
  mutate(V2 = ifelse(V2 == 1, "Koskenpää",
                     ifelse(V2 == 2, "Kummunsuo",
                            ifelse(V2 == 3, "Lauttasuo",
                                   ifelse(V2 == 4, "Lehtosuo", 
                                          ifelse(V2 == 5, "Nyrölä",
                                                 ifelse(V2 == 6, "Palosuo",
                                                        ifelse(V2 == 7, "Pihtissuo", 
                                                               ifelse(V2== 8, "Pirttilampi",
                                                                      ifelse(V2 == 9, "Pirttisuo",
                                                                             ifelse(V2 == 10, "Saarisuo",
                                                                                    ifelse(V2 == 11, "Teerisuo", "Utusuo"))))))))))))
pops <- data.frame(pops[seq(1, nrow(pops), 2)])
colnames(pops) <- "location"
pops$location <- as.character(pops$location)

plotQ(qlist=(slist)[c(1101, 201)], splab=spnames[c(1101, 201)], imgoutput = "join", grplab=pops,ordergrp = TRUE,sortind = "all",
      clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E",
                   "#56445D", "#1E3888", "#DDF9C1", "#772D8B", "#780116", "#E9FFF9"),
      outputfilename = "K9_K11_males", exportpath = "data/plots/", sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, titlelab = "a) Structure Barplot for males for K = 9 and K = 11",
      grplabangle = 25, panelspacer = 0.04, grplabsize = 2, splabsize = 6, showtitle = T, titlesize = 8, titlespacer = 2)

## females 
path_to_structure_out <- "data/structure/Results_stru_females/Run_files/"
path_to_struc_file <- "data/cleandata/Microsat.females.noLOCUS1+13.forstructure.stru"
all_files <- list.files(path_to_structure_out, pattern = "^results")
struc_out_paths <- paste0(path_to_structure_out, all_files)
slist <- readQ(files=struc_out_paths, filetype = "structure")
fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(slist,fn1))
ggthemr(palette = "solarized", layout = "clean",
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
struc_file <- path_to_struc_file
pops <- fread(struc_file) %>%
  select(V2) %>%
  mutate(V2 = ifelse(V2 == 1, "Koskenpää",
                     ifelse(V2 == 2, "Kummunsuo",
                            ifelse(V2 == 3, "Lauttasuo",
                                   ifelse(V2 == 4, "Lehtosuo", 
                                          ifelse(V2 == 5, "Nyrölä",
                                                 ifelse(V2 == 6, "Palosuo",
                                                        ifelse(V2 == 7, "Pihtissuo", 
                                                               ifelse(V2== 8, "Pirttilampi",
                                                                      ifelse(V2 == 9, "Pirttisuo",
                                                                             ifelse(V2 == 10, "Saarisuo",
                                                                                    ifelse(V2 == 11, "Teerisuo", "Utusuo"))))))))))))
pops <- data.frame(pops[seq(1, nrow(pops), 2)])
colnames(pops) <- "location"
pops$location <- as.character(pops$location)

## Females: best log likelihood K = 1, highest delta K = 2
plotQ(qlist=(slist)[c(1, 401)], splab=spnames[c(1, 401)], imgoutput = "join", grplab=pops,ordergrp = TRUE,
      #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
      outputfilename = "K1_K2", exportpath = getwd())

plotQ(qlist=(slist)[c(401)], splab=spnames[c(401)], imgoutput = "sep", grplab=pops,ordergrp = TRUE,
      clustercol=c("#1B9E77", "#D95F02"), sortind = "all",
      outputfilename = "K2_1_females", grplabangle = 25, panelspacer = 0.04, grplabsize = 2, exportpath = "data/plots/",
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, titlelab = "b) Structure Barplot for females for K = 2", 
      splabsize = 6, showtitle = T, titlesize = 8, titlespacer = 2)

#### chicks
#system("mkdir data/structure/Results_stru_chicks/Run_files/Pops-merged")
#system("cp data/structure/Results_stru_chicks/Run_files/outputclumpp/pop_K*/pop_*.txt data/structure/Results_stru_chicks/Run_files/Pops-merged")

path_to_structure_out <- "data/structure/Results_stru_chicks/Run_files/"
path_to_struc_file <- "data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru"
all_files <- list.files(path_to_structure_out, pattern = "^results")
struc_out_paths <- paste0(path_to_structure_out, all_files)
slist <- readQ(files=struc_out_paths, filetype = "structure", indlabfromfile=T)

clumpp_files<- list.files(paste0(path_to_structure_out,"Pops-merged/"))
clumpp_out_paths <- paste0(path_to_structure_out,"Pops-merged/", clumpp_files)
clist <- readQ(files=clumpp_out_paths, filetype = "clumpp")

fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(slist,fn1))

ggthemr(palette = "solarized", layout = "clean",
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()

struc_file <- path_to_struc_file
pops <- fread(struc_file) %>%
  select(V2) %>%
  mutate(V2 = ifelse(V2 == 1, "Koskenpää",
                     ifelse(V2 == 2, "Kummunsuo",
                            ifelse(V2 == 3, "Lauttasuo",
                                   ifelse(V2 == 4, "Lehtosuo", 
                                          ifelse(V2 == 5, "Nyrölä",
                                                 ifelse(V2 == 6, "Palosuo",
                                                        ifelse(V2 == 7, "Pihtissuo", 
                                                               ifelse(V2== 8, "Pirttilampi",
                                                                      ifelse(V2 == 9, "Pirttisuo",
                                                                             ifelse(V2 == 10, "Saarisuo",
                                                                                    ifelse(V2 == 11, "Teerisuo", "Utusuo"))))))))))))
pops <- data.frame(pops[seq(1, nrow(pops), 2)])
colnames(pops) <- "location"
pops$location <- as.character(pops$location)

## Chicks: best log likelihood K = 27, highest delta K = 17

plotQ(qlist=(slist)[c(81, 191)], splab=spnames[c(81, 191)], imgoutput = "join", grplab=pops,ordergrp = TRUE,
      sortind = "all", sharedindlab = FALSE,clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E",
                                                         "#56445D", "#1E3888", "#DDF9C1", "#772D8B", "#780116", "#E9FFF9",
                                                         "#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E",
                                                         "#56445D", "#1E3888", "#DDF9C1", "#772D8B", "#780116", "#E9FFF9",
                                                         "#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"), 
      outputfilename = "K17_27_chicks", grplabangle = 25, panelspacer = 0.04, grplabsize = 2, exportpath = "data/plots/",
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, titlelab = "c) Structure Barplot for chicks for K = 17 and K = 27",
      splabsize = 6, showtitle = T, titlesize = 8, titlespacer = 2)

