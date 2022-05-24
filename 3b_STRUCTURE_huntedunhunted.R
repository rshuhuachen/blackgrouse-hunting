#STRUCTURE WITH HUNTED/UNHUNTED DATA

library(data.table); library(ParallelStructure)
library(dplyr); library(tidyr); library(stringr)
library(pophelper)

#### NB: the raw results from structre are not uploaded, but can be fully recreated
#### using the workflow below

##### Running structure for males #####
### hunted
males_hunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.males.hunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_males_hunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.males.hunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_males_hunted")
outpath_males_hunted <- "data/structure/Results_stru_males_hunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("1,3,6,7,9,12") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_males_hunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_males_hunted.txt', 
                                      n_cpu=45, infile=infile_males_hunted,
                                      outpath=outpath_males_hunted,numinds = nrow(males_hunted)/2,
                                      numloci=ncol(males_hunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)


### unhunted
males_unhunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.males.unhunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_males_unhunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.males.unhunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_males_unhunted")
outpath_males_unhunted <- "data/structure/Results_stru_males_unhunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("2,4,5,8,10,11") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_males_unhunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_males_unhunted.txt', 
                                      n_cpu=45, infile=infile_males_unhunted,
                                      outpath=outpath_males_unhunted,numinds = nrow(males_unhunted)/2,
                                      numloci=ncol(males_unhunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)

##### Running structure for females #####
### hunted
females_hunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.females.hunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_females_hunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.females.hunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_females_hunted")
outpath_females_hunted <- "data/structure/Results_stru_females_hunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("1,3,6,7,9,12") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_females_hunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_females_hunted.txt', 
                                      n_cpu=45, infile=infile_females_hunted,
                                      outpath=outpath_females_hunted,numinds = nrow(females_hunted)/2,
                                      numloci=ncol(females_hunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)


### unhunted
females_unhunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.females.unhunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_females_unhunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.females.unhunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_females_unhunted")
outpath_females_unhunted <- "data/structure/Results_stru_females_unhunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("2,4,5,8,10,11") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_females_unhunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_females_unhunted.txt', 
                                      n_cpu=45, infile=infile_females_unhunted,
                                      outpath=outpath_females_unhunted,numinds = nrow(females_unhunted)/2,
                                      numloci=ncol(females_unhunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)


##### Running structure for chicks #####
### hunted
chicks_hunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.chicks.hunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_chicks_hunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.chicks.hunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_chicks_hunted")
outpath_chicks_hunted <- "data/structure/Results_stru_chicks_hunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 200

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("1,12") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_chicks_hunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_chicks_hunted.txt', 
                                      n_cpu=45, infile=infile_chicks_hunted,
                                      outpath=outpath_chicks_hunted,numinds = nrow(chicks_hunted)/2,
                                      numloci=ncol(chicks_hunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)


### unhunted
chicks_unhunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.chicks.unhunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_chicks_unhunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.chicks.unhunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_chicks_unhunted")
outpath_chicks_unhunted <- "data/structure/Results_stru_chicks_unhunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 200

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("2,4,5,10,11") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_chicks_unhunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_chicks_unhunted.txt', 
                                      n_cpu=45, infile=infile_chicks_unhunted,
                                      outpath=outpath_chicks_unhunted,numinds = nrow(chicks_unhunted)/2,
                                      numloci=ncol(chicks_unhunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)

##### Running structure for unrelated.chicks #####
### hunted
unrelated.chicks_hunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.unrelated.chicks.hunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_unrelated.chicks_hunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.unrelated.chicks.hunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_unrelated.chicks_hunted")
outpath_unrelated.chicks_hunted <- "data/structure/Results_stru_unrelated.chicks_hunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- () #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_unrelated.chicks_hunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_unrelated.chicks_hunted.txt', 
                                      n_cpu=45, infile=infile_unrelated.chicks_hunted,
                                      outpath=outpath_unrelated.chicks_hunted,numinds = nrow(unrelated.chicks_hunted)/2,
                                      numloci=ncol(unrelated.chicks_hunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)


### unhunted
unrelated.chicks_unhunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.unrelated.chicks.unhunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory
infile_unrelated.chicks_unhunted <- "data/cleandata/SplitHuntedUnhunted/Microsat.unrelated.chicks.unhunted.noLOCUS1+13.forstructure.stru"
system("mkdir data/structure/Results_stru_unrelated.chicks_unhunted")
outpath_unrelated.chicks_unhunted <- "data/structure/Results_stru_unrelated.chicks_unhunted/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- () #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_unrelated.chicks_unhunted.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_unrelated.chicks_unhunted.txt', 
                                      n_cpu=45, infile=infile_unrelated.chicks_unhunted,
                                      outpath=outpath_unrelated.chicks_unhunted,numinds = nrow(unrelated.chicks_unhunted)/2,
                                      numloci=ncol(unrelated.chicks_unhunted)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)



###### ANALYSING RESULTS ######
source("3c_STRUCTURE_functions.R")

## males hunted
make_directories("data/structure/Results_stru_males_hunted/")
load_and_clumpp("data/structure/Results_stru_males_hunted/Run_files/")
load_and_K("data/structure/Results_stru_males_hunted/Run_files/")
ks_ad_male_hunted <- load_and_K("data/structure/Results_stru_males_hunted/Run_files/")
save(ks_ad_male_hunted, file = "data/structure/ks_ad_male_hunted.R")
optimal_k(ks_ad_male_hunted) #8
both_k(ks_ad_male_hunted) #5, 8

## males unhunted
make_directories("data/structure/Results_stru_males_unhunted/")
load_and_clumpp("data/structure/Results_stru_males_unhunted/Run_files/")
load_and_K("data/structure/Results_stru_males_unhunted/Run_files/")
ks_ad_male_unhunted <- load_and_K("data/structure/Results_stru_males_unhunted/Run_files/")
save(ks_ad_male_unhunted, file = "data/structure/ks_ad_male_unhunted.R")
optimal_k(ks_ad_male_unhunted) #8
both_k(ks_ad_male_unhunted) #5, 8

## females hunted
make_directories("data/structure/Results_stru_females_hunted/")
load_and_clumpp("data/structure/Results_stru_females_hunted/Run_files/")
load_and_K("data/structure/Results_stru_females_hunted/Run_files/")
ks_ad_male_hunted <- load_and_K("data/structure/Results_stru_females_hunted/Run_files/")
save(ks_ad_male_hunted, file = "data/structure/ks_ad_male_hunted.R")
optimal_k(ks_ad_male_hunted) #8
both_k(ks_ad_male_hunted) #5, 8

## females unhunted
make_directories("data/structure/Results_stru_females_unhunted/")
load_and_clumpp("data/structure/Results_stru_females_unhunted/Run_files/")
load_and_K("data/structure/Results_stru_females_unhunted/Run_files/")
ks_ad_male_unhunted <- load_and_K("data/structure/Results_stru_females_unhunted/Run_files/")
save(ks_ad_male_unhunted, file = "data/structure/ks_ad_male_unhunted.R")
optimal_k(ks_ad_male_unhunted) #8
both_k(ks_ad_male_unhunted) #5, 8

## chicks hunted
make_directories("data/structure/Results_stru_chicks_hunted/")
load_and_clumpp("data/structure/Results_stru_chicks_hunted/Run_files/")
load_and_K("data/structure/Results_stru_chicks_hunted/Run_files/")
ks_ad_male_hunted <- load_and_K("data/structure/Results_stru_chicks_hunted/Run_files/")
save(ks_ad_male_hunted, file = "data/structure/ks_ad_male_hunted.R")
optimal_k(ks_ad_male_hunted) #8
both_k(ks_ad_male_hunted) #5, 8

## chicks unhunted
make_directories("data/structure/Results_stru_chicks_unhunted/")
load_and_clumpp("data/structure/Results_stru_chicks_unhunted/Run_files/")
load_and_K("data/structure/Results_stru_chicks_unhunted/Run_files/")
ks_ad_male_unhunted <- load_and_K("data/structure/Results_stru_chicks_unhunted/Run_files/")
save(ks_ad_male_unhunted, file = "data/structure/ks_ad_male_unhunted.R")
optimal_k(ks_ad_male_unhunted) #8
both_k(ks_ad_male_unhunted) #5, 8

## unrelated.chicks hunted
make_directories("data/structure/Results_stru_unrelated.chicks_hunted/")
load_and_clumpp("data/structure/Results_stru_unrelated.chicks_hunted/Run_files/")
load_and_K("data/structure/Results_stru_unrelated.chicks_hunted/Run_files/")
ks_ad_male_hunted <- load_and_K("data/structure/Results_stru_unrelated.chicks_hunted/Run_files/")
save(ks_ad_male_hunted, file = "data/structure/ks_ad_male_hunted.R")
optimal_k(ks_ad_male_hunted) #8
both_k(ks_ad_male_hunted) #5, 8

## unrelated.chicks unhunted
make_directories("data/structure/Results_stru_unrelated.chicks_unhunted/")
load_and_clumpp("data/structure/Results_stru_unrelated.chicks_unhunted/Run_files/")
load_and_K("data/structure/Results_stru_unrelated.chicks_unhunted/Run_files/")
ks_ad_male_unhunted <- load_and_K("data/structure/Results_stru_unrelated.chicks_unhunted/Run_files/")
save(ks_ad_male_unhunted, file = "data/structure/ks_ad_male_unhunted.R")
optimal_k(ks_ad_male_unhunted) #8
both_k(ks_ad_male_unhunted) #5, 8





