#STRUCTURE WITH HUNTED/UNHUNTED DATA

library(data.table); library(ParallelStructure)
library(dplyr); library(tidyr); library(stringr)
library(pophelper)

#### NB: the raw results from structre are not uploaded, but can be fully recreated
#### using the workflow below

##### Running structure for males #####

males_hunted <- fread("data/cleandata/SplitHuntedUnhunted/Microsat.males.hunted.noLOCUS1+13.forstructure.stru") #located in /data in github directory

infile <- "data/cleandata/SplitHuntedUnhunted/Microsat.males.hunted.noLOCUS1+13.forstructure.stru"
system("mkdir C:/Users/rchen2/Documents/Black Grouse PhD/GitHub/blackgrouse-hunting/data/structure/Results_stru_males_hunted")
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
