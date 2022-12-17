###### STRUCTURE ANALYSIS HUNTING MICROSATS ######

library(data.table); library(ParallelStructure)
library(dplyr); library(tidyr); library(stringr)
library(pophelper); library(ggthemr)

#### NB: the raw results from structure are not uploaded, but can be fully recreated
#### using the workflow below

##### Running structure for adults plus unrelated chicks #####
infile <- "/data/home/rchen/Hunting/blackgrouse-hunting/data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru"
outpath <- "/data/home/rchen/Hunting/blackgrouse-hunting/analyses/structure/results/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 12 #= number of sites

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- "1,2,3,4,5,6,7,8,9,10,11,12" #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "analyses/structure/structure_jobs.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='/data/home/rchen/Hunting/blackgrouse-hunting/analyses/structure/structure_jobs.txt', 
                                      n_cpu=45, infile=infile,
                                      outpath=outpath,
                                      numinds = nrow(all)/2,
                                      numloci=ncol(all)-2,
                                      noadmix = 0, alpha = 1.0, freqscorr=1, lambda = 1,
                                      printqhat=1, plot_output=0, onerowperind=0, locprior = 0)

# use admixture model (noadmix = 0), not using prior information regarding location (locprior =0) 
# not using prior population info (usepopinfo = 0)

#### Analysing output ####

source("analyses/R/STRUCTURE_functions.R") #contains functions to analyse structure output
make_directories("analyses/structure/results/")
load_and_clumpp("analyses/structure/results/Run_files/")
load_and_K("analyses/structure/results/Run_files/")
both_k(read.csv("analyses/structure/results/Run_files/output.csv")) # lnk = 4, deltak = 2

#### Plotting barcharts ####

## All: best log likelihood K = , highest delta K = 
path_to_structure_out <- "analyses/structure/results/"
path_to_struc_file <- "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru"
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

#choose K = 4 and K = 2
plotQ(qlist=(slist)[c(61, 41)], splab=spnames[c(61, 41)], imgoutput = "join", grplab=pops,ordergrp = TRUE,
      sortind = "all",
      clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3"),
      outputfilename = "K4_K2_all", exportpath = "figures", sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, 
      titlelab = "a) Structure Barplot for K = 4 and K = 2",
      grplabangle = 25, panelspacer = 0.04, grplabsize = 2, splabsize = 6, showtitle = T, titlesize = 8, titlespacer = 2)

