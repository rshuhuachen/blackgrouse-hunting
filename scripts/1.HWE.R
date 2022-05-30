#________________________________________________________________________
####### POPULATION STRUCTURE IN HUNTED VS UNHUNTED LEKKING SITES #######
#________________________________________________________________________

##### In this script, we will import the structure data files for adults and chicks 
##### separately and plot some basic distribution and statistics on the microsatellite data,        
##### plot PCA's, test for HWE 

### Load packages ###
library(dplyr);library(tibble)
library(adegenet); library(pegas); library(data.table)

adults <- read.structure("data/rawdata/Microsat.adults.forstructure.stru", n.ind = 1878, n.loc = 14, onerowperind = F,
                         col.lab = 1, col.pop = 2, col.others = NULL,
                         row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                         ask = F, quiet = FALSE) #just need this for HWE


#### Summary ####
summary(adults)

#### Testing for Hardy-Weinberg equilibrium ####

# First all together
adultHWE.all <- pegas::hw.test(adults, B = 1000) #B = 1000 for 1000 Monte Carlo permutations 
adultHWE.all 

# Then per population using a for loop
adultpop <- seppop(adults)
# Run loop
adultHWE = NULL
for(i in 1:length(adultpop)) {
  hwt <- pegas::hw.test(adultpop[[i]], B=1000)
  smry <- summary(adultpop[[i]])
  
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  pexact <- hwt[,4] #hw.test does chi2 test and exact test. We use p-values of exact test which are given in 4th col
  qval.FDR <- p.adjust(pexact, method = "fdr")
  qval.bon <- p.adjust(pexact, method = "bonferroni")
  adultHWE <- as.data.frame(cbind(adultHWE, Hobs, Hexp, pexact, qval.FDR, qval.bon))
  
}

sites<-rep(names(adultpop[1:length(adultpop)]),each=5)
adultHWE <- rbind(adultHWE, sites)
adultHWE <- adultHWE[c(nrow(adultHWE),1:(nrow(adultHWE)-1)),]
rownames(adultHWE)[1] <- "Site"

adultHWE.t <- as.data.frame(t(adultHWE))
nums <- c(2:15)
adultHWE.t[nums] <- lapply(adultHWE.t[nums], as.numeric)

adultHWE.t[,c(2:15)]<- round(adultHWE.t[,c(2:15)], 2)
head(adultHWE.t)

#exclude locus 1, and 13 as FDR-corrected q-values indicate violation
#of assumption of HWE in over 70% of the sites
