###### Here we redo all genetic analyses but taking only one chick per brood ######

library(readxl); library(tidyverse); library(data.table);library(hierfstat); 
library(plot.matrix); library(lme4); library(adegenet); library(forcats)
library(ape); library(cowplot); library(gridExtra)

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

### Running actual STRUCTURE analysis can be found in split 2b_STRUCTURE_huntedunhunted.R

### FST values:

#### Calculate Fst ####

## unrchicks 
#convert to hfstat object
unrchicks.stru <- read.structure("data/cleandata/Microsat.unrelated.chicks.noLOCUS1+13+14.forstructure.stru", n.ind = 200, n.loc = 11, onerowperind = F,
                              col.lab = 1, col.pop = 2, col.others = NULL,
                              row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                              ask = F, quiet = FALSE)

unrchicks.hfstat <- genind2hierfstat(unrchicks.stru)
#calculate stats
basicstat.unrchicks <- basic.stats(unrchicks.stru, diploid = TRUE, digits = 2) # Hobs, mean gene diversities, Fis and Fst

# per locus
fst.unrchicks.perlocus <- basicstat.unrchicks$perloc$Fst
fst.unrchicks.perlocus <- data.frame(Locus = seq(from = 1, to = 11), Fst = fst.unrchicks.perlocus)

# Pairwise Fst
fst.unrchicks <- pairwise.neifst(unrchicks.hfstat)
head(fst.unrchicks)

# Fst per population
boxplot(fst.unrchicks, col=funky(nPop(unrchicks.stru)), las=3,
        xlab="Population", ylab="Fst", main = "Pairwise Fst values per population only unrelated chicks") #pop 6 only has 1 sample

#Bootstrap
boot.fst.unrchicks <- boot.ppfst(unrchicks.hfstat, nboot = 1000) #confidence interval for Fst, upper triangle is upper limit, lower triangle lower limit

#create a long dataframe for pairwise Fst
boot.fst.unrchicks.UL <- boot.fst.unrchicks$ul
boot.fst.unrchicks.LL <- boot.fst.unrchicks$ll

flat.matrix <- function(d){
  data.frame(i=rep(row.names(d),ncol(d)),
             j=rep(colnames(d),each=nrow(d)),
             score=as.vector(d))
}

fst.unrchicks.flat <- flat.matrix(fst.unrchicks)
names(fst.unrchicks.flat) <- c("site.x", "site.y", "Fst")

boot.fst.unrchicks.LL.flat <- flat.matrix(boot.fst.unrchicks.LL)
names(boot.fst.unrchicks.LL.flat) <- c("site.x", "site.y", "LL")

boot.fst.unrchicks.UL.flat <- flat.matrix(boot.fst.unrchicks.UL)
names(boot.fst.unrchicks.UL.flat) <- c("site.x", "site.y", "UL")

pairwise.fst.unrchicks <- left_join(fst.unrchicks.flat, boot.fst.unrchicks.LL.flat, by = c("site.x", "site.y"))
pairwise.fst.unrchicks <- left_join(pairwise.fst.unrchicks, boot.fst.unrchicks.UL.flat, by = c("site.x", "site.y"))

pairwise.fst.unrchicks <- subset(pairwise.fst.unrchicks, site.x != "-9" & site.y != "-9")
pairwise.fst.unrchicks <- subset(pairwise.fst.unrchicks, !is.na(Fst) & !is.na(UL))
pairwise.fst.unrchicks <- pairwise.fst.unrchicks %>% mutate(Significance = case_when(
  UL > 0 & LL > 0 ~ "significant",
  UL > 0 & LL < 0 ~ "insignificant",
  UL < 0 & LL < 0 ~ "significant" ))

pairwise.fst.unrchicks <- pairwise.fst.unrchicks %>% mutate(Sig = case_when(
  UL > 0 & LL > 0 ~ "",
  UL > 0 & LL < 0 ~ "NS",
  UL < 0 & LL < 0 ~ "" ))

#### Make figure FST heatplot ####

# add abbreviation to data
sitenames <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")
# add abbreviations to sitenames
sitenames[1]
sitenames$abb <- c("KOS", "KUM", "LAU", "LEH", "NYR", "PAL", "PIH", "PIL", "PIS", "SAA", "TEE", "UTU")
str(sitenames)

pairwise.fst.unrchicks$site.x <- as.integer(pairwise.fst.unrchicks$site.x)
pairwise.fst.unrchicks$site.y <- as.integer(pairwise.fst.unrchicks$site.y)
pairwise.fst.unrchicks <- left_join(pairwise.fst.unrchicks, sitenames[,c(1,6,2)], by = c("site.x" = "pop_num"))
pairwise.fst.unrchicks <- left_join(pairwise.fst.unrchicks, sitenames[,c(1,6,2)], by = c("site.y" = "pop_num"))

### Unrelated chicks - Pairwise FST
mypalette3 <- c("#EDEDFD","#C2C1EC","#9795DB","#6C69C9","#413DB8","#1611A7") #purples from structure
theme_set(theme_classic())
unrchicks.fst <- ggplot(pairwise.fst.unrchicks, aes(abb.x, abb.y, fill = Fst)) + geom_tile() + theme_classic() + 
  scale_fill_gradientn(colors = mypalette3, limits = c(0,0.05)) + 
  geom_text(aes(label = Sig), size = 7)+ 
  theme(text = element_text(family = "Helvetica", size = 22),
        axis.text.x = element_text(angle = 90,
                                   face = c("bold", "plain", "plain", 
                                            "plain", "plain", "plain")), 
        axis.text.y = element_text(face = c("plain", "plain", "plain", "plain",
                                            "plain", "bold")),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(),
        legend.box.spacing = unit(0.5, 'cm'),
        legend.position = c(0.9, 0.36)) +
  ggtitle("(a)  ")

unrchicks.fst

##### Running structure for unrelated.chicks #####
### total
unrelated.chicks <- fread("data/cleandata/Microsat.unrelated.chicks.noLOCUS1+13+14.forstructure.stru") #located in /data in github directory
infile_unrelated.chicks <- "data/cleandata/Microsat.unrelated.chicks.noLOCUS1+13+14.forstructure.stru"
system("mkdir data/structure/Results_stru_unrelated.chicks")
outpath_unrelated.chicks <- "data/structure/Results_stru_unrelated.chicks/"

# job matrix and write to job file
nrep <- 10
burnin <- 10000
niter <- 10000
up_to_k <- 7

# job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# make the job matrix
pop <- ("1,2,4,5,10,11,12") #number of pops in the file

hunt_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(hunt_jobs), ncol = length(hunt_jobs[1,]), file = "data/structure/hunt_jobs_unrelated.chicks.txt")

# file path to structure

STR_path='/usr/local/bin/'

# Run Parallel Structure

# system("mkdir Results_adults")

# Run structure (from terminal, do not run this last part in Rstudio)


ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist='data/structure/hunt_jobs_unrelated.chicks.txt', 
                                      n_cpu=45, infile=infile_unrelated.chicks,
                                      outpath=outpath_unrelated.chicks,numinds = nrow(unrelated.chicks)/2,
                                      numloci=ncol(unrelated.chicks)-2,noadmix = 0, alpha = 1.0,freqscorr=1,lambda = 1,
                                      printqhat=1,plot_output=0,onerowperind=0, locprior = 0)


#### Analysing results STRUCTURE unrelated.chicks ####
source("scripts/STRUCTURE_functions.R")
make_directories("data/structure/Results_stru_unrelated.chicks/")
load_and_clumpp("data/structure/Results_stru_unrelated.chicks/Run_files/")
load_and_K("data/structure/Results_stru_unrelated.chicks/Run_files/")
ks_ad_unrchick <- load_and_K("data/structure/Results_stru_unrelated.chicks/Run_files/")
save(ks_ad_unrchick, file = "data/structure/ks_ad_unrchick.R")
optimal_k(ks_ad_unrchick) #
both_k(ks_ad_unrchick) #

#### Creating STRUCTURE plots #####
load("data/structure/ks_ad_unrchick.R")

unrchicks.ks.lnpr <- ggplot(ks_ad_unrchick, aes(x = k, y = elpdmean)) +
  geom_point() + geom_line()+ ylab("LnPr(X|k)") +
  geom_errorbar(aes(ymin = elpdmin, ymax = elpdmax), width = 0.1, col = "gray45") +
  scale_y_continuous(breaks = c(-4000000, -2000000, 0),
                     labels = c("-4M", "-2M", "0"))+
  theme_classic() +
  theme(text = element_text(size = 22, family = "Helvetica"),
        axis.title = element_text(),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(),
        axis.text = element_text(),
        plot.title = element_text())+ 
  labs(title = "(d)")

unrchicks.ks.delta <- ggplot(ks_ad_unrchick, aes(x = k, y = deltaK  )) +
  geom_point() + geom_line()+ ylab(expression(Delta*"k")) +
  theme_classic() +
  theme(text = element_text(size = 22, family = "Helvetica"),
        axis.title = element_text(),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(),
        axis.text = element_text(),
        plot.title = element_text())+ 
  labs(title = "(e)")

#### Correlogram ####
spatial.additional <- read_excel("data/tables/SpatialAutocor_Supplementary.xlsx", sheet = "ForR")


### spatial

unr.spatial <- spatial.additional %>% filter(Data  == "Unrelated_chicks" & Hunt == "All") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + 
  geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, 
                size=1, color="black", stat = "identity")+
  ylim(-0.05, 0.040) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+ labs(title = "(b)")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Helvetica", size = 22),
        plot.title = element_text(),
        axis.text.x = element_text(margin = margin(b = 10)), 
        axis.text.y = element_text( margin = margin(l = 10)),
        axis.title.x = element_text(), 
        axis.title.y = element_text()) 


### Make complete figure

png("data/figures/K_unrelatedchicks.png", 
    width = 1000, height = 700, units = "px")
gA <- ggplotGrob(unrchicks.fst)
gB <- ggplotGrob(unr.spatial)
gD <- ggplotGrob(unrchicks.ks.lnpr)
gE <- ggplotGrob(unrchicks.ks.delta)
grid::grid.newpage()
grid::grid.draw(cbind(rbind(gA, gD), rbind(gB, gE)))

dev.off()

ggsave(plot = unrchicks.fst, "data/figures/UnrChicks.FST.png",
       width = 30, height = 20, units = c("cm"))
