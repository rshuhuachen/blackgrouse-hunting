############# 6. Creating plots ##############

## In this script, you will find the code used to create the figures of the manuscript

### load packages ###

library(RColorBrewer); library(ggplot2); library(extrafont); 
library(dplyr); library(readxl); library(devtools); library(pophelper)
library(gridExtra); library(adegenet)

#palette
mypalette3 <- c("#EDEDFD","#C2C1EC","#9795DB","#6C69C9","#413DB8","#1611A7") #purples from structure
theme_set(theme_classic())
#### Figure 2: heatmaps pairwise Fst + correlograms spatial autocorrelation ####

### First: heatmap
pairwise.fst.males <- read.csv("data/Pairwise_Fst_males.csv")
pairwise.fst.females <- read.csv("data/Pairwise_Fst_females.csv")
pairwise.fst.chicks <- read.csv("data/Pairwise_Fst_chicks.csv")

#change the negative values to 0 in females
pairwise.fst.females$Fst[which (pairwise.fst.females$Fst < 0)] <- 0

# add abbreviation to chick data
sitenames <- read.csv("data/Codes.pops.both.filtered_withcoord.csv")
# add abbreviations to sitenames
sitenames[1]
sitenames$abb <- c("KOS", "KUM", "LAU", "LEH", "NYR", "PAL", "PIH", "PIL", "PIS", "SAA", "TEE", "UTU")
str(sitenames)

pairwise.fst.chicks <- left_join(pairwise.fst.chicks, sitenames[,c(1,6,2)], by = c("site.x" = "pop_num"))
pairwise.fst.chicks <- left_join(pairwise.fst.chicks, sitenames[,c(1,6,2)], by = c("site.y" = "pop_num"))

### Males - pairwise Fst heatmap
ggplot(pairwise.fst.males, aes(abb.x, abb.y, fill = Fst)) + geom_tile() + theme_classic() + 
  scale_fill_gradientn(colors = mypalette3, limits = c(0,0.05)) + 
  geom_text(aes(label = Sig), size = 8)+
  theme(text = element_text(family = "Arial", size = 22),
        axis.text.x = element_text(angle = 90, size = 22,
                                   face = c("bold", "plain", "bold", "plain", "plain",
                                            "bold", "bold", "plain", "bold",
                                            "plain", "plain")), 
        axis.text.y = element_text(size = 22, 
                                   face = c("plain", "bold", "plain", "plain",
                                            "bold", "bold", "plain", "bold",
                                            "plain", "plain", "bold")),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size = 38),
        legend.position = c(0.8, 0.3)) +
  ggtitle('(a) Males') 


### Females - pairwise Fst heatmap
ggplot(pairwise.fst.females, aes(abb.x, abb.y, fill = Fst)) + geom_tile() + theme_classic() + 
  scale_fill_gradientn(colors = mypalette3, limits = c(0,0.05)) + 
  geom_text(aes(label = Sig), size = 8)+ 
  theme(text = element_text(family = "Arial", size = 22),
        axis.text.x = element_text(angle = 90, size = 22,
                                   face = c("bold", "plain", "bold", "plain", "plain",
                                            "bold", "bold", "plain", "bold",
                                            "plain", "plain")), 
        axis.text.y = element_text(size = 22, 
                                   face = c("plain", "bold", "plain", "plain",
                                            "bold", "bold", "plain", "bold",
                                            "plain", "plain", "bold")),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size = 38),
        legend.position = c(0.8, 0.3)) +
  ggtitle('(b) Females') 

### Chicks - Pairwise FST
chicks.fst <- ggplot(pairwise.fst.chicks, aes(abb.x, abb.y, fill = Fst)) + geom_tile() + theme_classic() + 
  scale_fill_gradientn(colors = mypalette3, limits = c(0,0.05)) + 
  geom_text(aes(label = Sig), size = 8)+ 
  theme(text = element_text(family = "Helvetica", size = 22),
        axis.text.x = element_text(angle = 90, size = 22,
                                   face = c("bold", "plain", "plain", 
                                            "plain", "plain", "plain")), 
        axis.text.y = element_text(size = 20, 
                                   face = c("plain", "plain", "plain", "plain",
                                            "plain", "bold")),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size = 38),
        legend.position = c(0.8, 0.3)) +
  ggtitle("(c) Chicks ")

### Second - correlograms ###
spatial <- read_excel("data/SpatialAutocor_4.4.22.xlsx", sheet = "ForR")


### Males - spatial
spatial %>% filter(Who == "Male") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + 
  geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, 
                size=1.5, color="black", stat = "identity")+
  ylim(-0.03, 0.035) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_text(size = 38),
        axis.text.x = element_text(size = 26, margin = margin(b = 10)), 
        axis.text.y = element_text(size = 26, margin = margin(l = 10)),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30)) 

### Females - spatial

spatial %>% filter(Who == "Females") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, size=1.5, color="black", stat = "identity")+
  ylim(-0.03, 0.035) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Helvetica", size = 14),
        plot.title = element_text(size = 38),
        axis.text.x = element_text(size = 26, margin = margin(b = 10)), 
        axis.text.y = element_text(size = 26, margin = margin(l = 10)),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30))


### Chicks - spatial
spatial %>% filter(Who == "Chicks") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, size=1.5, color="black", stat = "identity")+
  ylim(-0.03, 0.035) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_text(size = 38),
        axis.text.x = element_text(size = 26, margin = margin(b = 10)), 
        axis.text.y = element_text(size = 26, margin = margin(l = 10)),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30)) 

### Figure 3: boxplots migration rates ####
male_run5_clean <- read.csv("data/run5_males_clean.csv")
female_run5_clean <- read.csv("data/run5_females_clean.csv")

names(female_run5_clean) == names(male_run5_clean)
female_run5_clean$sex <- "Female"
male_run5_clean$sex <- "Male"

#combine in one df
migration_both <- rbind(male_run5_clean, female_run5_clean)

#change levels
migration_both$hunt_in <- relevel(as.factor(migration_both$hunt_in), ref = "unhunted")
migration_both$hunt_out <- relevel(as.factor(migration_both$hunt_out), ref = "unhunted")

#first, exclude the 'non-migration rates' which are those where pop in = pop out
migration_both <- subset(migration_both, m_in != m_out)

# plot emigration rates 
emigration.plot <- ggplot(migration_both, aes(x = hunt_out, y = migration_ESSc, fill = sex)) + 
  geom_boxplot(outlier.shape = NA, aes(middle = mean(migration_ESSc))) + ylim(0, 0.03)+
  labs(title = "(a) Emigration rates") + 
  ylab("Migration rate")+
  theme(text = element_text(family = "Arial", size = 26),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size = 38),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_fill_manual(values = c("Male" = "#57939a",
                               "Female" = "#be4d5a"),
                    labels=c("Male", "Female"))+
  guides(fill = guide_legend("Sex")) 

#plot immigration rates

immmigration.plot <- ggplot(migration_both, aes(x = hunt_in, y = migration_ESSc, fill = sex)) + 
  geom_boxplot(outlier.shape = NA, aes(middle = mean(migration_ESSc))) + ylim(0, 0.03)+
  labs(title = "(b) Immigration rates") + 
  ylab("Migration rate")+
  theme(text = element_text(family = "Arial", size = 26),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size = 38),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_fill_manual(values = c("Male" = "#57939a",
                               "Female" = "#be4d5a"),
                    labels=c("Male", "Female"))+
  guides(fill = guide_legend("Sex"))


ggsave(plot = emigration.plot, "data/Emigration.png",
       width = 30, height = 20, units = c("cm"))

ggsave(plot = immmigration.plot, "data/Immigration.png",
       width = 30, height = 20, units = c("cm"))

#### Supplementary Figure 1: log likelihood per K (STRUCTURE) ##### 

load("data/ks_ad_fem.R")
load("data/ks_chick_2.R")
load("data/ks_ad_male.R")

#### plotting K vs mean log likelihood

ks_fem <- ggplot(ks_ad_fem, aes(x = k, y = elpdmean)) +
  geom_point() + geom_line()+ ylab("LnPr(X|k)") +
  scale_x_continuous(name = "k", 
                     breaks = c(1:12),
                     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")) +
  scale_y_continuous(breaks = c(-40000, -36000, -32000),
                     labels = c("-40K", "-36K", "-32K"))+
  geom_errorbar(aes(ymin = elpdmin, ymax = elpdmax), width = 0.1, col = "gray45") +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 22),
        plot.title = element_text(size = 22)) + 
  labs(title = "(b)") 

ks_male <- ggplot(ks_ad_male, aes(x = k, y = elpdmean)) +
  geom_point() + geom_line()+ ylab("LnPr(X|k)") +
  scale_x_continuous(name = "k", 
                     breaks = c(1:12),
                     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")) +
  geom_errorbar(aes(ymin = elpdmin, ymax = elpdmax), width = 0.1, col = "gray45") +
  scale_y_continuous(breaks = c(-43000, -41000, -39000),
                     labels = c("-43K", "-41K", "-39K"))+
  theme_classic()+
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 22),
        plot.title = element_text(size = 22))+ 
  labs(title = "(a)")

ks_chick <- ggplot(ks_chick, aes(x = k, y = elpdmean)) +
  geom_point() + geom_line()+ ylab("LnPr(X|k)") +
  geom_errorbar(aes(ymin = elpdmin, ymax = elpdmax), width = 0.1, col = "gray45") +
  scale_y_continuous(breaks = c(-4000000, -2000000, 0),
                     labels = c("-4M", "-2M", "0"))+
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 22),
        plot.title = element_text(size = 22))+ 
  labs(title = "(c)")

library(cowplot); library(gridExtra)
png("data/K_all.png", 
    width = 1000, height = 700, units = "px")
grid.arrange(ks_male, ks_fem, ks_chick)
dev.off()

##### Supplementary Figures 2 and 3: barplots STRUCTURE #####

##### First plot the bargraphs for 10 different K's an example of 10 different runs

## males 
males_path_to_structure_out <- "data/Results_stru/males/Run_files/"
males_path_to_struc_file <- "data/Microsat.males.noLOCUS1+13.forstructure.stru"
males_all_files <- list.files(males_path_to_structure_out, pattern = "^results")
males_struc_out_paths <- paste0(males_path_to_structure_out, males_all_files)
males_slist <- readQ(files=males_struc_out_paths, filetype = "structure")
fn1 <- function(x) attr(x,"k")
males_spnames <- paste0("K=",sapply(males_slist,fn1))
ggthemr(palette = "solarized", layout = "clean",
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
males_struc_file <- males_path_to_struc_file
males_pops <- fread(males_struc_file) %>%
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
males_pops <- data.frame(males_pops[seq(1, nrow(males_pops), 2)])
colnames(males_pops) <- "location"
males_pops$location <- as.character(males_pops$location)

## Males: best log likelihood and also highest delta K = 4

## just highest log likelihood K = 4 (run 601) and split up per group
plotQ(qlist=(males_slist)[c(601)],imgoutput = "sep", grplab=males_pops,ordergrp = TRUE,sortind = "all",
      clustercol=palette1,
      subsetgrp = c("Koskenpää","Kummunsuo","Lauttasuo","Lehtosuo","Nyrölä"),
      outputfilename = "K11_males_simplecolours_a", exportpath = getwd(), 
      sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = F, 
      grplabangle = 25, panelspacer = 0.02, grplabsize = 6, showtitle = F, titlesize = 14, titlespacer = 2,
      height = 8, width = 20, splab = "", splabcol = "white", grplabpos = 0.75, linepos = 1)

plotQ(qlist=(males_slist)[c(601)],imgoutput = "sep", grplab=males_pops,ordergrp = TRUE,sortind = "all",
      clustercol=palette1,
      subsetgrp = c("Palosuo", "Pihtissuo", "Pirttilampi", "Pirttisuo", "Saarisuo", "Teerisuo", "Utusuo"),
      outputfilename = "K11_males_simplecolours_b", exportpath = getwd(), 
      sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = F, 
      grplabangle = 25, panelspacer = 0.02, grplabsize = 6, showtitle = F, titlesize = 14, titlespacer = 2,
      height = 8, width = 20, splab = "", splabcol = "white", grplabpos = 0.75, linepos = 1)


## females 
females_path_to_structure_out <- "data/Results_stru_females/Run_files/"
females_path_to_struc_file <- "data/Microsat.females.noLOCUS1+13.forstructure.stru"
females_all_files <- list.files(females_path_to_structure_out, pattern = "^results")
females_struc_out_paths <- paste0(females_path_to_structure_out, females_all_files)
females_slist <- readQ(files=females_struc_out_paths, filetype = "structure")
females_spnames <- paste0("K=",sapply(females_slist,fn1))
females_struc_file <- females_path_to_struc_file
females_pops <- fread(females_struc_file) %>%
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
females_pops <- data.frame(females_pops[seq(1, nrow(females_pops), 2)])
colnames(females_pops) <- "location"
females_pops$location <- as.character(females_pops$location)

#system("mkdir /data/home/rchen/Hunting/Final_analysis/Results_females_11.1/Run_files/Plots/")
#highest delta k = 2 (run 401)
plotQ(qlist=(females_slist)[c(401)], splab="", imgoutput = "sep", grplab=females_pops,ordergrp = TRUE,
      clustercol=c("#BE4D5A", "#8B95C9"), 
      subsetgrp = c("Koskenpää","Kummunsuo","Lauttasuo","Lehtosuo","Nyrölä"),sortind = "all",
      outputfilename = "K2_females_simplecolours_a", grplabangle = 25, panelspacer = 0.04, grplabsize = 6, 
      exportpath = getwd(),
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4,  
      splabsize = 6, showtitle = F, 
      height = 8, width = 20,  splabcol = "white", grplabpos = 0.75, linepos = 1)

plotQ(qlist=(females_slist)[c(401)], splab="", imgoutput = "sep", grplab=females_pops,ordergrp = TRUE,
      clustercol=c("#BE4D5A", "#8B95C9"), 
      subsetgrp = c("Palosuo", "Pihtissuo", "Pirttilampi", "Pirttisuo", "Saarisuo", "Teerisuo", "Utusuo"),sortind = "all",
      outputfilename = "K2_females_simplecolours_b", grplabangle = 25, panelspacer = 0.04, grplabsize = 6, 
      exportpath = getwd(),
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4,  
      splabsize = 6, showtitle = F, 
      height = 8, width = 20,  splabcol = "white", grplabpos = 0.75, linepos = 1)


#### chicks

chicks_path_to_structure_out <- "data/Results_stru_chicks/Run_files/"
chicks_path_to_struc_file <- "data/Microsat.chicks.noLOCUS1+13+14.forstructure.stru"
chicks_all_files <- list.files(chicks_path_to_structure_out, pattern = "^results")
chicks_struc_out_paths <- paste0(chicks_path_to_structure_out, chicks_all_files)
chicks_slist <- readQ(files=chicks_struc_out_paths, filetype = "structure", indlabfromfile=T)

chicks_spnames <- paste0("K=",sapply(chicks_slist,fn1))

chicks_struc_file <- chicks_path_to_struc_file
chicks_pops <- fread(chicks_struc_file) %>%
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
chicks_pops <- data.frame(chicks_pops[seq(1, nrow(chicks_pops), 2)])
colnames(chicks_pops) <- "location"
chicks_pops$location <- as.character(chicks_pops$location)

palette1.chicks = c("#be4d5a","#474b4e","#8b95c9","#f7f6f0","#57939a","#d7d6d5","#4d9de0","#785964","#918ef4","#519e8a",
                    "#df2a3f","#224f71","#7184e5","#fbf8ea","#2ab6c5","#ebd6c1","#3d9ff0","#a62c56","#8e8afa","#28c89d",
                    "#ff0a27","#005494","#5773ff","#fffae5","#00d8f0","#ffd6ad","#2ea1ff","#d10049","#8985ff","#00f0b0",
                    "#be4d5a","#474b4e","#8b95c9","#f7f6f0","#57939a","#d7d6d5","#4d9de0","#785964","#918ef4","#519e8a",
                    "#df2a3f","#224f71","#7184e5","#fbf8ea","#2ab6c5","#ebd6c1","#3d9ff0","#a62c56","#8e8afa","#28c89d",
                    "#ff0a27","#005494","#5773ff","#fffae5","#00d8f0","#ffd6ad","#2ea1ff","#d10049","#8985ff","#00f0b0",
                    "#be4d5a","#474b4e","#8b95c9","#f7f6f0","#57939a","#d7d6d5","#4d9de0","#785964","#918ef4","#519e8a",
                    "#df2a3f","#224f71","#7184e5","#fbf8ea","#2ab6c5","#ebd6c1","#3d9ff0","#a62c56","#8e8afa","#28c89d",
                    "#ff0a27","#005494","#5773ff","#fffae5","#00d8f0","#ffd6ad","#2ea1ff","#d10049","#8985ff","#00f0b0",
                    "#be4d5a","#474b4e","#8b95c9","#f7f6f0","#57939a","#d7d6d5","#4d9de0","#785964","#918ef4","#519e8a")
#repeating colours

## Chicks: highest delta K = 20 only (run 1122)
plotQ(qlist=(chicks_slist)[c(1122)], imgoutput = "sep", grplab=palette1.chicks,ordergrp = TRUE, showsp = FALSE, showtitle = F,
      sortind = "all", sharedindlab = FALSE,clustercol=palette1, 
      subsetgrp = c("Koskenpää","Kummunsuo","Lehtosuo"),
      outputfilename = "K17_chicks_simplecolours_a", grplabangle = 25, panelspacer = 0.02, grplabsize = 6, 
      exportpath = getwd(),
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, 
      height = 8, width = 20,  splabcol = "white", grplabpos = 0.75, linepos = 1)

plotQ(qlist=(chicks_slist)[c(1122)], imgoutput = "sep", grplab=palette1.chicks,ordergrp = TRUE, showsp = FALSE, showtitle = F,
      sortind = "all", sharedindlab = FALSE,clustercol=palette1, 
      subsetgrp = c("Nyrölä","Saarisuo","Teerisuo", "Utusuo"),
      outputfilename = "K17_chicks_simplecolours_b", grplabangle = 25, panelspacer = 0.02, grplabsize = 6, 
      exportpath = getwd(),
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, 
      height = 8, width = 20,  splabcol = "white", grplabpos = 0.75, linepos = 1)

## Chicks: best log likelihood K = 27 only (run 1971)
plotQ(qlist=(chicks_slist)[c(1971)], imgoutput = "sep", grplab=palette1.chicks,ordergrp = TRUE,
      sortind = "all", sharedindlab = FALSE,clustercol=palette1, 
      subsetgrp = c("Koskenpää","Kummunsuo","Lehtosuo"),
      outputfilename = "K27_chicks_simplecolours_a", grplabangle = 25, panelspacer = 0.02, grplabsize = 6, 
      exportpath = getwd(),
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = FALSE, showtitle = F,
      height = 8, width = 20,  splabcol = "white", grplabpos = 0.75, linepos = 1)

plotQ(qlist=(chicks_slist)[c(1971)], imgoutput = "sep", grplab=palette1.chicks,ordergrp = TRUE,
      sortind = "all", sharedindlab = FALSE,clustercol=palette1, 
      subsetgrp = c("Nyrölä","Saarisuo","Teerisuo", "Utusuo"),
      outputfilename = "K27_chicks_simplecolours_b", grplabangle = 25, panelspacer = 0.02, grplabsize = 6, 
      exportpath = getwd(),
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = FALSE, showtitle = F,
      height = 8, width = 20,  splabcol = "white", grplabpos = 0.75, linepos = 1)

