############# 6. Creating plots ##############

## In this script, you will find the code used to create the figures of the manuscript

### load packages ###

library(RColorBrewer); library(ggplot2); library(extrafont); 
library(dplyr); library(readxl); library(devtools); library(pophelper)
library(gridExtra); library(adegenet); library(ggthemr)

#palette
palette1 = c("#be4d5a","#474b4e","#8b95c9","#f7f6f0","#57939a","#d7d6d5","#4d9de0","#785964","#918ef4","#519e8a",
                    "#df2a3f","#224f71","#7184e5","#fbf8ea","#2ab6c5","#ebd6c1","#3d9ff0","#a62c56","#8e8afa","#28c89d",
                    "#ff0a27","#005494","#5773ff","#fffae5","#00d8f0","#ffd6ad","#2ea1ff","#d10049","#8985ff","#00f0b0")
mypalette3 <- c("#EDEDFD","#C2C1EC","#9795DB","#6C69C9","#413DB8","#1611A7") #purples from structure
theme_set(theme_classic())
#### Figure 2: heatmaps pairwise Fst + correlograms spatial autocorrelation ####

### First: heatmap
pairwise.fst <- read.csv("tables/Pairwise_Fst_all.csv")

#change the negative values to 0 
pairwise.fst$Fst[which (pairwise.fst$Fst < 0)] <- 0

# add abbreviation to data
sitenames <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")
# add abbreviations to sitenames
sitenames[1]
sitenames$abb <- c("KOS", "KUM", "LAU", "LEH", "NYR", "PAL", "PIH", "PIL", "PIS", "SAA", "TEE", "UTU")
str(sitenames)

pairwise.fst <- left_join(pairwise.fst, sitenames[,c(1,6,2)], by = c("site.x" = "pop_num"))
pairwise.fst <- left_join(pairwise.fst, sitenames[,c(1,6,2)], by = c("site.y" = "pop_num"))

### All - pairwise Fst heatmap
pairwise.fst.fig <- ggplot(pairwise.fst, aes(abb.x, abb.y, fill = Fst)) + geom_tile() + theme_classic() + 
  scale_fill_gradientn(colors = mypalette3, limits = c(0,0.02)) + 
  geom_text(aes(label = Sig), size = 8)+
  theme(text = element_text(family = "Arial", size = 22),
        axis.text.x = element_text(angle = 90, size = 26,
                                   face = c("bold", "plain", "bold", "plain", "plain",
                                            "bold", "bold", "plain", "bold",
                                            "plain", "plain")), 
        axis.text.y = element_text(size = 26, 
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
  ggtitle('(a) Fst heatmap') 

ggsave(pairwise.fst.fig, filename="figures/Fst_all.png")

### Second - correlograms ###
spatial <- read_excel("tables/SpatialAutocor_04.10.22.xlsx", sheet = "ForR")

### Males - spatial
males.spatial <- spatial %>% filter(Who == "Male") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + 
  geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, 
                size=1.5, color="black", stat = "identity")+
  ylim(-0.05, 0.04) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_text(size = 38),
        axis.text.x = element_text(size = 26, margin = margin(b = 10)), 
        axis.text.y = element_text(size = 26, margin = margin(l = 10)),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30)) +
  ggtitle('(b) Correlogram adult males') 

ggsave(males.spatial, file="figures/Correlogram_males.png")

### Females - spatial

females.spatial <- spatial %>% filter(Who == "Females") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, size=1.5, color="black", stat = "identity")+
  ylim(-0.05, 0.04) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Helvetica", size = 14),
        plot.title = element_text(size = 38),
        axis.text.x = element_text(size = 26, margin = margin(b = 10)), 
        axis.text.y = element_text(size = 26, margin = margin(l = 10)),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30))+
  ggtitle('(c) Correlogram adult females') 

ggsave(females.spatial, file="figures/Correlogram_females.png")

### Chicks - spatial
chicks.spatial <- spatial %>% filter(Who == "Unrelated_chicks") %>% ggplot(aes(x = What)) + 
  geom_line(aes(y = r), col = "#8989D0", size = 1.5) + theme_classic() + geom_hline(yintercept = 0, col = "black")+
  geom_ribbon(aes(ymax = U, ymin = L), fill = "#8989D0", alpha = 0.5)+
  geom_errorbar(aes(y = r, ymin = r-Le, ymax = r+Ue), width=1, size=1.5, color="black", stat = "identity")+
  ylim(-0.05, 0.04) +
  xlab("Distance class") + ylab("Autocorrelation coefficient r")+
  scale_x_continuous(breaks = c(seq(0, 60, by = 10)), limits=c(4,61))+
  theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_text(size = 38),
        axis.text.x = element_text(size = 26, margin = margin(b = 10)), 
        axis.text.y = element_text(size = 26, margin = margin(l = 10)),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30)) +
  ggtitle('(d) Correlogram unrelated chicks') 

ggsave(chicks.spatial, file="figures/Correlogram_chicks.png")

png("figures/Fst_Correlograms.png", 
    width = 1500, height = 1000, units = "px")
gA <- ggplotGrob(pairwise.fst.fig)
gB <- ggplotGrob(males.spatial)
gD <- ggplotGrob(females.spatial)
gE <- ggplotGrob(chicks.spatial)
grid::grid.newpage()
grid::grid.draw(cbind(rbind(gA, gD), rbind(gB, gE)))
dev.off()

### Figure 3: boxplots migration rates ####
male_run5_clean <- read.csv("data/migrationanalysis/run5_males_clean.csv")
female_run5_clean <- read.csv("data/migrationanalysis/run5_females_clean.csv")

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
  scale_fill_manual(values = c("Male" = "#57939a", #44A1A0
                               "Female" = "#be4d5a"),#AF9BB6
                    labels=c("Male", "Female"))+
  guides(fill = guide_legend("Sex")) 

library(ggsignif)

emigration.plot.poster <- ggplot(migration_both, aes(x = hunt_out, y = migration_ESSc, fill = sex)) + 
  geom_boxplot(outlier.shape = NA, aes(middle = mean(migration_ESSc))) + ylim(0, 0.03)+
  ylab("Migration rate")+
  theme(text = element_text(family = "Arial", size = 36),
        legend.position = "none",
        plot.title = element_text(size = 38),
        axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_fill_manual(values = c("Male" = "#44A1A0", #44A1A0
                               "Female" = "#AF9BB6"),#AF9BB6
                    labels=c("Male", "Female"))+ 
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank() #remove minor gridlines
  ) + geom_text(x = 1.5, y = 0.028, label = "*", size = 10) + geom_segment(x = 1, xend = 2, y = 0.025, yend = 0.025)+
  geom_segment(x=1, xend = 1, y = 0.025, yend = 0.024)+geom_segment(x=2, xend = 2, y = 0.025, yend = 0.024)

emigration.plot
emigration.plot.poster

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

immmigration.plot.poster <- ggplot(migration_both, aes(x = hunt_in, y = migration_ESSc, fill = sex)) + 
  geom_boxplot(outlier.shape = NA, aes(middle = mean(migration_ESSc))) + ylim(0, 0.03)+
  ylab("Migration rate")+
  theme(text = element_text(family = "Arial", size = 36),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        axis.text = element_text(colour = "black"),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size = 38),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_fill_manual(values = c("Male" = "#44A1A0",
                               "Female" = "#AF9BB6"),
                    labels=c("Male", "Female"))+
  guides(fill = guide_legend("Sex"))+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'))+ #transparent legend bg )
  geom_text(x = 1.5, y = 0.028, label = "*", size = 10) + geom_segment(x = 1, xend = 2, y = 0.025, yend = 0.025)+
  geom_segment(x=1, xend = 1, y = 0.025, yend = 0.024)+geom_segment(x=2, xend = 2, y = 0.025, yend = 0.024)

immmigration.plot
immmigration.plot.poster

ggsave(plot = emigration.plot, "data/figures/Emigration.png",
       width = 30, height = 20, units = c("cm"))

ggsave(plot = immmigration.plot, "data/figures/Immigration.png",
       width = 30, height = 20, units = c("cm"))

ggsave(plot = emigration.plot.poster, "data/figures/Emigration_forposter.png",
       width = 30, height = 20, units = c("cm"))

ggsave(plot = immmigration.plot.poster, "data/figures/Immigration_forposter.png",
       width = 30, height = 20, units = c("cm"))

#### Supplementary Figure 1: log likelihood per K (STRUCTURE) ##### 

ks_all <- read.csv("analyses/structure/results/Run_files/output.csv")

#### plotting K vs mean log likelihood

ks_all_plot <- ggplot(ks_all, aes(x = k, y = elpdmean)) +
  geom_point(col="black") + geom_line(col="black")+ ylab("LnPr(X|k)") +
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
        plot.title = element_text(size = 22)) 

ks_all_plot

library(cowplot); library(gridExtra)
ggsave(plot=ks_all_plot, "figures/K_combined.png", device = "png",
       width = 3000, height = 1000, unit="px")

##### Supplementary Figure 2: barplots STRUCTURE #####

##### First plot the bargraphs for 10 different K's an example of 10 different runs
all_path_to_structure_out <- "analyses/structure/results/Run_files/"
all_path_to_struc_file <- "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru"
all_all_files <- list.files(all_path_to_structure_out, pattern = "^results")
all_struc_out_paths <- paste0(all_path_to_structure_out, all_all_files)
all_slist <- readQ(files=all_struc_out_paths, filetype = "structure")
fn1 <- function(x) attr(x,"k")
all_spnames <- paste0("K=",sapply(all_slist,fn1))
ggthemr(palette = "solarized", layout = "clean",
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
all_struc_file <- all_path_to_struc_file
all_pops <- fread(all_struc_file) %>%
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
all_pops <- data.frame(all_pops[seq(1, nrow(all_pops), 2)])
colnames(all_pops) <- "location"
all_pops$location <- as.character(all_pops$location)

## all: best log likelihood and also highest delta K = 4

## just highest log likelihood K = 4 (run 61) and split up per group
plotQ(qlist=(all_slist)[c(61)],imgoutput = "sep", 
      grplab=all_pops,ordergrp = TRUE,sortind = "all",
      clustercol=palette1,
      subsetgrp = c("Koskenpää","Kummunsuo","Lauttasuo","Lehtosuo","Nyrölä","Palosuo"),
      outputfilename = "K4_combined_a", exportpath = paste0(getwd(), "/figures"), 
      sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = F, 
      grplabangle = 90, panelspacer = 0.02, grplabsize = 5, showtitle = F, titlesize = 14, titlespacer = 2,
      height = 8, width = 20, splab = "", splabcol = "white", grplabpos = 0.54, linepos = 1)

plotQ(qlist=(all_slist)[c(61)],imgoutput = "sep", grplab=all_pops,ordergrp = TRUE,sortind = "all",
      clustercol=palette1,
      subsetgrp = c("Pihtissuo", "Pirttilampi", "Pirttisuo", "Saarisuo", "Teerisuo", "Utusuo"),
      outputfilename = "K4_combined_b", exportpath = paste0(getwd(), "/figures"), 
      sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = F, 
      grplabangle = 90, panelspacer = 0.02, grplabsize = 5, showtitle = F, titlesize = 14, titlespacer = 2,
      height = 8, width = 20, splab = "", splabcol = "white", grplabpos = 0.65, linepos = 1)

## K = 2
plotQ(qlist=(all_slist)[c(41)],imgoutput = "sep", 
      grplab=all_pops,ordergrp = TRUE,sortind = "all",
      clustercol=palette1,
      subsetgrp = c("Koskenpää","Kummunsuo","Lauttasuo","Lehtosuo","Nyrölä","Palosuo"),
      outputfilename = "K2_combined_a", exportpath = paste0(getwd(), "/figures"), 
      sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = F, 
      grplabangle = 25, panelspacer = 0.02, grplabsize = 6, showtitle = F, titlesize = 14, titlespacer = 2,
      height = 8, width = 20, splab = "", splabcol = "white", grplabpos = 0.75, linepos = 1)

plotQ(qlist=(all_slist)[c(41)],imgoutput = "sep", grplab=all_pops,ordergrp = TRUE,sortind = "all",
      clustercol=palette1,
      subsetgrp = c("Pihtissuo", "Pirttilampi", "Pirttisuo", "Saarisuo", "Teerisuo", "Utusuo"),
      outputfilename = "K2_combined_b", exportpath = paste0(getwd(), "/figures"), 
      sharedindlab = F,
      panelratio = c(1.8,1.2), font = "Arial", grplabjust = 0.4, showsp = F, 
      grplabangle = 25, panelspacer = 0.02, grplabsize = 6, showtitle = F, titlesize = 14, titlespacer = 2,
      height = 8, width = 20, splab = "", splabcol = "white", grplabpos = 0.75, linepos = 1)



