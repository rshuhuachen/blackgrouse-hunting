### Second - correlograms ###
spatial <- read_excel("data/tables/SpatialAutocor_Supplementary", sheet = "ForR")

make_corr <- function(data){
  plot <- ggplot(data, aes(x = What)) + 
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
  ggsave(plot = plot, paste0("data/figures/", data, ".png"),
         width = 30, height = 20, units = c("cm"))
}

males_hunted <- subset(spatial, data == "male" & hunt == "hunted")
males_unhunted <- subset(spatial, data == "male" & hunt == "unhunted")
females_hunted <- subset(spatial, data == "female" & hunt == "hunted")
females_unhunted <- subset(spatial, data == "female" & hunt == "unhunted")
chicks_hunted <- subset(spatial, data == "chick" & hunt == "hunted")
chicks_unhunted <- subset(spatial, data == "chick" & hunt == "unhunted")
unrelated_chick_hunted <- subset(spatial, data == "unrelated_chick" & hunt == "hunted")
unrelated_chick_unhunted <- subset(spatial, data == "unrelated_chick" & hunt == "unhunted")

datafiles <- list(males_hunted, males_unhunted, females_hunted, females_unhunted, 
     chicks_hunted, chicks_unhunted, unrelated_chick_hunted, unrelated_chick_unhunted)

lapply(datafiles, make_corr())