
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 2, 2024
# Code purpose: analyze and visualize serology data for cattle sacrifice experiment

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
#library(gridExtra) # for grid.arrange function in doing multiple plots for jpeg
#library(tidyverse) # for piping %>%


###############################################
#  Serology Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
serodat <- read.csv("data/cattle_sac/sero/serodat_clean_cattlesacrifice_final.csv", header = TRUE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
serodat$eartag <- as.factor(serodat$eartag)
serodat$barn <- as.factor(serodat$barn)

# Setting the line and symbol colors for each barn and status
barncol <- c("1" = "purple4", "3" = "purple4")
color2 <- c("purple4")


# Split data into each barn, and have combined barns
barn1_sero  <- subset(serodat, barn == 1) 
barn3_sero  <- subset(serodat, barn == 3) 
barn13_sero <- subset(serodat, barn == 1 | barn == 3) 


########################################
#  Saving figures to file (as jpeg) ----
########################################

#####################
#  Barn 1 ----
#####################
# Barn 1: 8 Experimental Cattle

jpeg("output/cattle_sac/pprv_cattlesac_sero_barn1_todpi28.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn1_sero, aes(dpi, sn, color ="purple4")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28, by=2)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 1", 
                     labels = c("Goats (Inoculated)"),
                     values = "purple4") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())


#####################
#  Barn 3 ----
#####################
# Barn 3: 8 Experimental Cattle

jpeg("output/cattle_sac/pprv_cattlesac_sero_barn3_todpi28.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn3_sero, aes(dpi, sn, color="purple4")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28, by=2)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 3", 
                     labels = c("Goats (Inoculated)"),
                     values = "purple4") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())


######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats

jpeg("output/cattle_sac/pprv_cattlesac_sero_barns13_todpi28.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn13_sero, aes(dpi, sn, color=barn)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28, by=2)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "", 
                     labels = c("Barn 1", "Barn 3"),
                     values = barncol) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), axis.text=element_text(size=10))
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_sero_barns13c_todpi28.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn13_sero, aes(dpi, sn, color=species)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28, by=2)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "", 
                     #labels = c("Cattle"),
                     values = color2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), axis.text=element_text(size=10))
invisible(dev.off())

# Need to put the three plots above each into a named object and then put names in here in line 115 to plot a panel
# jpeg("output/goat_sac/pprv_cattlesac_clin_3panel_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(barn1_sero, barn3_sero, barn13_sero, ncol=3)
# invisible(dev.off())


