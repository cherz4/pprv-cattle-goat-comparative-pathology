
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 2, 2024
# Code purpose: analyze and visualize serology data for goat sacrifice experiment

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
serodat <- read.csv("data/goat_sac/sero/serodat_clean_goatsacrifice_final.csv", header = TRUE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
serodat$eartag <- as.factor(serodat$eartag)
serodat$barn <- as.factor(serodat$barn)

# Setting the line and symbol colors for each barn and status
barncol <- c("2" = "cyan3", "6" = "cyan4")
color2 <- c("cyan4")


# Split data into each barn, and have combined barns
barn2_sero  <- subset(serodat, barn == 2) 
barn6_sero  <- subset(serodat, barn == 6) 
barn26_sero <- subset(serodat, barn == 2 | barn == 6) 

# Note although no sampling conducted on 24 dpi, this was in earlier trials to kept this on x axis for merging later.


########################################
#  Saving figures to file (as jpeg) ----
########################################

#####################
#  Barn 2 ----
#####################
# Barn 2: 8 Experimental Goats

jpeg("output/goat_sac/pprv_goatsac_sero_barn2_todpi8.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn2_sero, aes(dpi, sn, color ="cyan3")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,8), breaks = seq(0,8, by=1)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 2", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())


#####################
#  Barn 6 ----
#####################
# Barn 6: 8 Experimental Goats

jpeg("output/goat_sac/pprv_goatsac_sero_barn6_todpi8.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn6_sero, aes(dpi, sn, color="cyan3")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,8), breaks = seq(0,8, by=1)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 6", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
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

jpeg("output/goat_sac/pprv_goatsac_sero_barns26_todpi8.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn26_sero, aes(dpi, sn, color=barn)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,8), breaks = seq(0,8, by=1)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "", 
                     labels = c("Barn 2", "Barn 6"),
                     values = barncol) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), axis.text=element_text(size=10))
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_sero_barns26c_todpi8.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn26_sero, aes(dpi, sn, color=species)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,8), breaks = seq(0,8, by=1)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "", 
                     #labels = c("Goats"),
                     values = color2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), axis.text=element_text(size=10))
invisible(dev.off())


# Need to put the three plots above each into a named object and then put names in here in line 117 to plot a panel
# jpeg("output/goat_sac/pprv_goatsac_clin_3panel_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(barn2_sero, barn6_sero, barn26_sero, ncol=3)
# invisible(dev.off())


