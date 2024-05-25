
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 2, 2024
# Code purpose: analyze and visualize clinical data for cattle sacrifice experiment

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
# library(mgcv)    # need if doing a GAM (generalized additive model) for the smooth
library(gridExtra) # for grid.arrange function in doing multiple plots for jpeg
library(tidyverse) # for piping %>%


###############################################
#  Clinical Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
clindat <- read.csv("data/cattle_sac/clin/clindat_clean_cattlesac.csv", header = TRUE) # use csv file name on your computer here

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
clindat$eartag <- as.factor(clindat$eartag)
clindat$barn <- as.factor(clindat$barn)
clindat$tempc <- as.numeric(clindat$tempc)
clindat <- subset(clindat, !dpi <0) # making sure there are no dpi prior to experiment start in the data
#clindat$dpi <- as.integer(clindat$dpi)

# Create clinical score without temperature - this is used throughout the file
clindat <- clindat %>% mutate(clinscore_mod = clindat$general + clindat$feces + clindat$odnd + clindat$headmucos +
                                                  clindat$resp)

# Setting the line and symbol colors for each barn and status
barncol <- c("1" = "purple4", "3" = "purple4")
color2 <- c("purple4")


# Split data into each barn, and have combined barns
barn1 <- subset(clindat, barn == 1) 
barn3 <- subset(clindat, barn == 3)
barns13 <- subset(clindat, barn == 1 | barn == 3)


#####################
#  Barn 1 ----
#####################
# Barn 1: 8 Experimental Cattle
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn1i1 <- barn1[barn1$tempc>=36 & barn1$tempc<=41,]
barn1i1_fit <- loess(tempc ~ dpi, barn1i1)
barn1i1_nd <- data.frame(dpi=seq(min(barn1i1$dpi), max(barn1i1$dpi), length=100))
barn1i1_nd$fit <- predict(barn1i1_fit, newdata=barn1i1_nd)
barn1i1tmax <- barn1i1_nd$dpi[which.max(barn1i1_nd$fit)]
#plot(barn1i1$tempc~barn1i1$dpi)
#lines(barn1i1_nd$fit~barn1i1_nd$dpi, col="grey")

# clinical score
barn1i1c <- barn1
barn1i1c_fit <- loess(clinscore_mod ~ dpi, barn1i1c)
barn1i1c_nd <- data.frame(dpi=seq(min(barn1i1c$dpi), max(barn1i1c$dpi), length=100))
barn1i1c_nd$fit <- predict(barn1i1c_fit, newdata=barn1i1c_nd)
barn1i1cmax <- barn1i1c_nd$dpi[which.max(barn1i1c_nd$fit)]
# plot(barn1i1c$clinscore_mod~barn1i1c$dpi)
# lines(barn1i1c_nd$fit~barn1i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn1 <- ggplot(barn1, aes(dpi, tempc, color="purple4")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)"),
                     values = "purple4") +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired  NOTE: had to shift backward to see peak label in second line
  geom_vline(aes(xintercept = barn1i1tmax), color = "purple4") +
  annotate(geom = "text", x = (barn1i1tmax - 1.5), y = 41, label = paste(round(barn1i1tmax,1)), color = "purple4") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn1 <- ggplot(barn1, aes(dpi, clinscore_mod, color="purple4")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)"),
                     values = "purple4") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn1i1cmax), color = "purple4") +
  annotate(geom = "text", x = (barn1i1cmax + 1), y = 10, label = paste(round(barn1i1cmax,1)), color = "purple4") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



#####################
#  Barn 3 ----
#####################
# Barn 3: 8 Experimental Cattle
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn3i1 <- barn3[barn3$tempc>=36 & barn3$tempc<=41,]
barn3i1_fit <- loess(tempc ~ dpi, barn3i1)
barn3i1_nd <- data.frame(dpi=seq(min(barn3i1$dpi), max(barn3i1$dpi), length=100))
barn3i1_nd$fit <- predict(barn3i1_fit, newdata=barn3i1_nd)
barn3i1tmax <- barn3i1_nd$dpi[which.max(barn3i1_nd$fit)]
# plot(barn3i1$tempc~barn3i1$dpi)
# lines(barn3i1_nd$fit~barn3i1_nd$dpi, col="grey")

#clinical score
barn3i1c <- barn3
barn3i1c_fit <- loess(clinscore_mod ~ dpi, barn3i1c)
barn3i1c_nd <- data.frame(dpi=seq(min(barn3i1c$dpi), max(barn3i1c$dpi), length=100))
barn3i1c_nd$fit <- predict(barn3i1c_fit, newdata=barn3i1c_nd)
barn3i1cmax <- barn3i1c_nd$dpi[which.max(barn3i1c_nd$fit)]
# plot(barn3i1c$clinscore_mod~barn3i1c$dpi)
# lines(barn3i1c_nd$fit~barn3i1c_nd$dpi, col="grey")

# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn3 <- ggplot(barn3, aes(dpi, tempc, color="purple4")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)"),
                     values = "purple4") +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn3i1tmax), color = "purple4") +
  annotate(geom = "text", x = (barn3i1tmax + 1), y = 41, label = paste(round(barn3i1tmax,1)), color = "purple4") + 
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn3 <- ggplot(barn3, aes(dpi, clinscore_mod, color="purple4")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)"),
                     values = "purple4") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn3i1cmax), color = "purple4") +
  annotate(geom = "text", x = (barn3i1cmax + 1), y = 10, label = paste(round(barn3i1cmax,1)), color = "purple4") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))




######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Cattle
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barns13i1 <- barns13[barns13$tempc>=36 & barns13$tempc<=41,]
barns13i1_fit <- loess(tempc ~ dpi, barns13i1)
barns13i1_nd <- data.frame(dpi=seq(min(barns13i1$dpi), max(barns13i1$dpi), length=100))
barns13i1_nd$fit <- predict(barns13i1_fit, newdata=barns13i1_nd)
barns13i1tmax <- barns13i1_nd$dpi[which.max(barns13i1_nd$fit)]
# plot(barns13i1$tempc~barns13i1$dpi)
# lines(barns13i1_nd$fit~barns13i1_nd$dpi, col="grey")

#clinical score
barns13i1c <- barns13
barns13i1c_fit <- loess(clinscore_mod ~ dpi, barns13i1c)
barns13i1c_nd <- data.frame(dpi=seq(min(barns13i1c$dpi), max(barns13i1c$dpi), length=100))
barns13i1c_nd$fit <- predict(barns13i1c_fit, newdata=barns13i1c_nd)
barns13i1cmax <- barns13i1c_nd$dpi[which.max(barns13i1c_nd$fit)]
# plot(barns13i1c$clinscore_mod~barns13i1c$dpi)
# lines(barns13i1c_nd$fit~barns13i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barns13 <- ggplot(barns13, aes(dpi, tempc, color=barn, linetype = barn)) +
  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_linetype_manual(name = "", values = c("1" = "solid", "3" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Cattle", "Cattle"),
                     values = barncol) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired  NOTE: shifted value farther to right due to 3 digits
  geom_vline(aes(xintercept = barns13i1tmax), color = "grey50") +
  annotate(geom = "text", x = (barns13i1tmax + 1.5), y = 41, label = paste(round(barns13i1tmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barns13 <- ggplot(barns13, aes(dpi, clinscore_mod, color=barn, linetype = barn)) +
  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_linetype_manual(name = "", values = c("1" = "solid", "3" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Cattle", "Cattle"),
                     values = barncol) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barns13i1cmax), color = "grey50") +
  annotate(geom = "text", x = (barns13i1cmax + 1), y = 10, label = paste(round(barns13i1cmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barns13c <- ggplot(barns13, aes(dpi, tempc, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_linetype_manual(name = "", values = c("1" = "solid", "3" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Cattle"),
                     values = color2) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired  NOTE: shifted value farther to right due to 3 digits
  geom_vline(aes(xintercept = barns13i1tmax), color = "grey50") +
  annotate(geom = "text", x = (barns13i1tmax + 1.5), y = 41, label = paste(round(barns13i1tmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barns13c <- ggplot(barns13, aes(dpi, clinscore_mod, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, linewidth =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), linewidth =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_linetype_manual(name = "", values = c("1" = "solid", "3" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Cattle"),
                     values = color2) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barns13i1cmax), color = "grey50") +
  annotate(geom = "text", x = (barns13i1cmax + 1), y = 10, label = paste(round(barns13i1cmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/cattle_sac/pprv_cattlesac_clin_barn1_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn1, cs_barn1, ncol=2)
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_clin_barn3_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn3, cs_barn3, ncol=2)
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_clin_barns13_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barns13, cs_barns13, ncol=2)
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_clin_barns13c_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barns13c, cs_barns13c, ncol=2)
invisible(dev.off())

# If need symbols by unique eartag, sample code for this is in PPRV Trial Stage 1.1A code
