
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: April 8, 2024
# Code purpose: analyze and visualize AgELISA antigen data for goat sacrifice experiment

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
library(gridExtra) # for grid.arrange function in doing multiple plots for jpeg
library(tidyverse) # for piping %>%



###############################################
#  Molecular Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
molecdat <- read.csv("data/goat_sac/mol/molecdat_clean_goatsac.csv", header = TRUE, stringsAsFactors = FALSE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
molecdat_noc <- molecdat %>% select(eartag, expt, species, dpi, barn, od, sp, status, 
                                    od_f, sp_f, status_f, notes1, notes2)

molecdat_noc$eartag <- as.factor(molecdat_noc$eartag)
molecdat_noc$expt <- as.factor(molecdat_noc$expt)
molecdat_noc$species <- as.factor(molecdat_noc$species)
molecdat_noc$barn <- as.factor(molecdat_noc$barn)
molecdat_noc$status <- as.factor(molecdat_noc$status)
molecdat_noc$status_f <- as.factor(molecdat_noc$status_f)

# Setting the line and symbol colors for each barn and status
barncol <- c("2" = "cyan3", "6" = "cyan4")
color2 <- c("cyan4")


# Split data into each barn, and have combined barns
barn2_ag <- subset(molecdat_noc, barn == 2) 
barn6_ag <- subset(molecdat_noc, barn == 6)
barns26_ag <- subset(molecdat_noc, barn == 2 | barn == 6)


#####################
#  Barn 2 ----
#####################
# Barn 2: 8 Experimental Goats
# Smooth: LOESS 

ag_n_barn2 <- ggplot(barn2_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 2: Nasal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barn2 <- ggplot(barn2_ag, aes(dpi, sp_f, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 2: Rectal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 6 ----
#####################
# Barn 6: 8 Experimental Goats
# Smooth: LOESS 

ag_n_barn6 <- ggplot(barn6_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 6: Nasal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") + 
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


ag_r_barn6 <- ggplot(barn6_ag, aes(dpi, sp_f, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 6: Rectal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() +
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))




######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats
# Smooth: LOESS 

ag_n_barns26 <- ggplot(barns26_ag, aes(dpi, sp, color = barn, linetype = barn)) +
                geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
                scale_y_continuous(limits = c(-5, 515)) +
                labs(title = "Nasal", x = "Days Post Infection", y = "") +
                scale_color_manual(#name = "Nasal Status", 
                                   #labels = c("Goats", "Goats"),
                                   values = barncol) +  
                geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barns26 <- ggplot(barns26_ag, aes(dpi, sp_f, color = barn, linetype = barn)) +
                geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
                scale_y_continuous(limits = c(-5, 515)) +
                labs(title = "Rectal", x = "Days Post Infection", y = "") +
                scale_color_manual(#name = "Rectal Status", 
                                   #labels = c("Goats", "Goats"),
                                   values = barncol) +  
                geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


ag_n_barns26c <- ggplot(barns26_ag, aes(dpi, sp, color = species)) +
                 geom_line(aes(group = eartag), alpha = 0.4) +
                 geom_point(size = 1, alpha = 0.4) +
                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                 scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
                 scale_y_continuous(limits = c(-5, 515)) +
                 labs(title = "Nasal", x = "", y = "") +
                 scale_color_manual(#name = "Nasal Status",
                   #labels = c("Goats"),
                   values = color2) +
                 geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                 theme_minimal() +
                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barns26c <- ggplot(barns26_ag, aes(dpi, sp_f, color = species)) +
                 geom_line(aes(group = eartag), alpha = 0.4) +
                 geom_point(size = 1, alpha = 0.4) +
                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                 scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
                 scale_y_continuous(limits = c(-5, 515)) +
                 labs(title = "Rectal", x = "", y = "") +
                 scale_color_manual(#name = "Rectal Status", 
                   #labels = c("Goats"),
                   values = color2) +  
                 geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                 theme_minimal() +
                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/goat_sac/pprv_goatsac_agelisa_barn2_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn2, ag_r_barn2, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_agelisa_barn6_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn6, ag_r_barn6, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_agelisa_barns26_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barns26, ag_r_barns26, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_agelisa_barns26c_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barns26c, ag_r_barns26c, ncol=2, left = "S/P Ratio")
invisible(dev.off())
