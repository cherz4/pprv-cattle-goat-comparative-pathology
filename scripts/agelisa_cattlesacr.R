
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: April 8, 2024
# Code purpose: analyze and visualize AgELISA antigen data for cattle sacrifice experiment

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
molecdat <- read.csv("data/cattle_sac/mol/molecdat_clean_cattlesac.csv", header = TRUE, stringsAsFactors = FALSE)

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
barncol <- c("1" = "purple4", "3" = "purple4")
color2 <- c("purple4")


# Split data into each barn, and have combined barns
barn1_ag <- subset(molecdat_noc, barn == 1) 
barn3_ag <- subset(molecdat_noc, barn == 3)
barns13_ag <- subset(molecdat_noc, barn == 1 | barn == 3)


#####################
#  Barn 1 ----
#####################
# Barn 1: 8 Experimental Cattle
# Smooth: LOESS 

ag_n_barn1 <- ggplot(barn1_ag, aes(dpi, sp, color = "purple4")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 1: Nasal Status", 
                                 #labels = c("Cattle", "Cattle"),
                                 values = "purple4") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barn1 <- ggplot(barn1_ag, aes(dpi, sp_f, color = "purple4")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 1: Rectal Status", 
                                 #labels = c("Cattle", "Cattle"),
                                 values = "purple4") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 3 ----
#####################
# Barn 3: 8 Experimental Cattle
# Smooth: LOESS 

ag_n_barn3 <- ggplot(barn3_ag, aes(dpi, sp, color = "purple4")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 3: Nasal Status", 
                                 #labels = c("Cattle", "Cattle"),
                                 values = "purple4") + 
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barn3 <- ggplot(barn3_ag, aes(dpi, sp_f, color = "purple4")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 3: Rectal Status", 
                                 #labels = c("Cattle", "Cattle"),
                                 values = "purple4") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() +
              theme(legend.position = "none", axis.title.y = element_text(size = 8),
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))




######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Cattle
# Smooth: LOESS 

ag_n_barns13 <- ggplot(barns13_ag, aes(dpi, sp, color = barn, linetype = barn)) +
                geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
                scale_y_continuous(limits = c(-5, 515)) +
                labs(title = "Nasal", x = "Days Post Infection", y = "") +
                scale_color_manual(#name = "Nasal Status", 
                                   #labels = c("Cattle", "Cattle"),
                                   values = barncol) +  
                geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barns13 <- ggplot(barns13_ag, aes(dpi, sp_f, color = barn, linetype = barn)) +
                geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
                scale_y_continuous(limits = c(-5, 515)) +
                labs(title = "Rectal", x = "Days Post Infection", y = "") +
                scale_color_manual(#name = "Rectal Status", 
                                   #labels = c("Cattle", "Cattle"),
                                   values = barncol) +  
                geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_n_barns13c <- ggplot(barns13_ag, aes(dpi, sp, color = species)) +
                 geom_line(aes(group = eartag), alpha = 0.4) +
                 geom_point(size = 1, alpha = 0.4) +
                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                 scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
                 scale_y_continuous(limits = c(-5, 515)) +
                #labs(title = "Nasal", x = "Days Post Infection", y = "") +
                 labs(title = "", x = "", y = "") +
                 scale_color_manual(#name = "Nasal Status", 
                   #labels = c("Cattle"),
                   values = color2) +  
                 geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                 theme_minimal() +
                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barns13c <- ggplot(barns13_ag, aes(dpi, sp_f, color = species)) +
                 geom_line(aes(group = eartag), alpha = 0.4) +
                 geom_point(size = 1, alpha = 0.4) +
                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                 scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
                 scale_y_continuous(limits = c(-5, 515)) +
                 #labs(title = "Rectal", x = "Days Post Infection", y = "") +
                 labs(title = "", x = "", y = "") +
                 scale_color_manual(#name = "Rectal Status", 
                   #labels = c("Cattle"),
                   values = color2) +  
                 geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                 theme_minimal() +
                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/cattle_sac/pprv_cattlesac_agelisa_barn1_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn1, ag_r_barn1, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_agelisa_barn3_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn3, ag_r_barn3, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_agelisa_barns13_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barns13, ag_r_barns13, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_agelisa_barns13c_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barns13c, ag_r_barns13c, ncol=2, left = "S/P Ratio", bottom = "Days Post Infection")
invisible(dev.off())
