
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: April 8, 2024
# Code purpose: analyze and visualize qRT-PCR antigen data for cattle sacrifice experiment

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
molecdat$eartag <- as.factor(molecdat$eartag)
molecdat$expt <- as.factor(molecdat$expt)
molecdat$species <- as.factor(molecdat$species)
molecdat$barn <- as.factor(molecdat$barn)
molecdat$status <- as.factor(molecdat$status)
molecdat$status_f <- as.factor(molecdat$status_f)

molecdat$o1[molecdat$o1 == "U"] <- 35
molecdat$o2[molecdat$o2 == "U"] <- 35
molecdat$n1[molecdat$n1 == "U"] <- 35
molecdat$n2[molecdat$n2 == "U"] <- 35
molecdat$r1[molecdat$r1 == "U"] <- 35
molecdat$r2[molecdat$r2 == "U"] <- 35
molecdat$omean[molecdat$omean == 0] <- "35"
molecdat$nmean[molecdat$nmean == 0] <- "35"
molecdat$rmean[molecdat$rmean == 0] <- "35"
molecdat$o1 <- as.numeric(molecdat$o1)
molecdat$o2 <- as.numeric(molecdat$o2)
molecdat$n1 <- as.numeric(molecdat$n1)
molecdat$n2 <- as.numeric(molecdat$n2)
molecdat$r1 <- as.numeric(molecdat$r1)
molecdat$r2 <- as.numeric(molecdat$r2)
molecdat$omean <- as.numeric(molecdat$omean)
molecdat$nmean <- as.numeric(molecdat$nmean)
molecdat$rmean <- as.numeric(molecdat$rmean)


# Setting the line and symbol colors for each barn and status
barncol <- c("1" = "purple4", "3" = "purple4")
color2 <- c("purple4")


# Split data into each barn, and have combined barns
barn1_pcr <- subset(molecdat, barn == 1) 
barn3_pcr <- subset(molecdat, barn == 3)
barns13_pcr <- subset(molecdat, barn == 1 | barn == 3)


#####################
#  Barn 1 ----
#####################
# Barn 1: 8 Experimental Cattle
# Smooth: LOESS 

# ocular
go_barn1 <- ggplot(barn1_pcr, aes(dpi, 35-omean, color = "purple4")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Ocular", x = "", y = "") +
            scale_color_manual(name = "Barn 1", 
                               #labels = c("Goats", "Goats"),
                               values = "purple4") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
gn_barn1 <- ggplot(barn1_pcr, aes(dpi, 35-nmean, color = "purple4")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Nasal", x = "", y = "") +
            scale_color_manual(name = "Barn 1", 
                               #labels = c("Goats", "Goats"),
                               values = "purple4") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
gr_barn1 <- ggplot(barn1_pcr, aes(dpi, 35-rmean, color = "purple4")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Rectal", x = "", y = "") +
            scale_color_manual(name = "Barn 1", 
                               #labels = c("Goats", "Goats"),
                               values = "purple4") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



#####################
#  Barn 3 ----
#####################
# Barn 3: 8 Experimental Cattle
# Smooth: LOESS 

# ocular
go_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-omean, color = "purple4")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Ocular", x = "", y = "") +
            scale_color_manual(name = "Barn 3", 
                               #labels = c("Goats", "Goats"),
                               values = "purple4") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
gn_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-nmean, color = "purple4")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Nasal", x = "", y = "") +
            scale_color_manual(name = "Barn 3", 
                               #labels = c("Goats", "Goats"),
                               values = "purple4") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
gr_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-rmean, color = "purple4")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Rectal", x = "", y = "") +
            scale_color_manual(name = "Barn 3", 
                               #labels = c("Goats", "Goats"),
                               values = "purple4") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Cattle
# Smooth: LOESS 

# ocular
go_barns13 <- ggplot(barns13_pcr, aes(dpi, 35-omean, color = barn, linetype = barn)) +
              geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(0,20)) +
              labs(title = "Ocular", x = "", y = "") +
              scale_color_manual(#name = "Barn", 
                                 #labels = c("Goats", "Goats"),
                                 values = barncol) +
              geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
gn_barns13 <- ggplot(barns13_pcr, aes(dpi, 35-nmean, color = barn, linetype = barn)) +
              geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(0,20)) +
              labs(title = "Nasal", x = "", y = "") +
              scale_color_manual(#name = "Barn", 
                                 #labels = c("Goats", "Goats"),
                                 values = barncol) +
              geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
gr_barns13 <- ggplot(barns13_pcr, aes(dpi, 35-rmean, color = barn, linetype = barn)) +
              geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
              scale_y_continuous(limits = c(0,20)) +
              labs(title = "Rectal", x = "", y = "") +
              scale_color_manual(# name = "Barn", 
                                 #labels = c("Goats", "Goats"),
                                 values = barncol) +
              geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# ocular barns combined
go_barns13c <- ggplot(barns13_pcr, aes(dpi, 35-omean, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "", x = "", y = "") +
  scale_color_manual(#name = "Barn", 
    #labels = c("Goats"),
    values = color2) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal barns combined
gn_barns13c <- ggplot(barns13_pcr, aes(dpi, 35-nmean, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "", x = "", y = "") +
  scale_color_manual(#name = "Barn", 
    #labels = c("Goats"),
    values = color2) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal barns combined
gr_barns13c <- ggplot(barns13_pcr, aes(dpi, 35-rmean, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "", x = "", y = "") +
  scale_color_manual(# name = "Barn", 
    #labels = c("Goats"),
    values = color2) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/cattle_sac/pprv_cattlesac_qRTPCR_barn1_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barn1, gn_barn1, gr_barn1, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_qRTPCR_barn3_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barn3, gn_barn3, gr_barn3, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_qRTPCR_barns13_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barns13, gn_barns13, gr_barns13, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/cattle_sac/pprv_cattlesac_qRTPCR_barns13c_todpi28.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barns13c, gn_barns13c, gr_barns13c, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())
