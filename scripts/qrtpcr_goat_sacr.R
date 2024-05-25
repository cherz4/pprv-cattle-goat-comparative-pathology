
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: April 8, 2024
# Code purpose: analyze and visualize qRT-PCR antigen data for goat sacrifice experiment

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
barncol <- c("2" = "cyan3", "6" = "cyan4")
color2 <- c("cyan4")


# Split data into each barn, and have combined barns
barn2_pcr <- subset(molecdat, barn == 2) 
barn6_pcr <- subset(molecdat, barn == 6)
barns26_pcr <- subset(molecdat, barn == 2 | barn == 6)


#####################
#  Barn 2 ----
#####################
# Barn 2: 8 Experimental Goats
# Smooth: LOESS 

# ocular
go_barn2 <- ggplot(barn2_pcr, aes(dpi, 35-omean, color = "cyan3")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Ocular", x = "", y = "") +
            scale_color_manual(name = "Barn 2", 
                               #labels = c("Goats", "Goats"),
                               values = "cyan3") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
gn_barn2 <- ggplot(barn2_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Nasal", x = "", y = "") +
            scale_color_manual(name = "Barn 2", 
                               #labels = c("Goats", "Goats"),
                               values = "cyan3") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
gr_barn2 <- ggplot(barn2_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Rectal", x = "", y = "") +
            scale_color_manual(name = "Barn 2", 
                               #labels = c("Goats", "Goats"),
                               values = "cyan3") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



#####################
#  Barn 6 ----
#####################
# Barn 6: 8 Experimental Goats
# Smooth: LOESS 

# ocular
go_barn6 <- ggplot(barn6_pcr, aes(dpi, 35-omean, color = "cyan3")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Ocular", x = "", y = "") +
            scale_color_manual(name = "Barn 6", 
                               #labels = c("Goats", "Goats"),
                               values = "cyan3") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
gn_barn6 <- ggplot(barn6_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Nasal", x = "", y = "") +
            scale_color_manual(name = "Barn 6", 
                               #labels = c("Goats", "Goats"),
                               values = "cyan3") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
gr_barn6 <- ggplot(barn6_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
            geom_line(aes(group = eartag), alpha = 0.4) +
            geom_point(size = 1, alpha = 0.4) +
            stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
            scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
            scale_y_continuous(limits = c(0,20)) +
            labs(title = "Rectal", x = "", y = "") +
            scale_color_manual(name = "Barn 6", 
                               #labels = c("Goats", "Goats"),
                               values = "cyan3") +
            geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
            theme_minimal() + 
            theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                  axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats
# Smooth: LOESS 

# ocular
go_barns26 <- ggplot(barns26_pcr, aes(dpi, 35-omean, color = barn)) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(0,20)) +
              labs(title = "Ocular", x = "", y = "") +
              scale_color_manual(name = "Barn", 
                                 #labels = c("Goats", "Goats"),
                                 values = barncol) +
              geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
gn_barns26 <- ggplot(barns26_pcr, aes(dpi, 35-nmean, color = barn)) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(0,20)) +
              labs(title = "Nasal", x = "", y = "") +
              scale_color_manual(name = "Barn", 
                                 #labels = c("Goats", "Goats"),
                                 values = barncol) +
              geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
gr_barns26 <- ggplot(barns26_pcr, aes(dpi, 35-rmean, color = barn)) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
              scale_y_continuous(limits = c(0,20)) +
              labs(title = "Rectal", x = "", y = "") +
              scale_color_manual(name = "Barn", 
                                 #labels = c("Goats", "Goats"),
                                 values = barncol) +
              geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# ocular barns combined
go_barns26c <- ggplot(barns26_pcr, aes(dpi, 35-omean, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Ocular", x = "", y = "") +
  scale_color_manual(name = "Barn", 
                     #labels = c("Goats"),
                     values = color2) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal barns combined
gn_barns26c <- ggplot(barns26_pcr, aes(dpi, 35-nmean, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn", 
                     #labels = c("Goats"),
                     values = color2) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal barns combined
gr_barns26c <- ggplot(barns26_pcr, aes(dpi, 35-rmean, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,8), breaks = seq(0:8)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "Barn", 
                     #labels = c("Goats"),
                     values = color2) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/goat_sac/pprv_goatsac_qRTPCR_barn2_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barn2, gn_barn2, gr_barn2, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_qRTPCR_barn6_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barn6, gn_barn6, gr_barn6, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_qRTPCR_barns26_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barns26, gn_barns26, gr_barns26, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/goat_sac/pprv_goatsac_qRTPCR_barns26c_todpi8.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_barns26c, gn_barns26c, gr_barns26c, ncol=3, left = "qRT-PCR 35 - mean Ct")
invisible(dev.off())

