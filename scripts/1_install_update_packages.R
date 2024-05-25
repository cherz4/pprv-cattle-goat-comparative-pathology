
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 23, 2023
# Code purpose: run this code to ensure you are on the most recent R and RStudio interface, and that
#               all packages are installed and updated that are needed throughout the PPRV analysis code

#############################################################################################################

# Check your version of R
R.version
sessionInfo()

# Update R
library(installr)
updateR() # look for a popup here. 
# It is best to open R directly and run from there instead of from RStudio. You can say 'yes' to 
# moving your packages over from the old R installation to the new R installation and to deleting
# them from their old spot. Also yes to updating them.

# On your PC go to the search bar by the Windows logo in bottom left. Search Add or Remove Programs.
# Remove all old versions of R with the Uninstall button here. This will force RStudio to use the latest R version.

# Install latest RStudio (now called Posit)
# In RStudio go to Help -> About RStudio. Check to see if the version is the same or older than at this website
# https://posit.co/download/rstudio-desktop/
# If older, download the new, free Rstudio .exe and install.
# When you launch new RStudio, if it does not automatically load the correct, newest version of R, it should give a 
# pop up and you should confirm that it uses the correct version of R.


# Check to make sure all the packages used in the next PPRV code files are installed, and if they are not, install them
list.of.packages <- c("pkgbuild","devtools","installr", "ggplot2", "mgcv", "gridExtra", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages>0)) install.packages(new.packages)

# For packages you already have installed, see which packages have an update available
old.packages()

# Update all your installed packages (this may take a while!)
# It's ok to install to a personal library (if this popup comes up), as if the first path for the library
# does not work it is likely this is because you are not the computer admin and R is trying to install
# packages to the system path for the R library. You would have to run R and the code as an admin
# Instead of this code you could navigate to the packages tab in the packages pane and click Update
update.packages() 
#update.packages(ask = FALSE)  # this does the install without pop prompts for permissions but may mask install issues

# Run again - if 'Binaries will be installed' message appeared for some packages, they need to be installed from the source
# First you will have to install/update your RTools installation

# Checking your RTools installation
# https://search.r-project.org/CRAN/refmans/pkgbuild/html/has_rtools.html
library(pkgbuild)
find_rtools()  # Should be TRUE if it is installed
check_rtools() # Checks if version of Rtools matches version of R - it may not, giving errors in package installation above
rtools_path()  # Check PATH of where Rtools is installed

# If need to install R tools, do from here - pick the correct RTools for the version of R you have:
# https://cran.r-project.org/bin/windows/Rtools/  # look for the Rtools##installer blue hyperlink to get the .exe
# Download correct version of .exe for your R and install

# Attempt to install packages again from source
old.packages()
old <- old.packages()
install.packages(old, type = "source")

