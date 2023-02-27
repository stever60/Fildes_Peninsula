
# BChron KGI version 2.1 - for Fildes and Potter papers 

#Load packages
library(Bchron)
library(rcarbon)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parameters)
library(mclust)
library(NbClust)
library(cluster)

#clear previous console
remove (list = ls())
#clear plot window
dev.off()
# Clear plots
if(!is.null(dev.list())) dev.off()

# BChron - Combined probability density analysis --------------------------

# Background info - copied from package notes ----------------------------------

#The BchronDensity:Non-parametric phase model function runs a non-parametric phase model on 14C 
#and non-14C ages via Gaussian Mixture density estimation

#This model places a Gaussian mixture prior distribution on the calibrated ages and so estimates the density of 
#the overall set of radiocarbon ages. It is designed to be a probabilistic version of the Oxcal SUM command which 
#takes calibrated ages and sums the probability distributions with the aim of estimating activity through age as 
#a proxy.

#An object of class BchronDensityRun with the following elements is produced:
#theta = The posterior samples of the restricted ages
#p = Posterior samples of the mixture proportions
#mu = Values of the means of each Gaussian mixture
#calAges = The calibrated ages from BchronCalibrate
#G = The number of mixture components. Equal to numMix
#age_grid = A grid of ages used for the final density estimate
#density = The density estimate based on the above age grid


# SECTION 1: import datasets ------------------------------

#set working directory on mac
setwd("/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data")
#check working directory
getwd()

# Read the SSI database & rearrange in base R --------
SSI_df <- read.csv("SSI_db.csv")

# Datasets for density phase analysis from SSI_df- run first time only -------------------------------------------------------

# x2 - Fildes paper - new Artigas Beach C14 data - this study
ART_Terr20_0 <- subset(SSI_df, SiteName == "Artigas_Beach_Moraine" & calCurves == "shcal20", 
                   select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Terr20_1 <- cbind(ART_Terr20_0, thickness = 1)
names(ART_Terr20_1)[names(ART_Terr20_1) == "FeatureID"] <- "id"
names(ART_Terr20_1)[names(ART_Terr20_1) == "Altitude_masl"] <- "position"
colnames(ART_Terr20_1)
BChron_column_order <- c("id", "ages", "ageSds", "position", "thickness", "calCurves", "LabID")
ART_Terr20 <- ART_Terr20_1[, BChron_column_order]
ART_Terr20
write.csv(ART_Terr20,"Inputs/x2_ART_Terr20.csv", row.names = FALSE)
x2 <- ART_Terr20

# x2.1 - Fildes paper - all Shetland I Moraine terrestrial C14 data - this study & Hall 2007
ART_Hall_Terr20_0 <- subset(SSI_df, SiteName == "Artigas_Beach_Moraine" & 
                              calCurves == "shcal20",
                            select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Hall_Terr20_1 <- subset(SSI_df, SiteName == "Valle_Norte_Moraines" &
                            calCurves == "shcal20",
                            select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Hall_Terr20_2 <- subset(SSI_df, SiteName == "Valle_Klotz_Moraine" &
                            calCurves == "shcal20",
                            select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Hall_Terr20_3 <- rbind(ART_Hall_Terr20_0, ART_Hall_Terr20_1, ART_Hall_Terr20_2)
ART_Hall_Terr20_4 <- cbind(ART_Hall_Terr20_3, thickness = 1)
names(ART_Hall_Terr20_4)[names(ART_Hall_Terr20_4) == "FeatureID"] <- "id"
names(ART_Hall_Terr20_4)[names(ART_Hall_Terr20_4) == "Altitude_masl"] <- "position"
colnames(ART_Hall_Terr20_4)
BChron_column_order <- c("id", "ages", "ageSds", "position", "thickness", "calCurves", "LabID")
ART_Hall_Terr20 <- ART_Hall_Terr20_4[, BChron_column_order]
ART_Hall_Terr20
write.csv(ART_Hall_Terr20,"Inputs/x2.1_ART_Hall_Terr20.csv", row.names = FALSE)
x2.1 <- ART_Hall_Terr20

# x4 - Fildes paper - all Shetland I Moraine marine C14 data - this study and Hall 2007 <20ka
ART_Hall_Mar20_0 <- subset(SSI_df, SiteName == "Artigas_Beach_Moraine" & 
                              calCurves == "marine20666" & ages <20000,
                            select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Hall_Mar20_1 <- subset(SSI_df, SiteName == "Valle_Norte_Moraines" &
                              calCurves == "marine20666" & ages <20000,
                            select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Hall_Mar20_2 <- subset(SSI_df, SiteName == "Valle_Klotz_Moraine" &
                              calCurves == "marine20666" & ages <20000,
                            select=c(FeatureID, ages, ageSds, Altitude_masl, calCurves, LabID))
ART_Hall_Mar20_3 <- rbind(ART_Hall_Mar20_0, ART_Hall_Mar20_1, ART_Hall_Mar20_2)
ART_Hall_Mar20_4 <- cbind(ART_Hall_Mar20_3, thickness = 1)
names(ART_Hall_Mar20_4)[names(ART_Hall_Mar20_4) == "FeatureID"] <- "id"
names(ART_Hall_Mar20_4)[names(ART_Hall_Mar20_4) == "Altitude_masl"] <- "position"
colnames(ART_Hall_Mar20_4)
BChron_column_order <- c("id", "ages", "ageSds", "position", "thickness", "calCurves", "LabID")
ART_Hall_Mar20666 <- ART_Hall_Mar20_4[, BChron_column_order]
ART_Hall_Mar20666
write.csv(ART_Hall_Mar20666,"Inputs/x4_ART_Hall_Mar20666.csv", row.names = FALSE)
x4 <- ART_Hall_Mar20666


# Read in csv input files - fast start  -----------------------------------------

# Datasets for desnity phase analysis - fast start

# Fildes paper - new Artigas Beach C14 data - this study
x2 <- read.csv("Inputs/x2_ART_Terr20.csv")

# Fildes paper - all Shetland I Moraine terrestrial C14 data - this study & Hall 2007
x2.1 <- read.csv("Inputs/x2.1_ART_Hall_Terr20.csv")

# Fildes paper - all Shetland I Moraine marine C14 data - this study and Hall 2007 <20ka
x4 <- read.csv("Inputs/x4_ART_Mar20666.csv")

# Advance - max age constraints
x6 <- read.csv("Inputs/x6_Advance_Potter.csv")
x7 <- read.csv("Inputs/x7_Advance_Fildes.csv")
x8 <- read.csv("Inputs/x8_Advance_KGI.csv")

# Retreat - min age constraints
x9 <- read.csv("Inputs/x9_Retreat_KGI.csv")
x10 <- read.csv("Inputs/x10_Retreat_SSI.csv")
x15 <- read.csv("Inputs/x15_Lakes_basal_Fildes.csv")
x9.4 <- read.csv("Inputs/x9.4_Retreat_KGI_earlyHolocene_cosmo.csv")
x9.5 <- read.csv("Inputs/x9.5_Retreat_KGI_earlyHolocene_lakes.csv")

# Retreat - cosmogenic min age constraints
x13 <- read.csv("Inputs/x13_Cosmo_Fildes.csv")
v <- c('ages')  ## create vector of column names to recaluculate as cal ka BP
x13[v] <- x13[v] - (2010-1950)  # convert to BP by subtracting 1950 from samples collection date and assign back
x13.1 <- read.csv("Inputs/x13.1_Cosmo_Potter.csv")
x13.1[v] <- x13.1[v] - (2010-1950) # convert to BP
x14 <- read.csv("Inputs/x14_Cosmo_KGI.csv")
x14[v] <- x14[v] - (2010-1950) # convert to BP

# Warm conditions similar to present - Sub-aquatic moss & GDGT Temp - & tephra layers
x5 <- read.csv("Inputs/x5_Aq_moss_Fildes.csv") #without post-bomb ages
x5 <- subset(x5, calCurves=="shcal20", 
             select=c(id:altitude)) # remove post-bomb (normal) ages for density phase analysis as prob density not comparable to sh20 data & post-bomb data can't be calibrated in BChron and RCarbon
x11 <- read.csv("Inputs/x11_Aq_Moss_KGI.csv")
x12 <- read.csv("Inputs/x12_Aq_Moss_SSI.csv")
x16 <- read.csv("Inputs/x16_Aq_Moss_Fildes.csv") #includes post-bomb ages
x17 <- read.csv("Inputs/x17_Yanou_GDGT.csv")

# Published AP -regional - global data - density phase analysis
x20 <- read.csv("Inputs/x20_Kaplan2020_JRI_Advance.csv")
x20[v] <- x20[v] - (2010-1950) #convert to BP from sampling date (nearest 10 years) and assign back
x21 <- read.csv("Inputs/x21_Kaplan2020_JRI_Retreat.csv")
x21[v] <- x21[v] - (2010-1950) #convert to BP from sampling date (nearest 10 years) and assign back

# Calibration and density plotting --------------------------------------------------

# Density function runs a non-parametric phase model on 14C and non-14C ages via Gaussian Mixture density estimation
# numMix - number of mixture components - set to number of datapoints rounded to nearest 10 to avoid overfitting

# Fildes paper - new and published Shetland I Moraine data
xAges2 <- with (x2, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges2.1 <- with (x2.1, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges4 <- with (x4, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))

xDens2 <- with(ART_Terr20, BchronDensity(ages = ages,ageSds = ageSds,calCurves = calCurves, numMix = 10)) #final plot
xDens2.1 <- with(x2.1, BchronDensity(ages = ages,ageSds = ageSds,calCurves = calCurves, numMix = 100)) #final plot
xDens4 <- with(x4, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 5))

# Advance - max age constraints
xAges6 <- with (x6, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges7 <- with (x7, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges8 <- with (x8, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))

xDens6 <- with(x6, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 30)) #final plot
xDens7 <- with(x7, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 50)) #final plot
xDens8 <- with(x8, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 80)) #final plot

# Retreat - C14 min age constraints
xAges9 <- with (x9, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges10 <- with (x10, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges15 <- with (x15, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))

xDens9 <- with(x9, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 30)) #final table
xDens10 <- with(x10, BchronDensity(ages = ages, ageSds = ageSds, calCurves = calCurves, numMix = 40))
xDens15 <- with(x15, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 20)) #final plot

# Retreat - cosmogenic min age constraints
xAges13 <- with (x13, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges13.1 <- with (x13b, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges14 <- with (x14, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))

xDens13 <- with(x13, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 5)) #final plot
xDens13.1 <- with(x13b, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 5)) #final plot
xDens14 <- with(x14, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 10)) #final plot

# Warm conditions & tephra layers
xAges5 <- with (x5, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges11 <- with (x11, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges12 <- with (x12, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges16 <- with (x16, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges17 <- with (x17, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))

xDens5 <- with(x5, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 50)) #final plot
xDens11 <- with(x11, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 20)) #final plot
xDens12 <- with(x12, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 80))
xDens16 <- with(x16, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 50)) #final table
xDens17 <- with(x17a, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 10)) #final table

# Published AP data
xAges20 <- with (x20, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xAges21 <- with (x21, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))

xDens20 <- with(x20, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 30))
xDens21 <- with(x21, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 20))

# Plot age distribution for each group ------------------------------------
plot(xAges2, withPositions = TRUE)
plot(xAges2.1, withPositions = TRUE)
plot(xAges4, withPositions = TRUE)
plot(xAges6, withPositions = TRUE)
plot(xAges7, withPositions = TRUE)
plot(xAges8, withPositions = TRUE)
plot(xAges9, withPositions = TRUE)
plot(xAges10, withPositions = TRUE)
plot(xAges11, withPositions = TRUE)
plot(xAges12, withPositions = TRUE)
plot(xAges13, withPositions = TRUE)
plot(xAges13.1, withPositions = TRUE)
plot(xAges14, withPositions = TRUE)
plot(xAges15, withPositions = TRUE)
plot(xAges16, withPositions = TRUE)
plot(xAges17, withPositions = TRUE)
plot(xAges20, withPositions = TRUE)
plot(xAges21, withPositions = TRUE)
# Summary of highest density region for all input data --------------------
xAges2_sum <- summary (ART_Terr20_Ages, prob = 95.4)
xAges2.1_sum <- summary (xAges2.1, prob = 95.4)
xAges4_sum <- summary (xAges4, prob = 95.4)
xAges6_sum <- summary (xAges6, prob = 95.4)
xAges7_sum <- summary (xAges7, prob = 95.4)
xAges8_sum <- summary (xAges8, prob = 95.4)
xAges9_sum <- summary (xAges9, prob = 95.4)
xAges10_sum <- summary (xAges10, prob = 95.4)
xAges11_sum <- summary (xAges11, prob = 95.4)
xAges12_sum <- summary (xAges12, prob = 95.4)
xAges13_sum <- summary (xAges13, prob = 95.4)
xAges13.1_sum <- summary (xAges13b, prob = 95.4)
xAges14_sum <- summary (xAges14, prob = 95.4)
xAges15_sum <- summary (xAges15, prob = 95.4)
xAges16_sum <- summary (xAges16, prob = 95.4)
xAges17_sum <- summary (xAges17, prob = 95.4)
xAges20_sum <- summary (xAges20, prob = 95.4)
xAges21_sum <- summary (xAges21, prob = 95.4)
# Bind density data and write to file - takes ages! -----------------------

#xDens2_output <- merge(xDens2$ageGrid,xDens2$densities)
#xDens2.1_output <- merge(xDens2.1$ageGrid,xDens2.1$densities)
#xDens4_output <- merge(xDens4$ageGrid,xDens4$densities)
#xDens6_output <- merge(xDens6$ageGrid,xDens6$densities)
#xDens7_output <- merge(xDens7$ageGrid,xDens7$densities)
#xDens8_output <- merge(xDens8$ageGrid,xDens8$densities)
#xDens9_output <- merge(xDens9$ageGrid,xDens9$densities)
#xDens10_output <- merge(xDens10$ageGrid,xDens10$densities)
#xDens11_output <- merge(xDens11$ageGrid,xDens11$densities)
#xDens12_output <- merge(xDens12$ageGrid,xDens12$densities)
#xDens13_output <- merge(xDens13$ageGrid,xDens13$densities)
#xDens13.1_output <- merge(xDens13$ageGrid,xDens13$densities)
#xDens14_output <- merge(xDens14$ageGrid,xDens14$densities)
#xDens15_output <- merge(xDens15$ageGrid,xDens15$densities)
#xDens16_output <- merge(xDens16$ageGrid,xDens16$densities)
#xDens17_output <- merge(xDens17$ageGrid,xDens17$densities)

# Write to file for plotting in other programs
#write.csv(ART_Terr20_Dens_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/ART_Terr20_densities.csv", row.names = FALSE)
#write.csv(xDens2.1_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens2.1_densities.csv", row.names = FALSE)
#write.csv(xDens4_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens4_densities.csv", row.names = FALSE)
#write.csv(xDens6_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens6_densities.csv", row.names = FALSE)
#write.csv(xDens7_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens7_densities.csv", row.names = FALSE)
#write.csv(xDens8_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens8_densities.csv", row.names = FALSE)
#write.csv(xDens9_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens9_densities.csv", row.names = FALSE)
#write.csv(xDens10_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens10_densities.csv", row.names = FALSE)
#write.csv(xDens11_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens11_densities.csv", row.names = FALSE)
#write.csv(xDens12_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens12_densities.csv", row.names = FALSE)
#write.csv(xDens13_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens13_densities.csv", row.names = FALSE)
#write.csv(xDens13.1_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens13_densities.csv", row.names = FALSE)
#write.csv(xDens14_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens14_densities.csv", row.names = FALSE)
#write.csv(xDens15_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens15_densities.csv", row.names = FALSE)
#write.csv(xDens16_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens16_densities.csv", row.names = FALSE)
#write.csv(xDens17_output,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputs/xDens17_densities.csv", row.names = FALSE)

# Summary of phases at prob = 0.95 and 0.68 ----------------------------------------
summary(xDens2, type = "outliers", prob = 0.95) # Look at outlier probabilities 0.95 is the default
summary(xDens2.1, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens4, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens6, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens7, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens8, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens9, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens10, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens11, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens12, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens13, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens13.1, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens14, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens15, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens16, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens17, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens20, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens21, type = "outliers", prob = 0.95) # Look at outlier probabilities

summary(xDens6, prob = 0.95) # Look at outlier probabilities
summary(xDens7, prob = 0.95) # Look at outlier probabilities
summary(xDens8, prob = 0.95) # Look at outlier probabilities

#Summary of phases at prob = 0.68
summary(xDens2, type = "outliers", prob = 0.68) # Look at outlier probabilities 0.95 is the default
summary(xDens2.1, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens4, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens6, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens7, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens8, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens9, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens10, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens11, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens12, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens13, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens13.1, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens14, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens15, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens16, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens17, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens20, type = "outliers", prob = 0.68) # Look at outlier probabilities
summary(xDens21, type = "outliers", prob = 0.68) # Look at outlier probabilities




# SECTION 2 - RCarbon data and plots -------------------------------------------------

# Fildes Advance constraints - RCarbon for SPD and KDE plots to compare with B --------

r7 <- read.csv("Inputs/r7_Advance_KGI_RCarbon.csv")
head(r7)
#pooled dates from the same event prior to calibration into a weighted mean and std err
#to check for internal consistency before calibration
r7_C14 = r7$CRA
errors_r7 = r7$Error
id_r7 = r7$LabID
pool_r7 <- poolDates(r7_C14,errors_r7,id_r7)
pool_r7

#Calibration
r7_cal <- calibrate(r7$CRA, r7$Error, normalised=TRUE, calCurves=r7$calCurves,
                    resOffsets=r7$resOffsets,resErrors=r7$resErrors,ids=r7$LabID)
head(r7_cal)

#Output summary stats of prob density distributions for each
r7_sum <- summary(r7_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
r7_sum <- rename(r7_sum, LabID = DateID)
r7_cal_sum <- left_join(r7, r7_sum, by = c("LabID"))
r7_cal_sum

#Write summary stats for CRA and median Cal BP to file
r7_C14summary <- describe_distribution(r7_cal_sum$CRA)
r7_Calsummary  <- describe_distribution(r7_cal_sum$MedianBP)
r7_C14summary
r7_Calsummary 
write.csv(r7_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r7_C14summary.csv", row.names = FALSE)
write.csv(r7_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r7_Calsummary.csv", row.names = FALSE)
write.csv(r7_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r7_cal_sum.csv", row.names = FALSE)

#++++++ Investigate subsets +++++++
# Subsets of all dates that have a probability above 0.5 and are <6000 BP
r7_subset1 = which.CalDates(r7_cal,BP<6000,p=0.5)
r7_subset1
# Subsets of all dates that have a probability above 0.5 and are between 1000 and 2000 BP
r7_subset2 = which.CalDates(r7_cal,BP>1000&BP<2000,p=0.5)
r7_subset2

## ++++++++ Combine calibrated ages into a SPD - Summed Prob. Dist +++++++++
# Plot distributions based on Site Name column as the group 
# Use class 'CalDates' - first create a bin based on site name to group by later on 
# first need to use cutree(tree, k = NULL, h = NULL) to look at and define groups according to height (h)
# Arguments tree - a cluster tree/dendrogram produced by hclust - allows you to look and see if 
# similar ages come from the same site and if there us similarity between the sites 
# cutree() only expects a list with components merge, height, and labels, of appropriate content each.k = an integer scalar or vector with the desired number of groups
# h = numeric scalar or vector with heights where the tree should be cut.

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using hierachical cluster analysis 
#This clustering method defines the cluster distance between two clusters to
# be the maximum distance between their individual components
r7_dd <- dist(scale(r7$CRA), method = "euclidean")
r7_hc <- hclust(r7_dd, method = "ward.D2")
cutree(r7_hc, h = 100)
as.dendrogram(r7_hc)
plot(r7_hc, labels=r7$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(r7_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(r7$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
r7_bins <- binPrep(r7$SiteName, r7$CRA, h=100)
summary(r7_cal, prob = 0.954, calendar = "BP")
r7_bins

#Compute median date for each bin
r7_bm <- binMed(x=r7_cal,bins=r7_bins)
r7_bm
#Compute median date for each date
r7_dm <- medCal(r7_cal)
r7_dm

#set up Age scale 
r7_timeRange <- c(12000,-100)
r7_revtimeRange <- c(-100,12000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
r7_spd <- spd(r7_cal, bins=r7_bins, runm=50, timeRange=r7_timeRange) #runm is running mean 
plot(r7_spd,runm=50, xlim = r7_revtimeRange) 
#medians for each date
barCodes(r7_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(r7_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
r7_res = stackspd(x=r7_cal,timeRange=r7_timeRange,bins=r7_bins,group=r7$SiteName)
plot(r7_res,type='lines', xlim = r7_revtimeRange)
#medians for each bin
barCodes(r7_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot
plot(r7_res,type='stacked', xlim = r7_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion') #not that useful
plot(r7_res,type='multipanel', xlim = r7_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(r7_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = r7_revtimeRange)
multiplot(r7_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = r7_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages
r7_s = sampleDates(r7_cal,bins=r7_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(r7_s,timeRange=r7_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',interval = 0.95, xlim = r7_revtimeRange, 
     main = "Advance KGI: Composite Kernel Density Estimate (Quantile Int. = 0.95)")
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.01,0.0001), labels=rep("",101), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.01,0.0001), labels=rep("",101), tck=0.04))
par(new = FALSE) #remove hold on the plot frame to add median values

#set up for plotting later on and write to csv for use in other programs
r7_p1 <- plot(r7_spd,runm=50, xlim = r7_revtimeRange)
write.csv(r7_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/r7_spd.csv", row.names = FALSE)

# Aquatic moss Fildes - RCarbon - for SPD and KDE -------------------------

r5 <- read.csv("Inputs/r5_Aq_moss_Fildes_RCarbon.csv")
head(r5)
#pooled dates from the same event prior to calibration into a weighted mean and std err
#to check for internal consistency before calibration
r5_C14 = r5$CRA
errors_r5 = r5$Error
id_r5 = r5$LabID
pool_r5 <- poolDates(r5_C14,errors_r5,id_r5)
pool_r5

#Calibration
r5_cal <- calibrate(r5$CRA, r5$Error, normalised=TRUE, calCurves=r5$calCurves,
                    resOffsets=r5$resOffsets,resErrors=r5$resErrors,ids=r5$LabID)
head(r5_cal)

#Output summary stats of prob density distributions for each
r5_sum <- summary(r5_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
r5_sum <- rename(r5_sum, LabID = DateID)
r5_cal_sum <- left_join(r5, r5_sum, by = c("LabID"))
r5_cal_sum

#Write summary stats for CRA and median Cal BP to file
r5_C14summary <- describe_distribution(r5_cal_sum$CRA)
r5_Calsummary  <- describe_distribution(r5_cal_sum$MedianBP)
r5_C14summary
r5_Calsummary 
write.csv(r5_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data/Outputs/r5_C14summary.csv", row.names = FALSE)
write.csv(r5_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data/Outputs/r5_Calsummary.csv", row.names = FALSE)
write.csv(r5_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data/Outputs/r5_cal_sum.csv", row.names = FALSE)

#++++++ Investigate subsets +++++++
# Subsets of all dates that have a probability above 0.5 and are <6000 BP
r5_subset1 = which.CalDates(r5_cal,BP<6000,p=0.5)
r5_subset1
# Subsets of all dates that have a probability above 0.5 and are between 1000 and 2000 BP
r5_subset2 = which.CalDates(r5_cal,BP>1000&BP<2000,p=0.5)
r5_subset2

## ++++++++ Combine calibrated ages into a SPD - Summed Prob. Dist +++++++++
# Plot distributions based on Site Name column as the group 
# Use class 'CalDates' - first create a bin based on site name to group by later on 
# first need to use cutree(tree, k = NULL, h = NULL) to look at and define groups according to height (h)
# Arguments tree - a cluster tree/dendrogram produced by hclust - allows you to look and see if 
# similar ages come from the same site and if there us similarity between the sites 
# cutree() only expects a list with components merge, height, and labels, of appropriate content each.k = an integer scalar or vector with the desired number of groups
# h = numeric scalar or vector with heights where the tree should be cut.

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using hierachical cluster analysis 
#This clustering method defines the cluster distance between two clusters to
# be the maximum distance between their individual components
r5_dd <- dist(scale(r5$CRA), method = "euclidean")
r5_hc <- hclust(r5_dd, method = "ward.D2")
cutree(r5_hc, h = 100)
as.dendrogram(r5_hc)
plot(r5_hc, labels=r5$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(r5_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(r5$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
r5_bins <- binPrep(r5$SiteName, r5$CRA, h=100)
summary(r5_cal, prob = 0.954, calendar = "BP")
r5_bins

#Compute median date for each bin
r5_bm <- binMed(x=r5_cal,bins=r5_bins)
r5_bm
#Compute median date for each date
r5_dm <- medCal(r5_cal)
r5_dm

#set up Age scale 
r5_timeRange <- c(12000,-100)
r5_revtimeRange <- c(-100,12000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
r5_spd <- spd(r5_cal, bins=r5_bins, runm=50, timeRange=r5_timeRange) #runm is running mean 
plot(r5_spd,runm=50, xlim = r5_revtimeRange) 
#medians for each date
barCodes(r5_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(r5_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
r5_res = stackspd(x=r5_cal,timeRange=r5_timeRange,bins=r5_bins,group=r5$SiteName)
plot(r5_res,type='lines', xlim = r5_revtimeRange)
#medians for each bin
barCodes(r5_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot
plot(r5_res,type='stacked', xlim = r5_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion') #not that useful
plot(r5_res,type='multipanel', xlim = r5_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 

# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(r5_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = r5_revtimeRange)
multiplot(r5_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = r5_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages
r5_s = sampleDates(r5_cal,bins=r5_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(r5_s,timeRange=r5_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',interval = 0.95, xlim = r5_revtimeRange, 
     main = "Advance KGI: Composite Kernel Density Estimate (Quantile Int. = 0.95)")
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.01,0.0001), labels=rep("",101), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.01,0.0001), labels=rep("",101), tck=0.04))
par(new = FALSE) #remove hold on the plot frame to add median values

#set up for plotting later on and write to csv for use in other programs
r5_p1 <- plot(r5_spd,runm=50, xlim = r5_revtimeRange)
write.csv(r5_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/r5_spd.csv", row.names = FALSE)

# Aquatic moss KGI - RCarbon - for SPD and KDE -------------------------

r11 <- read.csv("Inputs/r11_Aq_moss_KGI_RCarbon.csv")
head(r11)
tail(r11)
#pooled dates from the same event prior to calibration into a weighted mean and std err
#to check for internal consistency before calibration
r11_C14 = r11$CRA
errors_r11 = r11$Error
id_r11 = r11$LabID
pool_r11 <- poolDates(r11_C14,errors_r11,id_r11)
pool_r11

#Calibration
r11_cal <- calibrate(r11$CRA, r11$Error, normalised=TRUE, calCurves=r11$calCurves,
                    resOffsets=r11$resOffsets,resErrors=r11$resErrors,ids=r11$LabID)
head(r11_cal)

#Output summary stats of prob density distributions for each
r11_sum <- summary(r11_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
r11_sum <- rename(r11_sum, LabID = DateID)
r11_cal_sum <- left_join(r11, r11_sum, by = c("LabID"))
r11_cal_sum

#Write summary stats for CRA and median Cal BP to file
r11_C14summary <- describe_distribution(r11_cal_sum$CRA)
r11_Calsummary  <- describe_distribution(r11_cal_sum$MedianBP)
r11_C14summary
r11_Calsummary 
write.csv(r11_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r11_C14summary.csv", row.names = FALSE)
write.csv(r11_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r11_Calsummary.csv", row.names = FALSE)
write.csv(r11_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r11_cal_sum.csv", row.names = FALSE)

#++++++ Investigate subsets +++++++
# Subsets of all dates that have a probability above 0.5 and are <6000 BP
r11_subset1 = which.CalDates(r11_cal,BP<6000,p=0.5)
r11_subset1
# Subsets of all dates that have a probability above 0.5 and are between 1000 and 2000 BP
r11_subset2 = which.CalDates(r11_cal,BP>1000&BP<2000,p=0.5)
r11_subset2

## ++++++++ Combine calibrated ages into a SPD - Summed Prob. Dist +++++++++
# Plot distributions based on Site Name column as the group 
# Use class 'CalDates' - first create a bin based on site name to group by later on 
# first need to use cutree(tree, k = NULL, h = NULL) to look at and define groups according to height (h)
# Arguments tree - a cluster tree/dendrogram produced by hclust - allows you to look and see if 
# similar ages come from the same site and if there us similarity between the sites 
# cutree() only expects a list with components merge, height, and labels, of appropriate content each.k = an integer scalar or vector with the desired number of groups
# h = numeric scalar or vector with heights where the tree should be cut.

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using hierachical cluster analysis 
#This clustering method defines the cluster distance between two clusters to
# be the maximum distance between their individual components
r11_dd <- dist(scale(r11$CRA), method = "euclidean")
r11_hc <- hclust(r11_dd, method = "ward.D2")
cutree(r11_hc, h = 100)
as.dendrogram(r11_hc)
plot(r11_hc, labels=r11$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(r11_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(r11$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
r11_bins <- binPrep(r11$SiteName, r11$CRA, h=100)
summary(r11_cal, prob = 0.954, calendar = "BP")
r11_bins

#Compute median date for each bin
r11_bm <- binMed(x=r11_cal,bins=r11_bins)
r11_bm
#Compute median date for each date
r11_dm <- medCal(r11_cal)
r11_dm

#set up Age scale 
r11_timeRange <- c(12000,-100)
r11_revtimeRange <- c(-100,12000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
r11_spd <- spd(r11_cal, bins=r11_bins, runm=50, timeRange=r11_timeRange) #runm is running mean 
plot(r11_spd,runm=50, xlim = r11_revtimeRange) 
#medians for each date
barCodes(r11_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(r11_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
r11_res = stackspd(x=r11_cal,timeRange=r11_timeRange,bins=r11_bins,group=r11$SiteName)
plot(r11_res,type='lines', xlim = r11_revtimeRange)
#medians for each bin
barCodes(r11_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot
plot(r11_res,type='stacked', xlim = r11_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion') #not that useful
plot(r11_res,type='multipanel', xlim = r11_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 

# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(r11_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = r11_revtimeRange)
multiplot(r11_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = r11_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages
r11_s = sampleDates(r11_cal,bins=r11_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(r11_s,timeRange=r11_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',interval = 0.95, xlim = r11_revtimeRange, 
     main = "Advance KGI: Composite Kernel Density Estimate (Quantile Int. = 0.95)")
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.01,0.0001), labels=rep("",101), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.01,0.0001), labels=rep("",101), tck=0.04))
par(new = FALSE) #remove hold on the plot frame to add median values

#set up for plotting later on and write to csv for use in other programs
r11_p1 <- plot(r11_spd,runm=50, xlim = r11_revtimeRange)
write.csv(r11_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/r11_spd.csv", row.names = FALSE)

# Fildes Lakes basal ages - RCarbon - for SPD and KDE  --------------------

r15 <- read.csv("Inputs/r15_Fildes_lakes_RCarbon.csv")
head(r15)
#pooled dates from the same event prior to calibration into a weighted mean and std err
#to check for internal consistency before calibration
r15_C14 = r15$CRA
errors_r15 = r15$Error
id_r15 = r15$LabID
pool_r15 <- poolDates(r15_C14,errors_r15,id_r15)
pool_r15

#Calibration
r15_cal <- calibrate(r15$CRA, r15$Error, normalised=TRUE, calCurves=r15$calCurves,
                    resOffsets=r15$resOffsets,resErrors=r15$resErrors,ids=r15$LabID)
head(r15_cal)

#Output summary stats of prob density distributions for each
r15_sum <- summary(r15_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
r15_sum <- rename(r15_sum, LabID = DateID)
r15_cal_sum <- left_join(r15, r15_sum, by = c("LabID"))
r15_cal_sum

#Write summary stats for CRA and median Cal BP to file
r15_C14summary <- describe_distribution(r15_cal_sum$CRA)
r15_Calsummary  <- describe_distribution(r15_cal_sum$MedianBP)
r15_C14summary
r15_Calsummary 
write.csv(r15_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r15_C14summary.csv", row.names = FALSE)
write.csv(r15_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r15_Calsummary.csv", row.names = FALSE)
write.csv(r15_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r15_cal_sum.csv", row.names = FALSE)

#++++++ Investigate subsets +++++++
# Subsets of all dates that have a probability above 0.5 and are <6000 BP
r15_subset1 = which.CalDates(r15_cal,BP<6000,p=0.5)
r15_subset1
# Subsets of all dates that have a probability above 0.5 and are between 1000 and 2000 BP
r15_subset2 = which.CalDates(r15_cal,BP>1000&BP<2000,p=0.5)
r15_subset2

## ++++++++ Combine calibrated ages into a SPD - Summed Prob. Dist +++++++++
# Plot distributions based on Site Name column as the group 
# Use class 'CalDates' - first create a bin based on site name to group by later on 
# first need to use cutree(tree, k = NULL, h = NULL) to look at and define groups according to height (h)
# Arguments tree - a cluster tree/dendrogram produced by hclust - allows you to look and see if 
# similar ages come from the same site and if there us similarity between the sites 
# cutree() only expects a list with components merge, height, and labels, of appropriate content each.k = an integer scalar or vector with the desired number of groups
# h = numeric scalar or vector with heights where the tree should be cut.

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using hierachical cluster analysis 
#This clustering method defines the cluster distance between two clusters to
# be the maximum distance between their individual components
r15_dd <- dist(scale(r15$CRA), method = "euclidean")
r15_hc <- hclust(r15_dd, method = "ward.D2")
cutree(r15_hc, h = 100)
as.dendrogram(r15_hc)
plot(r15_hc, labels=r15$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(r15_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(r15$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
r15_bins <- binPrep(r15$SiteName, r15$CRA, h=100)
summary(r15_cal, prob = 0.954, calendar = "BP")
r15_bins

#Compute median date for each bin
r15_bm <- binMed(x=r15_cal,bins=r15_bins)
r15_bm
#Compute median date for each date
r15_dm <- medCal(r15_cal)
r15_dm

#set up Age scale 
r15_timeRange <- c(12000,-100)
r15_revtimeRange <- c(-100,12000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
r15_spd <- spd(r15_cal, bins=r15_bins, runm=50, timeRange=r15_timeRange) #runm is running mean 
plot(r15_spd,runm=50, xlim = r15_revtimeRange) 
#medians for each date
barCodes(r15_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(r15_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
r15_res = stackspd(x=r15_cal,timeRange=r15_timeRange,bins=r15_bins,group=r15$SiteName)
plot(r15_res,type='lines', xlim = r15_revtimeRange)
#medians for each bin
barCodes(r15_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot
plot(r15_res,type='stacked', xlim = r15_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion') #not that useful
plot(r15_res,type='multipanel', xlim = r15_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 


# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(r15_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = r15_revtimeRange)
multiplot(r15_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = r15_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages
r15_s = sampleDates(r15_cal,bins=r15_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(r15_s,timeRange=r15_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',interval = 0.95, xlim = r15_revtimeRange, 
     main = "Advance KGI: Composite Kernel Density Estimate (Quantile Int. = 0.95)")
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.01,0.0001), labels=rep("",101), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.01,0.0001), labels=rep("",101), tck=0.04))
par(new = FALSE) #remove hold on the plot frame to add median values

#set up for plotting later on and write to csv for use in other programs
r15_p1 <- plot(r15_spd,runm=50, xlim = r15_revtimeRange)
write.csv(r15_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/r15_spd.csv", row.names = FALSE)




# Potter Advance constraints - RCarbon for SPD and KDE plots --------------------------

r6 <- read.csv("Inputs/r6_Potter_Terr_Holocene.csv")
head(r6)
#pooled dates from the same event prior to calibration into a weighted mean and std err
#to check for internal consistency before calibration
r6_C14 = r6$CRA
errors_r6 = r6$Error
id_r6 = r6$LabID
pool_r6 <- poolDates(r6_C14,errors_r6,id_r6)
pool_r6

#Calibration
r6_cal <- calibrate(r6$CRA, r6$Error, normalised=TRUE, calCurves=r6$calCurves,
                    resOffsets=r6$resOffsets,resErrors=r6$resErrors,ids=r6$LabID)
head(r6_cal)

#Output summary stats of prob density distributions for each
r6_sum <- summary(r6_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
r6_sum <- rename(r6_sum, LabID = DateID)
r6_cal_sum <- left_join(r6, r6_sum, by = c("LabID"))
r6_cal_sum

#Write summary stats for CRA and median Cal BP to file
r6_C14summary <- describe_distribution(r6_cal_sum$CRA)
r6_Calsummary  <- describe_distribution(r6_cal_sum$MedianBP)
r6_C14summary
r6_Calsummary 
write.csv(r6_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r6_C14summary.csv", row.names = FALSE)
write.csv(r6_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r6_Calsummary.csv", row.names = FALSE)
write.csv(r6_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/BChron/Data/Outputss/r6_cal_sum.csv", row.names = FALSE)

#++++++ Investigate subsets +++++++
# Subsets of all dates that have a probability above 0.5 and are <6000 BP
r6_subset1 = which.CalDates(r6_cal,BP<6000,p=0.5)
r6_subset1
# Subsets of all dates that have a probability above 0.5 and are between 1000 and 2000 BP
r6_subset2 = which.CalDates(r6_cal,BP>1000&BP<2000,p=0.5)
r6_subset2

## ++++++++ Combine calibrated ages into a SPD - Summed Prob. Dist +++++++++
# Plot distributions based on Site Name column as the group 
# Use class 'CalDates' - first create a bin based on site name to group by later on 
# first need to use cutree(tree, k = NULL, h = NULL) to look at and define groups according to height (h)
# Arguments tree - a cluster tree/dendrogram produced by hclust - allows you to look and see if 
# similar ages come from the same site and if there us similarity between the sites 
# cutree() only expects a list with components merge, height, and labels, of appropriate content each.k = an integer scalar or vector with the desired number of groups
# h = numeric scalar or vector with heights where the tree should be cut.

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using hierachical cluster analysis 
#This clustering method defines the cluster distance between two clusters to
# be the maximum distance between their individual components
r6_dd <- dist(scale(r6$CRA), method = "euclidean")
r6_hc <- hclust(r6_dd, method = "ward.D2")
cutree(r6_hc, h = 100)
as.dendrogram(r6_hc)
plot(r6_hc, labels=r6$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(r6_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(r6$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
r6_bins <- binPrep(r6$SiteName, r6$CRA, h=100)
summary(r6_cal, prob = 0.954, calendar = "BP")
r6_bins

#Compute median date for each bin
r6_bm <- binMed(x=r6_cal,bins=r6_bins)
r6_bm
#Compute median date for each date
r6_dm <- medCal(r6_cal)
r6_dm

#set up Age scale 
r6_timeRange <- c(12000,-100)
r6_revtimeRange <- c(-100,12000)
r6_timeRange_18ka <- c(18000,-100)
r6_revtimeRange_18ka <- c(-100,18000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
r6_spd <- spd(r6_cal, bins=r6_bins, runm=50, timeRange=r6_timeRange) #runm is running mean 
plot(r6_spd,runm=50, xlim = r6_revtimeRange) 
#medians for each date
barCodes(r6_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(r6_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
r6_res = stackspd(x=r6_cal,timeRange=r6_timeRange,bins=r6_bins,group=r6$SiteName)
plot(r6_res,type='lines', xlim = r6_revtimeRange)
#medians for each bin
barCodes(r6_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

# STACKED PLOTS
plot(r6_res,type='stacked', xlim = r6_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion')
plot(r6_res,type='multipanel', xlim = r6_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 

# Clear plots
if(!is.null(dev.list())) dev.off()
#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(r6_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = r6_revtimeRange)
multiplot(r6_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = r6_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages
r6_s = sampleDates(r6_cal,bins=r6_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(r6_s,timeRange=r6_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',interval = 0.95, xlim = r6_revtimeRange, 
     main = "Advance Potter: Composite Kernel Density Estimate (Quantile Int. = 0.95)")
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.01,0.0001), labels=rep("",101), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.01,0.0001), labels=rep("",101), tck=0.04))
par(new = FALSE) #remove hold on the plot frame to add median values

#set up for plotting later on and write to csv for use in other programs
r6_p1 <- plot(r6_spd,runm=50, xlim = r6_revtimeRange)
write.csv(r6_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/r6_spd.csv", row.names = FALSE)





# SECTION 3: Final plots for Fildes and Potter papers ---------------------------------------------------------------

# FILDES PAPER FIGURES ------------------------------------------------------------

# Figure 3 ----------------------------------------------------------------

# Clear plots
if(!is.null(dev.list())) dev.off()
## Translate 8 cm graph plot size *from* cm *to* inches:
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch <- 3.149608
# if the pdf A4 width is 8.27 inches the l and right margins are then b
# total pdf width - plot width then divide by 2 to give margin width on each size
# then -1 to for mai = 1 inch left and right
b <- ((8.27-plotinch)/2)-1
checksize <- b*2+plotinch+2
checksize #should = 8.27 inches A4 width

## FINAL PLOT TO PDF
#mai / omi = margin in inches or mar / oma = margin in lines - order is bottom, left, top, right 
#set up the page layout - A4 is 8.27 x 11.69 - closest to 
pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
layout(matrix(1:6, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.1,0.1), omi=c(0,0,0,0), pin=c(plotinch, plotinch/3), mgp=c(2,1,0), xaxs='i') # Set up internal margins
#xas = "i" tells plot device not to add 4% extra internal margin (which default plot style in R)

#  PLOT TO SCREEN
# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:5, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch, plotinch/3),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0)) 
# mai sets up internal margins around plot in inches
# pin sets up the size of the plot in inches - here its 8 cm wide by 8/3 cm high
# mgp sets: axis label distance from axis, tick labels from ticks, tick distance from the line
# xas = "i" tells plot device not to add 4% extra internal margin (which default plot style in R)
# omi sets up internal margin around output in inches
# example of math equation as text in graph
# text(7, 4, expression(bar(x) == sum(frac(x[i], n), i==1, n))) - or expression(hat(beta) == (X^t * X)^{-1} * X^t * y) 
# for text: font = 2 is bold, 3 is italic, 4 is bold italic

# Figure 3A) Artigas Terrestrial Moss in Moraines (shcal20) [this study]

# Advance - max age constraints
xAges2 <- with (x2, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, 
                                    calCurves = calCurves, allowOutside = TRUE))
xDens2 <- with(ART_Terr20, BchronDensity(ages = ages,ageSds = ageSds,calCurves = calCurves, 
                                         numMix = 10))
p2 <- plot(xDens2, main="Artigas Beach: terrestrial moss in moraines",
           #xlab="Age [a BP]", 
           xlim = c(-100, 4000), cex.main = 1,
           tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,4000,500), tck=-0.04), axis(1, seq(-100,4000,100),labels=rep("",42), tck=-0.02),
            axis(2, seq(0,0.012,0.002), labels=rep("",7), tck=-0.04), 
            axis(3, seq(0,4000,500), labels=rep("",9), tck=0.04),axis(3, seq(-100,4000,100),labels=rep("",42), tck=0.02),
            axis(4, seq(0,0.012,0.002), labels=rep("",7), tck=0.04))
text(100, 0.0065, "A", font = 2, cex = 1.5,  family = 'sans') 
par(new = FALSE) #remove hold on the plot frame to add median values

# B) Artigas-Valle Norte-Valle Klotz Terrestrial Moss in Moraines (shcal20) [this study & Hall, 2007]
xAges2.1 <- with (x2.1, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, 
                                        calCurves = calCurves, allowOutside = TRUE))
xDens2.1 <- with(x2.1, BchronDensity(ages = ages,ageSds = ageSds,calCurves = calCurves, 
                                     numMix = 100))
p2.1 <- plot(xDens2.1, main="Artigas-Valle Norte-Valle Klotz terrestrial moss in moraines",
             #xlab="Age [a BP]", 
             xlim = c(-100, 4000),  cex.main = 1,
             tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,4000,500), tck=-0.04), axis(1, seq(-100,4000,100),labels=rep("",42), tck=-0.02),
            axis(2, seq(0,0.012,0.002), labels=rep("",7), tck=-0.04), 
            axis(3, seq(0,4000,500), labels=rep("",9), tck=0.04),axis(3, seq(-100,4000,100),labels=rep("",42), tck=0.02),
            axis(4, seq(0,0.012,0.002), labels=rep("",7), tck=0.04))
text(100, 0.008, "B", font = 2, cex = 1.5,  family = 'sans') 
par(new = FALSE) #remove hold on the plot frame to add median values

# C) Fildes Peninsula: max. age constraints on glacier (re)advance  [this study, Hall, 2007, Watcham et al., 2011]
xAges7 <- with (x7, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, 
                                    calCurves = calCurves, allowOutside = TRUE))
xDens7 <- with(x7, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, 
                                 numMix = 50))
p7 <- plot(xDens7, main="Fildes Peninsula: constraints on glacier readvance",
           #xlab="Age [a BP]", 
           xlim = c(-100, 8000), cex.main = 1,
           tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,8000,1000), tck=-0.04), axis(1, seq(0,8000,500),labels=rep("",17), tck=-0.02),
            axis(2, seq(0,0.012,0.002), tck=-0.04),
            axis(3, seq(0,8000,1000), labels=rep("",9), tck=0.04), axis(3, seq(0,8000,500),labels=rep("",17), tck=0.02))
text(200, 0.009, "C", font = 2, cex = 1.5,  family = 'sans') 

#Adding new Fildes cosmo ages to plot 3  on secondary axis of altitude
par(new=TRUE)
x <- x13$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13$altitude
x.err <- x13$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "black", bg= "red", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,8000), ylim=c(0,100), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "black", bg= "red", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,8000), ylim=c(0,100), axes = FALSE)
axis(4,  seq(0,100,50), tck=-0.04, ylim=c(0,100), col = "red", col.lab = "red", col.axis = "red")
axis(4, seq(0,100,25),labels=rep("",5), tck=-0.02, ylim=c(0,100), col = "red", col.lab = "red", col.axis = "red")
mtext("Altitude (m a.s.l.)", side=4, line=2, col = "red", cex = 0.75)
par(new = FALSE) #remove hold on the plot frame to add median values

# D) KGI: max. age constraint phases for glacier readvance  [this study, HB et al., sub & Hall, 2007]
xAges8 <- with (x8, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, 
                                    calCurves = calCurves, allowOutside = TRUE))
xDens8 <- with(x8, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, 
                                 numMix = 80))
p8 <- plot(xDens8, main="KGI: constraints on glacier readvance",
           #xlab="Age [a BP]", 
           xlim = c(-100, 12000), cex.main = 1,
           tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
            axis(2, seq(0,0.012,0.002), tck=-0.04),
            axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02))
text(400, 0.009, "D", font = 2, cex = 1.5,  family = 'sans') 

#Add KGI cosmo ages to plot 5 on secondary axis of altitude
par(new=TRUE)
x <- x14$ages
#y <- c(0.004, 0.005, 0.006)
y <- x14$altitude
x.err <- x14$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 24, col= "black", bg= "red", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,12000), ylim=c(0,200), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 24, col= "black", bg= "red", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,12000), ylim=c(0,200), axes = FALSE)
axis(4,  seq(0,200,100), tck=-0.04, ylim=c(0,100), col = "red", col.lab = "red", col.axis = "red")
axis(4, seq(0,200,25),labels=rep("",9), tck=-0.02, ylim=c(0,100), col = "red", col.lab = "red", col.axis = "red")
mtext("", side=4, line=2, col = "red", cex = 0.75)

#Add Fildes cosmo ages to plot 5 on secondary axis of altitude
par(new=TRUE)
x <- x13$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13$altitude
x.err <- x13$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "black", bg= "red", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,12000), ylim=c(0,200), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "black", bg= "red", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,12000), ylim=c(0,200), axes = FALSE)
axis(4,  ylim=c(0,200), col = "red", col.lab = "red", col.axis = "red")
mtext("Altitude (m a.s.l.)", side=4, line=2, col = "red", cex = 0.75)
#par(new = FALSE) #remove hold on the plot frame to add median values

#par(mai=c(0.1,0.25,0.3,0.25), pin=c(plotinch, plotinch/3), mgp=c(3,1,0), xaxs='i') # Set up internal margins
#xas = "i" tells plot device not to add 4% extra internal margin (which default plot style in R)

# E1) Lakes Fildes - basal ages - deglaciation plot 
#par(pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i') 
xAges15 <- with (x15, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, 
                                      calCurves = calCurves, allowOutside = TRUE))
xDens15 <- with(x15, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, 
                                   numMix = 20))
p15 <- plot(xDens15, main="Fildes lakes: basal ages & aquatic moss layers ",
            xlab="Age [a BP]", xlim = c(-100, 12000), cex.main = 1,
            tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
            axis(2, seq(0,0.008,0.002), tck=-0.04),
            axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02))#,
            #axis(4, seq(0,0.008,0.002), labels=rep("",5), tck=0.04))
text(11500, 0.0065, "E", font = 2, cex = 1.5,  family = 'sans') 


# E2) RCarbon SPD Aquatic Moss stacked plot -----------------------------------------------------------

# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:2, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.25,0.25,0.25), pin=c(plotinch, plotinch),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 

r5 <- read.csv("Inputs/r5_Aq_moss_Fildes_RCarbon.csv")
r5_cal <- calibrate(r5$CRA, r5$Error, normalised=TRUE, calCurves=r5$calCurves,
                    resOffsets=r5$resOffsets,resErrors=r5$resErrors,ids=r5$LabID) #Calibration
r5_res = stackspd(x=r5_cal,timeRange=r5_timeRange,bins=r5_bins,group=r5$SiteName) #Stack SPD plot
plot(r5_res,type='stacked', xlab="",ylab="", xlim = r5_revtimeRange, legend = TRUE, xaxt = "n",yaxt = "n",
     legend.arg = NULL, axes = FALSE, main = "", tck=0)
p_axis <- c(axis(1, seq(0,12000,2000), tck=0), axis(1, seq(0,12000,2000),labels=rep("",7), tck=-0.02),
            axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.01, xlim=c(-100,12000), col = "black", col.lab = "black", col.axis = "black"),
            axis(4, seq(0,0.06,0.01), tck=-0.02, ylim=c(0,0.06), col = "black", col.lab = "black", col.axis = "black"),
            axis(4, seq(0,0.06,0.01),labels=rep("",7), tck=-0.01, ylim=c(0,0.06), col = "black", col.lab = "black", col.axis = "black"))
mtext("Summed probability", side=4, line=2, col = "black")

# Fildes aquatic moss 2sigma age-range bars at base 
x16b <- read.csv("Inputs/x16b_Aq_Moss_Fildes_Outer_calibrated.csv")
p16b <- plot(plot_y~MedianBP,data=x16b, pch = 3, col="#006400", 
             cex = 0.75, xlim = c(-100, 12000), ylim = c(2, -1),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE) #xaxt="n", yaxt="n" to remove x and y axis
par(new=TRUE)
#p_axis <- c(axis(1, labels = FALSE))
arrows(x0=x16b$MedianBP-x16b$TwoSigma_BP_lower,
       y0=x16b$plot_y,
       x1=x16b$MedianBP+x16b$TwoSigma_BP_upper,
       y1=x16b$plot_y,
       angle=90,
       code=3,
       length=0,
       col = "#006400",
       lwd=3)
par(new=TRUE)
p16b <- plot(plot_y~MedianBP,data=x16b, pch = 3, col="#006400", 
             cex = 0.75, xlim = c(-100, 12000), ylim = c(2, -1),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE) #xaxt="n", yaxt="n" to remove x and y axis

# Aq Moss Layers -  KGI - RCarbon plot -----------------------------------


#plot(r11_res,type='stacked', xlab="",ylab="", xlim = r11_revtimeRange, legend = TRUE, xaxt = "n",yaxt = "n",
#     legend.arg = NULL, axes = FALSE, main = "", tck=0)  #most useful for summaries of multiple sites / types
#p_axis <- c(axis(1, seq(0,12000,2000), tck=0), axis(1, seq(0,12000,500),labels=rep("",25), tck=0),
#            axis(2, seq(0,0.008,0.002), tck=0, labels=rep("",5)))
#par(new=TRUE)
#axis(4,  seq(0,0.06,0.02), tck=-0.04, ylim=c(0,0.06), col = "black", col.lab = "black", col.axis = "black")
#axis(4, seq(0,0.06,0.01),labels=rep("",7), tck=-0.02, ylim=c(0,0.06), col = "black", col.lab = "black", col.axis = "black")
#mtext("Summed probability", side=4, line=2, col = "black", cex = 0.75)



# or Fildes - aquatic moss - without post-bomb - desnity ph --------

#par(pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i') 
#p5 <- plot(xDens5, main="Fildes lakes: aquatic moss layers ",
#           xlab="Age [a BP]", xlim = c(-100, 12000), cex.main = 1,
#           tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
#par(new = TRUE) #hold the plot frame to add median values
#p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
#            axis(2, seq(0,0.008,0.002), tck=-0.04),
#            axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02))#,
#axis(4, seq(0,0.008,0.002), labels=rep("",5), tck=0.04))
#text(11500, 0.0065, "F", font = 2, cex = 1.5,  family = 'sans') 
#par(new = TRUE) #hold the plot frame to add median values
# add post bomb ages as horizontal black dashes
#par(new = TRUE)

# subset Fildes aquatic moss data to select only post-bomb ages 
#Fildes_AM_postbomb <- subset(x5, calCurves=="normal", 
#                                      select=c(id, ages, ageSds, LabID))
# or Fildes - aquatic moss - with post-bomb -----------------------------------------------

#par(pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i') 
#p16 <- plot(xDens16, main="Fildes Aquatic Moss",
#            xlab="Age [a BP]", xlim = c(-100, 12000), cex.main = 1,
#            tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
#par(new = TRUE) #hold the plot frame to add median values
#p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
#            #axis(2, seq(0,0.01,0.002), tck=-0.04),
#            axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02),
#            axis(4, seq(0,0.01,0.001), labels=rep("",11), tck=0.04))
#text(11250, 0.006, "F", font = 2, cex = 1.5,  family = 'sans') 
#par(new = FALSE) #remove hold on the plot frame to add median values

# or KGI aquatic moss ages -----------------------------------------------

#par(pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i') 
#p11 <- plot(xDens11, main="Fildes lakes: aquatic moss layers ",
#           xlab="Age [a BP]", xlim = c(-100, 12000), cex.main = 1,
#           tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
#par(new = TRUE) #hold the plot frame to add median values
#p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
#            axis(2, seq(0,0.008,0.002), tck=-0.04),
#            axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02))#,
#axis(4, seq(0,0.008,0.002), labels=rep("",5), tck=0.04))
#text(11500, 0.0065, "F", font = 2, cex = 1.5,  family = 'sans') 

# add post bomb ages as horizontal black dashes
#par(new = TRUE)

# subset Fildes aquatic moss data to select only post-bomb ages 
#KGI_AM_postbomb <- subset(x11, calCurves=="normal", 
#                           select=c(id, ages, ageSds, LabID))

# Run this when writing to pdf rather than screen
#dev.off() #need to include this to write to pdf file fully - will delete screen plot







# Figure 5 ----------------------------------------------------------------

library(Bchron)
#set working directory on mac
setwd("/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data")

# New Data
x8 <- read.csv("Inputs/x8_Advance_KGI.csv") # run density phase analysis once
# run once at start
xAges8 <- with (x8, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,
                                    positions = position, calCurves = calCurves, allowOutside = TRUE))
xDens8 <- with(x8, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, 
                                 numMix = 80))
x16a <- read.csv("Inputs/x16a_Aq_Moss_Fildes_calibrated.csv")
x16b <- read.csv("Inputs/x16b_Aq_Moss_Fildes_Outer_calibrated.csv")
x17 <- read.csv("Inputs/x17_Yanou_GDGT.csv")
x18 <- read.csv("Inputs/x18_Kite_diatoms.csv")
x19a <- read.csv("Inputs/x19a_YAN_A_tephra_ages.csv") #airfall only - with T3a lower error extending to T3c
x19b <- read.csv("Inputs/x19b_YAN_ARW_tephra_ages.csv") #airfall and reworked
x19c <- read.csv("Inputs/x19c_Ant2018_tephra_ages.csv") #airfall and reworked

# Published AP regional & global data
x20 <- read.csv("Inputs/x20_Kaplan2020_JRI_Advance.csv")
x20[v] <- x20[v] - (2010-1950) #convert to BP from sampling date (nearest 10 years) and assign back
xAges20 <- with (x20, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xDens20 <- with(x20, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 30))
x21 <- read.csv("Inputs/x21_Kaplan2020_JRI_Retreat.csv")
x21[v] <- x21[v] - (2010-1950) #convert to BP from sampling date (nearest 10 years) and assign back
xAges21 <- with (x21, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xDens21 <- with(x21, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 20))
x23 <- read.csv("Inputs/x23_JRI_ice_core.csv")
x24 <- read.csv("Inputs/x24_Marcott2013_30_90S.csv")
x25 <- read.csv("Inputs/x25_Kaufmann2020.csv")
x26 <- read.csv("Inputs/x26_Miliken2009.csv")
x35 <- read.csv("Inputs/x35_Vol_flux_WAIS.csv")

# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch <- 3.149608
b <- ((8.27-plotinch)/2)-1
checksize <- b*2+plotinch+2

## FINAL PLOT TO PDF
#pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
#layout(matrix(1:8, ncol=1)) # Set up layout
#par(mai=c(0.01,0.01,0.01,0.01), omi=c(0,0,0,0), pin=c(plotinch, plotinch/4), mgp=c(2,1,0), xaxs='i') # Set up internal margins

# PLOT TO SCREEN
# Clear plots
if(!is.null(dev.list())) dev.off()
layout(matrix(1:8, ncol=1)) # Set up layout
par(mai=c(0.01,0.01,0.01,0.01), pin=c(plotinch, plotinch/4),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0)) 

# A) KGI: max. age constraint phases for glacier readvance  [this study, HB et al., sub & Hall, 2007]
p8 <- plot(xDens8, xlab="", xlim = c(-100, 12000), cex.main = 1, las = 0, axes = FALSE,
           tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL), 
           xaxt="n", yaxt="n", ann=FALSE) #supress x-axis
par(new = TRUE) #hold the plot frame
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
            axis(2, seq(0,0.012,0.002), tck=-0.08, las = 1, line = 1)#,
            #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02)
            )
mtext("Density", side = 2, cex = 0.75, line = 5)
#par(new = FALSE) #remove hold

# B) Kiteschee Lake diatom plots
# DCCA plot
x18a <- na.omit(x18)
p18 <- plot(DCCA1~Age_Sh20, data = x18a, pch = 21, type = "l", col="black", bg="black", 
            cex= 1, xlim = c(-100, 12000), ylim = c(2, -0.8),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
par(new=TRUE)
plot(DCCA1~Age_Sh20, data = x18a, pch = 21, col="black", bg="black", 
     cex= 0.5, xlim = c(-100, 12000), ylim = c(2, -0.8), tck=-0.08,
     frame=FALSE, xlab="", ylab="", xaxt="n", ann=TRUE, las = 1, line = 1)
mtext("DCCA score", side = 2, cex = 0.75, line = 5)

# Add Moss ages - Fildes outer min-max calibrated ranges 
par(new=TRUE)
p16b <- plot(plot_y~MedianBP,data=x16b, pch = 3, col="#006400", 
             cex = 0.75, xlim = c(-100, 12000), ylim = c(2, -1), tck=-0.08,
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE) #xaxt="n", yaxt="n" to remove x and y axis
par(new=TRUE)
#p_axis <- c(axis(1, labels = FALSE))
arrows(x0=x16b$MedianBP-x16b$TwoSigma_BP_lower,
       y0=x16b$plot_y,
       x1=x16b$MedianBP+x16b$TwoSigma_BP_upper,
       y1=x16b$plot_y,
       angle=90,
       code=3,
       length=0,
       col = "#006400",
       lwd=3)
par(new=TRUE)
p16b <- plot(plot_y~MedianBP,data=x16b, pch = 3, col="#006400", 
             cex = 0.75, xlim = c(-100, 12000), ylim = c(2, -1),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE) #xaxt="n", yaxt="n" to remove x and y axis
# Gomph. sp. plot  
p18 <- plot(Gomphonema~Age_Sh20, data = x18, pch = 21, type = "l", col="black", bg="black", 
            cex= 1, xlim = c(-100, 12000), ylim = c(40, 0),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
par(new=TRUE)
plot(Gomphonema~Age_Sh20, data = x18, pch = 21, col="black", bg="black", 
     cex= 0.75, xlim = c(-100, 12000), ylim = c(40, 0),
     frame=FALSE, xlab="", ylab="", xaxt="n", ann=TRUE, las = 1, line = 1)
mtext("Gomph. sp.", side = 2, cex = 0.75, line = 5)

# C) Yanou GDGT and tephra plots
p17a <- plot(GDGT_Anomaly_FOSTER_2016~Age_Sh20, data = x17, pch = 21, type = "l", col="#D19C9C", bg="#D19C9C", 
             cex= 1, xlim = c(-100, 12000), ylim = c(-5, 16),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x17$mean_ANT_anomaly~age1, data = x17, 
      lty="dotted", col="#D19C9C", lwd=0.75) # add average median temp for 90-60S
par(new=TRUE)
plot(GDGT_Anomaly_FOSTER_2016~Age_Sh20, data = x17, pch = 21, col="#D19C9C", bg="#D19C9C",
     cex= 0.5, xlim = c(-100, 12000), ylim = c(-5, 16),
     frame=FALSE, ylab=NULL, xlab="",xaxt="n", yaxt="n", ann=FALSE)
lines(x17$RMSE_Max_anomaly~Age_Sh20, data = x17, 
      lty=1, col="grey", lwd=0.75) # add upper 95% CI error bounds
lines(x17$RMSE_min_anomaly~Age_Sh20, data = x17,
      lty=1, col="grey", lwd=0.75) # add lower 95% CI error bounds
par(new=TRUE)
p17b <- plot(GDGT_ANT_GLOBAL_anomaly~Age_Sh20, data = x17, pch = 21, type = "l", col="#800000", bg="#800000", 
             cex= 1, xlim = c(-100, 12000), ylim = c(-5, 16),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x17$mean_ANT_GLOBAL_anomaly~age1, data = x17, 
      lty="dotted", col="#800000", lwd=0.75) # add global data as a black line
par(new=TRUE)
plot(GDGT_ANT_GLOBAL_anomaly~Age_Sh20, data = x17, pch = 21, col="#800000", bg="#800000", 
     cex= 0.5, xlim = c(-100, 12000), ylim = c(-5, 16),
     frame=FALSE, ylab="Temp. anomaly", xaxt="n", yaxt="n", ann=FALSE)
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.08, line = 1), #axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.02),
  axis(2, seq(-5,16,5),  tck=-0.08, las = 1, line = 1, col = "#800000", col.axis = "#800000")
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(100,12000,500),labels=rep("",25), tck=0.02),
  #axis(4, seq(-5,16,5), tck=0.04), las = 1)
)
mtext("Temp. Anomaly 
      (deg. C)", side = 2, cex = 0.75, line = 4, adj=-0.1, col = "#C80000")
#Add Yanou tephra ages with  error bars
par(new=TRUE)
p19 <- plot(plot_y~Mean_age,data=x19c, pch = 21, col="black", bg="black", 
     cex= plot_cex, xlim = c(-100, 12000), ylim = c(-5, 16),
     frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE, xpd = 1) #xaxt="n", yaxt="n" to remove x and y axis
par(new=TRUE)
#p_axis <- c(axis(1, labels = FALSE))
arrows(x0=x19c$Mean_age-x19c$Lower_err,
       y0=x19c$plot_y,
       x1=x19c$Mean_age+x19c$Upper_err,
       y1=x19c$plot_y,
       angle=90,
       code=3,
       length=0,
       col = "grey",
       lwd=2)
par(new=TRUE)
plot(plot_y~Mean_age,data=x19a, pch = 19, col="black", 
     cex= plot_cex, xlim = c(-100, 12000), ylim = c(-5, 16), xaxt="n", 
     frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE, xpd = 1) # xpd = 1 plots over margins 
arrows(x0=x19a$Mean_age-x19a$Lower_err,
       y0=x19a$plot_y,
       x1=x19a$Mean_age+x19a$Upper_err,
       y1=x19a$plot_y,
       angle=90,
       code=3,
       length=0,
       col = "grey",
       lwd=2)
par(new=TRUE)
plot(plot_y~Mean_age,data=x19a, pch = 21, col="black", bg="white", 
     cex= plot_cex, xlim = c(-100, 12000), ylim = c(-5, 16),
     frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE, xpd = 1)
#par(new = FALSE) #remove hold

# C2) Miliken sea-ice record
p26 <- plot(TOC_pc~Age_BP, data = x26, pch = 21, type = "l", col="#B9C7F0", bg="#B9C7F0", 
            cex= 1, xlim = c(-100, 12000), ylim = rev(range(TOC_pc)),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x26$mean_6ka~age1, data = x26, 
      lty="dotted", col="black", lwd=0.75) # add average median temp for 90-60S
par(new=TRUE)
plot(TOC_pc~Age_BP, data = x26, pch = 21, col="#B9C7F0", bg="#B9C7F0",
     cex= 0.5, xlim = c(-100, 12000), ylim = c(0.6, 0),
     frame=FALSE, ylab=NULL, xlab="",xaxt="n", yaxt="n", ann=FALSE)
par(new = TRUE)
# set up and plot 100-yr  LOESS
# see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/loess for details
# locally weighted sum of squares, 2nd order polynomial, Gaussian fit 
# span α<1 determined by the time interval of the whole dataset and tri-cubic weighting proportional to (1-(distance/max. distance)3)3
# Ref: W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of Statistical Models in S eds J.M. Chambers and T.J. Hastie, Wadsworth & Brooks/Cole.
output_int_yrs = 1000 #enter the final output interval, e.g., 100 year LOESS output = 1000/12100
age_range = 10000-(-100) # enter the age range 
alpha26 <- output_int_yrs/age_range #calculate the span value
lw26 <- loess(TOC_pc~Age_BP, 
              data = x26,
              span = alpha26,
              degree = 2,
              family = "gaussian")
x26_count <- seq(from=-100, to=10000, by=1) # set up count fro 0-10 ka
x26_count_rev <- order(x26_count, decreasing=TRUE)
x26_pred <- predict(lw26, x26_count, se=TRUE) # Fit the LOESS
lines(x26_pred$fit, lty="solid", col="#00007D", lwd=1) # add it to graph
lines(x26_pred$fit-1.96*x26_pred$se.fit, 
      lty=1, col="grey", lwd=1) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x26_pred$fit+1.96*x26_pred$se.fit, 
      lty=1, col="grey", lwd=1) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x26_pred$fit, lty="solid", col="#00007D", lwd=1) # add it to graph

p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), #axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.02),
  #axis(2, seq(0.6, 05),  tck=-0.04, las = 1)
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(100,12000,500),labels=rep("",25), tck=0.02),
  axis(4, seq(0, 0.6, 0.3), tck=-0.08, las = 1, line = 1, col="#00007D", col.axis = "#00007D")
)
mtext("TOC (%)", side = 4, cex = 0.75, line = 4, col="#00007D")

# E) JRI Advance & retreat - Kaplan et al 2020
p20 <- plot(xDens20, main="JRI Advance (Kaplan et al. 2020)", cex.main = 0.75,
            xlab="Age [a BP]", xlim = c(-100, 12000), cex = 0.75, col = "blue", alpha = 0.5,
            tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL),
            xaxt="n", yaxt="n", ann=FALSE, axes = FALSE)
par(new = TRUE) #hold the plot frame
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
  #axis(2, seq(0,0.01,0.001), labels=rep("",11), tck=-0.04), 
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
  axis(4, seq(0,0.01,0.001), tck=-0.08, las = 1, line = 1))
mtext("Density", side = 4, cex = 0.75, line = 5)
#par(new = FALSE) #remove hold

p21 <- plot(xDens21, main="JRI Retreat (Kaplan et al. 2020)", cex.main = 0.75,
            xlab="Age [a BP]", xlim = c(-100, 12000), cex = 0.75, col = "red", alpha = 0.5,
            tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL),
            xaxt="n", yaxt="n", ann=FALSE, axes = FALSE)
par(new = TRUE)
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08, line=1), axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.04, line=1),
            #axis(2, seq(0,0.01,0.001), labels=rep("",11), tck=-0.04), 
            #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.01,0.001), tck=-0.08, las = 1, line = 1)
)
mtext("Density", side = 4, cex = 0.75, line = 5)
#par(new = FALSE) #remove hold

# F) JRI ice core record - Mulvaney et al 2012
p23 <- plot(T_anomaly~Age_BP, data = x23, pch = 21, type = "l", col="#800000", bg="#800000", 
             cex= 1, xlim = c(-100, 12000), ylim = c(-1.5, 2.5),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x23$T_upper_err~Age_BP, data = x23, 
      lty=1, col="grey", lwd=0.75) # add upper 95% CI error bounds
lines(x23$T_lower_err~Age_BP, data = x23,
      lty=1, col="grey", lwd=0.75) # add lower 95% CI error bounds
lines(x23$wmean~age1, data = x23, 
      lty=3, col="black", lwd=0.75) # add 12 ka mean
#lines(x23$upper~age1, data = x23, 
#      lty=3, col="grey", lwd=0.75)
#lines(x23$lower~age1, data = x23, 
#      lty=3, col="grey", lwd=0.75)
par(new = TRUE) 
p23 <- plot(T_anomaly~Age_BP, data = x23, pch = 21, type = "l", col="#800000", bg="#800000", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-1.5, 2.5),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(-100,12000,500),labels=rep("",122), tck=-0.02),
  axis(2, seq(-1.5, 2.5, 0.5), tck=-0.08,  las = 1, xpd = 1, line = 1, col.axis="#800000")
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
  #axis(4, seq(-1.5, 2.5, 0.5), tck=0.04, las = 1, xpd = 1)
  )
mtext("Temp. Anomaly (deg. C)", side = 2, cex = 0.75, line = 5, col="#800000")
#par(new = FALSE) #remove hold

# G) Kaufmann et al 2020 - 90-60S and global temp 
p25 <- plot(S90_S60_median~Age_Sh20, data = x25, pch = 21, type = "l", col="#800000", bg="#800000", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-1.5, 0.5),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
par(new = TRUE) 
p25 <- plot(S90_S60_median~Age_Sh20, data = x25, pch = 21, col="#800000", bg="#800000", 
            cex= 0.2, xlim = c(-100, 12000), ylim = c(-1.5, 0.5),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x25$global_median~Age_Sh20, data = x25, 
      lty=1, col="black", lwd=0.75) # add average median temp for 90-60S
lines(x25$global_median~age1, data = x25, 
      lty=3, col="grey", lwd=0.75) # add global data as a black line
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08), axis(1, seq(-100,12000,500),labels=rep("",122), tck=-0.02),
  #axis(2, seq(-1.5, 2.5, 0.5), tck=-0.04,  las = 1, xpd = 1),
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
  axis(4, seq(-1.5, 0.5, 0.5), tck=-0.08, las = 1, xpd = 1, line = 1)
  )
mtext("Temp. Anomaly (deg. C)", side = 4, cex = 0.75, line = 5)
mtext("Age (cal. ka BP)", side = 1, cex = 0.75, line = 3)
#par(new = FALSE) #remove hold

# WAIS Divide volcanic Event Flux -------------------------------------------------------------

#x35 <- read.csv("Inputs/x35_Vol_flux_WAIS.csv")
#par(pin=c(plotinch, plotinch/4),  mgp=c(2,0.5,0), xaxs='i') 
#p35 <- plot(x35, type = "l", lty = 1,  main="WAIS Divide volcanic Event Flux",
#            xlab="Age [a BP]", xlim = c(-100, 12000), cex.main = 1,
#            tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
#par(new = TRUE) #hold the plot frame to add median values
#p35 <- plot(x35, type ="p", pch =20, cex=0.5,  main="WASI Divide volcanic Event Flux",
#            xlab="Age [a BP]", xlim = c(-100, 12000), cex.main = 1,
#            tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
#p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.02),
#            axis(2, seq(0,200,50), tck=-0.04),
#            axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02))#,
#axis(4, seq(0,0.008,0.002), labels=rep("",5), tck=0.04))
#text(500, 160, "-", font = 2, cex = 1.5,  family = 'sans') 




#  Figure 7 - Mechanisms ---------------------------------------------------

library(Bchron)

#set working directory on mac
setwd("/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data")

# New Data
x8 <- read.csv("Inputs/x8_Advance_KGI.csv") # run density phase analysis once
xAges8 <- with (x8, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,
                                   positions = position, calCurves = calCurves, allowOutside = TRUE)) # run once at start
xDens8 <- with(x8, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, 
                                numMix = 80)) # run once at start

x9.4 <- read.csv("Inputs/x9.4_Retreat_KGI_earlyHolocene_cosmo.csv")
x9.5 <- read.csv("Inputs/x9.5_Retreat_KGI_earlyHolocene_lakes.csv")

x13 <- read.csv("Inputs/x13_Cosmo_Fildes.csv")
v <- c('ages')  ## create vector of column names to recaluculate as cal ka BP
x13[v] <- x13[v] - (2010-1950)  # convert to BP by subtracting 1950 from samples collection date and assign back

x13.1 <- read.csv("Inputs/x13.1_Cosmo_Potter.csv")
x13.1[v] <- x13.1[v] - (2010-1950) # convert to BP

# Published data
x30 <- read.csv("Inputs/x30_Laskar2004_InsAnomaly_62S.csv")
bp <- c('Age_ka')  # create vector of column names to recalculate as ka BP
x30[bp] <- x30[bp] - (2000-1950)  # convert to BP by subtracting 1950 from samples collection date and assign back

x31 <- read.csv("Inputs/x31_Baggenstos2019_Irradiance.csv")
x32 <- read.csv("Inputs/x32_Saunders2018_SPECIM.csv")
x33 <- read.csv("Inputs/x33_Saunders2018_DCond.csv") # updated age to Sh20
x34 <- read.csv("Inputs/x34_Moreno2018_SAM.csv")  # updated age to Sh20

# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch <- 3.149608
b <- ((8.27-plotinch)/2)-1
checksize <- b*2+plotinch+2
## FINAL PLOT TO PDF
pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
layout(matrix(1:5, ncol=1)) # Set up layout
par(mai=c(0.1,0.1,0.1,0.15), omi=c(0,0,0,0), pin=c(plotinch, plotinch/4), mgp=c(2,1,0), xaxs='i') # Set up internal margins
# PLOT TO SCREEN FIRST
# Clear plots
if(!is.null(dev.list())) dev.off()
layout(matrix(1:6, ncol=1)) # Set up layout
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch, plotinch/4),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0)) 

# A) Insolation and Irradiance
p30 <- plot(Spring_Summer_SONDJF~Age_ka, data = x30, pch = 21, type = "l", col="red", bg="red", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-8, 8),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
p30a <- plot(Winter_JJA~Age_ka, data = x30, pch = 21, type = "l", col="blue", bg="blue", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-8, 8),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE) 
p30b <- plot(Annual~Age_ka, data = x30, pch = 21, type = "l", lty = 3, col="black", bg="black", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-8, 8),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE) 
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), labels=rep("",37),
            axis(4, seq(-8,8,2), tck=-0.08, las = 1, xpd = 1, line = 1)
            )
mtext("62S Insolation 
      Anomaly (W m-2)", cex = 0.75, adj = 0.75, col = "black", xpd = TRUE, side=4, line=5)
par(new = TRUE)
p31 <- plot(LOESS_100_dTSI~LOESS_100yr, data = x31, pch = 21, type = "l", col="#800000", bg="#800000", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-0.5, 0.5),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
# add mean dTSI and error bar
p31 <- plot(mean~plot_age, data = x31, pch = 21, col="#800000", bg="white", 
            cex= 2, xlim = c(-100, 12000), ylim = c(-0.5, 0.5),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
#p_axis <- c(axis(1, labels = FALSE))
arrows(x0=x31$plot_age,
       y0=x31$mean-x31$error,
       x1=x31$plot_age,
       y1=x31$mean+x31$error,
       angle=180,
       code=3,
       length=0,
       col = "#800000",
       lwd=1)
par(new=TRUE)
p31 <- plot(mean~plot_age, data = x31, pch = 21, col="#800000", bg="white", 
            cex= 2, xlim = c(-100, 12000), ylim = c(-0.5, 0.5),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), labels=rep("",37),
            axis(2, seq(-0.5,0.5,0.5), tck=-0.08, las=1, col = "#800000", col.axis = "#800000", line = 1))
mtext("dTSI  W m-2", cex = 0.75, col = "#800000", side=2, line=5)

# B) Saunders et al 2018 - Macquarie SHW - get rest of data from Sigmaplot

p32 <- plot(R850_R900~Age_BP_Sh20_mean, data = x32, pch = 21, type = "l", col="#99CC99", bg="#99CC99", 
            cex= 1, xlim = c(-100, 12000), ylim = c(0.8,1.1),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
# set up the LOESS
output_int_yrs = 1000 #enter the final output interval, e.g., 100 year LOESS output = 1000/12100
age_range = 12000-(-100) # enter the age range 
alpha <- output_int_yrs/age_range #calculate the span value
lw32 <- loess(R850_R900~Age_BP_Sh20_mean, 
              data = x32,
              span = alpha)
x32_count <- seq(from=-100, to=12000, by=1) # set up count fro 0-12 ka
x32_count_rev <- order(x32_count, decreasing=TRUE)
x32_pred <- predict(lw32, x32_count, se=TRUE) # Fit the LOESS
lines(x32_pred$fit, lty="solid", col="#008000", lwd=0.75) # add it to graph
lines(x32_pred$fit-1.96*x32_pred$se.fit, 
      lty="dashed", col="grey", lwd=0.75) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x32_pred$fit+1.96*x32_pred$se.fit, 
      lty="dashed", col="grey", lwd=0.751) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x32_pred$fit, lty="solid", col="#008000", lwd=0.75) # add it to graph
#Create polygon error shading - doesn't work for negative values on x-axis
#x32_y_polygon <- c((x32_pred$fit+1.96*x32_pred$se.fit)[x32_count], 
#                   (x32_pred$fit-1.96*x32_pred$se.fit)[x32_count_rev]) 
#x32_x_polygon <- c(x32_count, x32_count_rev)
#polygon(x32_x_polygon, x32_y_polygon, col="#00009933", border=NA)
par(new = TRUE)
p32b <- plot(R850_R900_mean~age_plot, data = x32, pch = 21, type = "l", lty = 2, lwd = 1, col="#666666", bg="#666666", 
             cex= 1, xlim = c(-100, 12000), ylim = c(0.8,1.1),
             frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
p32c <- plot(R850_R900_mean+R850_R900_sd~age_plot, data = x32, pch = 21, type = "l", lty = 3, lwd = 1, col="#B2B2B2", bg="#B2B2B2", 
             cex= 1, xlim = c(-100, 12000), ylim = c(0.8,1.1),
             frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE)
p32d <- plot(R850_R900_mean-R850_R900_sd~age_plot, data = x32, pch = 21, type = "l", lty = 3, lwd = 1, col="#B2B2B2", bg="#B2B2B2", 
             cex= 1, xlim = c(-100, 12000), ylim = c(0.8,1.1),
             frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
p_axis <- c(#axis(1, seq(0,18000,2000), tck=-0.04), labels=rep("",37),
            axis(4, seq(0.8,1.1,0.1), tck=-0.08, col = "#008000", col.axis = "#008000", las=1, line=1))
mtext("R850/R900 index
      Minerogenic input", cex = 0.75, side=4, line=5, col = "#008000", adj = 0.8, xpd = TRUE)

# D) SAM plot
p34 <- plot(SAM_index_reconstruction~Age_Sh13, data = x34, pch = 21, type = "l", col="black", bg="black", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-2, 6),
            frame=FALSE, ylab = "", xlab="", yaxt="n", xaxt="n")
par(new = TRUE) #hold the plot frame to add median values
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), labels=rep("",37),
            axis(4, seq(-2,6,2), tck=-0.08, las=1, line=1))
mtext("Reconstructed
      SAM Index", cex = 0.75, side=4, line=5, col = "black", srt=180, adj = 1, xpd = TRUE)

# G) KGI: max. age constraint phases for glacier readvance  [this study, HB et al., sub & Hall, 2007]
p8 <- plot(xDens8, xlab="", xlim = c(-100, 12000), cex.main = 1, las = 0, axes = FALSE,
           tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL), 
           xaxt="n", yaxt="n", ann=FALSE) #supress x-axis
par(new = TRUE) #hold the plot frame
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08, line=1), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.04, line=1),
  axis(2, seq(0,0.012,0.002), tck=-0.08, las = 1, line = 1)#,
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02)
)
mtext("Age (cal. ka BP)", cex = 0.75, side=1, line=4, col = "black", xpd = TRUE)
mtext("Density", cex = 0.75, side=2, line=5, col = "black", xpd = TRUE)

# E) KGI Deglaciation - ** new plot 0-18 ka ** align in illustrator **

# 0-18 ka plot 12 cm wide
# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch1 <- 12 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch1 <- 3.149608 *(12/8)
b <- ((8.27-plotinch1)/2)-1
checksize <- b*2+plotinch1+2
## FINAL PLOT TO PDF
pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
layout(matrix(1:2, ncol=1)) # Set up layout
par(mai=c(0.25,0.25,0.25,0.25), omi=c(0,0,0,0), pin=c(plotinch1, plotinch1/1.99710788366826), mgp=c(2,0.5,0), xaxs='i') # Set up internal margins
# PLOT TO SCREEN FIRST
# Clear plots
if(!is.null(dev.list())) dev.off()
layout(matrix(1:2, ncol=1)) # Set up layout
par(mai=c(0.25,0.25,0.25,0.25), pin=c(plotinch1, plotinch1/1.99710788366826),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 

# Add new Fildes new cosmo ages to secondary axis of altitude
x <- x13$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13$altitude
x.err <- x13$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# KGI early Holocene deglaciation constraints - cosmo ages to secondary axis of altitude
par(new=TRUE)
x <- x9.4$ages
#y <- c(0.004, 0.005, 0.006)
y <- x9.4$altitude
x.err <- x9.4$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "#808080", bg= "#808080", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "#808080", bg= "#808080", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# KGI early Holocene deglaciation constraints - lake ages to secondary axis of altitude
par(new=TRUE)
x <- x9.5$ages
#y <- c(0.004, 0.005, 0.006)
y <- x9.5$altitude
x.err <- x9.5$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 22, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 22, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# Potter He-3 cosmo ages
par(new=TRUE)
x <- x13.1$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13.1$altitude
x.err <- x13.1$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "#C80000", bg= "#C80000", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "#C80000", bg= "#C80000", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
p_axis <- c(axis(1, seq(0,18000,2000), tck=-0.04, cex.axis = 0.75), labels=rep("",37), 
  axis(4, seq(0,300,50), cex.axis = 0.75, las=1))
mtext("Age (a BP)", side=1, line=2, col = "black", cex = 0.75)
mtext("Altitude (m a.s.l.)", side=4, line=2, col = "black", cex = 0.75)







# POTTER PAPER FIGURES  ---------------------------------------------------------

#set working directory on mac
setwd("/Users/Steve/Dropbox/BAS/Data/R/BChron/KGI/Data")

## FINAL PLOT TO PDF
#pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
#layout(matrix(1:6, ncol=1)) # Set up layout
#par(mai=c(0.1,0.1,0.1,0.1), omi=c(0,0,0,0), pin=c(plotinch, plotinch/4), mgp=c(2,0.5,0), xaxs='i') # Set up internal margins

#  Figure 6: Potter and KGI advance data ---------------------------------------------

# B) # 0-18 ka plot 12 cm wide
# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch1 <- 12 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch1 <- 3.149608 *(12/8)
b <- ((8.27-plotinch1)/2)-1
checksize <- b*2+plotinch1+2

# PLOT TO SCREEN FIRST
# Clear plots
if(!is.null(dev.list())) dev.off()
layout(matrix(1:2, ncol=1)) # Set up layout
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch1, plotinch1/1.99710788366826),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0)) 

# D) RCarbon Stacked Density Potter C14 ages prob plot 

# PLOT TO SCREEN FIRST
# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:2, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch1, plotinch1/3),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 

library(rcarbon)
r6 <- read.csv("Inputs/r6_Potter_Terr_Holocene.csv")
# define time ranges
r6_timeRange_18ka <- c(18000,-100)
r6_revtimeRange_18ka <- c(-100,18000)
r6_cal <- calibrate(r6$CRA, r6$Error, normalised=TRUE, calCurves=r6$calCurves,
                    resOffsets=r6$resOffsets,resErrors=r6$resErrors,ids=r6$LabID)
#create bins based using a suitable h split value based on maximum consensus in above  
r6_bins <- binPrep(r6$SiteName, r6$CRA, h=100)
r6_res = stackspd(x=r6_cal,timeRange=r6_timeRange_18ka,bins=r6_bins,group=r6$SiteName)
# make a stacked plot
p6 <- plot(r6_res,type='stacked', xlim = r6_revtimeRange_18ka, 
           legend = TRUE, legend.arg = NULL, ann = FALSE, axes=FALSE, xaxt="n", yaxt="n")  # most useful for summaries of multiple sites / types
axis(side=1, at=seq(0, 18000, by=2000))
axis(side=2, at=seq(0,0.030, by=0.005))
par(new = TRUE) # hold the plot frame

# E) KGI: max. age constraint phases for glacier readvance  [this study, HB et al., sub & Hall, 2007]
x8 <- read.csv("Inputs/x8_Advance_KGI.csv") # run density phase analysis once
#xAges8 <- with (x8, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,
#                                    positions = position, calCurves = calCurves, allowOutside = TRUE)) # run once at start
#xDens8 <- with(x8, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, 
#                                 numMix = 80)) # run once at start

p8 <- plot(xDens8, xlab="", xlim = c(-100, 12000), cex.main = 1, las = 0, axes = FALSE,
           tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL), 
           xaxt="n", yaxt="n", ann=FALSE) #supress x-axis
par(new = TRUE) #hold the plot frame
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08, line=1), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.04, line=1),
            axis(2, seq(0,0.012,0.002), tck=-0.08, las = 1, line = 1)#,
            #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04), axis(3, seq(0,12000,500),labels=rep("",25), tck=0.02)
)
mtext("Age (cal. ka BP)", cex = 0.75, side=1, line=4, col = "black", xpd = TRUE)
mtext("Density", cex = 0.75, side=2, line=5, col = "black", xpd = TRUE)

# Add new Fildes new cosmo ages to secondary axis of altitude
x <- x13$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13$altitude
x.err <- x13$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# KGI early Holocene deglaciation constraints - cosmo ages to secondary axis of altitude
par(new=TRUE)
x <- x9.4$ages
#y <- c(0.004, 0.005, 0.006)
y <- x9.4$altitude
x.err <- x9.4$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "#808080", bg= "#808080", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "#808080", bg= "#808080", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# KGI early Holocene deglaciation constraints - lake ages to secondary axis of altitude
par(new=TRUE)
x <- x9.5$ages
#y <- c(0.004, 0.005, 0.006)
y <- x9.5$altitude
x.err <- x9.5$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 22, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 22, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# Potter He-3 cosmo ages
par(new=TRUE)
x <- x13.1$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13.1$altitude
x.err <- x13.1$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "#C80000", bg= "#C80000", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "#C80000", bg= "#C80000", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
p_axis <- c(axis(1, seq(0,18000,2000), tck=-0.04, cex.axis = 0.75), labels=rep("",37), line=1, 
            axis(4, seq(0,300,50), tck=-0.04, cex.axis = 0.75, las=1, line=1))
mtext("Age (a BP)", side=1, line=2, col = "black", cex = 0.75)
mtext("Altitude (m a.s.l.)", side=4, line=3, col = "black", cex = 0.75)


#  Figure 9: Palaeo-records ---------------------------------------------

# Published datasets - import other datasets from Fildes plots above new datasets 
x17 <- read.csv("Inputs/x17_Yanou_GDGT.csv")
x20 <- read.csv("Inputs/x20_Kaplan2020_JRI_Advance.csv")
x20[v] <- x20[v] - (2010-1950) #convert to BP from sampling date (nearest 10 years) and assign back
xAges20 <- with (x20, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xDens20 <- with(x20, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 30))
x21 <- read.csv("Inputs/x21_Kaplan2020_JRI_Retreat.csv")
x21[v] <- x21[v] - (2010-1950) #convert to BP from sampling date (nearest 10 years) and assign back
xAges21 <- with (x21, BchronCalibrate(ids = id, ages = ages, ageSds = ageSds,positions = position, calCurves = calCurves, allowOutside = TRUE))
xDens21 <- with(x21, BchronDensity(ages = ages, ageSds = ageSds,calCurves = calCurves, numMix = 20))
x23 <- read.csv("Inputs/x23_JRI_ice_core.csv")
x26 <- read.csv("Inputs/x26_Miliken2009.csv")
x27 <- read.csv("Inputs/x27_Etourneaux2013_SST.csv")
x28 <- read.csv("Inputs/x28_Shevenell_TEX86.csv")
x29 <- read.csv("Inputs/x29_Roberts2017_WAP_diatoms.csv")

# LARGER HEIGHT PLOTS - plot to screen ----------------------------------

# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch <- 3.149608
b <- ((8.27-plotinch)/2)-1
checksize <- b*2+plotinch+2
# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:7, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch, plotinch/3),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0)) 

# E) Yanou GDGT and tephra plots
p17a <- plot(GDGT_Anomaly_FOSTER_2016~Age_Sh20, data = x17, pch = 21, type = "l", col="#D19C9C", bg="#D19C9C", 
             cex= 1, xlim = c(-100, 12000), ylim = c(-5, 16),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x17$mean_ANT_anomaly~age1, data = x17, 
      lty="dotted", col="#D19C9C", lwd=0.75) # add average median temp for 90-60S
par(new=TRUE)
plot(GDGT_Anomaly_FOSTER_2016~Age_Sh20, data = x17, pch = 21, col="#D19C9C", bg="#D19C9C",
     cex= 0.5, xlim = c(-100, 12000), ylim = c(-5, 16),
     frame=FALSE, ylab=NULL, xlab="",xaxt="n", yaxt="n", ann=FALSE)
lines(x17$RMSE_Max_anomaly~Age_Sh20, data = x17, 
      lty=1, col="grey", lwd=0.75) # add upper 95% CI error bounds
lines(x17$RMSE_min_anomaly~Age_Sh20, data = x17,
      lty=1, col="grey", lwd=0.75) # add lower 95% CI error bounds
par(new=TRUE)
p17b <- plot(GDGT_ANT_GLOBAL_anomaly~Age_Sh20, data = x17, pch = 21, type = "l", col="#800000", bg="#800000", 
             cex= 1, xlim = c(-100, 12000), ylim = c(-5, 16),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x17$mean_ANT_GLOBAL_anomaly~age1, data = x17, 
      lty="dotted", col="#800000", lwd=0.75) # add global data as a black line
par(new=TRUE)
plot(GDGT_ANT_GLOBAL_anomaly~Age_Sh20, data = x17, pch = 21, col="#800000", bg="#800000", 
     cex= 0.5, xlim = c(-100, 12000), ylim = c(-5, 16),
     frame=FALSE, ylab="Temp. anomaly", xaxt="n", yaxt="n", ann=FALSE)
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.08, line = 1), #axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.02),
  axis(2, seq(-5,16,5),  tck=-0.08, las = 1, line = 1, col = "#800000", col.axis = "#800000")
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(100,12000,500),labels=rep("",25), tck=0.02),
  #axis(4, seq(-5,16,5), tck=0.04), las = 1)
)
mtext("Temp. Anomaly 
      (deg. C)", side = 2, cex = 0.75, line = 4, adj=-0.1, col = "#C80000")

# J) JRI ice core record - Mulvaney et al 2012
p23 <- plot(T_anomaly~Age_BP, data = x23, pch = 21, type = "l", col="#800000", bg="#800000", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-2, 2.5),
            axes=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x23$T_upper_err~Age_BP, data = x23, 
      lty=1, col="grey", lwd=0.75) # add upper 95% CI error bounds
lines(x23$T_lower_err~Age_BP, data = x23,
      lty=1, col="grey", lwd=0.75) # add lower 95% CI error bounds
lines(x23$wmean~age1, data = x23, 
      lty=3, col="black", lwd=0.75) # add 12 ka mean
#lines(x23$upper~age1, data = x23, 
#      lty=3, col="grey", lwd=0.75)
#lines(x23$lower~age1, data = x23, 
#      lty=3, col="grey", lwd=0.75)
par(new = TRUE) 
p23 <- plot(T_anomaly~Age_BP, data = x23, pch = 21, type = "l", col="#800000", bg="#800000", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-2, 2.5),
            axes=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x23$warmphase_plot~RRR, data = x23, 
      lty=1, col="#FF0000", lwd=4) 
lines(x23$warmphase_plot~MCA, data = x23, 
      lty="dotted", col="#E26464", lwd=4) 
lines(x23$warmphase_plot~HH, data = x23, 
      lty=1, col="#950000", lwd=4) 
lines(x23$warmphase_plot~EHWP, data = x23, 
      lty=1, col="#950000", lwd=4) 
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08, line=1), axis(1, seq(0,12000,500),labels=rep("",25), tck=-0.04, line=1),
  axis(2, seq(-2, 2.5, 1), tck=-0.08,  las = 1, xpd = 1, line = 1, col = "#800000", col.axis = "#800000")
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
  #axis(4, seq(-1.5, 2.5, 0.5), tck=0.04, las = 1, xpd = 1)
)
mtext("Temp. Anomaly 
      (deg. C)", side = 2, cex = 0.75, adj = -0.1, line = 4, col="#800000")
#par(new = FALSE) #remove hold

# SMALL HEIGHT PLOTS - plot to screen ---------------------------------

# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch <- 3.149608
b <- ((8.27-plotinch)/2)-1
checksize <- b*2+plotinch+2
# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:7, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch, plotinch/5),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0))

# F) Miliken sea-ice record
p26 <- plot(TOC_pc~Age_BP, data = x26, pch = 21, type = "l", col="#B9C7F0", bg="#B9C7F0", 
             cex= 1, xlim = c(-100, 12000), ylim = rev(range(TOC_pc)),
             frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x26$mean_6ka~age1, data = x26, 
      lty="dotted", col="black", lwd=0.75) # add average median temp for 90-60S
par(new=TRUE)
plot(TOC_pc~Age_BP, data = x26, pch = 21, col="#B9C7F0", bg="#B9C7F0",
     cex= 0.5, xlim = c(-100, 12000), ylim = c(0.6, 0),
     frame=FALSE, ylab=NULL, xlab="",xaxt="n", yaxt="n", ann=FALSE)
par(new = TRUE)
# set up and plot 100-yr  LOESS
# see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/loess for details
# locally weighted sum of squares, 2nd order polynomial, Gaussian fit 
# span α<1 determined by the time interval of the whole dataset and tri-cubic weighting proportional to (1-(distance/max. distance)3)3
# Ref: W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of Statistical Models in S eds J.M. Chambers and T.J. Hastie, Wadsworth & Brooks/Cole.
output_int_yrs = 1000 #enter the final output interval, e.g., 100 year LOESS output = 1000/12100
age_range = 10000-(-100) # enter the age range 
alpha26 <- output_int_yrs/age_range #calculate the span value
lw26 <- loess(TOC_pc~Age_BP, 
              data = x26,
              span = alpha26,
              degree = 2,
              family = "gaussian")
x26_count <- seq(from=-100, to=10000, by=1) # set up count fro 0-10 ka
x26_count_rev <- order(x26_count, decreasing=TRUE)
x26_pred <- predict(lw26, x26_count, se=TRUE) # Fit the LOESS
lines(x26_pred$fit, lty="solid", col="#00007D", lwd=1) # add it to graph
lines(x26_pred$fit-1.96*x26_pred$se.fit, 
      lty=1, col="grey", lwd=1) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x26_pred$fit+1.96*x26_pred$se.fit, 
      lty=1, col="grey", lwd=1) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x26_pred$fit, lty="solid", col="#00007D", lwd=1) # add it to graph

p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), #axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.02),
  #axis(2, seq(0.6, 05),  tck=-0.04, las = 1)
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(100,12000,500),labels=rep("",25), tck=0.02),
  axis(4, seq(0, 0.6, 0.3), tck=-0.08, las = 1, line = 1, col="#00007D", col.axis = "#00007D")
  )
mtext("TOC (%)", side = 4, cex = 0.75, line = 4, col="#00007D")

# G) WAP sea-ice record - Roberts et al. (2017)
p29 <- plot(Fcurta_Fkerg~Age_Mar20, data = x29, pch = 21, type = "l", col="#808080", bg="#808080", 
            cex= 1, xlim = c(-100, 12000), ylim = rev(range(Fcurta_Fkerg)),
            frame=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x29$mean12ka~age1, data = x29, 
      lty="dotted", col="black", lwd=0.75) # mean
par(new=TRUE)
plot(Fcurta_Fkerg~Age_Mar20, data = x29, pch = 21, col="#808080", bg="#808080",
     cex= 0.5, xlim = c(-100, 12000), ylim = rev(range(Fcurta_Fkerg)),
     axes=FALSE, ylab=NULL, xlab="",xaxt="n", yaxt="n", ann=FALSE)
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), #axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.02),
  axis(2, seq(0, 4, 0.5),  tck=-0.08, las = 1, line=1, col="#808080", col.axis="#808080")
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(100,12000,500),labels=rep("",25), tck=0.02),
  #axis(4, seq(0, 3, 1), tck=-0.06, las = 1)
)
mtext("F. curta / 
F. Kerg.", side = 2, cex = 0.75, line = 4, adj = 0.5, col="#808080")

par(new = TRUE) # add open water pelagic plot on secondary axis

p29a <- plot(PelagicOW_pc~Age_Mar20, data = x29, pch = 21, type = "l", col="#0000C8", bg="#0000C8", 
            cex= 1, xlim = c(-100, 12000), ylim=c(0, 40),
            axes=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
par(new = TRUE)
plot(PelagicOW_pc~Age_Mar20, data = x29, pch = 21, col="#0000C8", bg="#0000C8",
     cex= 0.5, xlim = c(-100, 12000), ylim=c(0, 40),
     axes=FALSE, ylab="", xlab="",xaxt="n", yaxt="n", ann=FALSE)
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), #axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.02),
  axis(4, seq(0, 40, 10),  tck=-0.08, las = 1, line=1, col="#0000C8", col.axis = "#0000C8"), 
  axis(4, seq(0, 40, 10),  tck=-0.04, las = 1, line=1, col="#0000C8", col.axis = "#0000C8"))
mtext("Pelagic Open  
Water sp. (%)", side = 4, cex = 0.75, line = 5, adj = 0, col = "#0000C8")

# H) SST Temp - Etourneuax et al. (2013) 
p27 <- plot(SST_degC~Age_Mar20, data = x27, pch = 21, type = "l", col="#D19C9C", bg="#D19C9C", 
            cex= 1, xlim = c(-100, 12000), ylim = c(-1, 5),
            axes=FALSE, ylab="", xlab="", xaxt="n", yaxt="n", ann=FALSE)
lines(x27$mean_6ka~age1, data = x27, 
      lty="dotted", col="black", lwd=0.75) # mean
par(new=TRUE)
plot(SST_degC~Age_Mar20, data = x27, pch = 21, col="#D19C9C", bg="#D19C9C",
     cex= 0.5, xlim = c(-100, 12000), ylim = c(-1, 5),
     axes=FALSE, ylab=NULL, xlab="",xaxt="n", yaxt="n", ann=FALSE)
# set up and plot 100-yr  LOESS
output_int_yrs = 1000 #enter the final output interval, e.g., 100 year LOESS output = 1000/12100
age_range = 9000-(-100) # enter the age range 
alpha27 <- output_int_yrs/age_range #calculate the span value
lw27 <- loess(SST_degC~Age_Mar20, 
              data = x27,
              span = alpha27)
x27_count <- seq(from=-100, to=9000, by=1) # set up count fro 0-10 ka
x27_count_rev <- order(x27_count, decreasing=TRUE)
x27_pred <- predict(lw27, x27_count, se=TRUE) # Fit the LOESS
lines(x27_pred$fit, lty="solid", col="#800000", lwd=1) # add it to graph
lines(x27_pred$fit-1.96*x27_pred$se.fit, 
      lty=1, col="grey", lwd=1) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x27_pred$fit+1.96*x27_pred$se.fit, 
      lty=1, col="grey", lwd=1) #Fit 95% error lines - SE is very small here as there's so much data!
lines(x27_pred$fit, lty="solid", col="#800000", lwd=1) # add it to graph
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08, line=1), axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.04, line=1),
            axis(2, seq(-1, 5, 1),  tck=-0.08, las = 1, line=1, col="#800000", col.axis = "#800000"),
            axis(2, seq(-1, 5, 1),  tck=-0.04, las = 1, line=1, col="#800000", col.axis = "#800000")
)
mtext("SST Temp 
(deg. C)", side = 2, cex = 0.75, line = 4, col="#800000")

# SMALLER HEIGHT PLOTS - PLOT TO SCREEN ---------------------------------

# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch <- 3.149608
b <- ((8.27-plotinch)/2)-1
checksize <- b*2+plotinch+2
# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:7, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch, plotinch/10),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0))


# I) JRI Advance & retreat - Kaplan et al 2020
p20 <- plot(xDens20, main="JRI Advance (Kaplan et al. 2020)", cex.main = 0.75,
            xlab="Age [a BP]", xlim = c(-100, 12000), cex = 0.75, col = "blue", alpha = 0.5,
            tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL),
            xaxt="n", yaxt="n", ann=FALSE, axes = FALSE)
par(new = TRUE) #hold the plot frame
p_axis <- c(#axis(1, seq(0,12000,2000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
  #axis(2, seq(0,0.01,0.001), labels=rep("",11), tck=-0.04), 
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
  axis(4, seq(0,0.01,0.001), tck=-0.08, las = 1, line = 1))
mtext("Density", side = 4, cex = 0.75, line = 5)
#par(new = FALSE) #remove hold

p21 <- plot(xDens21, main="JRI Retreat (Kaplan et al. 2020)", cex.main = 0.75,
            xlab="Age [a BP]", xlim = c(-100, 12000), cex = 0.75, col = "red", alpha = 0.5,
            tck=-0.04, panel.first = grid(nx = NULL, col = NULL, lty = NULL),
            xaxt="n", yaxt="n", ann=FALSE, axes = FALSE)
par(new = TRUE)
p_axis <- c(axis(1, seq(0,12000,2000), tck=-0.08, line=1), axis(1, seq(0,12000,500), labels=rep("",25), tck=-0.04, line=1),
  #axis(2, seq(0,0.01,0.001), labels=rep("",11), tck=-0.04), 
  #axis(3, seq(0,12000,2000), labels=rep("",7), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
  axis(4, seq(0,0.01,0.001), tck=-0.08, las = 1, line = 1)
)
mtext("Density", side = 4, cex = 0.75, line = 5)
#par(new = FALSE) #remove hold


# Figure 10D  mechanisms plot to 18 ka -------------------------------------------------------------

# E) KGI Deglaciation - ** new plot 0-18 ka ** align in illustrator **

# 0-18 ka plot 12 cm wide
# Clear plots
if(!is.null(dev.list())) dev.off()
plotinch1 <- 12 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch1 <- 3.149608 *(12/8)
b <- ((8.27-plotinch1)/2)-1
checksize <- b*2+plotinch1+2

## FINAL PLOT TO PDF
pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
layout(matrix(1:2, ncol=1)) # Set up layout
par(mai=c(0.1,0.1,0.1,0.1), omi=c(0,0,0,0), pin=c(plotinch1, plotinch1/1.3), mgp=c(2,1,0), xaxs='i') # Set up internal margins
# PLOT TO SCREEN FIRST
# Clear plots
if(!is.null(dev.list())) dev.off()
layout(matrix(1:2, ncol=1)) # Set up layout
par(mai=c(0.1,0.1,0.1,0.1), pin=c(plotinch1, plotinch1/1.3),  mgp=c(2,1,0), xaxs='i', omi=c(0,0,0,0)) 

# Add new Fildes new cosmo ages to secondary axis of altitude
x <- x13$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13$altitude
x.err <- x13$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# Potter He-3 cosmo ages
par(new=TRUE)
x <- x13.1$ages
#y <- c(0.004, 0.005, 0.006)
y <- x13.1$altitude
x.err <- x13.1$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "#C80000", bg= "#C80000", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "#C80000", bg= "#C80000", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# KGI early Holocene deglaciation constraints - cosmo ages to secondary axis of altitude
par(new=TRUE)
x <- x9.4$ages
#y <- c(0.004, 0.005, 0.006)
y <- x9.4$altitude
x.err <- x9.4$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 21, col= "#808080", bg= "#808080", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 21, col= "#808080", bg= "#808080", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

# KGI early Holocene deglaciation constraints - lake ages to secondary axis of altitude
par(new=TRUE)
x <- x9.5$ages
#y <- c(0.004, 0.005, 0.006)
y <- x9.5$altitude
x.err <- x9.5$ageSds
#arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00)
#points(x, y, pch = 21, col= "black", bg = "red", cex=2)
plot(x, y, pch = 22, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)
arrows(x0=x-x.err, y0=y, x1=x+x.err, code=3, angle=90, length = 0.00, col= "darkgrey", lwd = 1)
par(new=TRUE)
plot(x, y, pch = 22, col= "black", bg= "white", cex = 1,
     xlab = "", ylab = "", xlim=c(-100,18000), ylim=c(0,300), axes = FALSE)

p_axis <- c(axis(1, seq(0,18000,2000), tck=-0.04), labels=rep("",37), line=1,
            axis(4, seq(0,300,100), las=1, line=1))
mtext("Age (a BP)", side=1, line=2, col = "black", cex = 1)
mtext("Altitude (m a.s.l.)", side=4, line=2, col = "black", cex = 1)

# BChron Potter C14 ages prob plot 
p6 <- plot(xDens6, main="Potter Peninsula: max. age constraints on glacier advance",
           xlab="Age [a BP]", xlim = c(-100, 18000), col.axis = "white",
           tck=-0.04, panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
axis(side=1, at=seq(0, 18000, by=2000))
axis(side=2, at=seq(0,0.015, by=0.005))
par(new = TRUE) #hold the plot frame
p_axis <- c(axis(1, seq(0,18000,2000), tck=-0.04), axis(1, seq(0,18000,500),labels=rep("",37), tck=-0.02),
            axis(2, seq(0,0.015,0.005), labels=rep("",16), tck=-0.02), 
            axis(3, seq(0,18000,2000), labels=rep("",10), tck=0.04),
            axis(4, seq(0,0.015,0.005), labels=rep("",16), tck=0.02))
par(new = FALSE) #remove hold on the plot frame to add median values

# Use Stacked SPD plot - note that S&J earlier distribution in 'stacked' plot shown in multipanel plot ends up plotting behind - fix this

# STACKED PLOTS
plot(r6_res,type='stacked', xlim = r6_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion')
plot(r6_res,type='multipanel', xlim = r6_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 

# Clear plots
if(!is.null(dev.list())) dev.off()
# plot XDens output together for comparison
layout(matrix(1:3, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.25,0.25,0.25), pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 

# RCarbon Stacked Density Potter C14 ages prob plot 
p6 <- plot(r6_res,type='stacked', xlim = r6_revtimeRange_18ka, 
           legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
axis(side=1, at=seq(0, 18000, by=2000))
axis(side=2, at=seq(0,0.015, by=0.005))
par(new = TRUE) #hold the plot frame
p_axis <- c(axis(1, seq(0,18000,2000), tck=-0.04), axis(1, seq(0,18000,500),labels=rep("",37), tck=-0.02),
            axis(2, seq(0,0.015,0.005), labels=rep("",16), tck=-0.02), 
            axis(3, seq(0,18000,2000), labels=rep("",10), tck=0.04),
            axis(4, seq(0,0.015,0.005), labels=rep("",16), tck=0.02))
par(new = FALSE) #remove hold on the plot frame to add median values



# ++++++++++++++++  END ++++++++++++++++  --------

# ++++++++++++++++ OTHER USEFUL CODE ++++++++++++++ --------

# Convert a to ka
x19a <- read.csv("Inputs/x19a_YAN_A_tephra_ages.csv") #airfall only - with T3a lower error extending to T3c
Uage <- c('Upper_95CI_age')
age <- c('Mean_age')  ## create vector of column names to recaluculate as cal ka BP
Uerr <- c('Upper_err')
Lerr <- c('Lower_err')
x19a[Uage]<- x19a[Uage]/10^3
x19a[age] <- x19a[age]/10^3
x19a[Uerr] <- x19a[Uerr]/10^3
x19a[Lerr] <- x19a[Lerr]/10^3
x19a
x19b <- read.csv("Inputs/x19b_YAN_ARW_tephra_ages.csv") #airfall and reworked
x19b[Uage]<- x19b[Uage]/10^3
x19b[age] <- x19b[age]/10^3
x19b[Uerr] <- x19b[Uerr]/10^3
x19b[Lerr] <- x19b[Lerr]/10^3
x19b
x19c <- read.csv("Inputs/x19c_Ant2018_tephra_ages.csv") #airfall and reworked
x19c[Uage]<- x19c[Uage]/10^3
x19c[age] <- x19c[age]/10^3
x19c[Uerr] <- x19c[Uerr]/10^3
x19c[Lerr] <- x19c[Lerr]/10^3
x19c

x19b <- read.csv("Inputs/x19b_YAN_ARW_tephra_ages.csv") #airfall and reworked
x19c <- read.csv("Inputs/x19c_Ant2018_tephra_ages.csv") #airfall and reworked

# STANDARD PLOTTING SIZE

# set to 8 cm width  = 12 ka to match Sigmaplot output for -100 a to 12 ka
# set to 12 cm width  = 18 ka to match Sigmaplot output for -100 a to 12 ka to match 8 cm = 12 ka

# Clear plots
if(!is.null(dev.list())) dev.off()
## Translate 8 cm graph plot size *from* cm *to* inches:
plotinch_8 <- 8 / cm(1) # -> 8 cm  is  3.149608 inches but R converts to 3.149606 - 
plotinch_12 <- 12 / cm(1) 
plotinch_9 <- 9 / cm(1) 
plotinch_7 <- 7 / cm(1)
## FINAL PLOT TO PDF
#mai / omi = margin in inches or mar / oma = margin in lines - order is bottom, left, top, right 
#set up the page layout - A4 is 8.27 x 11.69 - closest to 
#pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
layout(matrix(1:2, ncol=1)) # Set up layout
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.25,0.25,0.25), pin=c(plotinch_8, plotinch_8/2),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 
#dev.off() #need to include this to write to pdf file fully - will delete screen plot

#Summary of phases at prob = 0.95
summary(xDens1, type = "outliers", prob = 0.95) # Look at outlier probabilities 0.95 is the default
summary(xDens2, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens2.1, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens3, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens4, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens5, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens6, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens7, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens8, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens9, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens10, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens11, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens12, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens13, type = "outliers", prob = 0.95) # Look at outlier probabilities
summary(xDens14, type = "outliers", prob = 0.95) # Look at outlier probabilities

# Use Stacked SPD plot - note that S&J earlier distribution in 'stacked' plot shown in multipanel plot ends up plotting behind - fix this

# STACKED PLOTS
r6_timeRange <- c(12000,-100)
r6_revtimeRange <- c(-100,12000)
plot(r6_res,type='stacked', xlim = r6_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion')
#plot(r6_res,type='multipanel', xlim = r6_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 

# BChron slow cluster plot 
cc(
  ages, 
  ageSds,
  calCurves,
  pathToCalCurves = system.file("data", package = "Bchron"),
  dfs = rep(100, length(ages)),
  numMix = 50,
  iterations = 10000,
  burn = 2000,
  thin = 8,
  updateAges = FALSE,
  store_density = TRUE
)

# plot it
plot(xDens)

# Run the faster model
xDensFast <- with(
  x,
  BchronDensityFast(
    ages = ages,
    ageSds = ageSds,
    calCurves = calCurves
  )
)

# plot it
plot(xDensFast)


## Make an age depth model - doesnt work for ART_terr data 
# Run in Bchronology - using shcal20 in calCurves column
# Runs the Compound Poisson-Gamma chronology model of Haslett and Parnell (2008)
# Fits a non-parametric chronology model to age/position data according to the Compound Poisson-Gamma model defined by Haslett and Parnell (2008) <DOI:10.1111/j.1467-9876.2008.00623.x>. This version uses a slightly modified Markov chain Monte Carlo fitting algorithm which aims to converge quicker and requires fewer iterations. It also a slightly modified procedure for identifying outliers
xOut <- with(x,
             Bchronology(
               ages = ages,
               ageSds = ageSds,
               calCurves = calCurves,
               positions = position,
               positionThicknesses = thickness,
               jitterPositions = TRUE,
               ids = id,
               predictPositions = seq(0, 100, by = 10)
             )
)

# Summarise it a few different ways
summary(xOut) # Default is for quantiles of ages at predictPosition values
summary(xOut, type = "convergence") # Check model convergence
summary(xOut, type = "outliers") # Look at outlier probabilities

legend("topright",
       legend = rownames(gfg_data),
       pch = 15,
       col = 1:nrow(gfg_data),cex=2.5)

