#sets working directory on imac
setwd("/Users/Steve/Dropbox/BAS/Data/R/RCarbon")
getwd()

#sets working directory on macbook
setwd("/Users/Steve/Dropbox/BAS/Data/R/RCarbon")
#check working directory
getwd()
library(cowplot) # for plotting
library(ggplot2)
library(rcarbon)
library (dplyr)
library (tidyr)
library (parameters)
library (mclust)
library (NbClust)
library (cluster)

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Clear plots
if(!is.null(dev.list())) dev.off()

# +++++++++++ Example for calibrating a single date in INTCAL - simplest output ++++++++++++++++++++==

x_INTCAL <- calibrate(x=4000, errors=30)
plot(x1)
# Example with a Marine Date, using a DeltaR of 300 and a DeltaR error of 30
x1 <- calibrate(x=c(4000), errors=c(30), calCurves='marine20', resOffsets=300, resErrors=30)
plot(x1)
x1
## Example of how to combine multiple CalDates Class Objects into one with different calCurves for each.
x1 = calibrate(c(2000,3400),c(20,20),ids=1:2)
x2 = calibrate(c(4000,3000),c(30,30),calCurves=c('intcal20','marine20'),
               resOffsets=c(0,30),resErrors=c(0,20),ids=3:4)
## Example of how to create a mixed calibration curve
mcurve <- mixCurves('intcal20','marine20',p=0.7,resOffsets=300,resErrors=20)
x3 = calibrate(5300,20,calCurves=mcurve,ids=5)
x4 = combine(x1,x2,x3)
plot(x4)

## Example of how to create a marine20 calibration curve with DR+/-err to import/use in e.g., BChron
# p is the proportion of the first curve in the output - here p = 0 here means marine20 second curve with the offset 666 is generated
# if p = 1, which is the default, means marine20 without an offset is generated 
marine20_666 <- mixCurves('marine20', 'marine20', p=0, resOffsets=666,resErrors=76) %>%
  as_tibble()
marine20_666

marine20_0 <- mixCurves('marine20') %>%
  as_tibble()
marine20_0

plot(marine20_666$CALBP, marine20_666$C14BP, type="l", col="red", xlab="Cal BP", ylab = "C14", 
     xlim = c(0,40000), ylim = c(0,40000))
par (new=TRUE)
plot(marine20_0$CALBP, marine20_0$C14BP, axes= FALSE, xlab="", ylab = "",
     type="l", col="black", lty = 2) #alpha (rgb(0,0,0), 0.5) - black transparent 50%
legend("topleft", legend=c("Marine20 ΔR=666±76", "Marine20"),
       col=c("red", "black"), lty=1:2, cex=0.8)
#Create a csv file to import into RChron without row and column names 
write.table(marine20_666, file="/Users/Steve/Dropbox/BAS/Data/R/RCarbon/marine20666.csv", 
            sep=",", row.names = FALSE, col.names=FALSE)
#Switch to BChron for instructions on how to use that within BChron

## Example of how to create a marine13 calibration curve with DR+/-err to import/use in e.g., BChron
marine13_791 <- mixCurves('marine13', 'marine13', p=0, resOffsets=791,resErrors=121) %>%
  as_tibble()
marine13_791

marine13_0 <- mixCurves('marine13') %>%
  as_tibble()
marine13_0

plot(marine13_791$CALBP, marine13_791$C14BP, type="l", col="red", xlab="Cal BP", ylab = "C14", 
     xlim = c(0,10000), ylim = c(0,10000))
par (new=TRUE)
plot(marine13_0$CALBP, marine13_0$C14BP, axes= FALSE, xlab="", ylab = "",
     type="l", col="black", lty = 2) #alpha (rgb(0,0,0), 0.5) - black transparent 50%
legend("topleft", legend=c("Marine13 ΔR=791±121", "Marine13"),
       col=c("red", "black"), lty=1:2, cex=0.8)
#Create a csv file to import into RChron without row and column names 
write.table(marine13_791, file="/Users/Steve/Dropbox/BAS/Data/R/RCarbon/marine13791.csv", 
            sep=",", row.names = FALSE, col.names=FALSE)
#Switch to BChron for instructions on how to use that within BChron


#++++++++++++++++++++++++++ SINGLE AGE CALIBRATION ++++++++++++++++++++++++++=
#single age calibration fro eg ARTIGAS 
x2 = calibrate(c(7045),c(85),calCurves=c('marine20'),
               resOffsets=c(666),resErrors=c(76),ids=1)
x2
summary(x2, calendar = "BP")
plot(x2)


#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Clear plots
if(!is.null(dev.list())) dev.off()



#+++++++++++++++++++++++++++ ARTIGAS MARINE AND TERRESTRIAL HOLOCENE DATA FOR FILDES PAPER ++++++++++++++++++++
#clear previous console
remove (list = ls())
#clear plot window
dev.off()

#Example using Artigas terrestrial then  marine and fildes peninsula advance (FPAD) Holocene radiocarbon ages 
#change this to x2 and x3 - rearrange into x1, x2, x3 start to finish process
x1 <- read.csv("ART_Terr_Holocene.csv")
x2 <- read.csv("ART_Marine_Holocene.csv")
x3 <- read.csv("Advance_Fildes_Holocene.csv")
x4 <- read.csv("Potter_Terr_Holocene.csv")
x5 <- read.csv("Fildes_Aq_Moss.csv")
x6 <- read.csv("Fildes_lakes.csv")
x7 <- read.csv("Advance_KGI_Holocene.csv")
head(x1)
head(x2)
head(x3)
head(x4)
head(x5)
head(x6)
head(x7)

#pooled dates from the same event prior to calibration into a weighted mean and std err
#to check for internal consistency before calibration
#*** no need to run this to carry on further ****

x1_C14 = x1$CRA
errors_x1 = x1$Error
id_x1 = x1$LabID
pool_x1 <- poolDates(x1_C14,errors_x1,id_x1)
pool_x1

x2_C14 = x2$CRA
errors_x2 = x2$Error
id_x2 = x2$LabID
pool_x2 <- poolDates(x2_C14,errors_x2,id_x2)
pool_x2

x3_C14 = x3$CRA
errors_x3 = x3$Error
id_x3 = x3$LabID
pool_x3 <- poolDates(x3_C14,errors_x3,id_x3)
pool_x3

x4_C14 = x4$CRA
errors_x4 = x4$Error
id_x4 = x4$LabID
pool_x4 <- poolDates(x4_C14,errors_x4,id_x4)
pool_x4

x5_C14 = x5$CRA
errors_x5 = x5$Error
id_x5 = x5$LabID
pool_x5 <- poolDates(x5_C14,errors_x5,id_x5)
pool_x5

x6_C14 = x6$CRA
errors_x6 = x6$Error
id_x6 = x6$LabID
pool_x6 <- poolDates(x6_C14,errors_x6,id_x6)
pool_x6

#Calibration
x1_cal <- calibrate(x1$CRA, x1$Error, normalised=TRUE, calCurves=x1$calCurves,
                      resOffsets=x1$resOffsets,resErrors=x1$resErrors,ids=x1$LabID)

x2_cal <- calibrate(x2$CRA, x2$Error, normalised=TRUE, calCurves=x2$calCurves,
                      resOffsets=x2$resOffsets,resErrors=x2$resErrors, ids=x2$LabID)

x3_cal <- calibrate(x3$CRA, x3$Error, normalised=TRUE, calCurves=x3$calCurves,
                    resOffsets=x3$resOffsets,resErrors=x3$resErrors,ids=x3$LabID)

x4_cal <- calibrate(x4$CRA, x4$Error, normalised=TRUE, calCurves=x4$calCurves,
                    resOffsets=x4$resOffsets,resErrors=x4$resErrors,ids=x4$LabID)

x5_cal <- calibrate(x5$CRA, x5$Error, normalised=TRUE, calCurves=x5$calCurves,
                    resOffsets=x5$resOffsets,resErrors=x5$resErrors,ids=x5$LabID)

x6_cal <- calibrate(x6$CRA, x6$Error, normalised=TRUE, calCurves=x6$calCurves,
                    resOffsets=x6$resOffsets,resErrors=x6$resErrors,ids=x6$LabID)

#check prob density distributions for each sample
head(x1_cal)
head(x2_cal)
head(x3_cal)
head(x4_cal)
head(x5_cal)
head(x6_cal)


#******* Output summary stats of prob density distributions for each ******
x1_sum <- summary(x1_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
x1_sum <- rename(x1_sum, LabID = DateID)
x1_cal_sum <- left_join(x1, x1_sum, by = c("LabID"))
x1_cal_sum

x2_sum <- summary(x2_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
x2_sum <- rename(x2_sum, LabID = DateID)
x2_cal_sum <- left_join(x2, x2_sum, by = c("LabID"))
x2_cal_sum

x3_sum <- summary(x3_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
x3_sum <- rename(x3_sum, LabID = DateID)
x3_cal_sum <- left_join(x3, x3_sum, by = c("LabID"))
x3_cal_sum

x4_sum <- summary(x4_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
x4_sum <- rename(x4_sum, LabID = DateID)
x4_cal_sum <- left_join(x4, x4_sum, by = c("LabID"))
x4_cal_sum

x5_sum <- summary(x5_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
x5_sum <- rename(x5_sum, LabID = DateID)
x5_cal_sum <- left_join(x5, x5_sum, by = c("LabID"))
x5_cal_sum

x6_sum <- summary(x6_cal, calendar = "BP") #prob = 0.954 #for just two-sigma range
x6_sum <- rename(x6_sum, LabID = DateID)
x6_cal_sum <- left_join(x6, x6_sum, by = c("LabID"))
x6_cal_sum

#++++++++ Write summary stats for CRA and median Cal BP to file ++++++++
x1_C14summary <- describe_distribution(x1_cal_sum$CRA)
x1_Calsummary  <- describe_distribution(x1_cal_sum$MedianBP)
x1_C14summary
x1_Calsummary 
write.csv(x1_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x1_C14summary.csv", row.names = FALSE)
write.csv(x1_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x1_Calsummary.csv", row.names = FALSE)
write.csv(x1_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x1_cal_sum.csv", row.names = FALSE)

x2_C14summary <- describe_distribution(x2_cal_sum$CRA)
x2_Calsummary  <- describe_distribution(x2_cal_sum$MedianBP)
x2_C14summary
x2_Calsummary 
write.csv(x2_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x2_C14summary.csv", row.names = FALSE)
write.csv(x2_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x2_Calsummary.csv", row.names = FALSE)
write.csv(x2_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x2_cal_sum.csv", row.names = FALSE)

x3_C14summary <- describe_distribution(x3_cal_sum$CRA)
x3_Calsummary  <- describe_distribution(x3_cal_sum$MedianBP)
x3_C14summary
x3_Calsummary 
write.csv(x3_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x3_C14summary.csv", row.names = FALSE)
write.csv(x3_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x3_Calsummary.csv", row.names = FALSE)
write.csv(x3_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x3_cal_sum.csv", row.names = FALSE)

x4_C14summary <- describe_distribution(x4_cal_sum$CRA)
x4_Calsummary  <- describe_distribution(x4_cal_sum$MedianBP)
x4_C14summary
x4_Calsummary 
write.csv(x4_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x4_C14summary.csv", row.names = FALSE)
write.csv(x4_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x4_Calsummary.csv", row.names = FALSE)
write.csv(x4_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x4_cal_sum.csv", row.names = FALSE)

x5_C14summary <- describe_distribution(x5_cal_sum$CRA)
x5_Calsummary  <- describe_distribution(x5_cal_sum$MedianBP)
x5_C14summary
x5_Calsummary 
write.csv(x5_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x5_C14summary.csv", row.names = FALSE)
write.csv(x5_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x5_C14summary.csv", row.names = FALSE)
write.csv(x5_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x5_cal_sum.csv", row.names = FALSE)

x6_C14summary <- describe_distribution(x6_cal_sum$CRA)
x6_Calsummary  <- describe_distribution(x6_cal_sum$MedianBP)
x6_C14summary
x6_Calsummary 
write.csv(x6_C14summary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x6_C14summary.csv", row.names = FALSE)
write.csv(x6_Calsummary,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x6_Calsummary.csv", row.names = FALSE)
write.csv(x6_cal_sum,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x6_cal_sum.csv", row.names = FALSE)

#++++++ Investigate subsets of each dataset +++++++
# Example for x1 only - *** no need to run this to carry on further ****
# Subsets of all dates that have a probability mass above 0.8 <6000 BP
x1_subset1 = which.CalDates(x1_cal,BP<6000,p=0.5)
x1_subset1
# Subsets of all dates that have a probability mass above 0.5 between 1000 and 2000 BP
x1_subset2 = which.CalDates(x1_cal,BP>1000&BP<2000,p=0.5)
x1_subset2


## ++++++++ Combine calibrated ages into a SPD - Summed Prob. Dist +++++++++
# Plot distributions based on Site Name column as the group 

# S3 method for class 'CalDates' - first create a bin based on site name to group by later on 
# first need to use cutree(tree, k = NULL, h = NULL) to look at and define groups according to height (h)
# Arguments tree - a cluster tree/dendrogram produced by hclust - allows you to look and see if 
# similar ages come from the same site and if there us similarity between the sites 
# cutree() only expects a list with components merge, height, and labels, of appropriate content each.k = an integer scalar or vector with the desired number of groups
# h = numeric scalar or vector with heights where the tree should be cut.



#++++++++++ x1 - Artigas Terrestrial dates Summed prob density analysis / plotting ++++++++

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using hierachical cluster analysis 
#This clustering method defines the cluster distance between two clusters to
# be the maximum distance between their individual components
x1_dd <- dist(scale(x1$CRA), method = "euclidean")
x1_hc <- hclust(x1_dd, method = "ward.D2")
cutree(x1_hc, h = 100)
as.dendrogram(x1_hc)
plot(x1_hc, labels=x1$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")


## Compare the grouping:
group <- cutree(x1_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(x1$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
x1_bins <- binPrep(x1$SiteName, x1$CRA, h=100)
summary(x1_cal, prob = 0.954, calendar = "BP")
x1_bins

#Compute median date for each bin
x1_bm <- binMed(x=x1_cal,bins=x1_bins)
x1_bm
#Compute median date for each date
x1_dm <- medCal(x1_cal)
x1_dm

#set up Age scale 
x1_timeRange <- c(12000,-100)
x1_revtimeRange <- c(-100,12000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
x1_spd <- spd(x1_cal, bins=x1_bins, runm=50, timeRange=x1_timeRange) #runm is running mean 
plot(x1_spd,runm=50, xlim = x1_revtimeRange) 
#medians for each date
barCodes(x1_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(x1_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
x1_res = stackspd(x=x1_cal,timeRange=x1_timeRange,bins=x1_bins,group=x1$SiteName)
plot(x1_res,type='lines', xlim = x1_revtimeRange)
#medians for each bin
barCodes(x1_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot
plot(x1_res,type='stacked', xlim = x1_revtimeRange, legend = TRUE, legend.arg = NULL)  #most useful for summaries of multiple sites / types
#plot(res1,type='proportion') #not that useful
plot(x1_res,type='multipanel', xlim = x1_revtimeRange, legend = TRUE, legend.arg = NULL) #nice clear plot for one site 


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(x1_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = x1_revtimeRange)
multiplot(x1_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = x1_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages - this doesnt work for x1 - ART terrestrial dataset
x1_s = sampleDates(x1_cal,bins=x1_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(x1_s,timeRange=x1_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',xlim = x1_revtimeRange)


# Clear plots
if(!is.null(dev.list())) dev.off()

#set up for plotting later on and write to csv for use in other programs
p1 <- plot(x1_spd,runm=50, xlim = x1_revtimeRange)
write.csv(x1_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x1_spd.csv", row.names = FALSE)


#+++++++++ x2 - ARTIGAS Marine dates ++++++++++++++

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using cluster analysis 
x2_dd <- dist(scale(x2$CRA), method = "euclidean")
x2_hc <- hclust(x2_dd, method = "ward.D2")
cutree(x2_hc, h = 100)
as.dendrogram(x2_hc)
plot(x2_hc, labels=x2$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the 2 and 4 grouping:
group24 <- cutree(x2_hc, k = c(2,4))
table(grp2 = group24[,"2"], grp4 = group24[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(x2$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
x2_bins <- binPrep(x2$SiteName, x2$CRA, h=600)
summary(x2_cal, prob = 0.954, calendar = "BP")
x2_bins

#Compute median date for each bin
x2_bm <- binMed(x=x2_cal,bins=x2_bins)
x2_bm
#Compute median date for each date
x2_dm <- medCal(x2_cal)
x2_dm

#set up Age scale 
x2_timeRange <- c(12000,-100)
x2_revtimeRange <- c(-100,12000)

# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
x2_spd <- spd(x2_cal, bins=x2_bins, runm=50, timeRange=x2_timeRange)
plot(x2_spd,runm=50, xlim = x2_revtimeRange)
#median for ages
barCodes(x2_dm,yrng=c(0,0.0001), width = 25, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#median for bins
barCodes(x2_bm,yrng=c(0,0.0001), width = 25, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
x2_res = stackspd(x=x2_cal,timeRange=x2_timeRange,runm=50,bins=x2_bins,group=x2$SiteName)
plot(x2_res,type='lines', xlim = x2_revtimeRange)
barCodes(x2_bm,yrng=c(0,0.01), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 
#other types of plot - can flip the axes at this point
plot(x2_res,type='stacked', xlim = x2_revtimeRange) #nice summary for all sites in one
#plot(res2,type='proportion') #rubbish
plot(x2_res,type='multipanel', xlim = x2_revtimeRange) #nice plot for individual sites


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(x2_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = x2_revtimeRange)
multiplot(x2_cal,HPD=TRUE,decreasing=TRUE,label=TRUE,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = x2_revtimeRange, cex.id = 0.5)

#plot a kde for calibrated ages 
x2_s = sampleDates(x2_cal,bins=x2_bins,nsim=100,boot=TRUE)
ckdeNorm = ckde(x2_s,timeRange=x2_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP', xlim = x2_revtimeRange)


# Clear plots
if(!is.null(dev.list())) dev.off()

#set up for plotting later on and write to csv for use in other programs
p2 <- plot(x2_spd,runm=50, xlim = x2_revtimeRange)
write.csv(x2_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x2_spd.csv", row.names = FALSE)

#++++ x3 - Fildes advance max age constraints ++++++++

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using cluster analysis 
x3_dd <- dist(scale(x3$CRA), method = "euclidean")
x3_hc <- hclust(x3_dd, method = "ward.D2")
cutree(x3_hc, h = 100)
as.dendrogram(x3_hc)
plot(x3_hc, labels=x3$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")


## Compare the grouping:
group <- cutree(x3_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(x3$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
x3_bins <- binPrep(x3$SiteName, x3$CRA, h=100)
summary(x3_sum, prob = 0.954, calendar = "BP")
x3_bins

#Compute median date for each bin
x3_bm <- binMed(x=x3_cal,bins=x3_bins)
x3_bm
#Compute median date for each date
x3_dm <- medCal(x3_cal)
x3_dm

#set up Age range for x axis 
x3_timeRange <- c(12000,-100)
x3_revtimeRange <- c(-100,12000)


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
x3_spd <- spd(x3_cal, bins=x3_bins, runm=50, timeRange=x3_timeRange) #runm is running mean 
plot(x3_spd,runm=50, xlim = x3_revtimeRange)
#add median ages
barCodes(x3_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#or add median for bins (smaller number)
barCodes(x3_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
x3_res = stackspd(x=x3_cal,timeRange=x3_timeRange,runm=50, bins=x3_bins,group=x3$SiteName)
plot(x3_res,type='lines', xlim = x3_revtimeRange)
barCodes(x3_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot - can flip the axes at this point
plot(x3_res,type='stacked', xlim = x3_revtimeRange) #nice summary for all sites in one - reverses 
#plot(res3,type='proportion') #rubbish
plot(x3_res,type='multipanel', xlim = x3_revtimeRange) #nice plot for individual sites


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(x3_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = x3_revtimeRange)
multiplot(x3_cal,HPD=TRUE,decreasing=TRUE,label=TRUE, cex.id = 0.5, 
          gapFactor = 2, xlab = "Age [a cal BP]", xlim = x3_revtimeRange)

#plot a kde for calibrated ages - this works
x3_s = sampleDates(x3_cal,bins=x3_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(x3_s,timeRange=x3_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,xlim = x3_revtimeRange, type='multiline',calendar='BP')


# Clear plots
if(!is.null(dev.list())) dev.off()

#set up for plotting later on and write to csv for use in other programs
p3 <- plot(x3_spd,runm=50, xlim = x3_revtimeRange)
write.csv(x3_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x3_spd.csv", row.names = FALSE)
#++++ x4 - Potter Holocene terr advance max age constraints ++++++++

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using cluster analysis 
x4_dd <- dist(scale(x4$CRA), method = "euclidean")
x4_hc <- hclust(x4_dd, method = "ward.D2")
cutree(x4_hc, h = 100)
as.dendrogram(x4_hc)
plot(x4_hc, labels=x4$ID, main = "Cluster Dendrrogram", hang = 0.1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(x4_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(x4$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
x4_bins <- binPrep(x4$SiteName, x4$CRA, h=100)
summary(x4_sum, prob = 0.954, calendar = "BP")
x4_bins

#Compute median date for each bin
x4_bm <- binMed(x=x4_cal,bins=x4_bins)
x4_bm
#Compute median date for each date
x4_dm <- medCal(x4_cal)
x4_dm

#set up Age range for x axis 
x4_timeRange <- c(12000,-100)
x4_revtimeRange <- c(-100,12000)


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
x4_spd <- spd(x4_cal, bins=x4_bins, runm=50, timeRange=x4_timeRange) #runm is running mean 
plot(x4_spd,runm=50, xlim = x4_revtimeRange)
#add median ages
barCodes(x4_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#or add median for bins (smaller number)
barCodes(x4_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
x4_res = stackspd(x=x4_cal,timeRange=x4_timeRange,runm=50, bins=x4_bins,group=x4$SiteName)
plot(x4_res,type='lines', xlim = x4_revtimeRange)
barCodes(x4_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot - can flip the axes at this point
plot(x4_res,type='stacked', xlim = x4_revtimeRange) #nice summary for all sites in one - reverses 
#plot(x4_res,type='proportion') #rubbish
plot(x4_res,type='multipanel', xlim = x4_revtimeRange) #nice plot for individual sites


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(x4_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,gapFactor = 0.5, 
          xlim = x4_revtimeRange)
multiplot(x4_cal,HPD=TRUE,decreasing=TRUE,label=TRUE, cex.id = 0.5, gapFactor = 2, 
          xlab = "Age [a cal BP]", xlim = x4_revtimeRange)

#plot a kde for calibrated ages - this works
x4_s = sampleDates(x4_cal,bins=x4_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(x4_s,timeRange=x4_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,xlim = x4_revtimeRange, type='multiline',calendar='BP')


# Clear plots
if(!is.null(dev.list())) dev.off()

#set up for plotting later on and write to csv for use in other programs
p4 <- plot(x4_spd,runm=50, xlim = x4_revtimeRange)
write.csv(x4_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x4_spd.csv", row.names = FALSE)

#++++ x5 - Aquatic moss layers in Fildes Lakes ++++++++

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using cluster analysis 
x5_dd <- dist(scale(x5$CRA), method = "euclidean")
x5_hc <- hclust(x5_dd, method = "ward.D2")
cutree(x5_hc, h = 100)
as.dendrogram(x5_hc)
plot(x5_hc, labels=x5$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(x5_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(x5$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
x5_bins <- binPrep(x5$SiteName, x5$CRA, h=100)
summary(x5_sum, prob = 0.954, calendar = "BP")
x5_bins

#Compute median date for each bin
x5_bm <- binMed(x=x5_cal,bins=x5_bins)
x5_bm
#Compute median date for each date
x5_dm <- medCal(x5_cal)
x5_dm

#set up Age range for x axis 
x5_timeRange <- c(12000,-100)
x5_revtimeRange <- c(-100,12000)


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
x5_spd <- spd(x5_cal, bins=x5_bins, runm=50, timeRange=x5_timeRange) #runm is running mean 
plot(x5_spd,runm=50, xlim = x5_revtimeRange)
#add median ages
barCodes(x5_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#or add median for bins (smaller number)
barCodes(x5_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
x5_res = stackspd(x=x5_cal,timeRange=x5_timeRange,runm=50, bins=x5_bins,group=x5$SiteName)
plot(x5_res,type='lines', xlim = x5_revtimeRange)
barCodes(x5_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot - can flip the axes at this point
plot(x5_res,type='stacked', xlim = x5_revtimeRange) #nice summary for all sites in one - reverses 
#plot(x5_res,type='proportion') #rubbish
plot(x5_res,type='multipanel', xlim = x5_revtimeRange) #nice plot for individual sites


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(x5_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,gapFactor = 0.5, 
          xlim = x5_revtimeRange)
multiplot(x5_cal,HPD=TRUE,decreasing=TRUE,label=TRUE, cex.id = 0.5, gapFactor = 2, 
          xlab = "Age [a cal BP]", xlim = x5_revtimeRange)

#plot a kde for calibrated ages - this doesnt work
x5_s = sampleDates(x5_cal,bins=x5_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(x5_s,timeRange=x5_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,xlim = x5_revtimeRange, type='multiline',calendar='BP')

# Clear plots
if(!is.null(dev.list())) dev.off()

#set up for plotting later on and write to csv for use in other programs
p5 <- plot(x5_spd,runm=50, xlim = x4_revtimeRange)
write.csv(x5_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x4_spd.csv", row.names = FALSE)

#++++ x6 - Fildes peninsula Lake basal age constraints on deglaciation ++++++++

# split into inner and outer fildes groups - classify according to this

# Clear plots
if(!is.null(dev.list())) dev.off()

#examine groups using cluster analysis 
x6_dd <- dist(scale(x6$CRA), method = "euclidean")
x6_hc <- hclust(x6_dd, method = "ward.D2")
cutree(x6_hc, h = 100)
as.dendrogram(x6_hc)
plot(x6_hc, labels=x6$ID, main = "Cluster Dendrrogram", hang = -1, cex=0.8, xlab = NULL, ylab = "Height")

## Compare the grouping:
group <- cutree(x6_hc, k = c(3,4))
table(grp2 = group[,"3"], grp4 = group[,"4"])

# test how many clusters are present in the data - 
# returns the number of clusters based on the maximum consensus. 
# In case of ties, it will select the solution with fewer clusters.
if (require("mclust", quietly = TRUE) &&
    require("cluster", quietly = TRUE)) {
  n_clusters(x6$CRA, package = c("mclust", "cluster"))
}

#create bins based using a suitable h split value based on maximum consensus in above  
x6_bins <- binPrep(x6$SiteName, x6$CRA, h=670)
summary(x6_cal, prob = 0.954, calendar = "BP")
x6_bins

#Compute median date for each bin
x6_bm <- binMed(x=x6_cal,bins=x6_bins)
x6_bm
#Compute median date for each date
x6_dm <- medCal(x6_cal)
x6_dm

#set up Age scale 
x6_timeRange <- c(12000,-100)
x6_revtimeRange <- c(-100,12000)


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Create and plot SPD and barCodes of median dates for each date and then each bin
x6_spd <- spd(x6_cal, bins=x6_bins, runm=50, timeRange=x6_timeRange) #runm is running mean 
plot(x6_spd,runm=50, xlim = x6_revtimeRange) 
#medians for each date
barCodes(x6_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))
#medians for each bin
barCodes(x6_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))

#Stack SPD plot based groups defined as sitename in this case
x6_res = stackspd(x=x6_cal,timeRange=x6_timeRange,bins=x6_bins,group=x6$SiteName)
plot(x6_res,type='lines', xlim = x6_revtimeRange)
#medians for each bin
barCodes(x6_bm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255))


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:2, ncol=1)) # Set up layout
par(mfrow=c(2,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#other types of plot
plot(x6_res,type='stacked', xlim = x6_revtimeRange) #most useful for summaries of multiple sites / types
#plot(res1,type='proportion') #not that useful
plot(x6_res,type='multipanel', xlim = x6_revtimeRange) #nice clear plot for one site 


# Clear plots
if(!is.null(dev.list())) dev.off()

#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 

#Calibrated age ranges summaries
multiplot(x6_cal,type='b',calendar='BP',cex.id = 0.5,lwd=2,
          gapFactor = 0.5, xlim = x6_revtimeRange)
multiplot(x6_cal,HPD=TRUE,decreasing=TRUE,label=TRUE, cex.id = 0.5,
          gapFactor = 0.5, xlab = "Age [a cal BP]", xlim = x6_revtimeRange)

#plot a kde for calibrated ages - this doesnt work a lot of the time
x6_s = sampleDates(x6_cal,bins=x6_bins,nsim=100,boot=FALSE)
ckdeNorm = ckde(x6_s,timeRange=x6_timeRange, bw=100,normalised=TRUE)
plot(ckdeNorm,type='multiline',calendar='BP',xlim = x6_revtimeRange)


# Clear plots
if(!is.null(dev.list())) dev.off()

#set up for plotting later on and write to csv for use in other programs
p6 <- plot(x6_spd,runm=50, xlim = x6_revtimeRange)
write.csv(x6_spd$grid,"/Users/Steve/Dropbox/BAS/Data/R/RCarbon/x6_spd.csv", row.names = FALSE)


#+++++++++++ SET UP for FINAL OUTPUT PLOTS +++++++++++
#clear plot window
dev.off()
# plot combinations of above - for Fildes paper plots needed to be 8 cm wide and -100 to 12000 years to match Sigmplot
# inches to cm
cm(4)
cm(1)  # = 2.54

#++++++ set up to plot to A4 plots A-C to pdf to match existing Sigmaplot plot width for -100 to 12 ka ++++++
## Translate 8 cm graph plot size *from* cm *to* inches:
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch
# if the pdf A4 width is 8.27 inches the l and right margins are then b
# total pdf width - plot width then divide by 2 to give margin width on each size
# then -1 to for mai = 1 inch left and right
b <- ((8.27-plotinch)/2)-1
b
checksize <- b*2+plotinch+2
checksize #should = 8.27 inches A4 width
#mai / omi = margin in inches or mar / oma = margin in lines - order is bottom, left, top, right 
#set up the page layout - A4 is 8.27 x 11.69 - closest to 
pdf("output1.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t
layout(matrix(1:6, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mai=c(0.5,1,0.5,1),omi=c(0,b,0,b), xaxs='i') # Set up internal and external margins

#clear plot window - to check plotting OK on screen first - dont rung this when plotting to PDF
dev.off()
#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up internal margins
#xas = "i" tells plot device not to add 4% extra internal margin (which default plot style in R) 

#set up common y axis SPD limit greater than maximum in all graphs
y <- c(0,0.03)

p1 <- plot(x1_spd,runm=50, xlim = x1_revtimeRange, ylim=y,
           main="a) Artigas: SPD of terrestrial moss ages in moraines (shcal20) [this study]",
           panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.03,0.005), labels=rep("",7), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.03,0.005), labels=rep("",7), tck=0.04))
par(new = TRUE) #hold the plot frame to add median values
barCodes(x1_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255)) #red & 50% transparency

par(new = FALSE) #remove hold on the plot frame to add median values
p2 <- plot(x2_spd,runm=50, xlim = x2_revtimeRange, ylim=y,
           main="b) Artigas: SPD of Marine shell ages in moraines (shcal20) [this study & Hall, 2007]",
           panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted", p_axis))
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.03,0.005), labels=rep("",7), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.03,0.005), labels=rep("",7), tck=0.04))
par(new = TRUE) #hold the frame to add median values
barCodes(x2_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255)) #red & 50% transparency

par(new = FALSE) #remove hold on the plot frame to add median values
p3 <- plot(x3_spd,runm=50, xlim = x3_revtimeRange, ylim=y,
           main="c) Fildes Peninsula: SPD max. age constraints for glacier readvance  [this study & Hall, 2007]",
           panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
p_axis <- c(axis(1, seq(0,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, seq(0,0.03,0.005), labels=rep("",7), tck=-0.04), 
            axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.03,0.005), labels=rep("",7), tck=0.04))
par(new = TRUE) #hold the plot frame to add median values
barCodes(x3_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255)) #red & 50% transparency


#++++++ set up to plot to A4 plots D-F ++++++

#clear plot window - to check plotting OK on screen first - don't run this when plotting to PDF
dev.off()
#screen plot layout 
layout(matrix(1:3, ncol=1)) # Set up layout
par(mfrow=c(3,1))  # this sets up the graphics window to expect a 1x3 layout 
par(mar=c(4,4,2,2), oma=c(1,1,2,1),xaxs='i') # Set up margins, 
#xas = "i" tells device not to add 4% extra internal margin (which default plot style in R) 


p4 <- plot(x4_spd,runm=50, xlim = x4_revtimeRange, ylim=y,
           main="d) Potter Peninsula: Max. age constraints for glacier readvance  [this study & Hall, 2007]",
           panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
p_axis <- c(axis(1, seq(-100,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, tck=-0.04), axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),
            axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.03,0.005), labels=rep("",7), tck=0.04))
par(new = TRUE) #hold the frame to add median values
barCodes(x4_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255)) #red & 50% transparency

par(new = FALSE) #remove hold on the plot frame to add median values
p5 <- plot(x5_spd,runm=50, xlim = x5_revtimeRange, ylim=y,
           main="e) Fildes Lakes: Aquatic moss layer ages",
           panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
p_axis <- c(axis(1, seq(-100,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, tck=-0.04), axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),
            axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.03,0.005), labels=rep("",7), tck=0.04))
par(new = TRUE) #hold the frame to add median values
barCodes(x5_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255)) #red & 50% transparency

par(new = FALSE) #remove hold on the plot frame to add median values

p6 <- plot(x6_spd,runm=50, xlim = x6_revtimeRange, ylim=y,
           main="e) Fildes Lakes: deglaciation basal age constraints",
           panel.first = grid(nx = NULL, col = "lightgray", lty = "dotted"))
p_axis <- c(axis(1, seq(-100,12000,1000), tck=-0.04), axis(1, seq(-100,12000,100),labels=rep("",122), tck=-0.02),
            axis(2, tck=-0.04), axis(3, seq(0,12000,1000), labels=rep("",13), tck=0.04),
            axis(3, seq(-100,12000,100),labels=rep("",122), tck=0.02),
            axis(4, seq(0,0.03,0.005), labels=rep("",7), tck=0.04))
par(new = TRUE) #hold the frame to add median values
barCodes(x6_dm,yrng=c(0,0.001), width = 50, col = rgb(255, 0, 0, 50, maxColorValue = 255)) #red & 50% transparency

dev.off() #need to include this to write to pdf file fully - will delete screen plot


