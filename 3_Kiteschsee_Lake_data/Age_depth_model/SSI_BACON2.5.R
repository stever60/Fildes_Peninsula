#New BACON v2.5 INTCAL20 / SH20 / MARINE20

###To reinstall - copy existing BACON runs folder away first###
#install.packages('rbacon')
#only run this once 
#Bacon()
#this will load the default core (MSB2K), create a default Bacon_runs folder 
#in your top level directory and put data/outputs it in there
#old BACON location - setwd("/Users/Steve/Dropbox/BAS/Data/C14/MacBacon_2.2")

# Set working directories ----------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Mac
setwd("/Users/Steve/Dropbox/BAS/Data/R/RBacon")
#check working directory
getwd()
library(rbacon)

# Windows
setwd("C:/Users/sjro/Dropbox/BAS/Data/R/RBacon")
#check working directory
getwd()
library(rbacon)



# SSI ---------------------------------------------------------------------

### Kiteschee Lake

Bacon.cleanup()

## FINAL MODELS

# Clear plots
if(!is.null(dev.list())) dev.off()
# Translate cm size graph plot size *from* cm *to* inches to match eg illustrator/SP size:
plotwidthinch <- 6.5507 / cm(1) # -> 8 cm  is  3.149606 inches
plotheightinch <- 6.3083 / cm(1) # -> 8 cm  is  3.149606 inches
# total pdf width - plot width then divide by 2 to give margin width on each size
# then -1 to for mai = 1 inch left and right
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.25,0.25,0.25), pin=c(plotwidthinch,plotheightinch),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 
## FINAL PLOT TO PDF
#mai / omi = margin in inches or mar / oma = margin in lines - order is bottom, left, top, right 
#set up the page layout - A4 is 8.27 x 11.69 - closest to 
pdf("output.pdf", pointsize=12,  width = 8.27, height = 11.69) #will fit to t


## KITE-M4: lake reservoir age of 700±100 years applied to bulk sediment age data above 20 cm,  hiatus at 5 cm 320 years
# KITE4.1 same as above but with calibrated modelled tephra ages (light blue) included as test / for plotting 
# Run with: Y5a only Y5a age below T5 deposit at 62 cm included; then Y5ab Y5a age below T5 deposit at 62 cm included and Y5b above at 55 cm 
Bacon("KITE-M4Y5ab",depths.file=FALSE,thick=3,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50,rounded=2,
      mem.mean=0.5,mem.strength=10,yr.max=12000,
      slump=c(10.5, 12.1, 32, 34.3, 58.5, 63.5),
      d.min=0,d.max=80, hiatus.depths=5, hiatus.max=200,#boundary=5,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=15, plot.pdf = TRUE, title.location='topright',mgp=c(1.5, 0.7, 0))

dev.off() #need to include this to write to pdf file fully - will delete screen plot

## +++++++++ PmC data
Bacon.hist(50)
pMC.age(101.63, 0.38)
pMC.age(101.25, 0.38)


## INITIAL MODEL
## KITE-M1: no lake reservoir, no hiatus at 5 cm, no tephra correlation ages included
Bacon("KITE-M1",depths.file=FALSE,thick=3,rotate.axes=TRUE,cc=3,
      postbomb=4,acc.mean=100, 
      mem.mean=0.4,mem.strength=20,yr.max=12000,
      d.min=0,d.max=80,rounded=2,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=15, title.location='topright',mgp=c(1.5, 0.7, 0))

Bacon.cleanup()
# Clear plots
if(!is.null(dev.list())) dev.off()
## OTHER FINAL MODELS USED
## Potter Cove Pingfo and PC river section
#mjb samples set as 2.8 m depth - mean of Pingfo section depths
Bacon("PC",depths.file=FALSE,thick=5,rotate.axes=TRUE,cc=1,
      postbomb=4,acc.mean=0.5,
      mem.mean=0.4,mem.strength=20,yr.min=6500,yr.max=12000,
      d.min=0,d.max=500,rounded=2, dark=0.5,
      height=0.5, title.location='topright',mgp=c(1.5, 0.7, 0))

#+++++++Yanou Lake

# originally using DR=791+/-121 from Pablo paper instead of 700+/-100 from Roberts et al. (2017) for marine ages
# now DR = 666+76 recalculated for Marine20 - see Fildes paper Supp Material
# new age-model has slump added rather than correlation to Ardley Lake ages within T5 ash layer
Bacon("YANM5",depths.file=TRUE,thick=5,d.by=1,rotate.axes=TRUE, 
      postbomb=4,rounded=2,
      mem.mean=0.4,mem.strength=20,
      acc.mean=20,yr.max=13000,
      d.max=360, slump=c(40,189,252,256), #PCA>0.2 for terrestrial zone ash layers 
      res=1,dark=0.5,
      height=15, title.location='topright',mgp=c(1.5, 0.7, 0))

cc = copyCalibrationCurve(cc=3)

#+++++++Ardley Lake - SH20

# originally using DR=791+/-121 from Pablo paper instead of 700+/-100 from Roberts et al. (2017) for marine ages
# now DR = 666+76 recalculated for Marine20 - see Fildes paper Supp Material
# new age-model has slump added rather than correlation to Ardley Lake ages within T5 ash layer
Bacon("ARD-M6",depths.file=TRUE,thick=5,d.by=1,rotate.axes=TRUE, 
      postbomb=4,rounded=2,
      mem.mean=0.4,mem.strength=20,
      acc.mean=50,yr.max=9000,
      d.max=360, #PCA>0.2 for terrestrial zone ash layers 
      res=1,dark=0.3,
      height=150, title.location='topright',mgp=c(1.5, 0.7, 0))

# cc = copyCalibrationCurve(cc=3)


Bacon.cleanup()
#Other KITE models
## KITE-M6: all bulk ages with 700±100 years lake reservoir effect and T6 used instead of T7 age at base
Bacon("KITE-M6",depths.file=TRUE,thick=3,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=100,rounded=2, 
      mem.mean=0.4,mem.strength=20,yr.max=12000,
      d.min=0,d.max=77,hiatus.depths=5,hiatus.max=320,
      title.location='topright',dark=0.5)
## KITE-M5: without T5 correlation ages to test T5 ages produced by the model
Bacon("KITE-M5",depths.file=FALSE,thick=3,rotate.axes=TRUE,cc=3,
      postbomb=4,acc.mean=100,rounded=2,
      mem.mean=0.4,mem.strength=20,yr.max=12000,
      d.min=0,d.max=80,hiatus.depths=5,hiatus.max=320, #320
      title.location='bottomright',dark=0.5)
## KITE-M3: no lake reservoir, mean hiatus of 720±200 years at 5 cm
Bacon("KITE-M3",depths.file=FALSE,thick=3,rotate.axes=TRUE,cc=3,
      postbomb=4,acc.mean=100,rounded=2,
      mem.mean=0.4,mem.strength=20,yr.max=12000,
      d.min=0,d.max=80,hiatus.depths=c(5),hiatus.max=720,
      title.location='bottomright',dark=0.5)
## KITE-M2: lake reservoir age of 560±80 years applied to bulk sediment age data above 20 cm, no hiatus at 5 cm
Bacon("KITE-M2",depths.file=FALSE,thick=3,rotate.axes=TRUE,cc=3,
      postbomb=4,acc.mean=100, 
      mem.mean=0.4,mem.strength=20,yr.max=12000,
      d.min=0,d.max=80,rounded=2,
      title.location='bottomright',dark=0.5)

#++++++ Ardley Lake 
#ARD_M5c includes a marine reservoir effect of 700+/-50 years for collagen (bone), guano influenced sediments & and brackish sedimentary units
#see csv files for data input used in M1-M5 runs 
Bacon("ARD-M5c",title.location='topright',res=5,thick=5,depths.file=FALSE,
      rotate.axes=TRUE,cc=3, postbomb=4,acc.mean=50,rounded=2,
      mem.mean=0.4,mem.strength=20,yr.max=9000,d.max=360,dark=0.8)

#tephra layers and reworked tephra deposits greater than 3 cm thick max depth
#hiatus does the opposite of slump in CLAM
#hiatus mean is the age offset in years
#hiatus.depths=c(58,167,280,345),hiatus.mean=100


# Scotia Sea --------------------------------------------------------------

## PC460 - zero age at zero depth - DR = 2300+/-300 (surface age offset)
Bacon("PC460",depths.file=FALSE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=10,rounded=2,
      mem.mean=0.5,mem.strength=10,yr.max=15000,
      d.min=0,d.max=1000, 
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=5, title.location='topright',mgp=c(1.5, 0.7, 0))




