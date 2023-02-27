# Heredia_Barion_et_al_2020_Fildes_Holocene

install.packages("performance")
install.packages("see")
install.packages("tidyquant")

# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# setup workspace ----
library(tidyverse) # all core tidyverse packages
library(readr)
library(ggpubr)
library(GGally) # for correlation and Prob density matrix plotting
library(tidyquant)
citation("performance")

# Raw data folder is located here: 
# https://www.dropbox.com/sh/smvdhb4wyx4q9i6/AAD7pxA0RQiZAAZtEi8Q17cFa?dl=0

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Heredia_Barion_2020/Data/")
#check working directory
getwd() 

# Tidyquant has 6 types of moving average:

# simple moving averages (SMA)
# exponential moving averages (EMA)
# weighted moving averages (WMA)
# double exponential moving averages (DEMA)
# zero-lag exponential moving averages (ZLEMA)
# volume-weighted moving averages (VWMA)
# elastic, volume-weighted moving averages (EVMA)

# Precipitation

dp <- read_csv("ERA5_Fildes_precip_mmyr-1_1979_2021.csv")
tail(dp)

fildes_precip <- ggplot(dp, aes(x = Year_CE, y = Total_precipitation)) +
  geom_line(color = "darkgrey") +
  geom_point(color = "darkgrey") + 
  geom_ma(ma_fun = SMA, n = 5, color = "red", linetype = "solid", size = 1.5) +  # Plot 5-year SMA
  #coord_x_date(xlim = c("1999-01-01", "2013-08-01")) + # Zoom into certain section if have monthly dates
  ylab(bquote('Precipitation ('~ mm~yr^-1*')')) +
  labs(x = "Year CE", title = "Mean Annual Precipitation (1979-2021 CE)") +
  xlim(1969,2021) + 
  annotate("rect", fill = "grey", alpha = 0.25, xmin = c(1982,1995), xmax = c(1984,1997), ymin = -Inf, ymax = Inf) +
  theme_tq(base_size=10) +
  theme(
    plot.title = element_text(color="black", size=10, face="bold.italic"),
    axis.title.x = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=10))
fildes_precip

# Wind
dw <- read_csv("Bellingshausen_All_wind_speed_kmh-1_1969_2021.csv")
head (dw)
tail(dw)

fildes_wind <- ggplot(dw, aes(x = Year_CE, y = Annual_Mean)) +
  geom_line(color = "darkgrey") +
  geom_point(color = "darkgrey") + 
  geom_ma(ma_fun = SMA, n = 5, color = "red", linetype = "solid", size = 1.5) +  # Plot 5-year SMA
  #coord_x_date(xlim = c("1999-01-01", "2013-08-01")) + # Zoom into certain section if have monthly dates
  ylab(bquote('Velocity ('~ km~h^-1*')')) +
  labs(x = "Year CE", title = "Mean Annual Wind Velocity (1968-2021 CE)") +
  xlim(1969,2021) + 
  annotate("rect", fill = "grey", alpha = 0.25, xmin = c(1982,1995), xmax = c(1984,1997), ymin = -Inf, ymax = Inf) +
  theme_tq(base_size=10) +
  theme(
    plot.title = element_text(color="black", size=10, face="bold.italic"),
    axis.title.x = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=10))
fildes_wind 

# + scale_x_continuous(breaks=seq(1960,2030,10)) +
#  scale_x_continuous(breaks=seq(22,30,1))


# Temperature

dt <- read_csv("Bellingshausen_All_temperature_1969_2021.csv")
head (dt)
tail(dt)

fildes_temp <- ggplot(dt, aes(x = Year_CE, y = Annual_Mean)) +
  geom_line(color = "darkgrey") +
  geom_point(color = "darkgrey") + 
  geom_ma(ma_fun = SMA, n = 5, color = "red", linetype = "solid", size = 1.5) +  # Plot 5-year SMA
  #coord_x_date(xlim = c("1999-01-01", "2013-08-01")) + # Zoom into certain section if have monthly dates
  ylab(expression(paste("Temperature (",degree~C,")"))) +
  labs(x = "Year CE", title = "Mean Annual Temperature (1969-2021 CE)") +
  xlim(1969,2021) + 
  annotate("rect", fill = "grey", alpha = 0.25, xmin = c(1982,1995), xmax = c(1984,1997), ymin = -Inf, ymax = Inf) +
  theme_tq(base_size=10) +
  theme(
    plot.title = element_text(color="black", size=10, face="bold.italic"),
    axis.title.x = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=10))
fildes_temp

ggarrange(fildes_temp, fildes_wind, fildes_precip, nrow = 3, ncol = 1, align = "v", labels = c('A', 'B', 'C', 'D'))
ggsave("Figures/Fig 1G_temp_wind_vell.pdf", 
              height = c(16), width = c(10), dpi = 600, units = "cm")
