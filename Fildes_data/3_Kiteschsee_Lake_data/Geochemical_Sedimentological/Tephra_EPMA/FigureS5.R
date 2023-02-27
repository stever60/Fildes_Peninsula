# Set up ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/")
#check working directory
getwd()

library(tidyverse)
library(tidypaleo)
library(rioja)
library(repr)
library(patchwork)


# structure of datasets for comparison - from https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/ 
data("alta_lake_geochem")
alta_lake_geochem

# Set plotting parameters for base R  -------------

## set universal plot size
## Translate cm graph plot size *from* cm *to* inches:
plotinch_x <- 2 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch_y <- 7 / cm(1) # -> 8 cm  is  3.149606 inches
# or use aspect ratio
aspect_ratio <- 2.5
options(repr.plot.width=plotinch_x, repr.plot.height=plotinch_y)



# Figure S5 ----------------------------------------------------------------

# Import data  and data wrangling ------------------------------------------------------------

# import data, remove unwanted columns, rename and recalculate total shards by dividing by 1, 10, 100 1000 to enable easier plotting
# Tephra count csv datafile is here: https://www.dropbox.com/s/itx5xyaxwwzfpib/KITE_Tephra_Counts.csv?dl=0 
shards <- read_csv("/Users/Steve/Dropbox/BAS/Data/R/Papers/Heredia_Barion_2020/EPMA/Data/KITE_Tephra_Counts.csv") %>% 
  select(Lab_ID, strat_depth, Ab_Counts, Total_shards_gDM) %>% 
  mutate(Ab_Counts2 = Ab_Counts/1000) %>% 
  mutate(Total_shards_gDM2 = Total_shards_gDM/1000) %>%
  rename(depth = strat_depth, Shard_Counts = Ab_Counts2, Total_Shards = Total_shards_gDM2) %>% 
  select(Lab_ID, depth, Shard_Counts, Total_Shards)
shards

# Convert data to long format for ggplot and tidypaleo facet plotting --------------------------

# select out the parameters to plot,  and convert to long format - i.e., one value per row per variable 
shards_long <- select(shards, Lab_ID, depth, Shard_Counts, Total_Shards) %>%
  pivot_longer(c(`Shard_Counts`, `Total_Shards`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
shards_long

# Plot as point and line dataplot -----------------------------------------
theme_set(theme_paleo(base_size=12))

# Plotting as horizontal bars rather than line+points --------------------------------

#clear plot window
dev.off()

p3 <- ggplot(shards_long,
               aes(x = value, y = depth, ymin= 0, ymax = 80)) +
  geom_col_segsh(colour = "blue") +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  #facet_abundanceh(vars(param), grouping = NULL) +
  xlab(expression (Number~(x10^{3})))  +
  #xlab(expression (NULL)) +
  ylab(expression (Depth~(cm))) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(title = "Kiteschsee Lake") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("Shard_Counts" = "n", "Total_Shards" = "no. per g DM", "CONISS" = "SS")) +
  theme(text=element_text(size=16, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

p3

# add tephra zone shading for two main tephra layers
zone_data1 <- tibble(ymin = 55, ymax = 63, xmin = -Inf, xmax = Inf)
zone_data2 <- tibble(ymin = 31.5, ymax = 35, xmin = -Inf, xmax = Inf)

# add tephra zone shading for two main tephra layers
p4 <- p3 +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data1, 
    alpha = 0.2,
    fill = "grey",
    inherit.aes = FALSE
  ) +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data2, 
    alpha = 0.2,
    fill = "grey",
    inherit.aes = FALSE
  ) +
  # add dashed line to highlight peaks and two main tephra zones 
  geom_hline(yintercept = c(33, 58), col = "red", lty = 2, lwd = 1, alpha = 0.5) +
  theme(text=element_text(size=16, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="plain"),
        plot.margin = unit(c(1,0,1,1), "cm")
  )

p4

# adding dengrogram of layers with CONISS constrained cluster analysis from Rioja package
coniss <- shards_long %>%
  nested_data(qualifiers = c(depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

coniss_plot <- ggplot() +
  layer_dendrogram(coniss, aes(y = depth, ymin= 0, ymax = 80)) +
  layer_zone_boundaries(coniss, aes(y = depth)) +
  scale_y_reverse() +
  facet_geochem_gridh(vars("CONISS"),
                      units = c("CONISS" = "SS")) +
  xlab(expression (SSquares))  +
  #labs(x = NULL) +
  theme(text=element_text(size=16, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,0), "cm")
  )

# coniss_plot

# this wraps the previous plots + the CONISS plot - on the left hand side using patchwork package
wrap_plots(
  p4 + layer_zone_boundaries(coniss, aes(y = depth)) +
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  coniss_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(6, 1) # this is relative width scale 
)

## save plot in working directory
ggsave("FigureS5C_bars.pdf", height = c(8), width = c(8), dpi = 600) #units = "cm"


# Plot as points and lines  - did this first but swapped around as bar plot is better-----------------------------------------------


p1 <- shard_plot <- ggplot(shards_long, 
                           aes(x = value, y = depth, ymin = 0, ymax = 80)) +
  geom_lineh() +
  geom_point(size = 1) +
  # reverse y axis to show zero depth at top
  scale_y_reverse() +
  # facet wrap according to the number of parameters & add x and y axis titles 
  facet_geochem_gridh(vars(param)) +
  # xlab(expression (Number~(x10^{-3})))  +
  xlab(expression (NULL)) +
  ylab(expression (Depth~(cm))) +
  labs(title = "Kiteschsee Lake") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("Shard_Counts" = "n", "Total_Shards" = "no. per g DM", "CONISS" = "SS")
  ) +
  theme(text=element_text(size=20, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

p1

# axis.text=element_text(size=12),
# axis.title=element_text(size=12,face="plain"),
# legend.text=element_text(size=12), 
# legend.title=element_text(size=12),
# plot.margin = unit(c(1,1,1,1), "cm")

# add tephra zone shading for two main tephra layers
zone_data1 <- tibble(ymin = 55, ymax = 63, xmin = -Inf, xmax = Inf)
zone_data2 <- tibble(ymin = 32, ymax = 35, xmin = -Inf, xmax = Inf)

p2 <- p1 +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data1, 
    alpha = 0.2,
    fill = "grey",
    inherit.aes = FALSE
  ) +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data2, 
    alpha = 0.2,
    fill = "grey",
    inherit.aes = FALSE
  ) +
  # add dashed line to highlight peaks and two main tephra zones 
  geom_hline(yintercept = c(33, 58), col = "red", lty = 3, lwd = 1, alpha = 0.7) +
  theme(text=element_text(size=20, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

p2

# adding dengrogram of layers with CONISS constrained cluster analysis from Rioja package 
# with optimum number of clusters determined by broken stick analysis
coniss <- shards_long %>%
  nested_data(qualifiers = c(depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

p2.1 <- p2 + 
  layer_dendrogram(coniss, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss, aes(y = depth))


p2.1

## save plot in working directory
ggsave("Users/Steve/Dropbox/BAS/Data/R/Papers/Heredia_Barion_2020/EPMA/FigureS5C_line.pdf", height = c(12, 12, 12), width = c(8, 8, 4), dpi = 600) #units = "cm"






# Other stuff - for reference -------------------------------------------------------------

# this also plots CONISS but on the right hand side 
p5 <- p4 +
  layer_dendrogram(coniss, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss, aes(y = depth)) +
  facet_geochem_gridh(
    vars(param),
    units = c("Shard_Counts" = "n", "Total_Shards" = "no. per g DM", "CONISS" = "SS")) +
  theme(text=element_text(size=20, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
  )


p5

# Add secondary age scale  ------------------------------------------------

# import age-depth model 
adm_bacon <- read_csv("/Users/Steve/Dropbox/BAS/Data/R/Papers/Heredia_Barion_2020/EPMA/KITEM4.csv") 
adm_bacon

adm <- age_depth_model(
  adm_bacon, 
  depth = depth_cm,
  age = SH20_mean
)

# add a secondary age axis - this doesnt work for KITE - try to fix
p1 +
  scale_y_depth_age(
    adm,
    age_name = "Age (cal a BP)"
  )


