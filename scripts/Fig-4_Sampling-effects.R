
## Code to reproduce Fig. 4 of:
##
##  McFadden et al. Global plant-frugivore trait matching is shaped by
##  climate and biogeographic history. Ecology Letters (in revision).

## Note: An image editor was used to create the figure using the results
## from this script (SEM_summary_list, see below)


# Setup 
rm(list=ls())
library(sp); library(raster); library(tidyverse); library(patchwork)
load("data/all-data-SPDF-list_3spp-bc_ELE-R2.Rdata")
trait_clim_spdf_list <- trait_clim_spdf_list_3spp; rm(trait_clim_spdf_list_3spp)

# Convert gape width from mm to log cm 
for(i in 1:length(trait_clim_spdf_list)){
  trait_clim_spdf_list[[i]]$Average.Gape.Width_log_cm <- log(trait_clim_spdf_list[[i]]$Average.Gape.Width / 10)
  trait_clim_spdf_list[[i]]$Average.Gape.Width <- NULL
}

# 1. Fruit size range ~ palm richness
# Calculate R2 and p-value
R2_1_lm <- lm(trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log ~ log(trait_clim_spdf_list[["range"]]$palm_richness))
R2_1 <- summary(R2_1_lm) # lm is best model
R2_1$coefficients[2,4] # p is < 0.001, so using ***
R2_1 <- round(R2_1$r.squared, 2)
R2_1 <- as.expression(bquote(italic(R)^2 ~ "=" ~ .(R2_1) ~ "***"))

p1 <- ggplot(data=data.frame(x=log(trait_clim_spdf_list[["range"]]$palm_richness), y=trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log),
             aes(x=x, y=y)) +
  geom_point(size=0.75, color="#00938e") + 
  geom_smooth(method="lm", col="#00938e", fill="#00938e") +
  xlab("Palm richness (log)") +
  ylab("Fruit size range (log cm)") +
  annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, size=4.5,
           label=R2_1) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, color="black", vjust=-0.8)) +
  theme(axis.title.x=element_text(vjust=-2, size=12)) + 
  theme(axis.text.y=element_text(size=10, color="black", hjust=-3)) +
  theme(axis.title.y=element_text(vjust=3, hjust=0.15, size=12)) + 
  theme(axis.ticks=element_line(color="black")) + 
  theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

# 2. Gape size range ~ palm richness
# Calculate R2 and p-value
R2_2_poly2 <- lm(trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm ~ poly(log(trait_clim_spdf_list[["range"]]$palm_richness), 2))
R2_2 <- summary(R2_2_poly2) # 2nd degree polynomial is the best model 
R2_2$coefficients[2,4] # p is < 0.001, so using ***
R2_2 <- round(R2_2$r.squared, 2)
R2_2 <- as.expression(bquote(italic(R)^2 ~ "=" ~ .(R2_2) ~ "***"))

p2 <- ggplot(data=data.frame(x=log(trait_clim_spdf_list[["range"]]$palm_richness), y=trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm),
       aes(x=x, y=y)) +
  geom_point(size=0.75, color="#6c4b8d") + 
  geom_smooth(method="lm", col="#6c4b8d", fill="#6c4b8d", formula = y ~ poly(x, 2)) +
  xlab("Palm richness (log)") +
  ylab("Gape size range (log cm)") +
  annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
           label=R2_2, size=4.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, color="black", vjust=-0.8)) +
  theme(axis.title.x=element_text(vjust=-2, size=12)) + 
  theme(axis.text.y=element_text(size=10, color="black", hjust=0.5)) +
  theme(axis.title.y=element_text(vjust=3, hjust=0.15, size=12)) + 
  theme(axis.ticks=element_line(color="black")) + 
  theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

# 3. Gape size range ~ fruit size range 
# Calculate R2 and p-value
R2_3_poly2 <- lm(trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm ~ poly(trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log, 2))
R2_3 <- summary(R2_3_poly2) # 2nd degree polynomial is about equivalent, and the points have some curve so using this model
R2_3$coefficients[2,4] # p is < 0.001, so using ***
R2_3 <- round(R2_3$r.squared, 2)
R2_3 <- as.expression(bquote(italic(R)^2 ~ "=" ~ .(R2_3) ~ "***"))

p3 <- ggplot(data=data.frame(x=trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log, y=trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm),
             aes(x=x, y=y)) +
  geom_point(size=0.75, color="#6c4b8d") + 
  geom_smooth(method="lm", col="#6c4b8d", fill="#6c4b8d", formula = y ~ poly(x, 2)) +
  xlab("Fruit size range (log cm)") +
  ylab("Gape size range (log cm)") +
  annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
           label=R2_3, size=4.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, color="black", vjust=-0.8)) +
  theme(axis.title.x=element_text(vjust=-2, size=12)) + 
  theme(axis.text.y=element_text(size=10, color="black", hjust=0.5)) +
  theme(axis.title.y=element_text(vjust=3, hjust=0.15, size=12)) + 
  theme(axis.ticks=element_line(color="black")) + 
  theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Make the figure
pdf("tables_and_figures/Fig_4.pdf", 9*1.05, 3.2*1.05)
p1 + p2 + p3
dev.off()


## End of script
