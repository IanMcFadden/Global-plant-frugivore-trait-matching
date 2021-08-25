
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
r_1 <- cor.test(trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log, log(trait_clim_spdf_list[["range"]]$palm_richness))
r_1$p.value # p is < 0.001, so using ***
r_1 <- round(r_1$estimate, 2)
r_1 <- paste0("r = ", r_1, "***")
r_1

p1 <- ggplot(data=data.frame(x=log(trait_clim_spdf_list[["range"]]$palm_richness), y=trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log),
             aes(x=x, y=y)) +
  geom_point(size=0.75, color="#00938e") + 
  geom_smooth(method="lm", col="#00938e", fill="#00938e") +
  xlab("Palm richness") +
  ylab("Fruit size range") +
  annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
           label=r_1, size=4.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, color="black", vjust=-0.8)) +
  theme(axis.title.x=element_text(vjust=-2, size=12)) + 
  theme(axis.text.y=element_text(size=10, color="black", hjust=0.5)) +
  theme(axis.title.y=element_text(vjust=3, size=12)) + 
  theme(axis.ticks=element_line(color="black")) + 
  theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

# 2. Gape size range ~ palm richness
r_2 <- cor.test(trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm, log(trait_clim_spdf_list[["range"]]$palm_richness))
r_2$p.value # p is < 0.001, so using ***
r_2 <- round(r_2$estimate, 2)
r_2 <- paste0("r = ", r_2, "***")
r_2

p2 <- ggplot(data=data.frame(x=log(trait_clim_spdf_list[["range"]]$palm_richness), y=trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm),
       aes(x=x, y=y)) +
  geom_point(size=0.75, color="#6c4b8d") + 
  geom_smooth(method="lm", col="#6c4b8d", fill="#6c4b8d") +
  xlab("Palm richness") +
  ylab("Gape size range") +
  annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
           label=r_2, size=4.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, color="black", vjust=-0.8)) +
  theme(axis.title.x=element_text(vjust=-2, size=12)) + 
  theme(axis.text.y=element_text(size=10, color="black", hjust=0.5)) +
  theme(axis.title.y=element_text(vjust=3, size=12)) + 
  theme(axis.ticks=element_line(color="black")) + 
  theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

# 3. Gape size range ~ fruit size range 
r_3 <- cor.test(trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm, trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log)
r_3$p.value # p is < 0.001, so using ***
r_3 <- round(r_3$estimate, 2)
r_3 <- paste0("r = ", r_3, "***")
r_3

p3 <- ggplot(data=data.frame(x=trait_clim_spdf_list[["range"]]$AverageFruitSize_cm_w_allo_log, y=trait_clim_spdf_list[["range"]]$Average.Gape.Width_log_cm),
             aes(x=x, y=y)) +
  geom_point(size=0.75, color="#6c4b8d") + 
  geom_smooth(method="lm", col="#6c4b8d", fill="#6c4b8d") +
  xlab("Fruit size range") +
  ylab("Gape size range") +
  annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
           label=r_3, size=4.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, color="black", vjust=-0.8)) +
  theme(axis.title.x=element_text(vjust=-2, size=12)) + 
  theme(axis.text.y=element_text(size=10, color="black", hjust=0.5)) +
  theme(axis.title.y=element_text(vjust=3, size=12)) + 
  theme(axis.ticks=element_line(color="black")) + 
  theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Make the figure
pdf("tables_and_figures/Fig_4.pdf", 9, 3.2)
p1 + p2 + p3
dev.off()


## End of script


