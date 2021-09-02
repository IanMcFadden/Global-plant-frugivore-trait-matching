
## Code to reproduce Fig. 5 of:
##
##  McFadden et al. Global plant-frugivore trait matching is shaped by
##  climate and biogeographic history. Ecology Letters (in revision).

## Note: An image editor was used to create the final composite figure
## from the Fig. 5A, 5B and 5C PDF outputs from this script 


# Setup 
rm(list=ls())
library(viridis); library(patchwork); library(tidyverse); library(rgeos)
library(sp); library(raster); library(scales); library(colorspace)
load("data/all-data-SPDF-list_3spp-bc_ELE-R2.Rdata")

# Load the botanical country shapefiles
load("data/bc_shapeFiles_wherePalms_simplified.Rdata")

# Load in the Robinson projection map objects
load("data/sp_bbox_robinProj.Rdata")
load("data/sp_worldMap_robinProj.Rdata")

# Make shape ID dataframe and add ID to spatial polygons dataframe
id_df <- bc_shape_palms_simp[,c("LEVEL3_COD", "id")]
for(i in 1:length(trait_clim_spdf_list_3spp)){
  trait_clim_spdf_list_3spp[[i]] <- merge(trait_clim_spdf_list_3spp[[i]], id_df, by="LEVEL3_COD")  
}
lapply(trait_clim_spdf_list_3spp, head)

# Create spatial points dataframes for plotting 
trait_clim_sp_points_df_list_3spp <- trait_clim_spdf_list_3spp

for(i in 1:length(trait_clim_sp_points_df_list_3spp)){
  trait_clim_sp_points_df_list_3spp[[i]] <-
    SpatialPointsDataFrame(coords=cbind(trait_clim_sp_points_df_list_3spp[[i]]$centroid_long,
                                        trait_clim_sp_points_df_list_3spp[[i]]$centroid_lat),
                           data=as.data.frame(trait_clim_sp_points_df_list_3spp[[i]])) 
}
lapply(trait_clim_sp_points_df_list_3spp, class)

# Add projections back to the dataframes 
for(i in 1:length(trait_clim_sp_points_df_list_3spp)){
  proj4string(trait_clim_sp_points_df_list_3spp[[i]]) <- proj4string(trait_clim_spdf_list_3spp[[1]])
}
lapply(trait_clim_sp_points_df_list_3spp, proj4string)

# Project spatial points DF to Robinson projection
for(i in 1:length(trait_clim_sp_points_df_list_3spp)){
  trait_clim_sp_points_df_list_3spp[[i]]  <- spTransform(trait_clim_sp_points_df_list_3spp[[i]],
                                                    proj4string(sp_bbox_robinProj))
}
lapply(trait_clim_sp_points_df_list_3spp, proj4string) 

# Calculate new centroid positions 
for(i in 1:length(trait_clim_sp_points_df_list_3spp)){
    shape_centroid <- gCentroid(trait_clim_sp_points_df_list_3spp[[i]], byid = TRUE)
    trait_clim_sp_points_df_list_3spp[[i]]$x <- shape_centroid$x
    trait_clim_sp_points_df_list_3spp[[i]]$y <- shape_centroid$y
} 
lapply(trait_clim_sp_points_df_list_3spp, head) 

# Split out the maximum trait value dataframe and remove unnecessary columns
trait_points_max_3spp <- as.data.frame(trait_clim_sp_points_df_list_3spp[["maximum"]])
trait_points_max_3spp$coords.x1 <- NULL
trait_points_max_3spp$coords.x2 <- NULL
head(trait_points_max_3spp)

# Residual plot of gape width ~ fruit length (Fig. 5A)
# Colored by residual value, point size is fitted value 
gs_fs_mod_3spp <- lm(log(Average.Gape.Width/10) ~ AverageFruitSize_cm_w_allo_log, data=trait_points_max_3spp)
trait_points_max_3spp$Residual <- resid(gs_fs_mod_3spp)
trait_points_max_3spp$Fitted <- gs_fs_mod_3spp$fitted.values

# Plot residuals against absolute latitude (Fig. 5A)
mod_3spp <- lm(abs(trait_points_max_3spp$Residual)~abs(trait_points_max_3spp$centroid_lat))
lat_cor <- cor.test(abs(trait_points_max_3spp$Residual), abs(trait_points_max_3spp$centroid_lat))
lat_cor$p.value # p is < 0.001, so using ***
r_3spp <- round(lat_cor$estimate, 2)

# Make the plot
gs_fs_lat_3spp <- ggplot(trait_points_max_3spp, aes(x=abs(centroid_lat), y=abs(Residual))) +
  
    geom_point(aes(fill=Residual), size=2.2, pch=21) +
    scale_fill_gradient2(low="dodgerblue3", mid="white", high="brown3") + 
    theme_classic() +
    theme(legend.position="none") +
    xlab("Absolute latitude") +
    ylab("Absolute GS ~ FS residual") + 
    geom_smooth(method="lm", se=T, col="black", fill="grey70", lwd=0.9) +
    theme(plot.margin=unit(c(0,0.5,1.5,0.5), "cm")) +
    theme(axis.text.x = element_text(size=11, color="black", vjust=-0.8)) +
    theme(axis.text.y = element_text(size=11, color="black", hjust=0.5)) +
    theme(axis.title.x = element_text(vjust=-2, size=13)) + 
    theme(axis.title.y = element_text(vjust=3, size=13)) + 
    xlim(range(abs(trait_points_max_3spp$centroid_lat))[1]-4,
         range(abs(trait_points_max_3spp$centroid_lat))[2]+3) + 
    ylim(range(abs(trait_points_max_3spp$Residual))[1]-0.15,
         range(abs(trait_points_max_3spp$Residual))[2]+0.11) + 
    annotate("text", x=4.5, y=1.195, label=paste0("r = ", r_3spp, "***"), size=5) # CHANGE THIS!

dev.off()


# Plot residuals by Wallace realm (Fig. 5B)

# Calculate realm mean residual values and add to the dataframe 
trait_points_max_3spp$lm_resids_realmAve <- ave(trait_points_max_3spp$Residual, trait_points_max_3spp$wallace_realm)
head(trait_points_max_3spp)

# Reorder the dataset by increasing realm-level residual value
trait_points_max_3spp <- trait_points_max_3spp %>%
  mutate(wallace_realm=fct_reorder(wallace_realm, -lm_resids_realmAve))

# Make the plot
set.seed(30)
locs <- 1:10

# Make the plot
gs_fs_resids_by_region_BC_max_3spp <- ggplot(trait_points_max_3spp, aes(x=wallace_realm, y=Residual)) +

    geom_jitter(aes(color=Residual), size=3, alpha=1, width=0.1, pch=20) +
    scale_color_gradient2(low="dodgerblue3", mid="white", high="brown3", labels=NULL) +
    
    geom_vline(xintercept=1:132, color="black", size=0.4, lty="dotted") +
    geom_hline(yintercept=0, color="black", size=0.4, lty="dashed") +
    
    stat_summary(aes(fill=lm_resids_realmAve), fun=mean, geom="point", size=4, pch=21) +
    scale_fill_gradient2(low="dodgerblue3", mid="white", high="brown3", labels=NULL) +
    
    coord_flip() +
    scale_y_continuous(limits=range(trait_points_max_3spp$Residual), expand=c(0.04, 0.003)) +
    labs(x=NULL, y="Gape size ~ fruit size residuals") +
    
    theme_light(base_size=13) +
    theme(legend.position="none",
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          axis.title=element_text(size=13),
          plot.caption=element_text(size=9, color="black"),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.text=element_text(colour="black"),
          axis.ticks.x=element_line(colour="black"), 
          axis.text.x=element_text(vjust=-0.4),
          axis.title.x=element_text(vjust=-2),
          axis.ticks.y=element_line(colour="white"))

# Make multipanel regression figure (Fig. 5A and B)
pdf("tables_and_figures/Fig_5A_and_B.pdf", 7*1.1, 3.5*1.1)
gs_fs_lat_3spp + gs_fs_resids_by_region_BC_max_3spp
dev.off()

# Plot residuals in space (Fig. 5C)
pdf("tables_and_figures/Fig_5C.pdf", 7,4)
ggplot(trait_points_max_3spp, aes(x=x, y=y)) +
  
      # Plot the world map
      geom_polygon(data=sp_bbox_robinProj, aes(x=long, y=lat), colour="white", fill="white", size=0.25) +
      geom_polygon(data=sp_worldMap_robinProj, aes(long,lat, group=group), colour="grey70", fill="grey70", size=0.2) +
      geom_polygon(data=sp_bbox_robinProj, aes(x=long, y=lat), colour="grey10", fill="transparent", size=0.35) +
      coord_fixed(ratio = 1) +
      theme_void() +
      theme(legend.position=c(0.17,0.425),
            legend.key.size=unit(.3, "cm"),
            legend.title=element_text(size=9.5)) +
      guides(fill=guide_colourbar(ticks = FALSE)) +
      theme(plot.margin = unit(c(0,0,0,0), "cm")) +
      
      # Plot the points
      geom_point(aes(fill=Residual, size=Fitted), pch=21) + 
      guides(size=guide_legend("Fitted GS")) + 
      scale_size("Fitted", range=c(0.7, 3), breaks=c(0.6, 1, 1.6))  +
      scale_fill_gradient2(low="dodgerblue3", mid="white", high="brown3", labels=NULL) +
      annotate("text", x=-11800000, y=700000, label="GS > FS", size=3.3) +
      annotate("text", x=-11800000, y=-1300000, label="GS < FS", size=3.3)

dev.off()


# Determine the average residual value of tropical forest-dominated realms
# (Afrotropical, Indotropical, Neotropical etc.)
head(trait_points_max_3spp)
realm_ave_resids <- with(trait_points_max_3spp, tapply(lm_resids_realmAve, wallace_realm, unique)) 
mean(realm_ave_resids[names(realm_ave_resids) %in% c("Afrotropical", "Indotropical", "Neotropical")])


## End of script
