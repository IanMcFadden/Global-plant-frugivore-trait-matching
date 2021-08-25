
## Code to reproduce Fig. 3 of:
##
##  McFadden et al. Global plant-frugivore trait matching is shaped by
##  climate and biogeographic history. Ecology Letters (in revision).

## Note: An image editor was used to create the final composite figure
## from the Fig. 3A, 3B and 3C PDF outputs from this script 


# Setup 
rm(list=ls())
library(viridis); library(patchwork); library(tidyverse); library(rgeos)
library(sp); library(raster); library(scales); library(colorspace)
load("data/all-data-SPDF-list_3spp-bc_ELE-R2.Rdata")

# Load in the Robinson projection map objects
load("data/sp_bbox_robinProj.Rdata")
load("data/sp_worldMap_robinProj.Rdata")

# Load the botanical country shapefiles
load("data/bc_shapeFiles_wherePalms_simplified.Rdata")

# Make Wallace realm color dataframe
realm <- c("Neotropical", "Nearctic", "Indotropical", "Afrotropical", "Panamanian",
           "Oceanina", "Sino-Japanese", "Madagascan", "Australian", "Saharo-Arabian")
hex_color <- c("#4A79B1", "#77C2DF", "#F4D432", "#649F4F", "#759bd9", "#724B9B",
               "#BDD790", "#682D1F", "#BE3F91", "#7EB854")
wallace_realm_cols_df <- data.frame(realm=realm, hex_color=hex_color)

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

# Co-plot traits using dot-in-circle plot (Fig. 3A):
# Dot is plant trait, circle is bird trait 

# First split out maximum trait value dataframe and remove unnecessary columns
trait_points_max_3spp <- as.data.frame(trait_clim_sp_points_df_list_3spp[["maximum"]])
trait_points_max_3spp$coords.x1 <- NULL
trait_points_max_3spp$coords.x2 <- NULL
head(trait_points_max_3spp)

# Rescale the values between 0 and 1 for plotting
trait_points_max_3spp$AverageFruitSize_cm_w_allo_log_rescale <- rescale(trait_points_max_3spp$AverageFruitSize_cm_w_allo_log, newrange=c(0,1))
trait_points_max_3spp$Average.Gape.Width_cm_log_rescale <- rescale(log(trait_points_max_3spp$Average.Gape.Width / 10), newrange=c(0,1))

# Make the color scale
trait_cols_3spp <- viridis(nrow(trait_points_max_3spp), option="C")

# Make the plot for Fig. 3A
pdf("tables_and_figures/Fig_3A.pdf", 10, 4)

ggplot(trait_points_max_3spp, aes(x=x, y=y)) +
  
  # Plot the world map
  geom_polygon(data=sp_bbox_robinProj, aes(x=long, y=lat), colour="white", fill="lightskyblue", size=0.25) +
  geom_polygon(data=sp_worldMap_robinProj, aes(long,lat, group=group), colour="white", fill="white", size=0.2) +
  geom_polygon(data=sp_bbox_robinProj, aes(x=long, y=lat), colour="grey10", fill="transparent", size=0.35) +
  coord_fixed(ratio = 1) +
  theme_void() +

  # Plot the bird trait values as hollow circles 
  geom_point(fill="black", shape=21, size=1.8, stroke=1) +
  geom_point(aes(colour=Average.Gape.Width_cm_log_rescale), fill=NA, shape=21, size=1.5, stroke=1) +
  
  # Plot the palm trait values as dots 
  geom_point(fill="black", size=1.3, pch=20) + 
  geom_point(aes(colour=AverageFruitSize_cm_w_allo_log_rescale), size=1, pch=20) + 
  scale_color_viridis(option="C", name="") 

dev.off()

# Plot gape and fruit size regressions (Fig. 3B)

# Calculate the correlation and significance value 
cor_3spp <- cor.test(log(trait_points_max_3spp$Average.Gape.Width/10),
    trait_points_max_3spp$AverageFruitSize_cm_w_allo_log,
    method="pearson", use="pairwise.complete.obs")
r_val_3spp <- as.numeric(round(cor_3spp$estimate, 2))
cor_3spp$p.value # p is < 0.001, so using ***


# Reorder wallace_realm_cols_df to match trait_points_max and create color vector
wallace_realm_cols_df_3spp <- wallace_realm_cols_df[order(unique(trait_points_max_3spp$wallace_realm)),]
wallace_realm_cols_vec_3spp <- wallace_realm_cols_df_3spp$hex_color
names(wallace_realm_cols_vec_3spp) <- wallace_realm_cols_df_3spp$realm

# Plot gape and fruit size regression (Fig. 3B)
gs_fs_reg_3spp <- ggplot(data=trait_points_max_3spp,
                         aes(x=AverageFruitSize_cm_w_allo_log,
                             y=log(Average.Gape.Width/10),
                             fill=wallace_realm)) +
  
                 geom_abline(aes(intercept=0, slope=1), size=0.5, color="darkgrey", linetype="dashed") + 
                 geom_point(pch=21, size=2) +
                 scale_fill_manual(name="", values=c(lapply(as.list(wallace_realm_cols_vec_3spp), lighten, 0.25))) + 
                 geom_smooth(method="lm", se=T, col="black", fill="grey70", lwd=0.9) +
                 theme_classic() +
                 theme(legend.position="none") +
                 xlab("Fruit size (log cm)") +
                 ylab("Gape size (log cm)") +
                 theme(plot.margin=unit(c(0.5, 2.5, 1, 0.4), "cm")) +
                 theme(axis.text.x=element_text(size=11, color="black", vjust=-0.8)) +
                 theme(axis.title.x=element_text(vjust=-2, size=13)) + 
                 theme(axis.text.y=element_text(size=11, color="black", hjust=0.5)) +
                 theme(axis.title.y=element_text(vjust=3, size=13)) + 
                 annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
                          label=paste0("r = ", r_val_3spp, "***"), size=5.3) +
                 xlim(range(trait_points_max_3spp$AverageFruitSize_cm_w_allo_log)[1]-0.3,
                      range(trait_points_max_3spp$AverageFruitSize_cm_w_allo_log)[2]+0.3) + 
                 ylim(range(log(trait_points_max_3spp$Average.Gape.Width/10))[1]-0.3,
                      range(log(trait_points_max_3spp$Average.Gape.Width/10))[2]+0.3)

# Plot gape and fruit size richness residuals (Fig. 3C)
# Create the residual dataframe for plotting 
gs_fs_max_rich_resid_df <- trait_points_max_3spp
gs_fs_max_rich_resid_df <- gs_fs_max_rich_resid_df[,c("LEVEL3_NAM",
                                                      "wallace_realm",
                                                      "palm_richness",
                                                      "bird_richness",
                                                      "AverageFruitSize_cm_w_allo_log",
                                                      "Average.Gape.Width")]

gs_fs_max_rich_resid_df$Average.Gape.Width_cm <- gs_fs_max_rich_resid_df$Average.Gape.Width / 10
gs_fs_max_rich_resid_df$Average.Gape.Width <- NULL
gs_fs_max_rich_resid_df$fruit_size_rich_resids_log <- NA
gs_fs_max_rich_resid_df$gape_size_rich_resids_log <- NA
head(gs_fs_max_rich_resid_df)

# Calculate fruit size ~ palm richness residuals
fruit_size_rich_lm <- lm(gs_fs_max_rich_resid_df$AverageFruitSize_cm_w_allo_log~log(gs_fs_max_rich_resid_df$palm_richness))
gs_fs_max_rich_resid_df$fruit_size_rich_resids_log <- fruit_size_rich_lm$residuals

# Calculate gape size ~ bird richness residuals
gape_size_rich_lm <- lm(log(gs_fs_max_rich_resid_df$Average.Gape.Width_cm)~log(gs_fs_max_rich_resid_df$bird_richness))
gs_fs_max_rich_resid_df$gape_size_rich_resids_log <- gape_size_rich_lm$residuals
head(gs_fs_max_rich_resid_df)

# Extract the correlation between residuals
# Calculate the correlation and significance value 
resid_cor_3spp <- cor.test(gs_fs_max_rich_resid_df$gape_size_rich_resids_log,
                           gs_fs_max_rich_resid_df$fruit_size_rich_resids_log,
                           method="pearson", use="pairwise.complete.obs")
resid_r_val_3spp <- as.numeric(round(resid_cor_3spp$estimate, 2))
resid_cor_3spp$p.value # p is < 0.001, so using ***

# Plot gape and fruit size richness residuals (Fig. 3C)
gs_fs_reg_3spp_richResid <- ggplot(data=gs_fs_max_rich_resid_df, aes(x=fruit_size_rich_resids_log,
                                                                       y=gape_size_rich_resids_log, 
                                                                       fill=wallace_realm)) +
  
                 geom_abline(aes(intercept=0, slope=1), size=0.5, color="darkgrey", linetype="dashed") + 
                 geom_point(pch=21, size=2) +
                 scale_fill_manual(name="", values=c(lapply(as.list(wallace_realm_cols_vec_3spp), lighten, 0.25))) + 
                 geom_smooth(method="lm", se=T, col="black", fill="grey70", lwd=0.9) +
                 theme_classic() +
                 theme(legend.position = "none") +
                 xlab("Palm richness residuals") +
                 ylab("Bird richness residuals") +
                 theme(axis.text.x=element_text(size=11, color="black", vjust=-0.8)) +
                 theme(axis.title.x=element_text(vjust=-2, size=13)) + 
                 theme(axis.text.y=element_text(size=11, color="black", hjust=0.5)) +
                 theme(axis.title.y=element_text(vjust=3, size=13)) + 
                 annotate("text", x=Inf, y=-Inf, hjust=1, vjust=-1, 
                          label=paste0("r = ", resid_r_val_3spp, "***"), size=5.3) +
                 xlim(range(gs_fs_max_rich_resid_df$fruit_size_rich_resids_log)[1]-0.2,
                      range(gs_fs_max_rich_resid_df$fruit_size_rich_resids_log)[2]+0.2) + 
                 ylim(range(gs_fs_max_rich_resid_df$gape_size_rich_resids_log)[1]-0.2,
                      range(gs_fs_max_rich_resid_df$gape_size_rich_resids_log)[2]+0.2)

# Make the combined Fig. 3B and C plot
pdf("tables_and_figures/Fig_3B_and_C.pdf", 8.5, 3.5)
gs_fs_reg_3spp + gs_fs_reg_3spp_richResid
dev.off()


## End of script
