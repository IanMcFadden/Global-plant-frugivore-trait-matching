
## Code to reproduce Fig. 2, Table S1 and Fig. S6 of:
##
##  McFadden et al. Global plant-frugivore trait matching is shaped by
##  climate and biogeographic history. Ecology Letters (in press).

## Note: An image editor was used to create Figs. 2 and S6 using the results
## from this script (SEM_summary_list, see below)


# Setup 
rm(list=ls())
library(sp); library(lme4); library(raster); library(piecewiseSEM)
load("data/all-data-SPDF-list_3spp-bc_ELE-R2.Rdata")
trait_clim_spdf_list <- trait_clim_spdf_list_3spp; rm(trait_clim_spdf_list_3spp)

# Convert gape width from mm to log cm 
for(i in 1:length(trait_clim_spdf_list)){
  trait_clim_spdf_list[[i]]$Average.Gape.Width_log_cm <- log(trait_clim_spdf_list[[i]]$Average.Gape.Width / 10)
  trait_clim_spdf_list[[i]]$Average.Gape.Width <- NULL
}

# Function to prepare the data and then fit a piecewise SEM
prep_n_fit_pieceSEM <- function(df){
  
  # Subset to relevant columns
  df <- df[,c("wallace_realm", "LEVEL3_NAM", "MAT", "NPP",
              "palm_richness", "AverageFruitSize_cm_w_allo_log", "Average.Gape.Width_log_cm")]

  # Convert from a spatial dataframe to a regular dataframe
  df <- as.data.frame(df)
  
  # Make all values of temperature positive by adding the mean (so no NAs are produced when logging)
  df$MAT <- df$MAT + mean(df$MAT)
  
  # Log transform and scale continuous variables that were not logged
  rescale_cols <- c("palm_richness", "MAT", "NPP")
  df[, rescale_cols] <- scale(log(df[, rescale_cols]))
  
  # Scale the variables that were already logged
  df[, c("AverageFruitSize_cm_w_allo_log", "Average.Gape.Width_log_cm")] <- scale(df[, c("AverageFruitSize_cm_w_allo_log",
                                                                                         "Average.Gape.Width_log_cm")])
  
  # Mixed effect component model with bird gape size as response
  gape_size <- lmer(Average.Gape.Width_log_cm ~ AverageFruitSize_cm_w_allo_log + palm_richness + MAT + (1|wallace_realm), data=df)
  
  # Mixed effect component model with palm fruit size as response
  fruit_size <- lmer(AverageFruitSize_cm_w_allo_log ~ MAT + NPP + palm_richness + (1|wallace_realm), data=df)
  
  # Mixed effect component model with palm richness as response
  palm_rich <- lmer(palm_richness ~ MAT + NPP + (1|wallace_realm), data=df)
  
  # Make the full SEM and return it 
  full_model <- psem(gape_size, fruit_size, palm_rich)
  return(full_model)
  
}  

# Use the function to fit the piecewise SEMs for maximum and median botanical country values
SEM_mod_list <- lapply(trait_clim_spdf_list[c("maximum", "median")], prep_n_fit_pieceSEM)
SEM_summary_list <- lapply(SEM_mod_list, summary)

# Name the summaries
SEM_summary_list[["maximum"]]$name <- "maximum"
SEM_summary_list[["median"]]$name <- "median"

# Print the results
SEM_summary_list[["maximum"]] # Used as input for image editor make Fig. 2 and Table S1
SEM_summary_list[["median"]]  # Used as input for image editor make Fig. S6 and Table S1

# Save the results
save(SEM_summary_list, file="tables_and_figures/Fig_2_and_S6_and_Table_S1.Rdata")


## End of script
