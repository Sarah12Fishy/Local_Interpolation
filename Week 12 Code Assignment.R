#############################################
#     Week 12 Assignment 
#############################################

# Install Packages
install.packages("geoR")
install.packages("gstat")
install.packages("sp")
install.packages("raster")
install.packages("RColorBrewer")

# Load necessary libraries
library(geoR)          # Geostatistical analysis
library(gstat)         # Spatial interpolation
library(sp)            # Spatial data manipulation
library(raster)        # Raster data handling
library(RColorBrewer)  # Color palettes for plotting

### Step 1: Load and Inspect the Dataset

#----------------
data(soja98)     # Complete: Load the dataset
help(soja98)     # Complete: View documentation for the dataset
head(soja98)     # Complete: Inspect first few rows of the dataset


#----------------

### Step 2: Select Relevant Columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----------------
soja98 <- soja98[, c(1, 2, 3)]  # Complete: Select relevant columns (X, Y, P)

str(soja98)     # To view data


############ Question 1 on Assignment #############
# 1)	List the summary statistics of the soja98 phosphorus levels. Include mean (average), median, maximum, and minimum

# Calculate summary statistics for the P column
mean_P <- mean(soja98$P, na.rm = TRUE)  # Calculate mean
median_P <- median(soja98$P, na.rm = TRUE)  # Calculate median
max_P <- max(soja98$P, na.rm = TRUE)  # Calculate maximum
min_P <- min(soja98$P, na.rm = TRUE)  # Calculate minimum

# Display the results
summary_stats <- list(
  Mean = mean_P,
  Median = median_P,
  Maximum = max_P,
  Minimum = min_P
)

print(summary_stats)


#----------------

### Step 3: Convert to Spatial Points
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----------------
coordinates(soja98) <- ~X+Y  # Complete: Convert to spatial points using X and Y columns
#----------------

### Step 4: Split Data into Calibration and Validation Sets

set.seed(4887)                             # Ensures reproducibility in random sampling
IDX <- sample(1:nrow(soja98), 150)         # Complete: Randomly select indices for calibration data

cal <- soja98[IDX, ]                       # Calibration dataset (selected indices)
val <- soja98[-IDX, ]                      # Validation dataset (remaining indices)


###############################################
#           Inverse Distance Weighting (IDW)
###############################################

### Step 1: Create a Grid for Interpolation

cellSize <- 5  # Define grid cell size (5 meters)
range.all <- extent(soja98)  # Get spatial extent of the dataset

grd <- expand.grid(x=seq(range.all[1]-cellSize/2, range.all[2], by=cellSize), 
                   y=seq(range.all[3], range.all[4], by=cellSize)) # Complete: Generate Y coordinates

coordinates(grd) <- ~x + y   # Convert to spatial points object
gridded(grd) <- TRUE         # Convert to gridded structure


### Step 2: Perform IDW Interpolation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

P.idw <- idw(P~1, cal, grd)  # Complete: Specify formula, calibration data, and grid



### Step 3: Visualize IDW Results

cols <- brewer.pal(9, "BuPu")  # Generate color ramp for visualization

spplot(P.idw, "var1.pred", col.regions = cols, cuts = 8, main = 'Inverse Distance Weighting')


#plotting commands
cols <- brewer.pal(9, 'BuPu') #generates colorramp
spp.idw <- spplot(P.idw, "var1.pred", col.regions=cols, cuts = 8) #generates dummy map

cuts <- spp.idw$panel.args.common$at 
#gets coloramp bins from dummy map
lev <- cut(cal$P, breaks=cuts) #determines colors for points ftom cal
pts.1 <- list("sp.points", cal, pch = 24, cex=2, col='black', lwd = 2)
pts.2 <- list("sp.points", cal, pch = 17, cex=2, 
              col=cols[lev], lwd = 1)
#generates points to represent the calibration locations

spplot(P.idw, "var1.pred", col.regions=cols, cuts = 8,
       sp.layout=list(pts.2, pts.1),main = 'Inverse Distance')
#plots final map


### Step 4: Evaluate Accuracy with RMSE

r.idw <- raster(P.idw)               # Convert IDW output to a raster object

p.idw.cal <- extract(r.idw, cal)     # Extract interpolated values at calibration points
RMSE.idw.cal <- sqrt(mean((p.idw.cal-cal$P)^2, na.rm=T))  # Complete: Apply RMSE formula

p.idw.val <- extract(r.idw, val)     # Extract interpolated values at validation points
RMSE.idw.val <- sqrt(mean((p.idw.val-val$P)^2, na.rm=T))  # Complete: Apply RMSE formula

cat("IDW Calibration RMSE:", RMSE.idw.cal)
cat("IDW Validation RMSE:", RMSE.idw.val)



### Step 5: Use Yardstick Package for Alternative RMSE Calculation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

install.packages("yardstick")

library(yardstick)
library(dplyr)


rmse_df <- tibble(
  truth = val$P,
  estimate = p.idw.val
)
yardstick_rmse <- rmse(rmse_df, truth = truth, estimate = estimate)

cat("Validation RMSE using yardstick:", yardstick_rmse$.estimate)



###############################################
#           Effect of Sample Size on RMSE
###############################################

### Step 1: Initialize Sample Sizes and RMSE Array

n <- seq(5, 150, by = 5)  # Complete: Specify increment for sample sizes
RMSE <- array(NA, length(n))      # Array to store RMSE values


### Step 2: Simulation Loop


for (i in seq_along(n)) { #random sample of calibration 
  # Randomly sample calibration points
  samp <- sample(1:length(cal), n[i])  # Complete: Specify population and sample size
  
  # Perform IDW interpolation with sampled data
  P.idw <- idw(P~1, cal[samp ,], grd)   # Complete: Specify formula, data, and grid
  
  # Convert IDW output to raster and extract validation values
  r.idw <- raster(P.idw)
  p.idw.val <- extract(r.idw, val)
  
  # Calculate RMSE for validation data
  RMSE[i] <- sqrt(mean((p.idw.val-val$P)^2, na.rm=T))  # Complete: Apply RMSE formula
}


### Step 3: Plot Results


plot(n, RMSE, type = "b", 
     main = "Effect of Sample Size on RMSE",
     xlab = "Sample Size", 
     ylab = "RMSE")


###############################################
#           Ordinary Kriging (OK)
###############################################

### Step 1: Fit Variogram Models
#

variog <- variogram(P~1, cal) # Complete: Specify formula and data
sphr.fit <- fit.variogram(variog, vgm('Sph'))  # Complete: Fit spherical model
exp.fit <- fit.variogram(variog, vgm('Exp'))  # Complete: Fit exponential model



### Step 2: Plot Variogram Models

gamma <- variog$gamma                        # Experimental semivariance values
bins <- variog$dist                          # Lag distances

gamma.sphr <- variogramLine(sphr.fit, dist_vector = bins)$gamma   # Spherical model predictions
gamma.exp <- variogramLine(exp.fit, dist_vector = bins)$gamma     # Exponential model predictions

plot(bins, gamma, main = "Variogram Models", 
     xlab = "Lag Distance", ylab = "Semivariance")
lines(bins, gamma.sphr, col='blue')  # Complete: Add spherical model line (color = 'blue')
lines(bins, gamma.exp, col='red')  # Complete: Add exponential model line (color = 'red')
legend("topright", 
       legend = c("Spherical", "Exponential"), 
       col = c("blue", "red"), 
       lty = 1)



### Step 3: Perform Kriging
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

P.sph <- krige(P~1, cal, grd, model=sphr.fit) # Complete: Specify formula, data, grid, and model
P.exp <- krige(P~1, cal, grd, model=exp.fit)  # Complete: Specify formula, data, grid, and model



### Step 4: Visualize Kriging Results

spplot(P.sph, main = "Ordinary Kriging (Spherical)")
spplot(P.exp, main = "Ordinary Kriging (Exponential)")


### Step 5: Evaluate Kriging Accuracy
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Spherical model
r.sph <- raster(P.sph)
p.sph.cal <- extract(r.sph, cal)
RMSE.sph.cal <- sqrt(mean((p.sph.cal-cal$P)^2, na.rm=T))  # Complete: Apply RMSE formula

p.sph.val <- extract(r.sph, val)
RMSE.sph.val <- sqrt(mean((p.sph.val-val$P)^2, na.rm=T))  # Complete: Apply RMSE formula

# Exponential model
r.exp <- raster(P.exp)
p.exp.cal <- extract(r.exp, cal)
RMSE.exp.cal <- sqrt(mean((p.exp.cal-cal$P)^2, na.rm=T))  # Complete: Apply RMSE formula

p.exp.val <- extract(r.exp, val)
RMSE.exp.val <- sqrt(mean((p.exp.val-val$P)^2, na.rm=T))


cat("Spherical Model Calibration RMSE:", RMSE.sph.cal,
    "\nSpherical Model Validation RMSE:", RMSE.sph.val,
    "\nExponential Model Calibration RMSE:", RMSE.exp.cal,
    "\nExponential Model Validation RMSE:", RMSE.exp.val)
