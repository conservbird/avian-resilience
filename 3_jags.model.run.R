#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Running JAGS code - Bayesian Multi-scale Occupancy Model
# From "Quantifying avian resilience to habitat change to support conservation
# decision making"
# Code by Amanda L. Hayes-Puttfarcken
# Affiliated with Utah State University
# 10/17/2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Load libraries that will be needed
library(magrittr)
library(dplyr)
library(stringr)
library(jagsUI)
library(rjags)
library(coda)


## Load necessary data
load("jags.model.inputs.final3.RData") # model inputs at all transects and sites, scaled and centered

## List of species to iterate over
species_list <-  c("amro", "atfl", "bbma", "bchu", "bewr", "bggn", "bhco", "bhgr", "brbl", "brsp", "bthu", "btsp", "btyw",
                   "cafi", "chsp", "clnu", "deju", "dufl", "eucd", "grfl", "grvi", "gtto", "heth", "hofi", 
                   "hola", "hosp", "howr", "juti", "lasp", "lazb", "lego", "mgwa", "mobl", "moch", 
                   "modo", "nofl", "nomo", "ocwa", "pisi", "plvi", "rbnu", "rcki", "rowr", "rwbl", "sabs", "saph", "sath", "savs", "spto", "stja", "vesp", "viwa",
                   "wavi", "wcsp", "weki", "weme", "weta", "wewp", "wosj", "yewa", "yrwa") # species codes in the data

for (species in species_list) {
  
  # Load detection data for the current species
  # where $min is the minute category (ie. 1 is minutes 1-2, 2 is minutes 2-3, etc, and 4 is not detected) you detected it at each site-year
  species_data <- read.csv(paste0("jags.pa.matrices\\", species, ".uedep.csv")) # insert full path if necessary
  
  # Update JAGS data list with species-specific occurrence data
  JagsData <- list(
    nobservers = nobservers, 
    nPSUobs = nPSUobs, 
    nSSUobs = nSSUobs, 
    nsurveys = nsurveys,
    nWindows = nWindows,
    nyears = nyears,
    
    # PSU-related
    yearPSU = yearPSU,
    temp = temp_lin_PSU,
    precip = precip_lin_PSU, 
    temp.sq = temp_poly_PSU,
    precip.sq = precip_poly_PSU,
    VDEP.PSU = UEDEP_lin_PSU,
    VDEP.poly.PSU = UEDEP_poly_PSU,
    
    # SSU-related
    yearSSU = yearSSU, 
    VDEP = as.matrix(UEDEP_lin_SSU), 
    VDEP.sq = as.matrix(UEDEP_poly_SSU), 
    dist.road = as.matrix(droad_lin_SSU), 
    dist.road.sq = as.matrix(droad_poly_SSU), 
    c.cover = as.matrix(ccover.scale), 
    shrubland = as.matrix(shrubland_SSU[,3:7]),
    o.t.canopy = as.matrix(o.t.canopy_SSU[,3:7]), 
    herb = as.matrix(herb_SSU[,3:7]),
    PSUid = PSUid_SSU[,2],
    
    # Detection-related
    minutes = c(jags.obs.covs$min_scale), 
    julian = c(jags.obs.covs$jdate_scale), 
    ObsID = c(jags.obs.covs$observer.num), 
    year = c(jags.obs.covs$Year),
    SSUid = c(jags.obs.covs$SSU.num),
    
    # Time-removal detection data for the species
    Occ = c(species_data$min)
  )
  
  # Define initial values for JAGS
  jags_inits <- list(
    "pb_true" = matrix(rep(1, each = nPSUobs), nPSUobs, nyears),
    "ps_true" = matrix(rep(1, each = nSSUobs), nSSUobs, nyears)
  )
  
  # Initialize JAGS model
  jags_model <- rjags::jags.model(
    file = "multiscale_occurrence_yrs_IMBCR_timeRemoval.txt", # insert path if necessary
    data = JagsData,
    inits = jags_inits,
    n.chains = 1,
    n.adapt = 1000,  # Adjust for data-sparse species if needed
    quiet = FALSE
  )
  
  # Run JAGS model
  samples_jags <- rjags::jags.samples(
    jags_model,
    variable.names = c(
      "alpha.0","alpha.year", "alpha.meanprecip", "alpha.meantemp", "alpha.meanprecip.sq", "alpha.meantemp.sq", "alpha.VDEP", "alpha.VDEP.sq",
      "beta.0","beta.year","beta.VDEP","beta.dist.road","beta.c.cover", "beta.VDEP.sq", "beta.dist.road.sq", "beta.shrubland", "beta.o.t.canopy", "beta.herb",
      "gamma.0","gamma.min","gamma.jul",
      "u.PSU","u.observer",
      "pb_true",
      "ps_true",
      "detection_prob_in_one_temporal_window_conditional_on_occurrence","detection_prob_in_one_temporal_window","time_removal_detection_probability",
      "true_PSU_occurrence_prob","true_conditional_SSU_occurrence_prob"
    ),  
    n.iter = 20000,  # Increase iterations if necessary
    thin = 2
  )
  
  # Save the results for the species
  save(samples_jags, JagsData, jags_inits, jags_model, 
       file = paste0("jags_results_VDEP.", species, ".dat"))
}


#### figure creation #1-3 ####


### Figure 3 ###

load("real.occCurve.deriv.yewa.RData") # load the occupancy curve data points for the SSU level of species YEWA (see line 124 for PSU version)
par(mfrow = c(1,3), mar = numeric(4), oma = c(3, 3, .5, .5), # for 3 panels
    mgp = c(2, .6, 0), mai = c(.1, .1, .1, .1))

# Create spline fits up to x = 40
x_seq <- seq(min(xax), 40, by = 0.1)

spline_pred <- spline(xax, predOccCurve, xout = x_seq)
spline_lower <- spline(xax, predOccCurve.lower, xout = x_seq)
spline_upper <- spline(xax, predOccCurve.upper, xout = x_seq)

# Plot setup
plot(x_seq, spline_pred$y, type='n', las=1, ylim=range(0, 0.4), xlim=c(0,40),
     xlab="UEDEP", ylab="Occupancy probability", bty="n")

# Polygon for the credible interval
polygon(x = c(x_seq, rev(x_seq)),
        y = c(spline_lower$y, rev(spline_upper$y)),
        col = adjustcolor('skyblue', 0.5), border = NA)

# Main occupancy curve line
lines(spline_pred$x, spline_pred$y, lwd = 2, col = "blue")

# Add GMIP and RT based on Table S5
y0 <- approx(spline_pred$x, spline_pred$y, xout = 0)$y
y40 <- approx(spline_pred$x, spline_pred$y, xout = 40)$y
segments(0, 0, 0, y0, col = "#48D1CC", lty = 2, lwd = 2) # RT
segments(40, 0, 40, y40, col = "darkorange", lty = 2, lwd = 2) # GMIP

mtext("UEDEP", side = 1, outer = TRUE, line = 1.5)
mtext("Occupancy probability", side = 2, outer = TRUE, line = 1.5)

load("real.occCurve.deriv.atfl.RData") # panel 2 - species ATFL

# Create spline fits up to x = 40
x_seq <- seq(min(xax), 40, by = 0.1)

spline_pred <- spline(xax, predOccCurve, xout = x_seq)
spline_lower <- spline(xax, predOccCurve.lower, xout = x_seq)
spline_upper <- spline(xax, predOccCurve.upper, xout = x_seq)

# Plot setup
plot(x_seq, spline_pred$y, type='n', las=1, ylim=range(0, 0.4), xlim=c(0,40),
     xlab="UEDEP", ylab="Occupancy probability", bty="n")

# Polygon for the credible interval
polygon(x = c(x_seq, rev(x_seq)),
        y = c(spline_lower$y, rev(spline_upper$y)),
        col = adjustcolor('skyblue', 0.5), border = NA)

# Main occupancy curve line
lines(spline_pred$x, spline_pred$y, lwd = 2, col = "blue")

# Add GMIP and RT based on Table S5
y40a <- approx(spline_pred$x, spline_pred$y, xout = 40)$y
y40b <- approx(spline_pred$x, spline_pred$y, xout = 40)$y
segments(40, 0, 40, y40a, col = "#48D1CC", lty = 3, lwd = 2) # RT
segments(40, 0, 40, y40b, col = "darkorange", lty = 2, lwd = 2) # GMIP

load("real.occCurve.deriv.howr.RData") # panel 3 - species HOWR

# Create spline fits up to x = 40
x_seq <- seq(min(xax), 40, by = 0.1)

spline_pred <- spline(xax, predOccCurve, xout = x_seq)
spline_lower <- spline(xax, predOccCurve.lower, xout = x_seq)
spline_upper <- spline(xax, predOccCurve.upper, xout = x_seq)

# Plot setup
plot(x_seq, spline_pred$y, type='n', las=1, ylim=range(0, 0.4), xlim=c(0,40),
     xlab="UEDEP", ylab="Occupancy probability", bty="n")

# Polygon for the credible interval
polygon(x = c(x_seq, rev(x_seq)),
        y = c(spline_lower$y, rev(spline_upper$y)),
        col = adjustcolor('skyblue', 0.5), border = NA)

# Main occupancy curve line
lines(spline_pred$x, spline_pred$y, lwd = 2, col = "blue")

# Add GMIP and RT based on Table S5
y0a <- approx(spline_pred$x, spline_pred$y, xout = 0)$y
y0b <- approx(spline_pred$x, spline_pred$y, xout = 0)$y
segments(0, 0, 0, y0a, col = "#48D1CC", lty = 3, lwd = 2) # RT
segments(0, 0, 0, y0b, col = "darkorange", lty = 2, lwd = 2) # GMIP

legend("topright",                                 # location of the key
       legend = c("Occupancy", "RT", "GMIP"),  # labels
       col = c("blue", "#48D1CC", "darkorange"),     # colors for line segments
       lwd = c(2, 2, 2),                           # line widths
       lty = c(1, 2, 2),                           # line types: solid for curve, dashed for verticals
       bty = "n",                                  # remove box border
       cex = 0.9,                                  # text size
       y.intersp = 1.2)                            # spacing between lines


### Figure 2 ###

## creating simulated plots for resilient and less resilient species showing their GMIPs and RTs 
## plus a simulated plot for a positive species to show the meaning of the difference in RT calculation

### Less resilient species example

# Load necessary library
library(ggplot2)

# Set seed for reproducibility
set.seed(124)

# Generate predictor variable (x) from 0 to 100
x <- seq(0, 100, length.out = 100)

# Define parameters for the logistic function
a <- 0.1   # Steepness of the decline
b <- 50    # Inflection point (where occupancy ~ 0.5)

# Generate occupancy probability (y) using a logistic decline function
y_clean <- 1 / (1 + exp(a * (x - b)))  # Ensures values stay between 0 and 1

# Add slight noise but keep values constrained within [0,1]
noise <- runif(100, min=-0.02, max=0.02)  # Small noise between -0.02 and 0.02
y <- pmax(0, pmin(1, y_clean + noise))  # Keep within bounds
# Fit a logistic regression model
model <- glm(y ~ x, family = binomial(link = "logit"))
summary(model)

# Create data frame for plotting
plot_data <- data.frame(x, y)

# Create smooth prediction curve
pred_data <- data.frame(x = seq(0, 100, length.out = 100))
pred_data$y_pred <- predict(model, newdata = pred_data, type = "response")

# add derivatives
d1 <- diff(y)
d2 <- diff(d1)
if (y[1] > y[100]) {
  GMIP1 <- which.max(abs(d1)) + 1 
  resilience.theshold1 <- which.max(round(abs(d2), digits = 6)) + 1 # for negative species
} else {
  GMIP1 <- which.max(d1) + 1
  resilience.theshold1 <- which.max(d2) + 1 # for positive species
} 
GMIP1
resilience.theshold1

y_GMIP1 <- predict(model, newdata = data.frame(x = GMIP1), type = "response")
y_res.thresh1 <- predict(model, newdata = data.frame(x = resilience.theshold1), type = "response")

# Plot the logistic decline
spp1 <- ggplot(plot_data, aes(x, y)) +
  geom_line(data = pred_data, aes(x, y_pred), color = "black", size = 1.2) +  # Fitted curve
  geom_segment(aes(x = GMIP1, xend = GMIP1, y = 0, yend = y_GMIP1), 
               color = "blue", linetype = "dashed", size = 1) +  # Vertical line to regression line
  geom_segment(aes(x = resilience.theshold1, xend = resilience.theshold1, y = 0, yend = y_res.thresh1), 
               color = "darkorange", linetype = "dashed", size = 1) +  # Vertical line to regression line
  labs(x = "",
       y = "Occupancy Probability") +
  theme_classic()


### More resilient species example

# Set seed for reproducibility
set.seed(123)

# Generate predictor variable (x) from 0 to 100
x <- seq(0, 100, length.out = 100)

# Define parameters for the slow-then-fast decline
lambda <- 0.1  # Controls how fast the curve bends
b <- 75         # Inflection point (where decline speeds up)

# Generate occupancy probability (y) using an inverted exponential function
y_clean <- 1 / (1 + exp(lambda * (x - b)))  # Logistic function ensuring values in [0,1]

# Add slight noise but keep values within [0,1]
noise <- runif(100, min=-0.02, max=0.02)  # Small random variation
y <- pmax(0, pmin(1, y_clean + noise))  # Ensure values stay within bounds

# Fit a logistic regression model
model <- glm(y ~ x, family = binomial(link = "logit"))
summary(model)

# Create data frame for plotting
plot_data <- data.frame(x, y)

# Create smooth prediction curve
pred_data <- data.frame(x = seq(0, 100, length.out = 100))
pred_data$y_pred <- predict(model, newdata = pred_data, type = "response")

# add derivatives
d1 <- diff(y)
d2 <- diff(d1)
if (y[1] > y[100]) {
  GMIP <- which.max(abs(d1)) + 1 
  resilience.theshold <- which.max(round(abs(d2), digits = 6)) + 1 # for negative species
} else {
  GMIP <- which.max(d1) + 1
  resilience.theshold <- which.max(d2) + 1 # for positive species
} 
GMIP
resilience.theshold

y_GMIP <- predict(model, newdata = data.frame(x = GMIP), type = "response")
y_res.thresh <- predict(model, newdata = data.frame(x = resilience.theshold), type = "response")

# Plot the slow-then-fast decline
spp2 <- ggplot(plot_data, aes(x, y)) +
  geom_line(data = pred_data, aes(x, y_pred), color = "black", size = 1.2) +  # Fitted curve
  geom_segment(aes(x = GMIP, xend = GMIP, y = 0, yend = y_GMIP), 
               color = "blue", linetype = "dashed", size = 1) +  # Vertical line to regression line
  geom_segment(aes(x = resilience.theshold, xend = resilience.theshold, y = 0, yend = y_res.thresh), 
               color = "darkorange", linetype = "dashed", size = 1) +  # Vertical line to regression line
  ylim(0,1) +
  labs(x = "",
       y = "") +
  theme_classic()


### Positive species example 

# Set seed for reproducibility
set.seed(123)

# Predictor range
x <- seq(0, 100, length.out = 100)

# Quintic polynomial inside the logistic (to easily create the conditions of this relationship)
beta0 <- -6
beta1 <- 0.20
beta2 <- -0.004
beta3 <- 0.00006
beta4 <- -0.0000004
beta5 <- 0.000000001

eta <- beta0 + beta1*x + beta2*x^2 + beta3*x^3 + beta4*x^4 + beta5*x^5
y <- 1 / (1 + exp(-eta))   # logistic keeps y in [0,1]

# Numerical derivatives 
grad_central <- function(x, y) {
  n <- length(x)
  g <- numeric(n)
  g[1] <- (y[2] - y[1]) / (x[2] - x[1])
  g[n] <- (y[n] - y[n-1]) / (x[n] - x[n-1])
  for (i in 2:(n-1)) {
    g[i] <- (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
  }
  g
}

y_prime <- grad_central(x, y) # first derivative column
y_double <- grad_central(x, y_prime) # second derivative column

# combine in a data frame
df <- data.frame(x, y, y_prime, y_double) 

# calculate GMIP and RT
if (df$y[1] > df$y[100]) {
  GMIP2 <- which.max(abs(df$y_prime)) + 1 
  resilience.theshold2 <- which.max(round(abs(df$y_double), digits = 6)) + 1 # for negative species
} else {
  GMIP2 <- which.max(df$y_prime) + 1
  resilience.theshold2 <- which.max(df$y_double) + 1 # for positive species
} 
GMIP2
resilience.theshold2

y_GMIP2 <- df$y[GMIP2]
y_res.thresh2 <- df$y[resilience.theshold2]

# plot the positive trend example
spp3 <- ggplot(df, aes(x, y)) +
  geom_line(data = df, aes(x, y), color = "black", size = 1.2) +  # Fitted curve
  geom_segment(aes(x = GMIP2, xend = GMIP2, y = 0, yend = y_GMIP2), 
               color = "blue", linetype = "dashed", size = 1) +  # Vertical line to regression line
  geom_segment(aes(x = resilience.theshold2, xend = resilience.theshold2, y = 0, yend = y_res.thresh2), 
               color = "darkorange", linetype = "dashed", size = 1) +  # Vertical line to regression line
  labs(x = "",
       y = "") +
  theme_classic()


# plot them together
library(patchwork)
spp1 | spp2 | spp3


### Figure 1 ###

# Final UEDEP map

# Required libraries
library(terra)
library(sf)
library(raster)

# Load data
UEDEP <- rast("UT_UEDEP.tif") # the LANDFIRE UEDEP data
sampling.points <- vect("bird.spatial.pub.use.shp") # survey sites, in EPSG: 5070
crs(sampling.points) <- "epsg: 5070" # if it didn't load in
utah_utgis <- vect("utah.pub.use.shp") # state of Utah boundary shapefile, in EPSG: 5070
crs(utah_utgis) <- "epsg:5070" # if it didn't load in

# Reproject UEDEP and Utah boundary to a visually appealing EPSG (4269)
UEDEP.straight <- project(UEDEP, "epsg:4269")
utah.good <- project(utah_utgis, "epsg:4269")
UEDEP.straight <- mask(UEDEP.straight, utah.good)
sampling.points.good <- project(sampling.points, "epsg:4269")

# Define range for legend axis
r.range <- c(min(values(UEDEP.straight), na.rm = T), max(values(UEDEP.straight), na.rm = T))

# Plot UEDEP Raster with Legend
plot(UEDEP.straight, 
     #alpha = 0.7,
     legend.args = list(text = 'UEDEP Score', side = 4, font = 2, line = 2.5, cex = 0.8),
     axis.args = list(at = seq(r.range[1], r.range[2], length.out = 5),
                      labels = seq(r.range[1], r.range[2], length.out = 5), 
                      cex.axis = 0.6))

# Add Utah boundary and sampling points
plot(utah.good, add = TRUE)
plot(sampling.points.good, add = TRUE, 
     col = "black", 
     pch = 1,
     cex = 0.9)

# Add North Arrow and Scale Bar
north(type = 2, location = "topright") # Place north arrow in top right corner
sbar(100, xy = "bottomright", type = "bar", divs = 2, lonlat = TRUE, below = "Kilometers", labels = c(0, 50, 100))

# Manually added legend annotations in final version
