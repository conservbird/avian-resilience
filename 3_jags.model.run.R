#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Running JAGS code - Bayesian Multi-scale Occupancy Model
# From "Quantifying avian resilience to habitat change to support conservation
# decision making"
# Code by Amanda L. Hayes-Puttfarcken
# Affiliated with Utah State University
# 04/11/2025
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


## resilience metrics (GRIP and RT) are calculated as part of Figure 4 creation below (see lines 120-188 for an example)

#### figure creation ####

### Figure 4 ###

# make 2 dataframes, one for the SSU level and one for the PSU level, with GRIP and RT

## AMRO; one species example for the PSU level dataframe
load("jags_results_VDEP.amro.uedep.dat") # the model output
yrwa <- read.csv("amro.uedep.csv") # the detection input data
for.estimates2 <- yrwa 

colnames(for.estimates2)[2] <- "full.names" # check is correct column every time - should be pointnum_time_year
for.estimates2 <- for.estimates2 %>% 
  left_join(jags.site.covs %>% dplyr::select(lin_UEDEP, poly_UEDEP, UEDEP, lin_UEDEP_scale, poly_UEDEP_scale, full.names, PSU.num, pointnum), by = "full.names")
for.estimates2$pointnum <- str_extract(for.estimates2$pointnum, "^.*(?=(_))") 

psu.estimates <- for.estimates2 %>% 
  left_join(psu.extraction2 %>% dplyr::select(mean_UEDEP, pointnum), by = "pointnum")

psu.estimates <- psu.estimates %>% group_by(PSU.num) %>%
  dplyr::dplyr::summarise(detect = paste(min, collapse=", "),
                          UEDEP = mean(mean_UEDEP)) # grouping the PSUs together
psu.estimates$first <- rep(NA, nrow(psu.estimates)) # blank row

psu.ests <- strsplit(psu.estimates$detect, ", ") # splitting the conglomerated rows
letshope <- as.data.frame(do.call(rbind, psu.ests)) # re-binding the rows into a better form
psu.estimates$first <- apply(letshope, 1, FUN = min) # taking the minimum of the rows to select the first time period the species was spotted in in each PSU
psu.estimates <- psu.estimates[,-2] # removing the middle row with the list of detection info

psu_for_poly <- poly(psu.estimates$UEDEP, 2) # orthagonal polynomials
psu_for_df <- outer(psu.estimates$UEDEP, 0:2, `^`); colnames(psu_for_df) <- c("Intercept", "Linr", "Squ") #contains 'raw' polynomial values
psu.estimates$UEDEP <- psu_for_df[,2]
psu.estimates$UEDEP.sq <- psu_for_df[,3]
poly_lin <- lm(psu_for_poly[,1] ~ UEDEP, data = psu.estimates) #create a model that predicts the poly degree=1 value based on raw data
poly_qu <- lm(psu_for_poly[,2] ~ UEDEP + UEDEP.sq, data = psu.estimates) 


med.values <- c(median(samples_jags$alpha.0[2000:10000]), median(samples_jags$alpha.VDEP[2000:10000]), median(samples_jags$alpha.VDEP.sq[2000:10000])) # range here depends on the # of iterations and burn in the species model
hdi.values <- c(hdi(samples_jags$alpha.0[2000:10000]), hdi(samples_jags$alpha.VDEP[2000:10000]), hdi(samples_jags$alpha.VDEP.sq[2000:10000]))

xay <- seq(min(psu.estimates$UEDEP), max(psu.estimates$UEDEP), by = 1) 


predOccCurve.psu = plogis(med.values[1] + 
                            (med.values[2] * (((poly_lin$coefficients[["(Intercept)"]] + poly_lin$coefficients[["UEDEP"]] * xay)  - mean(psu_for_poly[,1])) / sd(psu_for_poly[,1])))  +
                            (med.values[3] * ((((poly_qu$coefficients[["(Intercept)"]] + poly_qu$coefficients[["UEDEP"]] * xay) + (poly_qu$coefficients[["UEDEP.sq"]] * xay^2)) - mean(psu_for_poly[,2])) / sd(psu_for_poly[,2])))) 
predOccCurve.psu.lower = plogis(hdi.values[1] + 
                                  (hdi.values[3] * (((poly_lin$coefficients[["(Intercept)"]] + poly_lin$coefficients[["UEDEP"]] * xay)  - mean(psu_for_poly[,1])) / sd(psu_for_poly[,1])))  +
                                  (hdi.values[5] * ((((poly_qu$coefficients[["(Intercept)"]] + poly_qu$coefficients[["UEDEP"]] * xay) + (poly_qu$coefficients[["UEDEP.sq"]] * xay^2)) - mean(psu_for_poly[,2])) / sd(psu_for_poly[,2])))) 
predOccCurve.psu.upper = plogis(hdi.values[2] + 
                                  (hdi.values[4] * (((poly_lin$coefficients[["(Intercept)"]] + poly_lin$coefficients[["UEDEP"]] * xay)  - mean(psu_for_poly[,1])) / sd(psu_for_poly[,1])))  +
                                  (hdi.values[6] * ((((poly_qu$coefficients[["(Intercept)"]] + poly_qu$coefficients[["UEDEP"]] * xay) + (poly_qu$coefficients[["UEDEP.sq"]] * xay^2)) - mean(psu_for_poly[,2])) / sd(psu_for_poly[,2]))))

d1 <- diff(predOccCurve.psu) # first derivative
d2 <- diff(d1) # second derivative

perc0.050 <- c(median(predOccCurve.psu) - (0.05 * median(predOccCurve.psu)), median(predOccCurve.psu) + (0.05 * median(predOccCurve.psu))) # decision rule of 5% change
if (between(min(predOccCurve.psu), perc0.050[1], perc0.050[2]) | between(max(predOccCurve.psu), perc0.050[1], perc0.050[2])) {
  GRIP <- NA
  resilience.theshold <- 37 # for linear species
  direct.1 <- "flat"
} else {
  if (predOccCurve[1] > predOccCurve[38]) {
    GRIP <- which.max(abs(d1)) + 1 
    resilience.theshold <- which.max(abs(d2)) + 1
    direct.1 <- "negative" # for negative species
  } else {
    GRIP <- which.max(d1) + 1
    resilience.theshold <- which.max(d2) + 1
    direct.1 <- "positive" # for positive species
  } 
}

derivatives.per.spp.psu[1, 2:4] <- c(GRIP, resilience.theshold, direct.1) # putting it into dataframe
save(xay, predOccCurve.psu, predOccCurve.psu.lower, predOccCurve.psu.upper, d1, d2, file = "real.occ.objs.psu.amro.uedep.RData") # will use this and the SSU version for occupancy probability plots 


## completed dataframes
derivatives.per.spp <- read.csv("ssu.derivs.w.rules.GRIPrange.csv") # derivatives at SSU level; same as "GRIP.RT.SSU.csv" but with 1 column for all GRIPs
derivatives.per.spp.psu <- read.csv("psu.derivs.w.rules.GRIPrange.csv") # derivatives at PSU level; same as "GRIP.RT.PSU.csv" but with 1 column for all GRIPs

# for spp with range of GRIP values, choose median
derivatives.per.spp.4grip <- derivatives.per.spp %>% drop_na("GRIP.range")
derivatives.per.spp.4grip$GRIP.plot <- c(38,36.5,38,3,38,2,33,38,29,9,2,2,16.5,38,38,38,2,38,2,2,38,37,37,4.5,38,38,38,38,38,2) # have to do manual average bc column is character

derivatives.per.spp.psu.4grip <- derivatives.per.spp.psu %>% drop_na("GRIP.psu")
derivatives.per.spp.psu.4grip$GRIP.plot <- c(38,2.5,5.5,25,3,20,3,2,2,2.5,2,5,21,38,3,5) # same as line 200 but for PSU


par(mfrow = c(2,2), mar = numeric(4), oma = c(4, 4, .5, .5), 
    mgp = c(2, .6, 0), mai = c(.1, .3, .3, .3), family = "arial")

hist(derivatives.per.spp.4grip$GRIP.plot, col = "#cc5500", main = "", # GRIP SSU values
     ylab = "", xlab = "", las = 1, ylim = c(0,20))
hist(derivatives.per.spp.psu.4grip$GRIP.plot, col = "#cc5500", main = "", # GRIP PSU values
     ylab = "", xlab = "", las = 1, ylim = c(0,20), breaks = 8)
mtext("GRIP", side = 1, outer = FALSE, line = 2) # ended up manually adding this
#mtext("# of Species", side = 2, outer = TRUE, line = 1)

hist(derivatives.per.spp$RT.newest, col = "#cc5500", main = "", # RT SSU values
     ylab = "", xlab = "", las = 1, ylim = c(0,50))
hist(derivatives.per.spp.psu$RT.psu, col = "#cc5500", main = "", # RT PSU values
     ylab = "", xlab = "", las = 1, ylim = c(0,50))
mtext("Resilience Threshold", side = 1, outer = FALSE, line = 2) # ended up manually adding this
mtext("# of Species", side = 2, outer = TRUE, line = 1)

# correlation between RT and GRIP @ PSU and SSU levels
d4cor.ssu <- derivatives.per.spp[, c(2,3)] # just GRIP and RT
d4cor.psu <- derivatives.per.spp.psu[, c(3,4)] # just GRIP and RT
cor(d4cor.psu,d4cor.ssu, method = "pearson", use = "pairwise.complete.obs")

# correlation between RT and GRIP at each level
cor(d4cor.ssu, method = "pearson", use = "pairwise.complete.obs")
cor(d4cor.psu, method = "pearson", use = "pairwise.complete.obs")


### Figure 3 ###

load("real.occCurve.deriv.yewa.RData") # load the occupancy curve data points for the SSU level of species YEWA (see line 124 for PSU version)
par(mfrow = c(1,3), mar = numeric(4), oma = c(3, 3, .5, .5), # for 3 panels
    mgp = c(2, .6, 0), mai = c(.1, .1, .1, .1))
plot(xax, predOccCurve, type='n', las=1, ylim=range(0, 0.4), xlim=range(0,40), bty="n",
     xlab="UEDEP", ylab="Occupancy probability")
polygon(x=c(xax, rev(xax)),
        y=c(predOccCurve.lower, rev(predOccCurve.upper)),
        col=adjustcolor('skyblue', 0.5), border=NA) # for the credible interval of occupancy probability
lines(predOccCurve~xax, lwd = 2, col = "blue") # occupancy probability
mtext("UEDEP", side = 1, outer = TRUE, line = 1.5)
mtext("Occupancy probability", side = 2, outer = TRUE, line = 1.5)

load("real.occCurve.deriv.atfl.RData") # panel 2 - species ATFL
plot(xax, predOccCurve, type='n', las=1, ylim=range(0, 0.4), xlim=range(0,40), bty="n",
     xlab="UEDEP", ylab="Occupancy probability")
polygon(x=c(xax, rev(xax)),
        y=c(predOccCurve.lower, rev(predOccCurve.upper)),
        col=adjustcolor('skyblue', 0.5), border=NA) # for the credible interval of occupancy probability
lines(predOccCurve~xax, lwd = 2, col = "blue") # occupancy probability

load("real.occCurve.deriv.howr.RData") # panel 3 - species HOWR
plot(xax, predOccCurve, type='n', las=1, ylim=range(0, 0.4), xlim=range(0,40), bty="n",
     xlab="UEDEP", ylab="Occupancy probability")
polygon(x=c(xax, rev(xax)),
        y=c(predOccCurve.lower, rev(predOccCurve.upper)),
        col=adjustcolor('skyblue', 0.5), border=NA) # for the credible interval of occupancy probability
lines(predOccCurve~xax, lwd = 2, col = "blue") # occupancy probability


### Figure 2 ###

## creating simulated plots for resilient and less resilient species showing their GRIPs and RTs 

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
  GRIP1 <- which.max(abs(d1)) + 1 
  resilience.theshold1 <- which.max(round(abs(d2), digits = 6)) + 1 # for negative species
} else {
  GRIP1 <- which.max(d1) + 1
  resilience.theshold1 <- which.max(d2) + 1 # for positive species
} 
GRIP1
resilience.theshold1

y_GRIP1 <- predict(model, newdata = data.frame(x = GRIP1), type = "response")
y_res.thresh1 <- predict(model, newdata = data.frame(x = resilience.theshold1), type = "response")

# Plot the logistic decline
spp1 <- ggplot(plot_data, aes(x, y)) +
  geom_line(data = pred_data, aes(x, y_pred), color = "black", size = 1.2) +  # Fitted curve
  geom_segment(aes(x = GRIP1, xend = GRIP1, y = 0, yend = y_GRIP1), 
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
  GRIP <- which.max(abs(d1)) + 1 
  resilience.theshold <- which.max(round(abs(d2), digits = 6)) + 1 # for negative species
} else {
  GRIP <- which.max(d1) + 1
  resilience.theshold <- which.max(d2) + 1 # for positive species
} 
GRIP
resilience.theshold

y_GRIP <- predict(model, newdata = data.frame(x = GRIP), type = "response")
y_res.thresh <- predict(model, newdata = data.frame(x = resilience.theshold), type = "response")

# Plot the slow-then-fast decline
spp2 <- ggplot(plot_data, aes(x, y)) +
  geom_line(data = pred_data, aes(x, y_pred), color = "black", size = 1.2) +  # Fitted curve
  geom_segment(aes(x = GRIP, xend = GRIP, y = 0, yend = y_GRIP), 
               color = "blue", linetype = "dashed", size = 1) +  # Vertical line to regression line
  geom_segment(aes(x = resilience.theshold, xend = resilience.theshold, y = 0, yend = y_res.thresh), 
               color = "darkorange", linetype = "dashed", size = 1) +  # Vertical line to regression line
  ylim(0,1) +
  labs(x = "",
       y = "") +
  theme_classic()


# plot both together
library(patchwork)
spp1 | spp2 


### Figure 1 ###

# Final UEDEP map

# Required libraries
library(terra)
library(sf)
library(raster)

# Load data
UEDEP <- rast("UT_UEDEP.tif") # the LANDFIRE UEDEP data
sampling.points <- vect("bird.spatial.pub.use.shp") # survey sites
utah_utgis <- vect("utah.pub.use.shp") # state of Utah boundary shapefile

# Reproject UEDEP and Utah boundary to a visually appealing EPSG (4269)
UEDEP.straight <- project(UEDEP, "epsg:4269")
utah.good <- project(utah_utgis, "epsg:4269")
UEDEP.straight <- mask(UEDEP.straight, utah.good)
sampling.points.good <- project(sampling.points, "epsg:4269")

# Define range for legend axis
r.range <- c(min(values(UEDEP.straight), na.rm = T), max(values(UEDEP.straight), na.rm = T))

# Plot UEDEP Raster with Legend
plot(UEDEP.straight, 
     legend.args = list(text = 'UEDEP Score', side = 4, font = 2, line = 2.5, cex = 0.8),
     axis.args = list(at = seq(r.range[1], r.range[2], length.out = 5),
                      labels = seq(r.range[1], r.range[2], length.out = 5), 
                      cex.axis = 0.6))

# Add Utah boundary and sampling points
plot(utah.good, add = TRUE)
plot(sampling.points.good, add = TRUE, col = "black", pch = 16)

# Add North Arrow and Scale Bar
north(type = 2, location = "topright") # Place north arrow in top right corner
sbar(100, xy = "bottomright", type = "bar", divs = 2, lonlat = TRUE, below = "Kilometers", labels = c(0, 50, 100))

# Manually added legend annotations in final version
