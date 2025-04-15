#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Formatting covatiates for JAGS 
# From "Quantifying avian resilience to habitat change to support conservation
# decision making"
# Code by Amanda L. Hayes-Puttfarcken
# Affiliated with Utah State University
# 04/11/2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Loading libraries
library(dplyr)
library(magrittr)
library(stringr)

# Loading base covariate data
jags.obs.covs <- read.csv("jags.obs.covs.final.forpub.csv") # loading detection level covariates
jags.site.covs <- read.csv("jags.site.covs.final.forpub.csv") # loading site level covariates

# with multiple years of data for PSU and SSU state its easier to write the jags code for covariates with wide data
# so, for each covariate, instead of a long dataset where year is a column
# we'll make a separate matrix for each covariate where each column contains values of the covariate
# for each year. We don't specify that detection probability changes between years, so we can keep that
# data in long format

## Distance to Roads
poly_droad <- as.matrix(poly(jags.site.covs$DisttoRoad, degree = 2)) # this line would take your raw droad data (your_data$droad) apply the poly function where degree = 2 means you want the square, and store the poly squared new data as a matrix in poly_droad
poly_droad <- as.data.frame(poly_droad) # linear and quadratic terms
jags.site.covs$ID_poly <- seq(1:nrow(jags.site.covs)) # only need to do this for the first poly cov
poly_droad$ID_poly <- seq(1:nrow(poly_droad))
jags.site.covs <- jags.site.covs %>%
  left_join(poly_droad, by = "ID_poly")
names(jags.site.covs)[names(jags.site.covs) == '1'] <- 'lin_droad'
names(jags.site.covs)[names(jags.site.covs) == '2'] <- 'poly_droad'
droad_lin_SSU <- jags.site.covs %>%
  dplyr::select(pointnum, Year, lin_droad, PSU.num) %>%
  tidyr::spread(Year, lin_droad)
droad_lin_SSU <- scale(droad_lin_SSU[,-c(1,2)]) # scaled
droad_poly_SSU <- jags.site.covs %>%
  dplyr::select(pointnum, Year, poly_droad, PSU.num) %>%
  tidyr::spread(Year, poly_droad)
droad_poly_SSU <- scale(droad_poly_SSU[,-c(1,2)]) # scaled

## UEDEP
poly_UEDEP <- as.matrix(poly(jags.site.covs$UEDEP, degree = 2)) # this line would take your raw VDEP data (your_data$UEDEP) apply the poly function where degree = 2 means you want the square, and store the poly squared new data as a matrix in poly_UEDEP
poly_UEDEP <- as.data.frame(poly_UEDEP)
poly_UEDEP$ID_poly <- seq(1:nrow(poly_UEDEP))
jags.site.covs <- jags.site.covs %>%
  left_join(poly_UEDEP, by = "ID_poly")
names(jags.site.covs)[names(jags.site.covs) == '1'] <- 'lin_UEDEP'
names(jags.site.covs)[names(jags.site.covs) == '2'] <- 'poly_UEDEP'
UEDEP_lin_SSU <- jags.site.covs %>%
  dplyr::select(pointnum, Year, lin_UEDEP, PSU.num) %>%
  tidyr::spread(Year, lin_UEDEP)
UEDEP_lin_SSU <- scale(UEDEP_lin_SSU[,-c(1,2)]) # scaled
UEDEP_poly_SSU <- jags.site.covs %>%
  dplyr::select(pointnum, Year, poly_UEDEP, PSU.num) %>%
  tidyr::spread(Year, poly_UEDEP)
UEDEP_poly_SSU <- scale(UEDEP_poly_SSU[,-c(1,2)]) # scaled

## Canopy Cover
ccover_SSU = jags.site.covs %>% 
  dplyr::select(pointnum, Year, lf_cc, PSU.num) %>%
  tidyr::spread(Year, lf_cc)
ccover.scale <- scale(ccover_SSU[,-c(1,2)])

## Index for which PSU the SSU is in
PSUid_SSU = jags.site.covs %>% 
  dplyr::select(pointnum, Year, PSU.num) %>%
  tidyr::spread(Year, PSU.num)

#year_SSU = jags.site.covs %>% 
#  dplyr::select(pointnum, Year, PSU.num) %>%
#  tidyr::spread(Year, Year)

## Percentage of shrubland habitat cover
shrubland_SSU = jags.site.covs %>%
  dplyr::select(pointnum, Year, PSU.num, shrubland) %>%
  tidyr::spread(Year, shrubland) 

## Percentage of open tree canopy habitat cover
o.t.canopy_SSU = jags.site.covs %>% 
  dplyr::select(pointnum, Year, PSU.num, open.tree.canopy) %>%
  tidyr::spread(Year, open.tree.canopy)

## Percentage of herbaceous habitat cover
herb_SSU = jags.site.covs %>% 
  dplyr::select(pointnum, Year, PSU.num, herbaceous) %>%
  tidyr::spread(Year, herbaceous)


#### PSU covariates: temp, precip, and UEDEP ----
psu.extraction2 <- read.csv("psu.extraction.wUEDEP.csv")
precip_PSU <- psu.extraction2[, c(3:7)]
temp_PSU <- psu.extraction2[, c(8:12)]

# Temperature
poly_temp2017 <- as.matrix(poly(temp_PSU$temp2017, degree = 2))
poly_temp2017 <- as.data.frame(poly_temp2017)
temp_PSU$ID_poly <- seq(1:nrow(temp_PSU))
poly_temp2017$ID_poly <- seq(1:nrow(poly_temp2017))
temp_PSU <- temp_PSU %>%
  left_join(poly_temp2017, by = "ID_poly")
names(temp_PSU)[names(temp_PSU) == '1'] <- 'lin_temp2017'
names(temp_PSU)[names(temp_PSU) == '2'] <- 'poly_temp2017'
poly_temp2018 <- as.matrix(poly(temp_PSU$temp2018, degree = 2))
poly_temp2018 <- as.data.frame(poly_temp2018)
temp_PSU$ID_poly <- seq(1:nrow(temp_PSU))
poly_temp2018$ID_poly <- seq(1:nrow(poly_temp2018))
temp_PSU <- temp_PSU %>%
  left_join(poly_temp2018, by = "ID_poly")
names(temp_PSU)[names(temp_PSU) == '1'] <- 'lin_temp2018'
names(temp_PSU)[names(temp_PSU) == '2'] <- 'poly_temp2018'
poly_temp2019 <- as.matrix(poly(temp_PSU$temp2019, degree = 2))
poly_temp2019 <- as.data.frame(poly_temp2019)
temp_PSU$ID_poly <- seq(1:nrow(temp_PSU))
poly_temp2019$ID_poly <- seq(1:nrow(poly_temp2019))
temp_PSU <- temp_PSU %>%
  left_join(poly_temp2019, by = "ID_poly")
names(temp_PSU)[names(temp_PSU) == '1'] <- 'lin_temp2019'
names(temp_PSU)[names(temp_PSU) == '2'] <- 'poly_temp2019'
poly_temp2020 <- as.matrix(poly(temp_PSU$temp2020, degree = 2))
poly_temp2020 <- as.data.frame(poly_temp2020)
temp_PSU$ID_poly <- seq(1:nrow(temp_PSU))
poly_temp2020$ID_poly <- seq(1:nrow(poly_temp2020))
temp_PSU <- temp_PSU %>%
  left_join(poly_temp2020, by = "ID_poly")
names(temp_PSU)[names(temp_PSU) == '1'] <- 'lin_temp2020'
names(temp_PSU)[names(temp_PSU) == '2'] <- 'poly_temp2020'
poly_temp2021 <- as.matrix(poly(temp_PSU$temp2021, degree = 2))
poly_temp2021 <- as.data.frame(poly_temp2021)
temp_PSU$ID_poly <- seq(1:nrow(temp_PSU))
poly_temp2021$ID_poly <- seq(1:nrow(poly_temp2021))
temp_PSU <- temp_PSU %>%
  left_join(poly_temp2021, by = "ID_poly")
names(temp_PSU)[names(temp_PSU) == '1'] <- 'lin_temp2021'
names(temp_PSU)[names(temp_PSU) == '2'] <- 'poly_temp2021'

temp_lin_PSU <- matrix(NA, nrow = length(temp_PSU$temp2017), ncol = 5)
temp_lin_PSU[, 1] <- temp_PSU[, 'lin_temp2017']
temp_lin_PSU[, 2] <- temp_PSU[, 'lin_temp2018']
temp_lin_PSU[, 3] <- temp_PSU[, 'lin_temp2019']
temp_lin_PSU[, 4] <- temp_PSU[, 'lin_temp2020']
temp_lin_PSU[, 5] <- temp_PSU[, 'lin_temp2021']
temp_lin_PSU <- scale(temp_lin_PSU) # scale

temp_poly_PSU <- matrix(NA, nrow = length(temp_PSU$temp2017), ncol = 5)
temp_poly_PSU[, 1] <- temp_PSU[, 'poly_temp2017']
temp_poly_PSU[, 2] <- temp_PSU[, 'poly_temp2018']
temp_poly_PSU[, 3] <- temp_PSU[, 'poly_temp2019']
temp_poly_PSU[, 4] <- temp_PSU[, 'poly_temp2020']
temp_poly_PSU[, 5] <- temp_PSU[, 'poly_temp2021']
temp_poly_PSU <- scale(temp_poly_PSU) # scale


# Precipitation
poly_precip2017 <- as.matrix(poly(precip_PSU$precip2017, degree = 2))
poly_precip2017 <- as.data.frame(poly_precip2017)
precip_PSU$ID_poly <- seq(1:nrow(precip_PSU))
poly_precip2017$ID_poly <- seq(1:nrow(poly_precip2017))
precip_PSU <- precip_PSU %>%
  left_join(poly_precip2017, by = "ID_poly")
names(precip_PSU)[names(precip_PSU) == '1'] <- 'lin_precip2017'
names(precip_PSU)[names(precip_PSU) == '2'] <- 'poly_precip2017'
poly_precip2018 <- as.matrix(poly(precip_PSU$precip2018, degree = 2))
poly_precip2018 <- as.data.frame(poly_precip2018)
precip_PSU$ID_poly <- seq(1:nrow(precip_PSU))
poly_precip2018$ID_poly <- seq(1:nrow(poly_precip2018))
precip_PSU <- precip_PSU %>%
  left_join(poly_precip2018, by = "ID_poly")
names(precip_PSU)[names(precip_PSU) == '1'] <- 'lin_precip2018'
names(precip_PSU)[names(precip_PSU) == '2'] <- 'poly_precip2018'
poly_precip2019 <- as.matrix(poly(precip_PSU$precip2019, degree = 2))
poly_precip2019 <- as.data.frame(poly_precip2019)
precip_PSU$ID_poly <- seq(1:nrow(precip_PSU))
poly_precip2019$ID_poly <- seq(1:nrow(poly_precip2019))
precip_PSU <- precip_PSU %>%
  left_join(poly_precip2019, by = "ID_poly")
names(precip_PSU)[names(precip_PSU) == '1'] <- 'lin_precip2019'
names(precip_PSU)[names(precip_PSU) == '2'] <- 'poly_precip2019'
poly_precip2020 <- as.matrix(poly(precip_PSU$precip2020, degree = 2))
poly_precip2020 <- as.data.frame(poly_precip2020)
precip_PSU$ID_poly <- seq(1:nrow(precip_PSU))
poly_precip2020$ID_poly <- seq(1:nrow(poly_precip2020))
precip_PSU <- precip_PSU %>%
  left_join(poly_precip2020, by = "ID_poly")
names(precip_PSU)[names(precip_PSU) == '1'] <- 'lin_precip2020'
names(precip_PSU)[names(precip_PSU) == '2'] <- 'poly_precip2020'
poly_precip2021 <- as.matrix(poly(precip_PSU$precip2021, degree = 2))
poly_precip2021 <- as.data.frame(poly_precip2021)
precip_PSU$ID_poly <- seq(1:nrow(precip_PSU))
poly_precip2021$ID_poly <- seq(1:nrow(poly_precip2021))
precip_PSU <- precip_PSU %>%
  left_join(poly_precip2021, by = "ID_poly")
names(precip_PSU)[names(precip_PSU) == '1'] <- 'lin_precip2021'
names(precip_PSU)[names(precip_PSU) == '2'] <- 'poly_precip2021'

precip_lin_PSU <- matrix(NA, nrow = length(precip_PSU$precip2017), ncol = 5)
precip_lin_PSU[, 1] <- precip_PSU[, 'lin_precip2017']
precip_lin_PSU[, 2] <- precip_PSU[, 'lin_precip2018']
precip_lin_PSU[, 3] <- precip_PSU[, 'lin_precip2019']
precip_lin_PSU[, 4] <- precip_PSU[, 'lin_precip2020']
precip_lin_PSU[, 5] <- precip_PSU[, 'lin_precip2021']
precip_lin_PSU <- scale(precip_lin_PSU) # scale

precip_poly_PSU <- matrix(NA, nrow = length(precip_PSU$precip2017), ncol = 5)
precip_poly_PSU[, 1] <- precip_PSU[, 'poly_precip2017']
precip_poly_PSU[, 2] <- precip_PSU[, 'poly_precip2018']
precip_poly_PSU[, 3] <- precip_PSU[, 'poly_precip2019']
precip_poly_PSU[, 4] <- precip_PSU[, 'poly_precip2020']
precip_poly_PSU[, 5] <- precip_PSU[, 'poly_precip2021']
precip_poly_PSU <- scale(precip_poly_PSU) # scale

# UEDEP polynomial in PSU vsn
poly_UEDEP2 <- as.matrix(poly(psu.extraction2$mean_UEDEP, degree = 2))
poly_UEDEP2 <- as.data.frame(poly_UEDEP2)
psu.extraction2$ID_poly <- seq(1:nrow(psu.extraction2))
poly_UEDEP2$ID_poly <- seq(1:nrow(poly_UEDEP2))
psu.extraction2 <- psu.extraction2 %>%
  left_join(poly_UEDEP2, by = "ID_poly")
names(psu.extraction2)[names(psu.extraction2) == '1'] <- 'lin_UEDEP'
names(psu.extraction2)[names(psu.extraction2) == '2'] <- 'poly_UEDEP'
UEDEP_lin_PSU <- matrix(NA, nrow = length(psu.extraction2$mean_UEDEP), ncol = 5)
UEDEP_lin_PSU[, c(1:5)] <- psu.extraction2[, "lin_UEDEP"]
UEDEP_lin_PSU <- scale(UEDEP_lin_PSU) # scale
UEDEP_poly_PSU <- matrix(NA, nrow = length(psu.extraction2$mean_UEDEP), ncol = 5)
UEDEP_poly_PSU[, c(1:5)] <- psu.extraction2[, "poly_UEDEP"]
UEDEP_poly_PSU <- scale(UEDEP_poly_PSU) # scale


#save(temp_lin_PSU, temp_poly_PSU, precip_lin_PSU, precip_poly_PSU, UEDEP_lin_PSU, 
#     UEDEP_poly_PSU, file = "psu.covs.final3.RData") # saving the final PSU covariates 

#### Detection covariates ----
# we want to match the SSUid from jags.site.covs and the matrices above, with SSUid in the jags.obs.covs data
# SSUid should be numeric

SSU.num = ccover_SSU %>%
  dplyr::select(pointnum) %>%
  mutate(SSU.num = as.numeric(1:nrow(ccover_SSU)))

jags.obs.covs = jags.obs.covs %>%
  left_join(SSU.num)

# so SSU.num is also represented by ID

# centering and scaling all continuous variables
jdate.scale <- scale(jags.obs.covs$jdate)
jags.obs.covs$jdate_scale <- jdate.scale[,1]

min.scale <- scale(jags.obs.covs$PointVisitStartTime)
jags.obs.covs$min_scale <- min.scale[,1]


#### last necessary numerical inputs ----
nyears = jags.obs.covs$year %>% unique() %>% length()
nPSUobs = jags.site.covs$PSU.num %>% unique() %>% length()
nSSUobs = jags.site.covs$pointnum %>% unique() %>% length()
nsurveys = jags.obs.covs %>% nrow()
nWindows = 3 # we had 3 (2min duration each) temporal survey windows
nobservers = jags.obs.covs$observer.num %>% unique() %>% length()

yearPSU = matrix(rep(1:nyears, each = nPSUobs), nPSUobs, nyears) #488 PSUs over 5 years sampling
yearSSU = matrix(rep(1:nyears, each = nSSUobs), nSSUobs, nyears) #6511 SSUs over 5 years sampling

#### saving everything into our final input .RData
#save(nobservers, nPSUobs, nSSUobs, nsurveys, nWindows, nyears,
#     yearPSU, temp_lin_PSU, precip_lin_PSU, precip_poly_PSU, temp_poly_PSU, UEDEP_lin_PSU, UEDEP_poly_PSU,
#     yearSSU, UEDEP_lin_SSU, UEDEP_poly_SSU, droad_lin_SSU, droad_poly_SSU, ccover.scale, shrubland_SSU, o.t.canopy_SSU, herb_SSU, PSUid_SSU,
#     jags.obs.covs, file = "jags.model.inputs.final3.RData")