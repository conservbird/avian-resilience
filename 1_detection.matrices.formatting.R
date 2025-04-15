#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Detection matrices for JAGS
# From "Quantifying avian resilience to habitat change to support conservation
# decision making"
# Code by Amanda L. Hayes-Puttfarcken
# Affiliated with Utah State University
# 04/11/2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#### Load libraries
library(magrittr)
library(dplyr)
library(stringr)

#### Load species detection data
bird.p <- read.csv("bird.presences.final.forpub.csv") # file to use of IMBCR bird detection data
get(load("all.site.names.RData")) # called names.vect

# List of BirdCodes to iterate over
bird_codes <- unique(bird.p$BirdCode)

# Loop through each BirdCode
for (code in bird_codes) {
  
  # Selecting one species for the matrix
  bird_detect <- subset(bird.p, BirdCode == code & year %in% c(2017:2021))
  
  # Creating the time categories for detection
  bird_detect$detect <- case_when(
    bird_detect$TimePeriod %in% c(1,2) ~ 1,
    bird_detect$TimePeriod %in% c(3,4) ~ 2,
    bird_detect$TimePeriod %in% c(5,6) ~ 3
  )
  
  # Setting up the blank output matrix
  bird_matrix <- matrix(4, nrow = length(names.vect), ncol = 2)
  bird_matrix[,1] <- names.vect
  colnames(bird_matrix) <- c("pointnum_time_year", "detect")
  bird_matrix <- as.data.frame(bird_matrix)
  
  # Convert ID fields to factors
  bird_matrix$pointnum_time_year <- as.factor(bird_matrix$pointnum_time_year)
  bird_detect$pointnum_time_year <- as.factor(bird_detect$pointnum_time_year)
  bird.p$pointnum_time_year <- as.factor(bird.p$pointnum_time_year)
  
  # Match IDs and update detection values
  bird_matrix$detect[match(bird_detect$pointnum_time_year, bird_matrix$pointnum_time_year)] <- bird_detect$detect
  
  # Assign NA to IDs not present in bird.p
  bird_matrix$detect[which(!(bird_matrix$pointnum_time_year %in% bird.p$pointnum_time_year))] <- NA
  
  # Process detection times
  str_sub(bird_matrix$pointnum_time_year, -7, -6) <- ""  # Rename to group same-site time periods
  bird_matrix <- bird_matrix %>%
    group_by(pointnum_time_year) %>%
    summarise(detect = paste(detect, collapse=", "))  # Group by site
  
  bird_matrix$min <- rep(NA, nrow(bird_matrix))  # Blank row for min value
  
  # Split detection records
  detect_split <- strsplit(bird_matrix$detect, ", ")  
  detect_reformatted <- as.data.frame(do.call(rbind, detect_split))  
  bird_matrix$min <- apply(detect_reformatted, 1, FUN = min)  # First time period seen
  
  # Clean up
  bird_matrix <- bird_matrix[, -2]  # Remove intermediate column
  bird_matrix$min <- ifelse(bird_matrix$min == "NA", NA, bird_matrix$min)
  
  # Assign inout values
  bird_matrix$inout <- rep(NA, nrow(bird_matrix))
  new_sites <- jags.obs.covs$pointnum_time_year
  for (i in unique(new_sites)) {
    bird_matrix$inout[bird_matrix$pointnum_time_year == i] <- 1
  }
  
  bird_matrix <- na.omit(bird_matrix)  # Remove NA rows
  
  # Save the CSV with a unique filename for each BirdCode
  filename <- paste0(code, ".uedep.csv")
  write.csv(bird_matrix, file = filename, row.names = FALSE)
}
