#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculating Resilience Metrics
# From "Quantifying avian resilience to habitat change to support conservation
# decision making"
# Code by Erica Stuber
# Affiliated with Utah State University and USGS
# 10/17/2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# R code to calculate species’ RT and GMIP values based on first and second derivatives
# of predicted occupancy curve, with special rule if predicted curve does not begin in a
# tail.

vdep_raw_sim = seq(0, 40 + 3 * .1, by=.1)
d1 <- apply(predicted_occupancy_matrix, 1, diff) # first derivative using differencing
d2 <- apply(d1, 2, diff) # second derivative using differencing
ind = 1:(ncol(predicted_occupancy_matrix)-3)
GRIP_i <- vdep_raw_sim[apply(abs(d1[ind,]),2, which.max)]
RT_i <- (abs(d2 [2:nrow(d2),]) < abs(d2 [1:(nrow(d2)-1),]))[ind,] |> apply(2,function(x){min(which(x))}) # first find left-most max d2
RT_i [sign(d1 [1,]) != sign(d2 [1,])] = 1 # then special rule if we don’t start in a tail; if we
# are already decelerating at 0, then rt=0
RT_i = ifelse(is.infinite(RT_i), RT_i, vdep_raw_sim[pmin(RT_i, max(ind))])