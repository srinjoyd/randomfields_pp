
# clean up existing workspace
graphics.off()
rm(list=ls())

# set arguments
args <- commandArgs(trailingOnly = TRUE)
gix <- as.numeric(args[1])
giy <- as.numeric(args[2])

# specify field dimensions
gl_field <- 101
gw_field <- 101

#### Run parameters
nruns <- 100

### Parallel computing, make the cores be 48 for point prediction
mc=40
run1 <- function(i) {
  source("../Yule_Walker_2d/betadiam.R")
  source("../Yule_Walker_2d/phin1n_b.R")
  source("../Yule_Walker_2d/dimN.R")
  source("../Yule_Walker_2d/dimS.R")
  source("../Yule_Walker_2d/R_0b.R")
  source("../Yule_Walker_2d/R_a0.R")
  source("../Yule_Walker_2d/R_ab.R")
  source("../Yule_Walker_2d/rho.R")
  source("../Yule_Walker_2d/rhovec.R")
  source("../Yule_Walker_2d/calc_sigma00.R")
  source("../Yule_Walker_2d/pred_YW_2d.R")
  source("../common_functions/exclude.R")
  source("../common_functions/repmatrix.R")
  source("../common_functions/rf_functions.R")
  library(snow)
  library(mvnfast)
  y <- as.matrix(read.table("../datafiles_halfdia_staterrors_sd_0pt1/Y_rf_staterrors_sd_0pt1_short.txt"))
  ix <- .GlobalEnv$gix
  iy <- .GlobalEnv$giy
  l_field <- .GlobalEnv$gl_field
  w_field <- .GlobalEnv$gw_field
  bans <- as.numeric(read.table("../bw_all_nruns_100/ban_opt_mb_ll.txt"))
  y1 <- y[i,]
  field <- matrix(y1, nrow=l_field, ncol=w_field)
  field_1s_with_NA <- exclude(field, ix, iy, colwise=FALSE)
  indices_rev <- which(!is.na(t(field_1s_with_NA)), arr.in=TRUE)
  temp <- indices_rev[,1]
  indices_rev[,1] <- indices_rev[,2]
  indices_rev[,2] <- temp
  field_1s_rev <- field_1s_with_NA[indices_rev]
  indices <- indices_rev[nrow(indices_rev):1,]
  field_1s <- rev(field_1s_rev)
  predict_mb.rf(field_1s, indices, ix, iy-1, l_field, w_field, dlength=500, ban=bans[i], pred=FALSE, lc=FALSE)
}

library(snow)
library(mvnfast)
source("../Yule_Walker_2d/betadiam.R")
source("../Yule_Walker_2d/phin1n_b.R")
source("../Yule_Walker_2d/dimN.R")
source("../Yule_Walker_2d/dimS.R")
source("../Yule_Walker_2d/R_0b.R")
source("../Yule_Walker_2d/R_a0.R")
source("../Yule_Walker_2d/R_ab.R")
source("../Yule_Walker_2d/rho.R")
source("../Yule_Walker_2d/rhovec.R")
source("../Yule_Walker_2d/calc_sigma00.R")
source("../Yule_Walker_2d/pred_YW_2d.R")
source("../common_functions/exclude.R")
source("../common_functions/repmatrix.R")
source("../common_functions/rf_functions.R")

print(paste(gix))
print(paste(giy))
# print(paste(gban))

cl <- makeCluster(mc, outfile="run.txt")
clusterExport(cl, c('predict_mb.rf','gix','giy','gl_field','gw_field'), envir=environment())
t2=system.time({
  test <- do.call(c, parLapply(cl, 1:nruns, run1) )
})
stopCluster(cl)

################################################
# collect the results for post-processing
################################################

# extract the predicted values
results <- t(matrix(test,ncol=nruns))

# extract the predicted values
pred_val <- as.numeric(results[,1])

# estimate the MSE of prediction
y <- as.matrix(read.table("../datafiles_halfdia_staterrors_sd_0pt1/Y_rf_staterrors_sd_0pt1_short.txt"))

mse_error_pred_cum <- 0
for (i in 1:nruns) {
  y1 <- y[i,]
  field <- matrix(y1, nrow=gl_field, ncol=gw_field)
  mse_error_pred_cum <- mse_error_pred_cum + ((field[gix, (giy-1)] - pred_val[i])^2)
}

mse_error_pred <- (mse_error_pred_cum)/nruns

print(paste("MSE of prediction"))
print(mse_error_pred)

# save.image("YW_2D")

quit()
