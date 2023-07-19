rho <- function(jj,kk,field){
#
if (jj < 0){
    jj <- -jj
    kk <- -kk
}
y <- 0
# a = abs(jj+kk);
# b = abs(jj-kk);
# if ( mod(a,2)==0 ) & ( mod(b,2)==0 )
# y = 0.25^(a/2)*0.20^(b/2);

########################################################################
# start translated R code
########################################################################

# get the dimensions of the subfield to be used in calculation of autocorrelation
T1 <- nrow(field)
S1 <- ncol(field)
nobs <- T1*S1

# calculate the mean of the field
fbar  <- mean(field[( (!is.na(field)) & (!is.infinite(field)) )])

# initialize the autocovariance
sigma <- 0

# calculate the number of NA positions
num_NA1 <- 0

# calculation of autocovariance
if (sign(jj*kk)==-1){
  jj <- abs(jj)
  kk <- abs(kk)
  for (t1 in 1:(T1-jj)){
   for (s1 in (kk+1):S1){
     temp <- (field[t1+jj,s1-kk]-fbar)*(field[t1,s1]-fbar)
     if ( (!is.na(temp)) & (!is.infinite(temp)) ) {
       sigma <- sigma + temp
     } else {
       num_NA1 <- num_NA1 + 1
     }
   }
  }
} else {
  jj <- abs(jj)
  kk <- abs(kk)
  for (t1 in 1:(T1-jj)){
   for (s1 in 1:(S1-kk)){
     temp <- (field[t1+jj,s1+kk]-fbar)*(field[t1,s1]-fbar)
     if ( (!is.na(temp)) & (!is.infinite(temp)) ) {
       sigma <- sigma + temp
     } else {
       num_NA1 <- num_NA1 + 1
     }
   }
  }
}

nobs1 <- nobs - num_NA1 
sigma <- sigma/nobs1

# initialize the variance
sigma_00 <- 0

# calculate the number of NA positions
num_NA2 <- 0

# calculation of variance
for (t1 in 1:T1){
 for (s1 in 1:S1){
     temp <- (field[t1,s1]-fbar)*(field[t1,s1]-fbar)
     if ( (!is.na(temp)) & (!is.infinite(temp)) ) {
     sigma_00 <- sigma_00 + temp
     } else {
     num_NA2 <- num_NA2 + 1
     }
 }
}

nobs2 <- nobs - num_NA2
sigma_00 <- sigma_00/nobs2

y <- sigma/sigma_00

return(y)

}
