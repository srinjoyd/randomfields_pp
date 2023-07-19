
calc_sigma00 <- function(field) {

# get dimensions of field 
sobs <- nrow(field)
tobs <- ncol(field)
nobs <- sobs*tobs

# initialize sigma00
sigma00 <- 0

# calculate the mean of the field
fbar  <- mean(field[( (!is.na(field)) & (!is.infinite(field)) )])

# calculate the number of NA positions
num_NA <- 0

# calculate sigma00
for (t1 in 1:tobs){
 for (s1 in 1:sobs){
     temp <- (field[t1,s1]-fbar)%*%(field[t1,s1]-fbar)
     if ( (!is.na(temp)) & (!is.infinite(temp)) ) {
       sigma00 <- sigma00 + temp
     } else {
       num_NA <- num_NA + 1
     }
 }
}

nobs <- nobs - num_NA
sigma00 <- sigma00/nobs

return(sigma00)

}

