
###########################################################
# calculate weights for Nadaraya-Watson and local linear
###########################################################

calc.weights.rf <- function(X, x, ban, lc=TRUE, negw=TRUE) {

# Notes:
# 1) X specifies the coordinates of all points on the random field 
# 2) x is the coordinate of a specific point
x1 <- rep.row(x, nrow(X))
a <- X-x1
w <- dmvn(a, mu=rep(0,2), sigma=ban*diag(2))

if (lc==TRUE) {
  w <- w/sum(w)
} else {
  knl <- w

  s_t1_1 <- sum(knl*a[,1])
  s_t2_1 <- sum(knl*a[,2])
  s_t1_2 <- sum(knl*(a[,1])^2)
  s_t2_2 <- sum(knl*(a[,2])^2)
  s_t1_t2_1 <- sum(knl*a[,1]*a[,2])

  w <- knl * ( (s_t1_2*s_t2_2 - (s_t1_t2_1)^2) - (a[,1]*(s_t1_1*s_t2_2 - s_t2_1*s_t1_t2_1)) + (a[,2]*(s_t1_1*s_t1_t2_1 - s_t1_2*s_t2_1)) )
  if (negw==TRUE) {
    w <- w/sum(w)
  } else {
    w[w<0] <- 0
    w <- w/sum(w)
  }

}

return(w)
}

#############################################################
# calculate mean using LC or LL
#############################################################

predict_mb.rf.uncorr <- function(Y, X, ix, iy, ban, lc=TRUE){

x <- c(ix,iy)

w <- calc.weights.rf(X, x, ban, lc=lc)
mu_next <- sum(Y*w)/(sum(w))

return(mu_next)

}

#############################################################
# calculate innovations for model-based case on random field
#############################################################

getW_mb.rf <- function(Y, X, x, ban, dlength, lc=TRUE, pred=FALSE, basic=FALSE) {

# get the total number of samples
nsamples <- length(Y)

W         <- numeric(nsamples)
mu_s      <- numeric(nsamples)
sigma_y_s <- numeric(nsamples)

# calculate the innovations
for (j in (nsamples-dlength):nsamples) {
  if (pred==TRUE) {
    Yp <- Y[1:(j-1)]
    Xp <- X[1:(j-1),]
    xp <- X[j,]
    w <- calc.weights.rf(Xp, xp, ban, lc=lc)
    mu <- sum(Yp*w)/(sum(w))

    if (basic==FALSE) {
      if (lc==TRUE) {
        M  <- sum(Yp^2*w)/(sum(w))
        sigma_y <- sqrt(M - mu^2)
        } else {
        w1 <- calc.weights.rf(Xp, xp, ban, lc=FALSE, negw=FALSE)
        mu1 <- sum(Yp*w1)/(sum(w1))
        M1  <- sum(Yp^2*w1)/(sum(w1))
        sigma_y <- sqrt(M1 - mu1^2)
        }
    } else {
      sigma_y <- 1
    }

  } else {

    Yp <- Y[1:j]
    Xp <- X[1:j,]
    xp <- X[j,]
    w <- calc.weights.rf(Xp, xp, ban, lc=lc)
    mu <- sum(Yp*w)/(sum(w))

    if (basic==FALSE) {
      if (lc==TRUE) {
        M  <- sum(Yp^2*w)/(sum(w))
        sigma_y <- sqrt(M - mu^2)
        } else {
        w1 <- calc.weights.rf(Xp, xp, ban, lc=FALSE, negw=FALSE)
        mu1 <- sum(Yp*w1)/(sum(w1))
        M1  <- sum(Yp^2*w1)/(sum(w1))
        sigma_y <- sqrt(M1 - mu1^2)
        }
    } else {
      sigma_y <- 1
    }
   }

  # calculate the error proxy
  W[j] <- (Y[j]-mu)/sigma_y

  # store the mu and sigma values
  mu_s[j] <- mu
  sigma_y_s[j] <- sigma_y
}

# calculate the predicted values of mu and sigma_y
w <- calc.weights.rf(X, x, ban, lc=lc)
mu_next <- sum(Y*w)/(sum(w))

if (basic==FALSE) {
  if (lc==TRUE) {
    M <- sum(Y^2*w)/(sum(w))
    sigma_y_next <- sqrt(M-mu_next^2)
  } else {
    w1 <- calc.weights.rf(X, x, ban, lc=FALSE, negw=FALSE)
    mu1 <- sum(Y*w1)/(sum(w1))
    M1 <- sum(Y^2*w1)/(sum(w1))
    sigma_y_next <- sqrt(M1 - mu1^2)
  }
} else {
  sigma_y_next <- 1
}

if (lc==TRUE){
  return(list(W, mu_next, sigma_y_next, mu_s, sigma_y_s, w)) 
} else{
  return(list(W, mu_next, sigma_y_next, mu_s, sigma_y_s, w1))
}


}

#########################################################################
# calculate point-prediction value for model-based case on random field
#########################################################################

predict_mb.rf <- function(Y, X, ix, iy, l_field, w_field, dlength, ban, lc=TRUE, pred=FALSE, boot=FALSE, basic=FALSE) {

# get the total number of samples
nsamples <- length(Y)

# set the data length for which processing is done
dlength <- dlength

# extract the innovations 
x <- c(ix,iy)
mb_out <- getW_mb.rf(Y, X, x, ban=ban, dlength=dlength, lc=lc, pred=pred, basic=basic)
W <- mb_out[[1]]
mu_next <- mb_out[[2]]
sigma_y_next <- mb_out[[3]]

# extract the mu and sigma values obtained during 1-sided estimation of innovations (used for bootstrap only)
mu <- mb_out[[4]]
sigma_y <- mb_out[[5]]

# recreate the random field with the extracted innovations
field_ext <- matrix(NA, nrow=l_field, ncol=w_field)
indices <- X[(nsamples-dlength):nsamples,]
field_ext[indices] <- W[(nsamples-dlength):nsamples]

# solve the Yule-Walker 2d equations to get the predictor
YW_out <- pred_YW_2d(field_ext, x[1], x[2])
W_est <- YW_out[[5]]

# calculate the estimated value
Y_est <- mu_next + (sigma_y_next*W_est)

# extract residuals if data is required for bootstrap
if (boot==TRUE) {
  return(list(mb_out, field_ext, indices, YW_out, Y_est))
} else if (boot==FALSE) {
  return(Y_est)
}

}

###########################################################
# force monotonicity
###########################################################

force_monotone <- function(z, dir="right", pava=FALSE) {

flag <- TRUE
iter <- 0

if(pava==FALSE) {
  if (dir=="right") {
    while (flag == TRUE) {
      z1 <- shift(z,1)
      z1[1] <- -Inf
      zs <- sign(z-z1)
      if (any(zs<0)==TRUE) {
       z <- (zs<0)*(shift(z,1)) + ((zs>=0)*z)
    #  iter <- iter+1
      }
      else {
       flag = FALSE
      }
    }
  } else {
    while (flag == TRUE) {
      z1 <- shift(z,1,dir="left")
      z1[length(z)] <- Inf
      zs <- sign(z-z1)
      if (any(zs>0)==TRUE) {
       z <- (zs>0)*(shift(z,1,dir="left")) + ((zs<=0)*z)
    #  iter <- iter+1
      }
      else {
       flag = FALSE
      }
    }
  }
} else {
  z <- pava(z)
}

return(z)
}

#########################################################################
# generate uniform samples
#########################################################################

uniformize.kernel <- function(x, x0, inverse=FALSE, w, forcem=FALSE) {
# APRIL 2012 
# New R function that does uniformize on the basis
# of a kernel smoothed edf (need kernel > 0 everywhere
# to ensure the resulting edf is invertible)
n <- length(x)
if (n <= 1) return(NA) # ("sample size too small!")
if(missing(w)) {
w <- c(1:n)*0+1/n }
if(any(w < 0)) warning("Negative weight encountered.")
if(sum(w) != 1) {
w <- w/sum(w)
warning("Weights do not sum to 1.  Renormalizing.")
}
if ( is.na(sum(x)) ) {
w<- w[!is.na(x)]
w <- w/sum(w)
x<- x[!is.na(x)]
warning("missing values in x") }
# first compute the smoothed edf, i.e., the inverse=FALSE case 
if (forcem == TRUE) {
  body(density.default)[[15]][[4]][[4]][[3]] = substitute(print(paste("Negative weight encountered.")))
}

ppdf<- density(x, weights = w ) # default uses gaussian kernel
if(sum(ppdf$y<0)>0){
  print(1)
}
npdf<- length(ppdf$x)
delta<- ppdf$x[2]-ppdf$x[1]
pcdf<- cumsum(ppdf$y)*delta

if (forcem == TRUE) {
# pcdf1 <- force_monotone(pcdf)
  pcdf  <- pcdf/max(pcdf)
}

if(inverse) {
FFinv <- approxfun(pcdf, ppdf$x)
return(FFinv(x0))
#In the inverse case x0 denotes the vector of y values to which FFinv is applied 
}
else {
FF <- approxfun(ppdf$x, pcdf)
return(FF(x0))
#Here x0 denotes the vector of x values to which FF is applied 
}}

#########################################################################
# calculate autocorrelation matrix using 2D AR recursion
#########################################################################

calc_armaacf_2d <- function(Y, X, x, l_field, w_field) {
  
  indices <- X
  len_data <- length(Y)
  
  # recreate the random field with the extracted innovations
  field_1s_with_NA <- matrix(NA, nrow=l_field, ncol=w_field)
  field_1s_with_NA[X] <- Y
  
  xy.acf <- numeric(len_data)
  
  #c1 <- 0.5
  for (i1 in 1:len_data) {
    lag_2d <- indices[i1,] - indices[1,] 
    xy.acf[i1] <- rho(lag_2d[1], lag_2d[2], field_1s_with_NA) 
  }
  
  # Use Yule-Walker 2d equations to fit AR model
  YW_out <- pred_YW_2d(field_1s_with_NA, x[1], x[2])
  bic_array <- unlist(YW_out[[1]])
  bicmin_idx <- which(bic_array==min(bic_array))
  p = q = bicmin_idx
  
  #if(abs((bic_array[bicmin_idx-1]-bic_array[bicmin_idx])/bic_array[bicmin_idx-1])<0.01
     #& length(which((X[,1]==X[1,][1]-p) & X[,2]==X[1,][2]))==0){
    #p=p-1
    #q=q-1
  #} else{
    #return(NA)
  #}
  if(length(which((X[,1]==X[1,][1]-p) & X[,2]==X[1,][2]))==0){
    return(NA)
  }
  temp <- which((X[,1]==X[1,][1]-p) & X[,2]==X[1,][2])
  xy.ar <- numeric(len_data)
  xy.ar[1:temp] <- xy.acf[1:temp]
  
  field_temp <- matrix(NA, nrow=l_field, ncol=w_field)
  field_temp[X] <- xy.acf
  # Find the following by iterating the difference equation of AR model
  field_temp <- rbind(matrix(0,nrow=p,ncol=w_field),field_temp)
  field_temp <- rbind(field_temp, matrix(0,nrow=p,ncol=w_field))
  field_temp <- cbind(field_temp, matrix(0,nrow=l_field+2*p,ncol=q))
  field_temp <- cbind(matrix(0,nrow=l_field+2*p,ncol=q),b=field_temp)
  for (i in (temp+1):len_data){
    ix <- X[i,][1]+p
    iy <- X[i,][2]+q
    gamma_next = pred_YW_2d(field_temp, ix, iy, YW_model=YW_out, boot=TRUE)[[5]]
    xy.ar[i] <- gamma_next
  }
  
  return(xy.ar)
}

###########################################################
# remove na values from list
###########################################################

na.rm<-function(x)
{       x[!is.na(x)]
}

#########################################################################
# calculate values from covariance matrix for model-free prediction
#########################################################################

mf.predict_vector.rf <- function(Y, X, x, l_field, w_field, taper=FALSE, lmf=FALSE, legacy = FALSE, trap=FALSE, thresh=1000) {

n <- length(Y)
Ybar <- mean(Y)
Y_sd <- sqrt(var(Y))
Y <- (Y-Ybar)/Y_sd

# calculate acf values for flat-top tapered and shrinkage applied covariance matrix
# xy.acf <- calc_autocorr_banded_2d(Y, X, l_field, w_field, taper=taper)
if (legacy==FALSE){
  xy.acf <- calc_armaacf_2d(Y, X, x, l_field, w_field)
} else{
  xy.acf <- calc_autocorr_banded_2d(Y, X, l_field, w_field, taper=taper)
}
if(is.na(xy.acf)==TRUE){
  flag <- FALSE
  
  pred <- NA
  en_out <- NA
  xy.acf <- NA
  ratio <- NA
  return(list(pred, en_out, xy.acf, ratio, flag))
}
# calculate the full autocorrelation matrix
Gamman <- toeplitz(xy.acf)
Gammanp1 <- toeplitz(c(xy.acf,0))

if (trap==TRUE){
  temp <- eigen(Gamman)$values
  temp1 <- eigen(Gammanp1)$values
  ratio <- max(temp)/min(temp)
  ratio1 <- max(temp1)/min(temp1)
}

if ((sum(temp <0) == 0) && (sum(temp1 <0) == 0) && (ratio < thresh) && (ratio1 < thresh)) {
  flag <- TRUE

  GammanChol <- chol(Gamman)
  Gammanp1Chol <- chol(Gammanp1)

  if (lmf==FALSE) {
     en <- forwardsolve(t(GammanChol), Y)

     constPart <- en %*% t(Gammanp1Chol)[n+1, -(n+1)]
     varPart   <- en*Gammanp1Chol[n+1, n+1]

     en_out    <- en

  } else {
     en_actual <- forwardsolve(t(GammanChol), Y)
     en_ideal  <- rnorm(length(Y))

     constPart <- en_actual %*% t(Gammanp1Chol)[n+1, -(n+1)]
     varPart   <- en_ideal*Gammanp1Chol[n+1, n+1]

     en_out    <- en_actual

  }

  pred <- as.vector((constPart + varPart)*Y_sd + Ybar)
} else {
  flag <- FALSE

  pred <- NA
  en_out <- NA
  xy.acf <- NA
}

if (trap==TRUE){
  return(list(pred, en_out, xy.acf, ratio, flag))
} else{
  return(list(pred, en_out, xy.acf))
}

}


#########################################################################
# calculate point-prediction values for model-free case
#########################################################################

predict_mf.rf <- function(Y, X, ix, iy, l_field, w_field, ban, dlength, lc=TRUE, pred=FALSE, forcem=FALSE, taper=FALSE, lmf=FALSE, legacy=FALSE, boot=FALSE, trap=FALSE) {

x <- c(ix,iy)

# set the correct weight switch for force monotone (for local linear only)
negw=TRUE
if (lc==FALSE) {
  if (forcem == TRUE) {
    negw=TRUE
  } else {
    negw=FALSE
  }
}

# get the total number of samples
nsamples <- length(Y)
u <- numeric(nsamples)

# set the data length for which processing is done
dlength <- dlength

# apply the uniformizing transformation
for (j in (nsamples-dlength):nsamples) {
  if (pred==TRUE) {
    Yp <- Y[1:(j-1)]
    Xp <- X[1:(j-1),]
    xp <- X[j,]
  } else {
    Yp <- Y[1:j]
    Xp <- X[(1:j),]
    xp <- X[j,]
  }
  w <- calc.weights.rf(Xp, xp, ban, lc=lc, negw=negw)
  u[j] <- uniformize.kernel(Yp, Y[j], inverse=FALSE, w=w, forcem=forcem)
}

u2 <- (u[(nsamples-dlength):nsamples])

u2_len <- length(u2)
# apply the normalizing transform
normsamples <- qnorm(u2)
normsamples <- replace(normsamples, is.infinite(normsamples), NA)
normsamples <- na.rm(normsamples)

# extract the coordinates corresponding to which normal samples were generated
#indices <- X[(nsamples-dlength):nsamples,]
indices <- X[(nsamples-length(normsamples)+1):nsamples,]

# build the reverse path to predict the (n+1) normal sample
normsamples_pred_pair <- mf.predict_vector.rf(normsamples, indices, x, l_field, w_field, taper=taper, lmf=lmf, legacy=legacy, trap=trap)
normsamples_pred <- normsamples_pred_pair[[1]]
en_out <- normsamples_pred_pair[[2]] # model-free residuals, uncorrelated standard normal
xy.acf <- normsamples_pred_pair[[3]]
if (trap==TRUE){
  ratio <- normsamples_pred_pair[[4]]
  flag <- normsamples_pred_pair[[5]]
}

# apply the reverse of the normalizing transform (outputs are correlated uniform)
if (flag == TRUE) {
  unifsamples <- pnorm(normsamples_pred)

  # predict the value at coordinate "x" 
  w <- calc.weights.rf(X, x, ban, lc=lc, negw=negw)
  u1 <- na.rm(uniformize.kernel(Y, unifsamples, inverse=TRUE, w=w, forcem=forcem))
  Y_est <- mean(u1)
} else {
  Y_est <- NA
}

if (trap==FALSE){
  if (boot==TRUE) {
    return(list(en_out, indices, xy.acf, normsamples, Y_est))
  } else if (boot==FALSE) {
    return(Y_est)
  }
} else{
  if (boot==TRUE) {
    return(list(en_out, indices, xy.acf, normsamples, Y_est, ratio, flag))
  } else if (boot==FALSE) {
    return(list(Y_est, ratio, flag))
  }
}

}
