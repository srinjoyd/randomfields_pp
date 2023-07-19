cv <- function(y, ix, iy, l_field, w_field, dlength, bans, frac, lc=TRUE, pred=TRUE, boot=FALSE, basic=FALSE, window = TRUE, mb=TRUE, forcem=FALSE, legacy=FALSE, taper=FALSE, trap=FALSE){
  
  refer_point <- w_field*(ix-1)+iy
  end_point <- refer_point+dlength*(1-frac)
  mses <- numeric(length(bans))
  
  for(j in 1:length(bans)){
    sse <- 0
    num_est <- 0
    ban <- bans[j]
    n <- dlength * frac
    ratios <- c()
    print(paste("Bandwidth"))
    print(paste(ban))
    for(i in end_point:(refer_point+1)){
      #     print(paste("Iteration number"))
      #     print(paste(i))
      ix <- ceiling(i/l_field)
      iy <- i %% l_field
      if ((ix+3>l_field) | (iy+3>w_field) | (iy-3<1)){
        sse <- sse
      } else {
        field <- matrix(y, nrow=l_field, ncol=w_field)
        field_1s_with_NA <- exclude(field, ix, iy, colwise=FALSE)
        indices_rev <- which(!is.na(t(field_1s_with_NA)), arr.in=TRUE)
        temp <- indices_rev[,1]
        indices_rev[,1] <- indices_rev[,2]
        indices_rev[,2] <- temp
        field_1s_rev <- field_1s_with_NA[indices_rev]
        indices <- indices_rev[nrow(indices_rev):1,]
        field_1s <- rev(field_1s_rev)
        if (mb==TRUE){
          y_hat <- predict_mb.rf(field_1s, indices, ix, iy-1, l_field, w_field, ban=ban, dlength = n, lc=lc, pred=pred) 
          flag <- TRUE
        } else{
          pred_res <- predict_mf.rf(field_1s, indices, ix, iy-1, l_field, w_field, ban=ban, dlength = n, lc=lc, pred=pred, forcem=forcem, legacy=legacy, taper=taper, trap=trap)
          if (trap==TRUE) {
            y_hat <- pred_res[[1]]
            eig_ratio <- pred_res[[2]]
            flag <- pred_res[[3]]
            ratios <- rbind(ratios,c(i,eig_ratio))
            #           print(ratios)
          } else {
            y_hat <- pred_res
          }
        }
        if (flag==TRUE) {
          sse <- sse + (y_hat-field[ix,iy-1])^2
          num_est <- num_est + 1
        }
      }  
      if(window == FALSE){
        n <- n+1
      }
    }
    mses[j] <- sse/num_est
    print(paste("Number of points used"))
    print(paste(num_est))
    #num_est <- 0
  }
  ban_opt <- bans[which.min(mses)]
  print('1')
  # print(paste(ban_opt))
  return(ban_opt)
}