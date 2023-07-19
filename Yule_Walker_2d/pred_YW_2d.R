
pred_YW_2d <- function(field, ix, iy, YW_model, debug=FALSE, debug_boot=FALSE, boot=FALSE) {

if (boot == FALSE){
  # [bic_array, beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9] = betadiam(1.0,field);
  betadiam_out <- betadiam(field, debug=debug)
  
  # determine which bic is minimum
  bic_array <- unlist(betadiam_out[[1]])
  bicmin_idx <- which(bic_array==min(bic_array))
  
  # extract the AR coefficients
  beta1 <- betadiam_out[[2]]
  beta2 <- betadiam_out[[3]]
  beta3 <- betadiam_out[[4]]
  if (debug==TRUE) {
    beta4 <- betadiam_out[[5]]
    beta5 <- betadiam_out[[6]]
    beta6 <- betadiam_out[[7]]
    beta7 <- betadiam_out[[8]]
    beta8 <- betadiam_out[[9]]
    beta9 <- betadiam_out[[10]]
  }
} else{
  YW_model <- YW_model
  bic_array <- unlist(YW_model[[1]])
  bicmin_idx <- which(bic_array==min(bic_array))
  beta1 <- YW_model[[2]]
  beta2 <- YW_model[[3]]
  beta3 <- YW_model[[4]]
}

# subtract the mean of the field
field_mean <- mean(field[( (!is.na(field)) & (!is.infinite(field)) )])
field <- field - field_mean

if(debug_boot == TRUE){
  pred_val1 <- 0.25*field[(ix+1),(iy+1)]
  pred_val2 <- 0.20*field[(ix+1),(iy-1)]
  pred_val3 <- (-0.05)*field[(ix+2),iy]
  
  pred_val <- pred_val1 + pred_val2 + pred_val3 + rnorm(1, sd=0.1)
} else{
  # print out the beta values corresponding to minimum BIC
  if (bicmin_idx==1) {
    c01 <- beta1[[1]][1]
    c10 <- beta1[[1]][2]
    pred_val <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
  } else if (bicmin_idx==2) {
    c01 <- beta2[[1]][1]
    c10 <- beta2[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta2[[2]][1]
    c11 <- beta2[[2]][2]
    c20 <- beta2[[2]][3]
    c1m1 <- beta2[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    pred_val <- pred_val1 + pred_val2
  } else if (bicmin_idx==3) {
    c01 <- beta3[[1]][1]
    c10 <- beta3[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta3[[2]][1]
    c11 <- beta3[[2]][2]
    c20 <- beta3[[2]][3]
    c1m1 <- beta3[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta3[[3]][1]
    c12 <- beta3[[3]][2]
    c21 <- beta3[[3]][3]
    c30 <- beta3[[3]][4]
    c2m1 <- beta3[[3]][5]
    c1m2 <- beta3[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3
  } else if (bicmin_idx==4) {
    c01 <- beta4[[1]][1]
    c10 <- beta4[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta4[[2]][1]
    c11 <- beta4[[2]][2]
    c20 <- beta4[[2]][3]
    c1m1 <- beta4[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta4[[3]][1]
    c12 <- beta4[[3]][2]
    c21 <- beta4[[3]][3]
    c30 <- beta4[[3]][4]
    c2m1 <- beta4[[3]][5]
    c1m2 <- beta4[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    c04 <- beta4[[4]][1]
    c13 <- beta4[[4]][2]
    c22 <- beta4[[4]][3]
    c31 <- beta4[[4]][4]
    c40 <- beta4[[4]][5]
    c3m1 <- beta4[[4]][6]
    c2m2 <- beta4[[4]][7]
    c1m3 <- beta4[[4]][8]
    pred_val4 <- (c04*field[ix,iy+4]) + (c13*field[ix+1,iy+3]) + (c22*field[ix+2,iy+2]) + (c31*field[ix+3,iy+1]) + (c40*field[ix+4,iy]) + (c3m1*field[ix+3,iy-1]) + (c2m2*field[ix+2,iy-2]) + (c1m3*field[ix+1,iy-3])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3 + pred_val4
  } else if (bicmin_idx==5) {
    c01 <- beta5[[1]][1]
    c10 <- beta5[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta5[[2]][1]
    c11 <- beta5[[2]][2]
    c20 <- beta5[[2]][3]
    c1m1 <- beta5[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta5[[3]][1]
    c12 <- beta5[[3]][2]
    c21 <- beta5[[3]][3]
    c30 <- beta5[[3]][4]
    c2m1 <- beta5[[3]][5]
    c1m2 <- beta5[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    c04 <- beta5[[4]][1]
    c13 <- beta5[[4]][2]
    c22 <- beta5[[4]][3]
    c31 <- beta5[[4]][4]
    c40 <- beta5[[4]][5]
    c3m1 <- beta5[[4]][6]
    c2m2 <- beta5[[4]][7]
    c1m3 <- beta5[[4]][8]
    pred_val4 <- (c04*field[ix,iy+4]) + (c13*field[ix+1,iy+3]) + (c22*field[ix+2,iy+2]) + (c31*field[ix+3,iy+1]) + (c40*field[ix+4,iy]) + (c3m1*field[ix+3,iy-1]) + (c2m2*field[ix+2,iy-2]) + (c1m3*field[ix+1,iy-3])
    
    c05 <- beta5[[5]][1]
    c14 <- beta5[[5]][2]
    c23 <- beta5[[5]][3]
    c32 <- beta5[[5]][4]
    c41 <- beta5[[5]][5]
    c50 <- beta5[[5]][6]
    c4m1 <- beta5[[5]][7]
    c3m2 <- beta5[[5]][8]
    c2m3 <- beta5[[5]][9]
    c1m4 <- beta5[[5]][10]
    pred_val5 <- (c05*field[ix,iy+5]) + (c14*field[ix+1,iy+4]) + (c23*field[ix+2,iy+3]) + (c32*field[ix+3,iy+2]) + (c41*field[ix+4,iy+1]) + (c50*field[ix+5,iy]) + (c4m1*field[ix+4,iy-1]) + (c3m2*field[ix+3,iy-2]) + (c2m3*field[ix+2,iy-3]) + (c1m4*field[ix+1,iy-4])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3 + pred_val4 + pred_val5
  } else if (bicmin_idx==6) {
    c01 <- beta6[[1]][1]
    c10 <- beta6[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta6[[2]][1]
    c11 <- beta6[[2]][2]
    c20 <- beta6[[2]][3]
    c1m1 <- beta6[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta6[[3]][1]
    c12 <- beta6[[3]][2]
    c21 <- beta6[[3]][3]
    c30 <- beta6[[3]][4]
    c2m1 <- beta6[[3]][5]
    c1m2 <- beta6[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    c04 <- beta6[[4]][1]
    c13 <- beta6[[4]][2]
    c22 <- beta6[[4]][3]
    c31 <- beta6[[4]][4]
    c40 <- beta6[[4]][5]
    c3m1 <- beta6[[4]][6]
    c2m2 <- beta6[[4]][7]
    c1m3 <- beta6[[4]][8]
    pred_val4 <- (c04*field[ix,iy+4]) + (c13*field[ix+1,iy+3]) + (c22*field[ix+2,iy+2]) + (c31*field[ix+3,iy+1]) + (c40*field[ix+4,iy]) + (c3m1*field[ix+3,iy-1]) + (c2m2*field[ix+2,iy-2]) + (c1m3*field[ix+1,iy-3])
    
    c05 <- beta6[[5]][1]
    c14 <- beta6[[5]][2]
    c23 <- beta6[[5]][3]
    c32 <- beta6[[5]][4]
    c41 <- beta6[[5]][5]
    c50 <- beta6[[5]][6]
    c4m1 <- beta6[[5]][7]
    c3m2 <- beta6[[5]][8]
    c2m3 <- beta6[[5]][9]
    c1m4 <- beta6[[5]][10]
    pred_val5 <- (c05*field[ix,iy+5]) + (c14*field[ix+1,iy+4]) + (c23*field[ix+2,iy+3]) + (c32*field[ix+3,iy+2]) + (c41*field[ix+4,iy+1]) + (c50*field[ix+5,iy]) + (c4m1*field[ix+4,iy-1]) + (c3m2*field[ix+3,iy-2]) + (c2m3*field[ix+2,iy-3]) + (c1m4*field[ix+1,iy-4])
    
    c06 <- beta6[[6]][1]
    c15 <- beta6[[6]][2]
    c24 <- beta6[[6]][3]
    c33 <- beta6[[6]][4]
    c42 <- beta6[[6]][5]
    c51 <- beta6[[6]][6]
    
    c60 <- beta6[[6]][7]
    c5m1 <- beta6[[6]][8]
    c4m2 <- beta6[[6]][9]
    c3m3 <- beta6[[6]][10]
    c2m4 <- beta6[[6]][11]
    c1m5 <- beta6[[6]][12]
    pred_val6 <- (c06*field[ix,iy+6]) + (c15*field[ix+1,iy+5]) + (c24*field[ix+2,iy+4]) + (c33*field[ix+3,iy+3]) + (c42*field[ix+4,iy+2]) + (c51*field[ix+5,iy+1]) + (c60*field[ix+6,iy]) + (c5m1*field[ix+5,iy-1]) + (c4m2*field[ix+4,iy-2]) + (c3m3*field[ix+3,iy+3]) + (c2m4*field[ix+2,iy-4]) + (c1m5*field[ix+1,iy-5])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3 + pred_val4 + pred_val5 + pred_val6
  } else if (bicmin_idx==7) {
    c01 <- beta7[[1]][1]
    c10 <- beta7[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta7[[2]][1]
    c11 <- beta7[[2]][2]
    c20 <- beta7[[2]][3]
    c1m1 <- beta7[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta7[[3]][1]
    c12 <- beta7[[3]][2]
    c21 <- beta7[[3]][3]
    c30 <- beta7[[3]][4]
    c2m1 <- beta7[[3]][5]
    c1m2 <- beta7[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    c04 <- beta7[[4]][1]
    c13 <- beta7[[4]][2]
    c22 <- beta7[[4]][3]
    c31 <- beta7[[4]][4]
    c40 <- beta7[[4]][5]
    c3m1 <- beta7[[4]][6]
    c2m2 <- beta7[[4]][7]
    c1m3 <- beta7[[4]][8]
    pred_val4 <- (c04*field[ix,iy+4]) + (c13*field[ix+1,iy+3]) + (c22*field[ix+2,iy+2]) + (c31*field[ix+3,iy+1]) + (c40*field[ix+4,iy]) + (c3m1*field[ix+3,iy-1]) + (c2m2*field[ix+2,iy-2]) + (c1m3*field[ix+1,iy-3])
    
    c05 <- beta7[[5]][1]
    c14 <- beta7[[5]][2]
    c23 <- beta7[[5]][3]
    c32 <- beta7[[5]][4]
    c41 <- beta7[[5]][5]
    c50 <- beta7[[5]][6]
    c4m1 <- beta7[[5]][7]
    c3m2 <- beta7[[5]][8]
    c2m3 <- beta7[[5]][9]
    c1m4 <- beta7[[5]][10]
    pred_val5 <- (c05*field[ix,iy+5]) + (c14*field[ix+1,iy+4]) + (c23*field[ix+2,iy+3]) + (c32*field[ix+3,iy+2]) + (c41*field[ix+4,iy+1]) + (c50*field[ix+5,iy]) + (c4m1*field[ix+4,iy-1]) + (c3m2*field[ix+3,iy-2]) + (c2m3*field[ix+2,iy-3]) + (c1m4*field[ix+1,iy-4])
    
    c06 <- beta7[[6]][1]
    c15 <- beta7[[6]][2]
    c24 <- beta7[[6]][3]
    c33 <- beta7[[6]][4]
    c42 <- beta7[[6]][5]
    c51 <- beta7[[6]][6]
    
    c60 <- beta7[[6]][7]
    c5m1 <- beta7[[6]][8]
    c4m2 <- beta7[[6]][9]
    c3m3 <- beta7[[6]][10]
    c2m4 <- beta7[[6]][11]
    c1m5 <- beta7[[6]][12]
    pred_val6 <- (c06*field[ix,iy+6]) + (c15*field[ix+1,iy+5]) + (c24*field[ix+2,iy+4]) + (c33*field[ix+3,iy+3]) + (c42*field[ix+4,iy+2]) + (c51*field[ix+5,iy+1]) + (c60*field[ix+6,iy]) + (c5m1*field[ix+5,iy-1]) + (c4m2*field[ix+4,iy-2]) + (c3m3*field[ix+3,iy+3]) + (c2m4*field[ix+2,iy-4]) + (c1m5*field[ix+1,iy-5])
    
    c07 <- beta7[[7]][1]
    c16 <- beta7[[7]][2]
    c25 <- beta7[[7]][3]
    c34 <- beta7[[7]][4]
    c43 <- beta7[[7]][5]
    c52 <- beta7[[7]][6]
    c61 <- beta7[[7]][7]
    
    c70 <- beta7[[7]][8]
    c6m1 <- beta7[[7]][9]
    c5m2 <- beta7[[7]][10]
    c4m3 <- beta7[[7]][11]
    c3m4 <- beta7[[7]][12]
    c2m5 <- beta7[[7]][13]
    c1m6 <- beta7[[7]][14]
    pred_val7 <- (c07*field[ix,iy+7]) + (c16*field[ix+1,iy+6]) + (c25*field[ix+2,iy+5]) + (c34*field[ix+3,iy+4]) + (c43*field[ix+4,iy+3]) + (c52*field[ix+5,iy+2]) + (c61*field[ix+6,iy+1]) + (c70*field[ix+7,iy]) + (c6m1*field[ix+6,iy-1]) +(c5m2*field[ix+5,iy-2]) + (c4m3*field[ix+4,iy-3]) + (c3m4*field[ix+3,iy-4]) + (c2m5*field[ix+2,iy-5]) + (c1m6*field[ix+1,iy-6])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3 + pred_val4 + pred_val5 + pred_val6 + pred_val7
  } else if (bicmin_idx==8) {
    c01 <- beta8[[1]][1]
    c10 <- beta8[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta8[[2]][1]
    c11 <- beta8[[2]][2]
    c20 <- beta8[[2]][3]
    c1m1 <- beta8[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta8[[3]][1]
    c12 <- beta8[[3]][2]
    c21 <- beta8[[3]][3]
    c30 <- beta8[[3]][4]
    c2m1 <- beta8[[3]][5]
    c1m2 <- beta8[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    c04 <- beta8[[4]][1]
    c13 <- beta8[[4]][2]
    c22 <- beta8[[4]][3]
    c31 <- beta8[[4]][4]
    c40 <- beta8[[4]][5]
    c3m1 <- beta8[[4]][6]
    c2m2 <- beta8[[4]][7]
    c1m3 <- beta8[[4]][8]
    pred_val4 <- (c04*field[ix,iy+4]) + (c13*field[ix+1,iy+3]) + (c22*field[ix+2,iy+2]) + (c31*field[ix+3,iy+1]) + (c40*field[ix+4,iy]) + (c3m1*field[ix+3,iy-1]) + (c2m2*field[ix+2,iy-2]) + (c1m3*field[ix+1,iy-3])
    
    c05 <- beta8[[5]][1]
    c14 <- beta8[[5]][2]
    c23 <- beta8[[5]][3]
    c32 <- beta8[[5]][4]
    c41 <- beta8[[5]][5]
    c50 <- beta8[[5]][6]
    c4m1 <- beta8[[5]][7]
    c3m2 <- beta8[[5]][8]
    c2m3 <- beta8[[5]][9]
    c1m4 <- beta8[[5]][10]
    pred_val5 <- (c05*field[ix,iy+5]) + (c14*field[ix+1,iy+4]) + (c23*field[ix+2,iy+3]) + (c32*field[ix+3,iy+2]) + (c41*field[ix+4,iy+1]) + (c50*field[ix+5,iy]) + (c4m1*field[ix+4,iy-1]) + (c3m2*field[ix+3,iy-2]) + (c2m3*field[ix+2,iy-3]) + (c1m4*field[ix+1,iy-4])
    
    c06 <- beta8[[6]][1]
    c15 <- beta8[[6]][2]
    c24 <- beta8[[6]][3]
    c33 <- beta8[[6]][4]
    c42 <- beta8[[6]][5]
    c51 <- beta8[[6]][6]
    
    c60 <- beta8[[6]][7]
    c5m1 <- beta8[[6]][8]
    c4m2 <- beta8[[6]][9]
    c3m3 <- beta8[[6]][10]
    c2m4 <- beta8[[6]][11]
    c1m5 <- beta8[[6]][12]
    pred_val6 <- (c06*field[ix,iy+6]) + (c15*field[ix+1,iy+5]) + (c24*field[ix+2,iy+4]) + (c33*field[ix+3,iy+3]) + (c42*field[ix+4,iy+2]) + (c51*field[ix+5,iy+1]) + (c60*field[ix+6,iy]) + (c5m1*field[ix+5,iy-1]) + (c4m2*field[ix+4,iy-2]) + (c3m3*field[ix+3,iy+3]) + (c2m4*field[ix+2,iy-4]) + (c1m5*field[ix+1,iy-5])
    
    c07 <- beta8[[7]][1]
    c16 <- beta8[[7]][2]
    c25 <- beta8[[7]][3]
    c34 <- beta8[[7]][4]
    c43 <- beta8[[7]][5]
    c52 <- beta8[[7]][6]
    c61 <- beta8[[7]][7]
    
    c70 <- beta8[[7]][8]
    c6m1 <- beta8[[7]][9]
    c5m2 <- beta8[[7]][10]
    c4m3 <- beta8[[7]][11]
    c3m4 <- beta8[[7]][12]
    c2m5 <- beta8[[7]][13]
    c1m6 <- beta8[[7]][14]
    pred_val7 <- (c07*field[ix,iy+7]) + (c16*field[ix+1,iy+6]) + (c25*field[ix+2,iy+5]) + (c34*field[ix+3,iy+4]) + (c43*field[ix+4,iy+3]) + (c52*field[ix+5,iy+2]) + (c61*field[ix+6,iy+1]) + (c70*field[ix+7,iy]) + (c6m1*field[ix+6,iy-1]) +(c5m2*field[ix+5,iy-2]) + (c4m3*field[ix+4,iy-3]) + (c3m4*field[ix+3,iy-4]) + (c2m5*field[ix+2,iy-5]) + (c1m6*field[ix+1,iy-6])
    
    c08 <- beta8[[8]][1]
    c17 <- beta8[[8]][2]
    c26 <- beta8[[8]][3]
    c35 <- beta8[[8]][4]
    c44 <- beta8[[8]][5]
    c53 <- beta8[[8]][6]
    c62 <- beta8[[8]][7]
    c71 <- beta8[[8]][8]
    
    c80 <- beta8[[8]][9]
    c7m1 <- beta8[[8]][10]
    c6m2 <- beta8[[8]][11]
    c5m3 <- beta8[[8]][12]
    c4m4 <- beta8[[8]][13]
    c3m5 <- beta8[[8]][14]
    c2m6 <- beta8[[8]][15]
    c1m7 <- beta8[[8]][16]
    pred_val8 <- (c08*field[ix,iy+8]) + (c17*field[ix+1,iy+7]) + (c26*field[ix+2,iy+6]) + (c35*field[ix+3,iy+5]) + (c44*field[ix+4,iy+4]) + (c53*field[ix+5,iy+3]) + (c62*field[ix+6,iy+2]) + (c71*field[ix+7,iy+1]) + (c80+field[ix+8,iy]) + (c7m1*field[ix+7,iy-1]) + (c6m2*field[ix+6,iy-2]) + (c5m3*field[ix+5,iy-3]) + (c4m4*field[ix+4,iy+4]) + (c3m5*field[ix+3,iy-5]) + (c2m6*field[ix+2,iy-6]) + (c1m7*field[ix+1,iy-7])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3 + pred_val4 + pred_val5 + pred_val6 + pred_val7 + pred_val8
  } else {
    c01 <- beta9[[1]][1]
    c10 <- beta9[[1]][2]
    pred_val1 <- (c01*field[ix,iy+1]) + (c10*field[ix+1,iy])
    
    c02 <- beta9[[2]][1]
    c11 <- beta9[[2]][2]
    c20 <- beta9[[2]][3]
    c1m1 <- beta9[[2]][4]
    pred_val2 <- (c02*field[ix,iy+2]) + (c11*field[ix+1,iy+1]) + (c20*field[ix+2,iy]) + (c1m1*field[ix+1,iy-1])
    
    c03 <- beta9[[3]][1]
    c12 <- beta9[[3]][2]
    c21 <- beta9[[3]][3]
    c30 <- beta9[[3]][4]
    c2m1 <- beta9[[3]][5]
    c1m2 <- beta9[[3]][6]
    pred_val3 <- (c03*field[ix,iy+3]) + (c12*field[ix+1,iy+2]) + (c21*field[ix+2,iy+1]) + (c30*field[ix+3,iy]) + (c2m1*field[ix+2,iy-1]) + (c1m2*field[ix+1,iy-2])
    
    c04 <- beta9[[4]][1]
    c13 <- beta9[[4]][2]
    c22 <- beta9[[4]][3]
    c31 <- beta9[[4]][4]
    c40 <- beta9[[4]][5]
    c3m1 <- beta9[[4]][6]
    c2m2 <- beta9[[4]][7]
    c1m3 <- beta9[[4]][8]
    pred_val4 <- (c04*field[ix,iy+4]) + (c13*field[ix+1,iy+3]) + (c22*field[ix+2,iy+2]) + (c31*field[ix+3,iy+1]) + (c40*field[ix+4,iy]) + (c3m1*field[ix+3,iy-1]) + (c2m2*field[ix+2,iy-2]) + (c1m3*field[ix+1,iy-3])
    
    c05 <- beta9[[5]][1]
    c14 <- beta9[[5]][2]
    c23 <- beta9[[5]][3]
    c32 <- beta9[[5]][4]
    c41 <- beta9[[5]][5]
    c50 <- beta9[[5]][6]
    c4m1 <- beta9[[5]][7]
    c3m2 <- beta9[[5]][8]
    c2m3 <- beta9[[5]][9]
    c1m4 <- beta9[[5]][10]
    pred_val5 <- (c05*field[ix,iy+5]) + (c14*field[ix+1,iy+4]) + (c23*field[ix+2,iy+3]) + (c32*field[ix+3,iy+2]) + (c41*field[ix+4,iy+1]) + (c50*field[ix+5,iy]) + (c4m1*field[ix+4,iy-1]) + (c3m2*field[ix+3,iy-2]) + (c2m3*field[ix+2,iy-3]) + (c1m4*field[ix+1,iy-4])
    
    c06 <- beta9[[6]][1]
    c15 <- beta9[[6]][2]
    c24 <- beta9[[6]][3]
    c33 <- beta9[[6]][4]
    c42 <- beta9[[6]][5]
    c51 <- beta9[[6]][6]
    
    c60 <- beta9[[6]][7]
    c5m1 <- beta9[[6]][8]
    c4m2 <- beta9[[6]][9]
    c3m3 <- beta9[[6]][10]
    c2m4 <- beta9[[6]][11]
    c1m5 <- beta9[[6]][12]
    pred_val6 <- (c06*field[ix,iy+6]) + (c15*field[ix+1,iy+5]) + (c24*field[ix+2,iy+4]) + (c33*field[ix+3,iy+3]) + (c42*field[ix+4,iy+2]) + (c51*field[ix+5,iy+1]) + (c60*field[ix+6,iy]) + (c5m1*field[ix+5,iy-1]) + (c4m2*field[ix+4,iy-2]) + (c3m3*field[ix+3,iy+3]) + (c2m4*field[ix+2,iy-4]) + (c1m5*field[ix+1,iy-5])
    
    c07 <- beta9[[7]][1]
    c16 <- beta9[[7]][2]
    c25 <- beta9[[7]][3]
    c34 <- beta9[[7]][4]
    c43 <- beta9[[7]][5]
    c52 <- beta9[[7]][6]
    c61 <- beta9[[7]][7]
    
    c70 <- beta9[[7]][8]
    c6m1 <- beta9[[7]][9]
    c5m2 <- beta9[[7]][10]
    c4m3 <- beta9[[7]][11]
    c3m4 <- beta9[[7]][12]
    c2m5 <- beta9[[7]][13]
    c1m6 <- beta9[[7]][14]
    pred_val7 <- (c07*field[ix,iy+7]) + (c16*field[ix+1,iy+6]) + (c25*field[ix+2,iy+5]) + (c34*field[ix+3,iy+4]) + (c43*field[ix+4,iy+3]) + (c52*field[ix+5,iy+2]) + (c61*field[ix+6,iy+1]) + (c70*field[ix+7,iy]) + (c6m1*field[ix+6,iy-1]) +(c5m2*field[ix+5,iy-2]) + (c4m3*field[ix+4,iy-3]) + (c3m4*field[ix+3,iy-4]) + (c2m5*field[ix+2,iy-5]) + (c1m6*field[ix+1,iy-6])
    
    c08 <- beta9[[8]][1]
    c17 <- beta9[[8]][2]
    c26 <- beta9[[8]][3]
    c35 <- beta9[[8]][4]
    c44 <- beta9[[8]][5]
    c53 <- beta9[[8]][6]
    c62 <- beta9[[8]][7]
    c71 <- beta9[[8]][8]
    
    c80 <- beta9[[8]][9]
    c7m1 <- beta9[[8]][10]
    c6m2 <- beta9[[8]][11]
    c5m3 <- beta9[[8]][12]
    c4m4 <- beta9[[8]][13]
    c3m5 <- beta9[[8]][14]
    c2m6 <- beta9[[8]][15]
    c1m7 <- beta9[[8]][16]
    pred_val8 <- (c08*field[ix,iy+8]) + (c17*field[ix+1,iy+7]) + (c26*field[ix+2,iy+6]) + (c35*field[ix+3,iy+5]) + (c44*field[ix+4,iy+4]) + (c53*field[ix+5,iy+3]) + (c62*field[ix+6,iy+2]) + (c71*field[ix+7,iy+1]) + (c80+field[ix+8,iy]) + (c7m1*field[ix+7,iy-1]) + (c6m2*field[ix+6,iy-2]) + (c5m3*field[ix+5,iy-3]) + (c4m4*field[ix+4,iy+4]) + (c3m5*field[ix+3,iy-5]) + (c2m6*field[ix+2,iy-6]) + (c1m7*field[ix+1,iy-7])
    
    c09 <- beta9[[9]][1]
    c18 <- beta9[[9]][2]
    c27 <- beta9[[9]][3]
    c36 <- beta9[[9]][4]
    c45 <- beta9[[9]][5]
    c54 <- beta9[[9]][6]
    c63 <- beta9[[9]][7]
    c72 <- beta9[[9]][8]
    c81 <- beta9[[9]][9]
    
    c90 <- beta9[[9]][10]
    c8m1 <- beta9[[9]][11]
    c7m2 <- beta9[[9]][12]
    c6m3 <- beta9[[9]][13]
    c5m4 <- beta9[[9]][14]
    c4m5 <- beta9[[9]][15]
    c3m6 <- beta9[[9]][16]
    c2m7 <- beta9[[9]][17]
    c1m8 <- beta9[[9]][18]
    pred_val9 <- (c09*field[ix,iy+9]) + (c18*field[ix+1,iy+8]) + (c27*field[ix+2,iy+7]) + (c36*field[ix+3,iy+6]) + (c45*field[ix+4,iy+5]) + (c54*field[ix+5,iy+4]) + (c63*field[ix+6,iy+3]) + (c72*field[ix+7,iy+2]) + (c81*field[ix+8,iy+1]) + (c90*field[ix+9,iy]) + (c8m1*field[ix+8,iy-1]) + (c7m2*field[ix+7,iy-2]) + (c6m3*field[ix+6,iy+3]) + (c5m4*field[ix+5,iy-4]) + (c4m5*field[ix+4,iy-5]) + (c3m6*field[ix+3,iy-6]) + (c2m7*field[ix+2,iy-7]) + (c1m8*field[ix+1,iy-8])
    
    pred_val <- pred_val1 + pred_val2 + pred_val3 + pred_val4 + pred_val5 + pred_val6 + pred_val7 + pred_val8
  } 
}

pred_val_final <- pred_val + field_mean

if (debug==TRUE) {
  return(list(bic_array, beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9, pred_val_final))
} else {
  return(list(bic_array, beta1, beta2, beta3, pred_val_final))
}

}
