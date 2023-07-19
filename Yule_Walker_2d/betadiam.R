betadiam <- function(field, debug=FALSE){

# sobs = 64
# tobs = 64
# nobs = sobs%*%tobs
sobs <- nrow(field)
tobs <- ncol(field)
nobs <- sobs%*%tobs

sigma00 <- calc_sigma00(field)

# specify arrays to store results
if (debug==TRUE) {
  bic_array <- numeric(9)
} else {
  bic_array <- numeric(3)
}

# Initial Stage :
m <- 1
# system("date")
# print(paste(m))
invR11 <- solve(R_ab(1,1,field))
beta11 <- invR11%*%R_a0(1,field)
phi21_0 <- invR11%*%R_ab(1,2,field)
lambda2 <- 1 - R_0b(1,field)%*%beta11
sigsig <- sigma00%*%lambda2
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs

bic_array[1] <- bic

##### m=2
m <- 2
# system("date")
# print(paste(m))
theta2 <- R_ab(2,2,field) - R_ab(2,1,field)%*%phi21_0
h2 <- R_a0(2,field) - R_ab(2,1,field)%*%beta11
beta22 <- solve(theta2)%*%h2
beta21 <- beta11 - phi21_0%*%beta22
lambda3 <- lambda2 - t(h2)%*%beta22
sigsig <- sigma00%*%lambda3
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs

phi21_1 <- invR11%*%R_ab(1,3,field)

bic_array[2] <- bic

n <- 2# b = m-n=0
b <- m-n
#
# PHIn_b1 =
#            | phi_n,1  (b+1) |
#            | phi_n,2  (b+1) |
#            |   --------     |
#            | phi_n,n-1(b+1) |
#
PHI2_1 <- phi21_1
phi32_0 <- phin1n_b(m,n,PHI2_1,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI2_0 <- phi21_0
#
# PHIn_0 =
#            | phi_n,1  (0) |
#            | phi_n,2  (0) |
#            |   --------   |
#            | phi_n,n-1(0) |
#
#
# PHIn_b1 =
#            | phi_n,1  (b+1) |
#            | phi_n,2  (b+1) |
#            |   --------     |
#            | phi_n,n-1(b+1) |
#
#
#             | PHIn_b1                |      | -PHIn_0  |
# PHIn_b1 =   |                        |  +   |          |
#             | dumz = O_N(n)XN(n+a+b) |      | I_N(n)   |
#
PHI3_0 <- rbind(PHI2_1,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_0 
phi31_0 <- PHI3_0[1:dimN(1),]

##### m=3
m <- 3
# system("date")
# print(paste(m))
theta3 <- R_ab(3,3,field) - R_ab(3,1,field)%*%phi31_0 - R_ab(3,2,field)%*%phi32_0 
h3 <- R_a0(3,field) - R_ab(3,1,field)%*%beta21 - R_ab(3,2,field)%*%beta22
beta33 <- solve(theta3)%*%h3
beta31 <- beta21 - phi31_0%*%beta33
beta32 <- beta22 - phi32_0%*%beta33
lambda4 <- lambda3 - t(h3)%*%beta33
sigsig <- sigma00%*%lambda4
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_2 <- invR11%*%R_ab(1,m+1,field)

bic_array[3] <- bic

n <- 2# b = m-n = 1
b <- m-n
PHI2_2 <-  phi21_2 
phi32_1 <- phin1n_b(m,n,PHI2_2,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI3_1 <- rbind(PHI2_2,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_1
phi31_1 <- PHI3_1[1:dimN(1),]

n <- 3# b = m-n = 0
b <- m-n
phi43_0 <- phin1n_b(m,n,PHI3_1,theta3,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI4_0 <- rbind(PHI3_1,dumz) + rbind(-PHI3_0,diag(dimN(n)))%*%phi43_0
phi41_0 <- PHI4_0[1:dimN(1),]
phi42_0 <- PHI4_0[ (dimS(1)+1):dimS(2),]

######## m=4
if (debug==TRUE) {
m <- 4
# system("date")
# print(paste(m))
theta4 <- R_ab(4,4,field) - R_ab(4,1,field)%*%phi41_0 - R_ab(4,2,field)%*%phi42_0 - R_ab(4,3,field)%*%phi43_0 
h4 <- R_a0(4,field) - R_ab(4,1,field)%*%beta31 - R_ab(4,2,field)%*%beta32 - R_ab(4,3,field)%*%beta33
beta44 <- solve(theta4)%*%h4
beta41 <- beta31 - phi41_0%*%beta44
beta42 <- beta32 - phi42_0%*%beta44
beta43 <- beta33 - phi43_0%*%beta44
lambda5 <- lambda4 - t(h4)%*%beta44
sigsig <- sigma00%*%lambda5
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_3 <- invR11%*%R_ab(1,5,field)

bic_array[4] <- bic

n <- 2# b = m-n = 2
b <- m-n
PHI2_3 <- phi21_3
phi32_2 <- phin1n_b(m,n,PHI2_3,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI3_2 <- rbind(PHI2_3,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_2
phi31_2 <- PHI3_2[1:dimN(1),]

n <- 3# b = m-n = 1
b <- m-n
phi43_1 <- phin1n_b(m,n,PHI3_2,theta3,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI4_1 <- rbind(PHI3_2,dumz) + rbind(-PHI3_0,diag(dimN(n)))%*%phi43_1
phi41_1 <- PHI4_1[1:dimN(1),]
phi42_1 <- PHI4_1[ (dimS(1)+1):dimS(2),]

n <- 4# b = m-n = 0
b <- m-n
phi54_0 <- phin1n_b(m,n,PHI4_1,theta4,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI5_0 <- rbind(PHI4_1,dumz) + rbind(-PHI4_0,diag(dimN(n)))%*%phi54_0
phi51_0 <- PHI5_0[1:dimN(1),]
phi52_0 <- PHI5_0[ (dimS(1)+1):dimS(2),]
phi53_0 <- PHI5_0[ (dimS(2)+1):dimS(3),]

########### m=5
m <- 5
# system("date")
# print(paste(m))
theta5 <- R_ab(5,5,field) - R_ab(5,1,field)%*%phi51_0 - R_ab(5,2,field)%*%phi52_0 - R_ab(5,3,field)%*%phi53_0 
theta5 <- theta5 - R_ab(5,4,field)%*%phi54_0 
h5 <- R_a0(5,field) - R_ab(5,1,field)%*%beta41 - R_ab(5,2,field)%*%beta42 - R_ab(5,3,field)%*%beta43
h5 <- h5 - R_ab(5,4,field)%*%beta44
beta55 <- solve(theta5)%*%h5
beta51 <- beta41 - phi51_0%*%beta55
beta52 <- beta42 - phi52_0%*%beta55
beta53 <- beta43 - phi53_0%*%beta55
beta54 <- beta44 - phi54_0%*%beta55
lambda6 <- lambda5 - t(h5)%*%beta55
sigsig <- sigma00%*%lambda6
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_4 <- invR11%*%R_ab(1,6,field)

bic_array[5] <- bic

n <- 2# b = m-n = 3
b <- m-n
PHI2_4 <- phi21_4
phi32_3 <- phin1n_b(m,n,PHI2_4,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI3_3 <- rbind(PHI2_4,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_3
phi31_3 <- PHI3_3[1:dimN(1),]

n <- 3# b = m-n = 2
b <- m-n
phi43_2 <- phin1n_b(m,n,PHI3_3,theta3,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI4_2 <- rbind(PHI3_3,dumz) + rbind(-PHI3_0,diag(dimN(n)))%*%phi43_2
phi41_2 <- PHI5_0[1:dimN(1),]
phi42_2 <- PHI5_0[ (dimS(1)+1):dimS(2),]

n <- 4# b = m-n = 1
b <- m-n
phi54_1 <- phin1n_b(m,n,PHI4_2,theta4,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI5_1 <- rbind(PHI4_2,dumz) + rbind(-PHI4_0,diag(dimN(n)))%*%phi54_1
phi51_1 <- PHI5_1[1:dimN(1),]
phi52_1 <- PHI5_1[ (dimS(1)+1):dimS(2),]
phi53_1 <- PHI5_1[ (dimS(2)+1):dimS(3),]

n <- 5# b = m-n = 0
b <- m-n
phi65_0 <- phin1n_b(m,n,PHI5_1,theta5,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI6_0 <- rbind(PHI5_1,dumz) + rbind(-PHI5_0,diag(dimN(n)))%*%phi65_0
phi61_0 <- PHI6_0[1:dimN(1),]
phi62_0 <- PHI6_0[ (dimS(1)+1):dimS(2),]
phi63_0 <- PHI6_0[ (dimS(2)+1):dimS(3),]
phi64_0 <- PHI6_0[ (dimS(3)+1):dimS(4),]


#########  m=6
m <- 6
# system("date")
# print(paste(m))
theta6 <- R_ab(6,6,field) - R_ab(6,1,field)%*%phi61_0 - R_ab(6,2,field)%*%phi62_0 - R_ab(6,3,field)%*%phi63_0 
theta6 <- theta6 - R_ab(6,4,field)%*%phi64_0 - R_ab(6,5,field)%*%phi65_0
h6 <- R_a0(6,field) - R_ab(6,1,field)%*%beta51 - R_ab(6,2,field)%*%beta52 - R_ab(6,3,field)%*%beta53
h6 <- h6 - R_ab(6,4,field)%*%beta54 - R_ab(6,5,field)%*%beta55
beta66 <- solve(theta6)%*%h6
beta61 <- beta51 - phi61_0%*%beta66
beta62 <- beta52 - phi62_0%*%beta66
beta63 <- beta53 - phi63_0%*%beta66
beta64 <- beta54 - phi64_0%*%beta66
beta65 <- beta55 - phi65_0%*%beta66
lambda7 <- lambda6 - t(h6)%*%beta66
sigsig <- sigma00%*%lambda7
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_5 <- invR11%*%R_ab(1,7,field)

bic_array[6] <- bic

n <- 2# b = m-n = 4
b <- m-n
PHI2_5 <- phi21_5
phi32_4 <- phin1n_b(m,n,PHI2_5,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI3_4 <- rbind(PHI2_5,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_4
phi31_4 <- PHI3_4[1:dimS(1),]

n <- 3# b = m-n = 3
b <- m-n
phi43_3 <- phin1n_b(m,n,PHI3_4,theta3,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI4_3 <- rbind(PHI3_4,dumz) + rbind(-PHI3_0,diag(dimN(n)))%*%phi43_3
phi41_3 <- PHI4_3[1:dimS(1),]
phi42_3 <- PHI4_3[ (dimS(1)+1):dimS(2),]

n <- 4# b = m-n = 2
b <- m-n
phi54_2 <- phin1n_b(m,n,PHI4_3,theta4,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI5_2 <- rbind(PHI4_3,dumz) + rbind(-PHI4_0,diag(dimN(n)))%*%phi54_2
phi51_2 <- PHI5_2[1:dimS(1),]
phi52_2 <- PHI5_2[ (dimS(1)+1):dimS(2),]
phi53_2 <- PHI5_2[ (dimS(2)+1):dimS(3),]

n <- 5# b = m-n = 1
b <- m-n
phi65_1 <- phin1n_b(m,n,PHI5_2,theta5,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI6_1 <- rbind(PHI5_2,dumz) + rbind(-PHI5_0,diag(dimN(n)))%*%phi65_1
phi61_1 <- PHI6_1[1:dimS(1),]
phi62_1 <- PHI6_1[ (dimS(1)+1):dimS(2),]
phi63_1 <- PHI6_1[ (dimS(2)+1):dimS(3),]
phi64_1 <- PHI6_1[ (dimS(3)+1):dimS(4),]

n <- 6# b = m-n = 0
b <- m-n
phi76_0 <- phin1n_b(m,n,PHI6_1,theta6,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI7_0 <- rbind(PHI6_1,dumz) + rbind(-PHI6_0,diag(dimN(n)))%*%phi76_0
phi71_0 <- PHI7_0[1:dimS(1),]
phi72_0 <- PHI7_0[ (dimS(1)+1):dimS(2),]
phi73_0 <- PHI7_0[ (dimS(2)+1):dimS(3),]
phi74_0 <- PHI7_0[ (dimS(3)+1):dimS(4),]
phi75_0 <- PHI7_0[ (dimS(4)+1):dimS(5),]

#########  m=7
m <- 7
# system("date")
# print(paste(m))
theta7 <- R_ab(7,7,field) - R_ab(7,1,field)%*%phi71_0 - R_ab(7,2,field)%*%phi72_0 - R_ab(7,3,field)%*%phi73_0 
theta7 <- theta7 - R_ab(7,4,field)%*%phi74_0 - R_ab(7,5,field)%*%phi75_0 - R_ab(7,6,field)%*%phi76_0
h7 <- R_a0(7,field) - R_ab(7,1,field)%*%beta61 - R_ab(7,2,field)%*%beta62 - R_ab(7,3,field)%*%beta63
h7 <- h7 - R_ab(7,4,field)%*%beta64 - R_ab(7,5,field)%*%beta65 - R_ab(7,6,field)%*%beta66
beta77 <- solve(theta7)%*%h7
beta71 <- beta61 - phi71_0%*%beta77
beta72 <- beta62 - phi72_0%*%beta77
beta73 <- beta63 - phi73_0%*%beta77
beta74 <- beta64 - phi74_0%*%beta77
beta75 <- beta65 - phi75_0%*%beta77
beta76 <- beta66 - phi76_0%*%beta77
lambda8 <- lambda7 - t(h7)%*%beta77
sigsig <- sigma00%*%lambda8
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_6 <- invR11%*%R_ab(1,8,field)

bic_array[7] <- bic

n <- 2# b = m-n = 5
b <- m-n
PHI2_6 <- phi21_6
phi32_5 <- phin1n_b(m,n,PHI2_6,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI3_5 <- rbind(PHI2_6,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_5
phi31_5 <- PHI3_5[1:dimS(1),]

n <- 3# b = m-n = 4
b <- m-n
phi43_4 <- phin1n_b(m,n,PHI3_5,theta3,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI4_4 <- rbind(PHI3_5,dumz) + rbind(-PHI3_0,diag(dimN(n)))%*%phi43_4
phi41_4 <- PHI4_4[1:dimS(1),]
phi42_4 <- PHI4_4[ (dimS(1)+1):dimS(2),]

n <- 4# b = m-n = 3
b <- m-n
phi54_3 <- phin1n_b(m,n,PHI4_4,theta4,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI5_3 <- rbind(PHI4_4,dumz) + rbind(-PHI4_0,diag(dimN(n)))%*%phi54_3
phi51_3 <- PHI5_3[1:dimS(1),]
phi52_3 <- PHI5_3[ (dimS(1)+1):dimS(2),]
phi53_3 <- PHI5_3[ (dimS(2)+1):dimS(3),]

n <- 5# b = m-n = 2
b <- m-n
phi65_2 <- phin1n_b(m,n,PHI5_3,theta5,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI6_2 <- rbind(PHI5_3,dumz) + rbind(-PHI5_0,diag(dimN(n)))%*%phi65_2
phi61_2 <- PHI6_2[1:dimS(1),]
phi62_2 <- PHI6_2[ (dimS(1)+1):dimS(2),]
phi63_2 <- PHI6_2[ (dimS(2)+1):dimS(3),]
phi64_2 <- PHI6_2[ (dimS(3)+1):dimS(4),]

n <- 6# b = m-n = 1
b <- m-n
phi76_1 <- phin1n_b(m,n,PHI6_2,theta6,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI7_1 <- rbind(PHI6_2,dumz) + rbind(-PHI6_0,diag(dimN(n)))%*%phi76_1
phi71_1 <- PHI7_1[1:dimS(1),]
phi72_1 <- PHI7_1[ (dimS(1)+1):dimS(2),]
phi73_1 <- PHI7_1[ (dimS(2)+1):dimS(3),]
phi74_1 <- PHI7_1[ (dimS(3)+1):dimS(4),]
phi75_1 <- PHI7_1[ (dimS(4)+1):dimS(5),]


n <- 7# b = m-n = 0
b <- m-n
phi87_0 <- phin1n_b(m,n,PHI7_1,theta7,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI8_0 <- rbind(PHI7_1,dumz) + rbind(-PHI7_0,diag(dimN(n)))%*%phi87_0
phi81_0 <- PHI8_0[1:dimS(1),]
phi82_0 <- PHI8_0[ (dimS(1)+1):dimS(2),]
phi83_0 <- PHI8_0[ (dimS(2)+1):dimS(3),]
phi84_0 <- PHI8_0[ (dimS(3)+1):dimS(4),]
phi85_0 <- PHI8_0[ (dimS(4)+1):dimS(5),]
phi86_0 <- PHI8_0[ (dimS(5)+1):dimS(6),]

#########  m=8
m <- 8
# system("date")
# print(paste(m))
theta8 <- R_ab(8,8,field) - R_ab(8,1,field)%*%phi81_0 - R_ab(8,2,field)%*%phi82_0 - R_ab(8,3,field)%*%phi83_0 
theta8 <- theta8 - R_ab(8,4,field)%*%phi84_0 - R_ab(8,5,field)%*%phi85_0 - R_ab(8,6,field)%*%phi86_0 - R_ab(8,7,field)%*%phi87_0
h8 <- R_a0(8,field) - R_ab(8,1,field)%*%beta71 - R_ab(8,2,field)%*%beta72 - R_ab(8,3,field)%*%beta73
h8 <- h8 - R_ab(8,4,field)%*%beta74 - R_ab(8,5,field)%*%beta75 - R_ab(8,6,field)%*%beta76 - R_ab(8,7,field)%*%beta77
beta88 <- solve(theta8)%*%h8
beta81 <- beta71 - phi81_0%*%beta88
beta82 <- beta72 - phi82_0%*%beta88
beta83 <- beta73 - phi83_0%*%beta88
beta84 <- beta74 - phi84_0%*%beta88
beta85 <- beta75 - phi85_0%*%beta88
beta86 <- beta76 - phi86_0%*%beta88
beta87 <- beta77 - phi87_0%*%beta88
lambda9 <- lambda8 - t(h8)%*%beta88
sigsig <- sigma00%*%lambda9
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_7 <- invR11%*%R_ab(1,9,field)

bic_array[8] <- bic

n <- 2# b = m-n = 6
b <- m-n
PHI2_7 <- phi21_7
phi32_6 <- phin1n_b(m,n,PHI2_7,theta2,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI3_6 <- rbind(PHI2_7,dumz) + rbind(-PHI2_0,diag(dimN(n)))%*%phi32_6
phi31_6 <- PHI3_6[1:dimS(1),]

n <- 3# b = m-n = 5
b <- m-n
phi43_5 <- phin1n_b(m,n,PHI3_6,theta3,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI4_5 <- rbind(PHI3_6,dumz) + rbind(-PHI3_0,diag(dimN(n)))%*%phi43_5
phi41_5 <- PHI4_5[1:dimS(1),]
phi42_5 <- PHI4_5[ (dimS(1)+1):dimS(2),]

n <- 4# b = m-n = 4
b <- m-n
phi54_4 <- phin1n_b(m,n,PHI4_5,theta4,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI5_4 <- rbind(PHI4_5,dumz) + rbind(-PHI4_0,diag(dimN(n)))%*%phi54_4
phi51_4 <- PHI5_4[1:dimS(1),]
phi52_4 <- PHI5_4[ (dimS(1)+1):dimS(2),]
phi53_4 <- PHI5_4[ (dimS(2)+1):dimS(3),]

n <- 5# b = m-n = 3
b <- m-n
phi65_3 <- phin1n_b(m,n,PHI5_4,theta5,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI6_3 <- rbind(PHI5_4,dumz) + rbind(-PHI5_0,diag(dimN(n)))%*%phi65_3
phi61_3 <- PHI6_3[1:dimS(1),]
phi62_3 <- PHI6_3[ (dimS(1)+1):dimS(2),]
phi63_3 <- PHI6_3[ (dimS(2)+1):dimS(3),]
phi64_3 <- PHI6_3[ (dimS(3)+1):dimS(4),]

n <- 6# b = m-n = 2
b <- m-n
phi76_2 <- phin1n_b(m,n,PHI6_3,theta6,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI7_2 <- rbind(PHI6_3,dumz) + rbind(-PHI6_0,diag(dimN(n)))%*%phi76_2
phi71_2 <- PHI7_2[1:dimS(1),]
phi72_2 <- PHI7_2[ (dimS(1)+1):dimS(2),]
phi73_2 <- PHI7_2[ (dimS(2)+1):dimS(3),]
phi74_2 <- PHI7_2[ (dimS(3)+1):dimS(4),]
phi75_2 <- PHI7_2[ (dimS(4)+1):dimS(5),]


n <- 7# b = m-n = 1
b <- m-n
phi87_1 <- phin1n_b(m,n,PHI7_2,theta7,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI8_1 <- rbind(PHI7_2,dumz) + rbind(-PHI7_0,diag(dimN(n)))%*%phi87_1
phi81_1 <- PHI8_1[1:dimS(1),]
phi82_1 <- PHI8_1[ (dimS(1)+1):dimS(2),]
phi83_1 <- PHI8_1[ (dimS(2)+1):dimS(3),]
phi84_1 <- PHI8_1[ (dimS(3)+1):dimS(4),]
phi85_1 <- PHI8_1[ (dimS(4)+1):dimS(5),]
phi86_1 <- PHI8_1[ (dimS(5)+1):dimS(6),]

n <- 8# b = m-n = 0
b <- m-n
phi98_0 <- phin1n_b(m,n,PHI8_1,theta8,field)
dumz <- matrix(0, nrow=dimN(n), ncol=dimN(n+1+b) )
PHI9_0 <- rbind(PHI8_1,dumz) + rbind(-PHI8_0,diag(dimN(n)))%*%phi98_0
phi91_0 <- PHI9_0[1:dimS(1),]
phi92_0 <- PHI9_0[ (dimS(1)+1):dimS(2),]
phi93_0 <- PHI9_0[ (dimS(2)+1):dimS(3),]
phi94_0 <- PHI9_0[ (dimS(3)+1):dimS(4),]
phi95_0 <- PHI9_0[ (dimS(4)+1):dimS(5),]
phi96_0 <- PHI9_0[ (dimS(5)+1):dimS(6),]
phi97_0 <- PHI9_0[ (dimS(6)+1):dimS(7),]

#########  m=9
m <- 9
# system("date")
# print(paste(m))
theta9 <- R_ab(9,9,field) - R_ab(9,1,field)%*%phi91_0 - R_ab(9,2,field)%*%phi92_0 - R_ab(9,3,field)%*%phi93_0 
theta9 <- theta9 - R_ab(9,4,field)%*%phi94_0 - R_ab(9,5,field)%*%phi95_0 - R_ab(9,6,field)%*%phi96_0 - R_ab(9,7,field)%*%phi97_0
theta9 <- theta9 - R_ab(9,8,field)%*%phi98_0
h9 <- R_a0(9,field) - R_ab(9,1,field)%*%beta81 - R_ab(9,2,field)%*%beta82 - R_ab(9,3,field)%*%beta83
h9 <- h9 - R_ab(9,4,field)%*%beta84 - R_ab(9,5,field)%*%beta85 - R_ab(9,6,field)%*%beta86 - R_ab(9,7,field)%*%beta87
h9 <- h9 - R_ab(9,8,field)%*%beta88
beta99 <- solve(theta9)%*%h9
beta91 <- beta81 - phi91_0%*%beta99
beta92 <- beta82 - phi92_0%*%beta99
beta93 <- beta83 - phi93_0%*%beta99
beta94 <- beta84 - phi94_0%*%beta99
beta95 <- beta85 - phi95_0%*%beta99
beta96 <- beta86 - phi96_0%*%beta99
beta97 <- beta87 - phi97_0%*%beta99
beta98 <- beta88 - phi98_0%*%beta99
lambda10 <- lambda9 - t(h9)%*%beta99
sigsig <- sigma00%*%lambda10
npara <-  m%*%m + m + 1
aic <- log(sigsig) + 2%*%npara/nobs
bic <- log(sigsig) + log(nobs)%*%npara/nobs
bmin <- log(log(nobs))/log(nobs)
eic <-  log(sigsig) + nobs^bmin%*%log(log(nobs))%*%npara/nobs
phi21_8 <- invR11%*%R_ab(1,10,field)

bic_array[9] <- bic

return(list(list(bic_array),list(beta11),list(beta21,beta22),list(beta31,beta32,beta33),list(beta41,beta42,beta43,beta44),list(beta51,beta52,beta53,beta54,beta55),list(beta61,beta62,beta63,beta64,beta65,beta66),list(beta71,beta72,beta73,beta74,beta75,beta76,beta77),list(beta81,beta82,beta83,beta84,beta85,beta86,beta87,beta88),list(beta91,beta92,beta93,beta94,beta95,beta96,beta97,beta98,beta99)))
} else {
return(list(list(bic_array),list(beta11),list(beta21,beta22),list(beta31,beta32,beta33)))
}

}
