phin1n_b <- function(m,n,PHIn_b1,thetam,field){
# Calculate phi_n+1,n (m-n)
#
# dumPHI =
#  | I_N(m+1)              |
#  | -Phi_n,1   (m-n+1)    |
#  | -Phi_n,2   (m-n+1)    |
#  |   -------------       |
#  | -Phi_n,n-1 (m-n+1)    |
#
# dumR = [ R_n,m+1 , R_n,1 , R_n,2 , .... , R_n,n-1 ]
#
b <- m-n
dumR <- R_ab(n,m+1,field)
for (a in 1:(n-1)){
   dumR <- cbind(dumR,R_ab(n,a,field))
}
dumPHI <- diag(dimN(m+1))
dumPHI <- rbind(dumPHI,-PHIn_b1)
phin1n_b <- solve(thetam)%*%dumR%*%dumPHI

return(phin1n_b)

}
