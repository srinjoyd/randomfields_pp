rhovec <- function(j,k,n,field){
#
# rhovec = [ rho(j-a_n1,k-b_n1)  rho(j-a_n2,k-b_n2) ...   rho(j-a_nN(n),k-b_nN(n)) ]
#
if (n==0){
  rhovec <- rho(j,k,field)
} else {
    rho1 <- NULL
    rho2 <- NULL
    for (duma in 1:n){
        rho1 <- cbind(rho1, rho(j-duma+1,k-n+duma-1,field))
        rho2 <- cbind(rho2, rho(j-n+duma-1,k+duma-1,field))
    }
    rhovec <- cbind(rho1,rho2)
}

return(rhovec)

}
