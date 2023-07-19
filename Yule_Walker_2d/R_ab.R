R_ab <- function(n,l, field){
#
# R =
#  | rhovec(a_n1 ,    b_n1 ; l)  |
#  | rhovec(a_n2 ,    b_n2 ; l)  |
#  |   -------------             |
#  | rhovec(a_nN(n),b_nN(n) ; l) |
#
R1 <- NULL
R2 <- NULL
for (ndum in 1:n){
   dumR1 <- rhovec(ndum-1,n-ndum+1,l, field)
   dumR2 <- rhovec(n-ndum+1,1-ndum,l, field)
   R1 <- rbind(R1,dumR1)
   R2 <- rbind(R2,dumR2)
}
R <- rbind(R1,R2)

return(R)
}
