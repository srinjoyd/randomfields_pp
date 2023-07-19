dimS <- function(n){
# Value of S(n) = N_n = sum_a=1^n N(a)
dimS <- 0
for (a in 1:n){
    dimS <- dimS + dimN(a)
}

return(dimS)
}
