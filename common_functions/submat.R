
submat <- function(b, offset_lr, offset_rr, offset_lc, offset_rc) {

# get number of rows and columns of original matrix
nrb <- nrow(b)
ncb <- ncol(b)

# initialize the matrix from which indices will be extracted
b1 <- matrix(0, nrow=nrow(b), ncol=ncol(b))

# form the matrix with only specified rows and columns populated
b1[((1+offset_lr):(nrb-offset_rr)),((1+offset_lc):(ncb-offset_rc))] <- b[((1+offset_lr):(nrb-offset_rr)),((1+offset_lc):(ncb-offset_rc))]

# diff the 2 matrices
diffb <- b - b1

# extract the indices corresponding to the original matrix where the submatrix was populated
indices <- which(diffb==0, arr.in=TRUE)

# extract the submatrix
b1_out <- b[((1+offset_lr):(nrb-offset_rr)),((1+offset_lc):(ncb-offset_rc))]

return(list(indices, b1_out))

}
