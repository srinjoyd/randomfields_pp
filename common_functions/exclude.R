
exclude <- function(b, ix, iy, colwise=TRUE) {

b <- as.matrix(b)

if (colwise==TRUE) {

# identify and exclude columns
dc <- iy-1
if (dc >= 1) {
  c_ex <- c(1:dc)
  b[,c_ex] =NA
}

# identify and exclude rows
dr <- ix-1
if (dr >= 1) {
  r_ex <- c(1:dr)
  b[r_ex,iy] <- NA
}

} else {

# identify and exclude rows
dr <- ix-1
if (dr >= 1) {
  r_ex <- c(1:dr)
  b[r_ex,] =NA
}

# identify and exclude columns
dc <- iy-1
if (dc >= 1) {
  c_ex <- c(1:dc)
  b[ix,c_ex] <- NA
}

}

return(b)

}
