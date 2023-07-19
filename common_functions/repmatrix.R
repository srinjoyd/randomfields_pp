
# repeat rows in a matrix (nrow=n)
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

# repeat columns in a matrix (ncol=n)
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
