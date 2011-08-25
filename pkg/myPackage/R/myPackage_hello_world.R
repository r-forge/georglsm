
####################################
#### Trivial
####################################
myPackage_hello_world <- function(){
	.Call( "myPackage_hello_world", PACKAGE = "myPackage" )
}


colNorm <- function(M) {
  stopifnot(is.matrix(M))
	res <- .Call("colNorm", M, package="myPackage")
  res
}
####################################
#### END
####################################