getSVD <- function(ii, tX){

	## Data
	myy <- tX[ii,]
	myX <- t(tX[-ii,])
	colnames(myX) <- paste("X",1:ncol(myX),sep="")

	# SVD
	mysvd <- svd(myX)
	myF <- mysvd$u %*% diag(mysvd$d) # F=UD; X=FV^{T}
	colnames(myF) <- paste("F",1:ncol(myF),sep="")
	then <- nrow(myF)
	thep <- ncol(myF)
	FTF <- crossprod(myF)

	list("u"=mysvd$u, "d"=mysvd$d, "v"=mysvd$v, "myF"=myF, "FTF"=FTF, "myy"=myy, "myX"=myX)
}

