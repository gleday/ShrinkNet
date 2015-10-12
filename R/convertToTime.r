# Internal function
# Convert seconds into a "HH:MM:SS" format

# Author: Gwenael G.R. Leday

.convertToTime <- function(x){
	h <- as.character(x%/%3600)
	m <- as.character((x%%3600)%/%60)
	s <- as.character(round((x%%3600)%%60))
	if(nchar(m)==1) m <- paste(0,m,sep="")
	if(nchar(s)==1) s <- paste(0,s,sep="")
	return(paste(h,m,s,sep=":"))
}