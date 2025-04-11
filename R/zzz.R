#----------------------#
# Internal R functions
#----------------------#

# Convert seconds into a "HH:MM:SS" format
# Author: Gwenael G.R. Leday
.convert2time <- function(x) {
  h <- as.character(x %/% 3600)
  m <- as.character((x %% 3600) %/% 60)
  s <- as.character(round((x %% 3600) %% 60))
  if (nchar(m) == 1) m <- paste(0, m, sep = "")
  if (nchar(s) == 1) s <- paste(0, s, sep = "")
  paste(h, m, s, sep = ":")
}
