source("https://bioconductor.org/biocLite.R")
biocLite("RNASeqPower")

### R code from vignette source 'samplesize.Rnw'

###################################################
### code chunk number 1: premable
###################################################
options(width=60, continue=" ")
options(SweaveHooks=list(fig=function() par(mar=c(5.1, 4.1, .3, 1.1))))
latex.list <- function(mat, digits=2, label= TRUE, ...) {
  dname <- dimnames(mat)
  if (length(dname[[2]]) & label) 
    cat('&', paste(dname[[2]], collapse=' &'), '\\\\\n')
  
  temp <- format(mat, digits=digits, trim=TRUE, ...)
  temp <- ifelse(temp=='NA', ' ', temp)
  rname <- dname[[1]]
  for (i in 1:nrow(mat)) {
    if (length(rname) && label)
      cat(rname[i,], ' & ', paste(temp[i], collapse=' & '), '\\\\\n')
    else  cat(paste(temp[i,], collapse=' & '), '\\\\\n')
  }
  invisible(NULL)
}


###################################################
### code chunk number 2: samplesize.Rnw:116-119
###################################################
library(RNASeqPower)
rnapower(depth=20, cv=.4, effect=c(1.25, 1.5, 1.75, 2), 
         alpha= .05, power=c(.8, .9))


###################################################
### code chunk number 3: samplesize.Rnw:127-135
###################################################
rnapower(depth=100, cv=.4, effect=c(1.25, 1.5, 1.75, 2), 
         alpha= .05, power=c(.8, .9))         

rnapower(depth=1000, cv=.4, effect=c(1.25, 1.5, 1.75, 2), 
         alpha= .05, power=c(.8, .9))         

rnapower(depth=20, cv=.3, effect=c(1.25, 1.5, 1.75, 2), 
         alpha= .05, power=c(.8, .9))         


###################################################
### code chunk number 4: samplesize.Rnw:148-150
###################################################
rnapower(depth=8, n=10, cv=0.1, effect=c(1.5, 1.75, 2),
         alpha=.05)

