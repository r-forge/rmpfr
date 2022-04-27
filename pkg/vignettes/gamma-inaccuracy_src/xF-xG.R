## i01 <- sort(outer(c(0:1,3,5,7,9), 10*(0:7), `+`))
## sixt <- (1:15)/16
## iu16 <- 1L+ 3L*(1:15)
## u01 <- i01/80
## if(!identical(sixt, u01[iu16])) {  # unneeded (always ?)
##     cat("fix the i/16 : "); u01[iu16] <- sixt ; cat("[Ok]\n") }
## if(FALSE)## alternative to the "tricky" above: has the same number of non-multiples of 1/16: 11168
##     u01 <- (0:39)/40
## if(FALSE)## even smaller: ==> length(xF) == 6803
    u01 <- (0:19)/20
nF <- -177:171
## For factorial(): {xF = xG-1  but (very slightly) more accurately}
## For gamma():     {xG = xF+1  but (very slightly) more accurately}
xG <- c(outer(u01, nF+1, `+`))
xF <- c(outer(u01, nF,   `+`))
xF <- xF[!(xF <  0 & xF == as.integer(xF))] # dropping -1, -2, ... where gamma() is known to be NaN
xG <- xG[!(xG <= 0 & xG == as.integer(xG))] # dropping 0, -1, -2, ...

