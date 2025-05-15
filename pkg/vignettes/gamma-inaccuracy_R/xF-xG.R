u01 <- (0:19)/20
nF <- -177:171
## For factorial(): {xF = xG-1  but (very slightly) more accurately}
## For gamma():     {xG = xF+1  but (very slightly) more accurately}
xG <- c(outer(u01, nF+1, `+`))
xF <- c(outer(u01, nF,   `+`))
xF <- xF[!(xF <  0 & xF == as.integer(xF))] # dropping -1, -2, ... where gamma() is known to be NaN
xG <- xG[!(xG <= 0 & xG == as.integer(xG))] # dropping 0, -1, -2, ...

