stopifnot(require("Rmpfr"))

x <- mpfr(0:7, 64)/7
mx <- x
dim(mx) <- c(4,2)
mx # "print"
m.x <- matrix((0:7)/7, 4,2)

y <- 7 * mpfr(1:12, 80)
my <- y
dim(my) <- 3:4
m.y <- matrix(7 * 1:12, 3,4)

stopifnot(my[2,2] == 35,
          my[,1] == 7*(1:3))

noDN <- function(.) { dimnames(.) <- NULL ; . }
allEQ <- function(x,y) all.equal(x,y, tol=1e-15)

stopifnot(allEQ(m.x, noDN(as(mx, "matrix"))),
          allEQ(m.y, noDN(as(my, "matrix"))),
          allEQ(noDN(as(my %*% mx,"matrix")), m.y %*% m.x),
          allEQ(noDN(as(crossprod(mx, t(my)),"matrix")), crossprod(m.x, t(m.y))),
          allEQ(noDN(as(tcrossprod(my, t(mx)),"matrix")),
                        tcrossprod(m.y, t(m.x))),

          identical(mx, t(t(mx))),
          identical(my, t(t(my))),

          identical(noDN(as(my %*% 1:4,"matrix")),
                         as(my,"matrix") %*% 1:4 )
          )

mx[3,1] <- Const("pi", 64)
stopifnot(allEQ(sum(mx[,1]), pi + 4/7))
m2 <- mx[c(1,4),]
stopifnot(dim(m2) == c(2,2), sum(m2) == 2)

## "mpfrArray" o "mpfr" :
Tmx <- array(TRUE, dim(mx), dimnames=dimnames(mx))
stopifnot(identical(Tmx, mx == (mx - mpfr(0, 10))),
	  identical(Tmx, mx - mpfr(1, 10) * mx == 0))

## %*% with vectors on LHS, ...
y <- t(2:4) # 1 x 3 matrix
m1 <-     (0:10)     %*% y
m2 <- mpfr(0:10, 50) %*% y
stopifnot((d <- m1 - m2) == 0, identical(dim(m1), dim(d)),
          m2 == m1, m1 == m2)

r <- 10*(0:4)
y <- t(2:6)
m1 <- 1:3 %*% y  %*% r
y. <- t(mpfr(2:6, 20))
m2 <- 1:3 %*% y. %*% r
stopifnot(m1 == m2, m1 - m2 == 0, identical(dim(m1), dim(m2)))

cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
