.First.lib <- function(lib, pkg) {
    library.dynam("Rmpfr",pkg,lib)
    if(FALSE) ## does not make sense: is *not* in package subsequently;  need assign() / Namespace ??
    Pi <- Const("pi")
    ## for *stable() -- for the moment -- but this *still* fails here.. methods are not yet "active"
    ## pi2 <- Const("pi", 256) / 2
}
