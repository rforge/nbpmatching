# add phantom values to a matrix
# return matrix with nphantoms more rows and nphantoms more columns

setGeneric("make.phantoms", function(x, nphantoms, name="phantom", maxval=2*10^8, ...) standardGeneric("make.phantoms"))
setMethod("make.phantoms", signature(x="matrix", nphantoms="numeric"), function(x, nphantoms, name="phantom", maxval=2*10^8, ...) {
    if(nphantoms < 1) return(x)
    if(missing(name)) {
        name <- "phantom"
    } else if(!is.character(name)) {
        stop("name argument is not character")
    }
    if(missing(maxval)) {
        maxval <- 2*10^8
    } else if(!is.numeric(maxval)) {
        stop("maxval argument is not numeric")
    }
    nr <- nrow(x)
    # preserve rownames
    mynames <- rownames(x)
    if(is.null(mynames)) {
        mynames <- seq(nr)
    }
    # phantom index
    p.index <- seq(from=nr+1, length.out=nphantoms)
    newvals <- rep(0, nphantoms)
    # add phantom columns
    m <- do.call("cbind", c(list(x), newvals))
    # add phantom rows
    m <- do.call("rbind", c(list(m), newvals))
    # distance between phantoms should be maxval
    m[p.index, p.index] <- maxval
    # create row names
    mynames <- c(mynames, sprintf("%s%s", name, p.index))
    dimnames(m) <- list(mynames, mynames)
    m
})
# x is data.frame instead of matrix
setMethod("make.phantoms", signature(x="data.frame", nphantoms="numeric"), function(x, nphantoms, name, maxval, ...) {
    # convert to matrix and back
    as.data.frame(make.phantoms(as.matrix(x), nphantoms, name, maxval, ...))
})
# don't do anything when nphantoms is missing
setMethod("make.phantoms", signature(nphantoms="missing"), function(x, nphantoms, name, maxval, ...) {
    return(x)
})
