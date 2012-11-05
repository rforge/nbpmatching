# extra utility functions

# list of names for created/fake elements found in matched dataset
created.names <- c("phantom", "ghost", "chameleon")

# get a factor variable out of a nonbimatch match
setGeneric("get.sets", function(matches, remove.unpaired=TRUE, ...) standardGeneric("get.sets"))
setMethod("get.sets", "data.frame", function(matches, remove.unpaired=TRUE, ...) {
    # thanks to Jake Bowers for providing this function
    sets <- matches[,grep("ID", names(matches))]
    f.sets <- apply(sets, MARGIN=1, FUN=function(x) paste(sort(x), collapse='-'))
    names(f.sets) <- sets[,1]
    if(remove.unpaired) f.sets <- f.sets[grep(paste(created.names, collapse="|"), f.sets, invert=TRUE)]
    factor(f.sets)
})

# calculate scalar distance
setGeneric("scalar.dist", function(x, ...) standardGeneric("scalar.dist"))
setMethod("scalar.dist", "vector", function(x, ...) {
    # thanks to Jake Bowers for providing this function
    if(!is.numeric(x)) stop("x should be numeric")
    outer(x, x, FUN=function(i,j) abs(i-j))
})
