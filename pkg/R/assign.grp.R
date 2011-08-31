# assign pairs to groups
# return original data.frame with added group column

setGeneric("assign.grp", function(matches, seed=68, ...) standardGeneric("assign.grp"))
setMethod("assign.grp", "data.frame", function(matches, seed=68, ...) {
    n <- nrow(matches)
    if(n%%2 == 1) {
        stop("There must be an even number of elements")
    }
    if(!all(sapply(matches[,c(2,4)], is.numeric))) {
        stop("matches must contain numeric values in columns two and four")
    }
    matches <- cbind(matches, treatment.grp=NA)
    pairs <- matches[matches[,2] < matches[,4], c(2,4)]
    restoreRandom <- FALSE
    if(exists(".Random.seed")) {
        save.seed <- .Random.seed
        restoreRandom <- TRUE
    }
    if(!is.numeric(seed)) seed <- 68
    set.seed(seed)
    choices <- sample(c(TRUE, FALSE), nrow(pairs), replace=TRUE)
    if(restoreRandom) set.seed(save.seed)
    else rm(.Random.seed, inherits=TRUE)

    for(i in seq_len(nrow(pairs))) {
        matches$treatment.grp[pairs[i, 1]] <- ifelse(choices[i], "A", "B")
        matches$treatment.grp[pairs[i, 2]] <- ifelse(choices[i], "B", "A")
    }
    matches
})
