# run quality of match diagnostics
# return data.frame with QoM information

setGeneric("qom", function(covariate, matches, iterations=10000, probs=NA, seed=101, all.vals=FALSE, ...) standardGeneric("qom"))
setMethod("qom", "data.frame", function(covariate, matches, iterations=10000, probs=NA, seed=101, all.vals=FALSE, ...) {
    n <- nrow(matches)
    if(n%%2 == 1) {
        stop("There must be an even number of elements")
    }
    if(!all(sapply(matches[,c(2,4)], is.numeric))) {
        stop("matches must contain numeric values in columns two and four")
    }
    pairs <- matches[matches[,2] < matches[,4], c(2,4)]
    if(all.vals != TRUE) {
        # obtain actual pairs (meaning not matched to phantoms)
        pairs <- pairs[pairs[,2] <= nrow(covariate),]
    }
    npairs <- nrow(pairs)
    ignorecols <- which(sapply(covariate, FUN=function(x) { length(setdiff(x, suppressWarnings(as.numeric(x)))) > 0 }))
    if(length(ignorecols) > 0) covariate <- covariate[, -ignorecols]

    if(is.na(iterations) || !is.numeric(iterations) || iterations < 2) iterations <- 10000
    restoreRandom <- FALSE
    if(exists(".Random.seed")) {
        save.seed <- .Random.seed
        restoreRandom <- TRUE
    }
    if(!is.numeric(seed)) seed <- 101
    set.seed(seed)
    choices <- matrix(sample(c(-1,1), npairs*iterations, replace=TRUE), ncol=iterations)
    if(restoreRandom) .Random.seed <<- save.seed
    else rm(.Random.seed, inherits=TRUE)

    group.one <- sapply(covariate[pairs[,1],], FUN=function(x) {
      colMeans(choices * x, na.rm=TRUE)
    })
    group.two <- sapply(covariate[pairs[,2],], FUN=function(x) {
      colMeans(choices * x, na.rm=TRUE)
    })
    pairdiff.sums.mat <- abs(group.one - group.two)
    if(is.na(probs) || !is.numeric(probs)) probs <- c(0,25,50,75,90,95,100)/100
    pairsumm <- round(t(apply(pairdiff.sums.mat, MARGIN=2, FUN=function(x) { quantile(x, probs=probs) })), 4)
    row.names(pairsumm) <- names(covariate)
    pairsumm
})
