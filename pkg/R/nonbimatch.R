# distance matrix method - nonbimatch
# returns a list of matched pairs
setGeneric("nonbimatch", function(mdm, threshold=NA, precision=6, ...) standardGeneric("nonbimatch"))
setMethod("nonbimatch", "distancematrix", function(mdm, threshold=NA, precision, ...) {
    if(any(is.na(as.integer(mdm)))) {
        stop("Elements of a distance matrix must be integers")
    }
    n <- nrow(mdm)
    if(n%%2 == 1) {
        stop("There must be an even number of elements")
    }
    if(is.null(rownames(mdm))) {
        rownames(mdm) <- seq.int(n)
    }
    ids <- rownames(mdm)
    if(!is.na(threshold) && threshold >= 0) {
        # extend the distance matrix for chameleons
        mdm2 <- data.frame(matrix(threshold, nrow=n*2, ncol=n*2))
        mdm2[seq_len(n), seq_len(n)] <- mdm
        mdm <- mdm2
        rm(mdm2)
        n <- n*2
    } else threshold <- NA

    wt <- as.vector(t(mdm))
    nmatch <- seq.int(n)
    numdigits <- floor(log10(max(wt))) + 1
    if(!is.numeric(precision) || precision < 1) {
        precision <- 6
        warning("Precision value is too small.  Setting precision to six.")
    }
    shift <- 10^(precision-numdigits)
    # the largest number will have at most [precision] digits (defaulting to six)
    # warning: if the vector is large and there are too many digits, the Fortran call will crash
    if(shift < 1) {
        wt <- wt*shift
        print(sprintf("Note: Distances scaled by %s to ensure all data can be handled", shift))
    }
    match <- .Fortran('mwrap', PACKAGE="nbpMatching", n=as.integer(n), wt=as.integer(wt), nmatch=as.vector(nmatch), prcn=precision)$nmatch

    # remove chameleon to chameleon matches
    if(!is.na(threshold)) {
        c2c <- which(match[seq(from=(n/2)+1, to=n)] > n/2)
        # number of chameleons that match elements
        ncham <- n/2 - length(c2c)
        if(length(c2c) > 0L) {
            match <- match[-(c2c + n/2)]
            index.to.replace <- which(match > n/2)
            # re-label chameleon matches from 1 to ncham
            if(length(index.to.replace) > 0L) {
                match[index.to.replace] <- seq(from=(n/2+1), to=(n/2+ncham))
                match[seq(from=(n/2+1), to=(n/2+ncham))] <- index.to.replace
            }
            n <- length(match)
            ids <- append(ids, sprintf('chameleon%s', seq_len(ncham)))
        }
    }

    result <- data.frame(match)
    matches <- data.frame("Group1.ID"=NA, "Group1.Row"=numeric(n), "Group2.ID"=NA, "Group2.Row"=numeric(n), "Distance"=numeric(n))

    for(i in seq.int(n)) {
        matches[i,c(1,3)] <- c(ids[i], ids[result[i, 'match']])
        matches[i,c(2,4,5)] <- c(i, result[i, 'match'], mdm[i, result[i,1]])
    }
    halves <- matches[matches[,2] < matches[,4],]
    distance <- sum(matches[,5])
    total <- distance/2
    mean <- distance/n
    list(matches=matches, halves=halves, total=total, mean=mean)
})

# runner starts with covariate matrix, generates distances, finds matches, and reports QoM
setGeneric("runner", function(covariate, seed=101, ..., mate.random=FALSE) standardGeneric("runner"))
setMethod("runner", "data.frame", function(covariate, seed=101, ..., mate.random=FALSE) {
    step1 <- gendistance(covariate, ...)
    step2 <- distancematrix(step1, ...)
    if(mate.random) {
        restoreRandom <- FALSE
        if(exists(".Random.seed")) {
            save.seed <- .Random.seed
            restoreRandom <- TRUE
        }
        if(!is.numeric(seed)) seed <- 101
        set.seed(seed)
        n <- nrow(step2)
        options <- seq_len(n)
        result <- numeric(n)
        for(i in seq_len(n/2)) {
            mates <- sample(options, 2)
            result[mates[1]] <- mates[2]
            result[mates[2]] <- mates[1]
            options <- options[-match(mates, options)]
        }
        result <- data.frame(match=result)
        if(restoreRandom) set.seed(save.seed)
        else rm(.Random.seed, inherits=TRUE)

        ids <- rownames(step2)
        matches <- data.frame("Group1.ID"=NA, "Group1.Row"=numeric(n), "Group2.ID"=NA, "Group2.Row"=numeric(n), "Distance"=numeric(n))
        for(i in seq.int(n)) {
            matches[i,c(1,3)] <- c(ids[i], ids[result[i, 'match']])
            matches[i,c(2,4,5)] <- c(i, result[i, 'match'], step2[i, result[i,1]])
        }
        halves <- matches[matches[,2] < matches[,4],]
        distance <- sum(matches[,5])
        total <- distance/2
        mean <- distance/n
        step3 <- list(matches=matches, halves=halves, total=total, mean=mean)
    } else {
        step3 <- nonbimatch(step2, ...)
    }
    step4 <- qom(step1$cov, step3$matches, ...)
    step5 <- assign.grp(step3$matches, ...)
    step6 <- cbind(step1$cov, step5[seq_len(nrow(step1$cov)),c(3,4,6)])
    return(list(setup=step1, mdm=step2, matches=step3, qom=step4, grps=step5, final=step6))
})

setGeneric("full.qom", function(covariate, iterations=NA, ...) standardGeneric("full.qom"))
setMethod("full.qom", "data.frame", function(covariate, iterations=NA, ...) {
    qom <- vector('list', ncol(covariate)+3)
    if(is.na(iterations) || !is.numeric(iterations) || iterations < 2) iterations <- 10000
    a <- runner(covariate, iterations=iterations, ...)
    ids <- rownames(a$setup$cov)
    cnames <- colnames(covariate)
    qom[[1]] <- a$qom
    qom[[2]] <- runner(covariate, iterations=iterations, mate.random=TRUE)$qom
    qom[[3]] <- runner(covariate, iterations=iterations, missingness=0)$qom
    for(i in seq_len(ncol(covariate))) {
        weights <- rep(0, ncol(covariate))
        weights[i] <- 1
        qom[[3+i]] <- tryCatch(runner(covariate, iterations=iterations, weights=weights, missingness=0)$qom, error=function(e) {})
    }
    names(qom) <- c("user.specified", "random.mates", "eq.weight", cnames)
    qom
})
