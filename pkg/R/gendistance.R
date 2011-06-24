# prepare a data.frame for distance matrix
# return a list of eight
# 1) distance matrix
# 2) covariates used to create distances
# 3) ignored covariates
# 4) column weights applied to create distances
# 5) column names used to prevent matches
# 6) column indeces, requested to use rank of values
# 7) missingness weight, applied to missingness indicator columns
# 8) number of phantoms created
# missing values are imputed, which creates new columns for missingness

setGeneric("gendistance", function(covariate, idcol=NULL, weights=NULL, prevent=NULL, rankcols=NULL, missing.weight=0.1, ndiscard=0, ...) standardGeneric("gendistance"))
setMethod("gendistance", "data.frame", function(covariate, idcol=NULL, weights=NULL, prevent=NULL, rankcols=NULL, missing.weight=0.1, ndiscard=0, ...) {
    nr<-nrow(covariate)
    nc<-ncol(covariate)
    myrownames <- character(0)
    mycolnames <- names(covariate)
    bad.data <- NULL

    # columns that aren't numeric should be marked bad
    badcol <- which(sapply(1:nc, FUN=function(x) suppressWarnings(!is.numeric(covariate[,x]))))
    # if all values in a column are the same, mark as bad column
    badcol <- union(badcol, which(sapply(1:nc, FUN=function(x) length(unique(covariate[,x])) == 1L)))

    if(length(idcol) == 1L && is.numeric(idcol) && idcol > 0 && idcol <= nc) {
        myrownames <- as.character(covariate[,idcol])
        row.names(covariate) <- myrownames
        badcol <- union(badcol, idcol)
    }

    if(is.null(weights)) {
        weights <- rep(1, nc)
    } else {
        weights <- as.numeric(weights)
        weights <- ifelse(is.na(weights) | weights < 0, 0, weights)
        weights <- c(weights, numeric(nc-length(weights)))
    }

    if(!is.null(prevent)) {
        prevent <- as.integer(prevent)
        prevent <- setdiff(na.omit(ifelse(prevent < 1 | prevent > nc, NA, prevent)), badcol)
        if(length(prevent) >= 1L) {
            badcol <- union(badcol, prevent)
            # convert column index to column name
            prevent <- mycolnames[prevent]
        }
    }
    badcol <- badcol[order(badcol)]

    # remove bad columns
    if(length(badcol) >= 1L) {
        weights <- weights[-badcol]
        bad.data <- covariate[,badcol]
        covariate <- covariate[,-badcol]
    }
    # validate rankcols
    if(!is.null(rankcols)) {
        rankcols <- as.integer(rankcols)
        rankcols <- setdiff(na.omit(ifelse(rankcols < 1 | rankcols > nc, NA, rankcols)), badcol)
        # bad columns have been removed, so if there are any rank columns, their column index has changed
        if(length(rankcols) >= 1L) {
            rankcols <- sapply(rankcols, FUN=function(y) { y - sum((y > badcol)*1) })
        }
    }
    # validate missing.weight and ndiscard
    if(!is.numeric(missing.weight) || missing.weight < 0) missing.weight <- 0.1
    if(!is.numeric(ndiscard) || ndiscard < 0 || (nr - ndiscard) < 2) ndiscard <- 0

    # impute any missing values
    if(any(is.na(covariate))) {
        orig.colnames <- names(covariate)
        covariate <- fill.missing(covariate)
        new.colnames <- names(covariate)
        # calculate new weights
        weight.lookups <- sapply(sub(".missing", "", new.colnames[(length(orig.colnames)+1):length(new.colnames)]), FUN=function(y) { which(orig.colnames == y)  })
        weights <- append(weights, weights[weight.lookups]*missing.weight)
        if(length(rankcols) >= 1L) {
            rankcols <- append(rankcols, (length(orig.colnames)+1):length(new.colnames))
        } else {
            rankcols <- (length(orig.colnames)+1):length(new.colnames)
        }
    }

    # Define your matrix of covariates covariate
    X <- as.matrix(covariate)
    if(length(rankcols) >= 1L) {
        for(i in rankcols) {
            X[,i] <- rank(covariate[,i])
        }
    }

    # Create the covariance matrix and invert it
    X.cov <- cov.wt(X)
    # use pseudo-inverse if matrix is singular
    Sinv <- tryCatch(solve(X.cov$cov), error=function(e) {})
    if(is.null(Sinv)) Sinv <- ginv(X.cov$cov)
    # prevent negative distances
    for(i in seq_len(nrow(Sinv))) {
        Sinv[i,] <- Sinv[i,]*weights[i]
        Sinv[,i] <- Sinv[,i]*weights[i]
    }

    # Define a function to create the distance matrix using Sinv and X
    mdistmaker <- function(row1, row2) { t(X[row1,]-X[row2,]) %*% Sinv %*% (X[row1,]-X[row2,]) }
    # Create the distance matrix mdists
    mdists <- sapply(seq_len(nr), FUN=function(x) mapply(mdistmaker, x, seq_len(nr)))

    # pick a big value for points that shouldn't match - like the diagonal
    numdigits<-floor(log10(max(mdists))) + 1
    shift<-10^(8-numdigits)
    maxval<-2*10^8
    mdists<-floor(mdists*shift)
    diag(mdists) <- maxval

    # penalize "prevent" columns with matching values by setting distance to maxval
    # prevent is a vector of column names found in bad.data
    if(length(prevent) >= 1L) {
        indeces <- match(prevent, names(bad.data))
        for(colnum in indeces) {
            for(rownum in seq_len(nr)) {
                # get the row number of all rows that have the same value in the given column
                eqrows <- setdiff(which(bad.data[rownum,colnum] == bad.data[,colnum]), colnum)
                if(length(eqrows) > 0L) {
                    mdists[rownum,eqrows] <- maxval
                    mdists[eqrows,rownum] <- maxval
                }
            }
        }
    }

    # need to add phantom rows/columns of max distance
    GROUPS <- 2
    nphantoms <- ndiscard + ((GROUPS - (nr - ndiscard) %% GROUPS) %% GROUPS)
    if(nphantoms) {
        m2 <- data.frame(matrix(0, nrow=nr+nphantoms, ncol=nr+nphantoms))
        m2[seq_len(nr), seq_len(nr)] <- mdists
        mdists <- m2
        names(mdists)[(nr+1):nrow(mdists)] <- sprintf("phantom%s", seq(nr+1, nrow(mdists)))
        row.names(mdists)[(nr+1):nrow(mdists)] <- sprintf("phantom%s", seq(nr+1, nrow(mdists)))
    } else {
        # convert matrix to data frame
        mdists <- data.frame(mdists)
    }
    # add back row names
    if(length(myrownames) > 0L) {
        names(mdists)[seq_len(nr)] <- myrownames
        row.names(mdists)[seq_len(nr)] <- myrownames
    } else {
        names(mdists)[seq_len(nr)] <- seq_len(nr)
        row.names(mdists)[seq_len(nr)] <- seq_len(nr)
    }

    list(dist=mdists, cov=covariate, ignored=bad.data, weights=weights, prevent=prevent, rankcols=rankcols, missing.weight=missing.weight, ndiscard=nphantoms)
})
