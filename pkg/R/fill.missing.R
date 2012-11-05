# data.frame method - fill.missing
# returns data with NAs replaced with imputed values
# appends "missing" columns with missingness indicator
setGeneric("fill.missing", function(x, seed=101, simplify=TRUE, idcol="id", ...) standardGeneric("fill.missing"))
setMethod("fill.missing", "data.frame", function(x, seed=101, simplify=TRUE, idcol="id", ...) {
    if(exists(".Random.seed", envir = .GlobalEnv)) {
        save.seed <- get(".Random.seed", envir= .GlobalEnv)
        on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    miss.cols <- which(apply(x, MARGIN=2, FUN=function(y) { xm<-is.na(y);(any(xm) & !all(xm))}))
    if(length(miss.cols) > 0) {
        missingness <- is.na(x[,miss.cols, drop=FALSE])+0
        dup <- duplicated(missingness, MARGIN=2)
        # simplify will remove duplicate "missing" columns
        if(simplify == TRUE && any(dup)) {
            missingness <- missingness[,!dup, drop=FALSE]
        }
        colnames(missingness) <- sprintf("%s.missing", colnames(missingness))
    } else {
      # no missing data
      return(x)
    }
    if(!is.numeric(seed)) seed <- 101
    set.seed(seed)
    cnames <- names(x)
    if(!is.null(idcol) && length(idcol) == 1) {
        if(is.character(idcol)) {
            idcol <- match(idcol, cnames)
        }
        if(!is.na(idcol) && is.numeric(idcol) && idcol > 0 && idcol <= ncol(x)) {
            cnames <- cnames[-idcol]
        }
    }
    fmla <- as.formula(paste("~", paste(cnames, collapse= "+")))
    # run until no errors
    fails <- 0
    while(fails < 10) {
        if(!is.null(tryCatch(imputer <- transcan(fmla, data=x, transformed=FALSE, pl=FALSE, pr=FALSE, imputed=TRUE, ...),  error = function(e){}))) break
        fails <- fails+1
    }
    if(fails >= 10) stop("cannot impute values")
    invisible(lapply(names(imputer$imputed), FUN=function(cname) {
        if(length(imputer$imputed[[cname]]) > 0L) {
            imputed.vals <- imputer$imputed[[cname]]
            x[names(imputed.vals), cname] <<- imputed.vals
        }
    }))
    cbind(x, missingness)
})
