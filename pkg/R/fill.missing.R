# data.frame method - fill.missing
# returns data with NAs replaced with imputed values
# appends "missing" columns with missingness indicator
setGeneric("fill.missing", function(x, seed=101, ...) standardGeneric("fill.missing"))
setMethod("fill.missing", "data.frame", function(x, seed=101, ...) {
    if(exists(".Random.seed", envir = .GlobalEnv)) {
        save.seed <- get(".Random.seed", envir= .GlobalEnv)
        on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    miss.cols <- which(apply(x, MARGIN=2, FUN=function(y) { xm<-is.na(y);(any(xm) & !all(xm))}))
    if(length(miss.cols) > 0) {
        missingness <- data.frame(is.na(x[,miss.cols])+0)
        names(missingness) <- sprintf("%s.missing", names(x)[miss.cols])
    }
    if(!is.numeric(seed)) seed <- 101
    set.seed(seed)
    fmla <- as.formula(paste("~", paste(names(x), collapse= "+")))
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
