# distance matrix S4 object
# a NxN matrix whose elements are Integers
setClass("distancematrix", contains="matrix")
setMethod("initialize", "distancematrix", function(.Object, ...) {
    .Object<-callNextMethod()
    nr<-nrow(.Object)
    nc<-ncol(.Object)
    if(abs(nr-nc) > 1)
        stop("Row and column lengths must be equal")
    if(nr != nc) {
        if(nr > nc) {
            mynames<-.Object[1,]
            .Object@.Data<-.Object@.Data[-1,]
            nr<-nr-1
        } else {
            mynames<-.Object[,1]
            .Object@.Data<-.Object@.Data[,-1]
            nc<-nc-1
        }
        colnames(.Object)<-mynames
        rownames(.Object)<-mynames
    }
    if(any(is.na(as.integer(.Object)))) {
        stop("Elements of a distance matrix must be integers")
    }
    if(nr%%2 == 1) {
        warning("There must be an even number of elements\nAdding a ghost value")
        .Object@.Data<-rbind(.Object@.Data, rep(0, nc))
        nr<-nr+1
        .Object@.Data<-cbind(.Object@.Data, rep(0, nr))
        nc<-nc+1
        colnames(.Object)[nr]<-'ghost'
        rownames(.Object)[nc]<-'ghost'
    }
    .Object@.Data<-matrix(as.integer(.Object@.Data), nrow=nr)
    # set the diagonal to zero
    diag(.Object@.Data)<-0
    if(any(.Object@.Data != t(.Object@.Data))) {
        stop("A distancematrix must be symmetric")
    }
    .Object
})
setGeneric("distancematrix", function(x, ...) standardGeneric("distancematrix"))
setMethod("distancematrix", "matrix", function(x, ...) {
    new('distancematrix', x, ...)
})
setMethod("distancematrix", "character", function(x, ...) {
    if(file.access(x) == -1) {
        stop("File not accessible")
    }
    if(length(grep("[.]csv$", x, ignore.case=TRUE)) != 1) {
        stop("Function expects a CSV file")
    }
    distancematrix(read.csv(x, ...))
})
setMethod("distancematrix", "data.frame", function(x, ...) {
    distancematrix(as.matrix(x), ...)
})
setMethod("distancematrix", "list", function(x, ...) {
    if("dist" %in% names(x)) {
        distancematrix(x$dist, ...)
    } else {
        stop("x should contain dist element")
    }
})
setMethod("[<-", "distancematrix", function(x, i, j, value) { stop("You may not re-assign elements of a distancematrix") })
setMethod("[[<-", "distancematrix", function(x, i, j, value) { stop("You may not re-assign elements of a distancematrix") })
