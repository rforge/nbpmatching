# distance matrix method - nonbimatch
# returns a list of matched pairs
setGeneric("nonbimatch", function(mdm, precision=6, ...) standardGeneric("nonbimatch"))
setMethod("nonbimatch", "distancematrix", function(mdm, precision, ...) {
    if(any(is.na(as.integer(mdm)))) {
        stop("Elements of a distance matrix must be integers")
    }
    wt<-as.vector(t(mdm))
    n<-nrow(mdm)
    if(n%%2 == 1) {
        stop("There must be an even number of elements")
    }
    nmatch<-1:n
    numdigits<-floor(log10(max(wt)))+1
    if(!is.numeric(precision) || precision < 1) {
        precision<-6
        warning("Precision value is too small.  Setting precision to six.")
    }
    shift<-10^(precision-numdigits)
    # the largest number will have at most [precision] digits (defaulting to six)
    # warning: if the vector is large and there are too many digits, the Fortran call will crash
    if(shift < 1) {
        wt<-wt*shift
        warning(paste("Values were too large!  Multiplied by ", shift, "to ensure all data can be handled"))
    }
    matched<-.Fortran('mwrap', PACKAGE="nbpMatching", n=as.integer(n), wt=as.integer(wt), nmatch=as.vector(nmatch), prcn=precision)$nmatch
    result<-data.frame(matched)
    matches<-data.frame("X"=numeric(n), "Xrow"=numeric(n), "Y"=numeric(n), "Yrow"=numeric(n), "Distance"=numeric(n))
    halves<-data.frame("Group1"=numeric(n/2), "Group1row"=numeric(n/2), "Group2"=numeric(n/2), "Group2row"=numeric(n/2), "Distance"=numeric(n/2))
    if(is.null(rownames(mdm))) {
        matches<-matches[,c(-1,-3)]
        halves<-halves[,c(-1,-3)]
    }
    count<-1
    for(i in 1:n) {
        dist<-mdm[i,result[i,1]]
        result[i,'dist']<-dist
        if(is.null(rownames(mdm))) {
            matches[i,]<-c(i, result[i, 'matched'], dist)
        } else {
            matches[i,]<-c(rownames(mdm)[i], i, rownames(mdm)[result[i, 'matched']], result[i, 'matched'], dist)
        }
        if(result[i,'matched'] > i) {
            halves[count,]<-matches[i,]
            count<-count+1
        }
    }
    distance<-sum(result[,2])
    total<-distance/2
    mean<-distance/n
    list(matches=matches, halves=halves, total=total, mean=mean)
})
