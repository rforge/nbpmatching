# distance matrix method - quantile
# returns quantiles for upper triangular portion of distance matrix
setMethod("quantile", "distancematrix", function(x, probs, ...) {
    if(missing(probs)) probs=seq(from=0, to=100)/100
    quantile(x[upper.tri(x)], probs=probs, ...)
})
