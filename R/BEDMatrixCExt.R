numSamples <- function(X) {
    .Call(C_numSamples, X@xptr)
}

numVariants <- function(X) {
    .Call(C_numVariants, X@xptr)
}

extractGenotypes <- function(X, i, j) {
    if (missing(j)) {
        .Call(C_extractGenotypeLinear, X@xptr, as.integer(i))
    } else {
        .Call(C_extractGenotypeCartesian, X@xptr, as.integer(i), as.integer(j))
    }
}
