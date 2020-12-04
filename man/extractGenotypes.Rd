\name{extractGenotypes}
\alias{extractGenotypes}
\title{Extract genotypes}
\description{
    This function extracts genotypes from a \code{BEDMatrix} object at the C
    level.
}
\usage{
extractGenotypes(X, i, j)
}
\arguments{
    \item{X}{
        A \code{BEDMatrix} instance.
    }
    \item{i}{
        The row indices to extract.
    }
    \item{j}{
        The column indices to extract.
    }
}
\value{
    A one-dimensional or two-dimensional integer vector of genotypes.
}
