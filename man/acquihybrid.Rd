\name{acquihybrid}
\alias{acquihybrid}
\title{
An hybrid CFA-EFA R package for controlling acquiescence in restricted solutions
}
\description{
Hybrid multi-stage factor analytic procedure for controlling the effects of ACQ in multidimensional questionnaires for which a restricted factorial solution (Ferrando & Lorenzo-Seva, 2000) can be specified. In the first stage, an ACQ factor is estimated and its effects are partialled-out from the inter-item correlations. In the second, a specified CFA solution is fitted to the cleaned matrix. This solution can be specified in two ways: (a) using the cleaned correlation matrix as input, or (b) specifying a full solution which has fixed loadings on the ACQ factor. The use of the program allows (a) the structural properties of the questionnaire to be assessed and (b) unbiased factor score estimates to be obtained for each respondent.
}
\usage{
acquihybrid(x, content_factors, target, corr = "Pearson", raw_data=TRUE,
   method = "fixed", display = TRUE)

}
\arguments{
  \item{x}{
Raw sample item scores or a covariance/correlation matrix, as a data.frame or a numerical matrix.
}
  \item{content_factors}{
The number of content factors to be retained. At least, each factor has to be defined by at least 3 items, so the maximum number of content factors will be the number of items / 3.
}
  \item{target}{
The pattern loading target matrix, which provides the signed dominant loading (higher in absolute values) of each item on its corresponding factor. The target is only used as a reference for determining which items have significant loadings on which factors. The specific loading estimates are not used.
}
  \item{corr}{
Determines the type of matrices to be used in the factor analysis. "Pearson": Computes Pearson correlation matrices (linear model); "Polychoric": Computes Polychoric/Tetrachoric correlation matrices (graded model). If the matrix is not positive definite, the smoothing procedure proposed by Bentler and Yuan (2011) is used.
}
 \item{raw_data}{
Logical argument, if TRUE, the entered data will be treated as raw scores (default). If FALSE, the entered data will be treated as a covariance/correlation matrix.
}
  \item{method}{
Two possibilities are provided: fixed, which use the ACQ loadings obtained in the first step to use it as input for CFA and resid, which use the ACQ free covariance matrix as input for the CFA.
}
  \item{display}{
Determines if the output will be displayed in the console, TRUE by default. If it is TRUE, the output is printed in console and if it is FALSE, the output is returned silently to the output variable.
}
}
\details{
Hybrid multi-stage factor analytic procedure for controlling the effects of ACQ in multidimensional questionnaires for which a restricted factorial solution (Ferrando & Lorenzo-Seva, 2000) can be specified.

The procedure is structurated in several steps:

1. Estimate the ACQ loadings obtaining the first MRFA factor of the inter-item correlation matrix, thereafter, obtaining the centroid of this matrix (using it as target) and rotating the MRFA factor to the position of maximal congruence with respect to target.

2. Obtain the corrected inter-item residual matrix, that is, free of ACQ. This matrix will be trated as a residual covariance matrix.

3. Specify the prescribed CFA solution. There are two alternatives. (a) using the cleaned correlation matrix as input, or (b) specifying a full solution which has fixed loadings on the ACQ factor. This step uses the \code{\link[lavaan]{cfa}} (Rosseel, 2012), with the model obtained in the previous steps.

4. Obtain the goodness fit indices and factor scores. The factor scores are obtained using \code{\link[lavaan]{fitmeasures}} (Rosseel, 2012).

}
\value{
\item{loadings}{Obtained loading matrix.}
\item{fit_indices}{Fit indices: Goodness of Fit Index (GFI), Root Mean Square Residuals (RMSR), Root Mean Square Error of Approximation (RMSEA) and Comparative Fit Index (CFA)}
\item{AC_variance}{The amount of variance explained by ACQ}
\item{pfactors}{factor scores. If corr="Pearson", the method is "regression". In the categorical case, with Polychoric correlation, the option is the empirical bayes Modal approach (EBM)}
}
\references{
Bentler, P. M., & Yuan, K. H. (2011). Positive definiteness via off-diagonal scaling of a symmetric indefinite matrix. Psychometrika, 76(1), 119-123. \doi{10.1007/s11336-010-9191-3}

Ferrando, P. J., & Seva, U. L. (2000). Unrestricted versus restricted factor analysis of multidimensional test items: Some aspects of the problem and some suggestions. Psicológica, 21(2), 301-323.

Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2).\doi{10.18637/jss.v048.i02}
}
\author{
David Navarro-Gonzalez

Pere J. Ferrando

Fabia Morales-Vives

Ana Hernandez-Dorado
}

\examples{
## This is an example using the psymas dataset, which contains 1309 responders to PSYMAS
## questionnaire. For this example, only 10 items are used, measuring two subscales.

psymas_target=cbind(c(-9,-9,0,0,0,9,0,0,9,0),c(0,0,-9,9,-9,0,9,-9,0,9))

acquihybrid(psymas,content_factors=2,psymas_target)


## For speeding purposes, some options have been disabled. Pearson correlation matrix is
## being used, but Polychoric matrix may be advisable.
}
