\name{kPPR-package}
\alias{kPPR-package}
\docType{package}
\title{
  Projection pursuit regression with local linear kernel smoother
}
\description{
  Obtains additive regression models from a basis consisting of ridge functions (i.e., those which only vary along one direction).
Individual functions are estimated using local linear kernel smoother.
}
\details{
  The function kPPR(X, y) fits a projection pursuit regression (ppr) model of y on X. For more details on the function, use help(kPPR). See Friedman (1984)
  for a description of the ppr model. The present implementation uses local linear kernel smoothing
  using the fast kernel method of Hofmeyr (2019).
}
\references{
Friedman, J., and Stuetzle, W. (1981) "Projection pursuit regression." \emph{Journal of the American statistical Association} \bold{76}.376.

Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\author{
David P. Hofmeyr
}
