\name{kPPR}
\alias{kPPR}
\title{Projection pursuit regression with local linear kernel smoother}
\description{
  Generates a projection pursuit regression model for covariates X and response y
}
\usage{
kPPR(X, y, nterms=1, hmult=1, betas=NULL, loss=NULL, dloss=NULL, initialisation="fancy")
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions) of covariates.}
  \item{y}{a numeric vector of responses.}
  \item{nterms}{The number of terms in the additive regression model. The default is a single term.}
  \item{betas}{The coefficients in the expression of the kernel. See Hofmeyr (2019) for details. The
  default is the smooth order one kernel described in the paper.}
  \item{hmult}{The bandwidth in the kernel smoother is set to sd(X%*%v0)/nrow(X)^.2*hmult, where v0
  is the initial projection given to the optimiser.}
  \item{loss}{The (additive) loss function to be used. Allows for generalised responses by specifying an appropriate likelihood/deviance function for the loss. Note: loss is to be minimised, so deviance or negative log-likelihood would be appropriate. If specified then must be a function of y and hy (the fitted values, yhat), returning a vector of "errors"" which are added as the total loss. The default is the squared error.}
  \item{dloss}{The derivative of the loss function. Also takes arguments y and hy, and returns the vector of partial derivatives of the loss w.r.t. hy.}
  \item{initialisation}{Method use to initialise the projection vectors. Must be one of "fancy", "lm" and "random". The default is "fancy", which fits a linear model to the (scaled) univariate predictions of y on each X[,i]. "lm" is just a simple set of linear model coefficients. In both cases the linear model uses a LASSO penalty with very small "lambda", to ensure a solution is found. "random" uses random initialisation.}
}
\value{
  A named list with class kPPR containing
  \item{$mu}{The estimated (global) mean of the response.}
  \item{$mu_X}{The vector of estimated means of the covariates.}
  \item{$y}{The responses given as argument to the function.}
  \item{$X}{The covariates given as argument to the function.}
  \item{$hs}{A vector of bandwidths used for each term in the model.}
  \item{$vs}{A matrix whose rows are the projection vectors.}
  \item{$fitted}{The fitted values on each projection. Note that these are based on the residuals used to fit that component of the model, and not the original y values. $fitted is used for prediction on test data.}
  \item{$bets}{The beta coefficients in the kernel formulation.}
}
\references{
  Friedman, J., and Stuetzle, W. (1981) "Projection pursuit regression." \emph{Journal of the American statistical Association} \bold{76}.376.

Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
set.seed(2)

### Generate a set of data

X = matrix(rnorm(10000), ncol = 10)

### Generate some "true" projection vectors

beta1 = (runif(10)>.5)*rnorm(10)
beta2 = (runif(10)>.5)*rnorm(10)

### Generate responses, dependent on X\%*\%beta1 and X\%*\%beta2

y = 1 + X\%*\%beta1 + ((X\%*\%beta2)>2)*(X\%*\%beta2-2)*10 + (X\%*\%beta1)*(X\%*\%beta2)/5 + rnorm(1000)

### Fit a PPR model with 2 terms on a sample of the data

train_ids = sample(1:1000, 500)

model = kPPR(X[train_ids,], y[train_ids], nterms = 2)

### Predict on left out data, and compute estimated coefficient of determination

yhat = predict(model, X[-train_ids,])

1-sum((yhat-y[-train_ids])^2)/sum((y[-train_ids]-mean(y[-train_ids]))^2)


#################################################

### Add some residuals in the training data and fit the model again, as well as one with an absolute loss

y[train_ids] = y[train_ids] + (runif(500)<.05)*(rnorm(500)*100)

model1 <- kPPR(X[train_ids,], y[train_ids], nterms = 2)

model2 <- kPPR(X[train_ids,], y[train_ids], nterms = 2, loss = function(y, hy) abs(y-hy), dloss = function(y, hy) sign(hy-y))

### Plot the resulting components in the model on the test data

par(mar = c(2, 2, 2, 2))
par(mfrow = c(2, 2))

plot(X[-train_ids,]\%*\%model1$vs[1,], y[-train_ids])
plot(X[-train_ids,]\%*\%model1$vs[2,], y[-train_ids])
plot(X[-train_ids,]\%*\%model2$vs[1,], y[-train_ids])
plot(X[-train_ids,]\%*\%model2$vs[2,], y[-train_ids])

### estimate comparative estimated coefficients of determination

1-sum((predict(model1, X[-train_ids,])-y[-train_ids])^2)/sum((y[-train_ids]-mean(y[-train_ids]))^2)
1-sum((predict(model2, X[-train_ids,])-y[-train_ids])^2)/sum((y[-train_ids]-mean(y[-train_ids]))^2)
}
\keyword{file}
