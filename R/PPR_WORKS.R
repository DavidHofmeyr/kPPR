f_ppr <- function(v, X, y, h, betas, loss = NULL, dloss = NULL){
  p <- X%*%v/sqrt(sum(v^2))
  n <- length(p)
  o <- order(p)
  sK <- ksum(p[o], numeric(n)+1, p[o], h, n, n, length(betas)-1, betas)-betas[1]
  sKx <- ksum(p[o], p[o], p[o], h, n, n, length(betas)-1, betas)-p[o]*betas[1]
  sKy <- ksum(p[o], y[o], p[o], h, n, n, length(betas)-1, betas)-y[o]*betas[1]
  sKx2 <- ksum(p[o], p[o]^2, p[o], h, n, n, length(betas)-1, betas)-p[o]^2*betas[1]
  sKxy <- ksum(p[o], p[o]*y[o], p[o], h, n, n, length(betas)-1, betas)-p[o]*y[o]*betas[1]
  hy <- ((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*p[o])/(sK*sKx2-sKx^2)
  if(is.null(loss)) sum((y[o]-hy)^2)
  else sum(loss(y[o]/100, hy/100)*100)
}


kLLreg <- function(x, y, h, betas){
  n <- length(x)
  o <- order(x)
  sK <- ksum(x[o], numeric(n)+1, x[o], h, n, n, length(betas)-1, betas)-betas[1]
  sKx <- ksum(x[o], x[o], x[o], h, n, n, length(betas)-1, betas)-x[o]*betas[1]
  sKy <- ksum(x[o], y[o], x[o], h, n, n, length(betas)-1, betas)-y[o]*betas[1]
  sKx2 <- ksum(x[o], x[o]^2, x[o], h, n, n, length(betas)-1, betas)-x[o]^2*betas[1]
  sKxy <- ksum(x[o], x[o]*y[o], x[o], h, n, n, length(betas)-1, betas)-x[o]*y[o]*betas[1]
  (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*x[o])/(sK*sKx2-sKx^2))[rank(x)]
}

df_ppr <- function(v, X, y, h, betas, loss = NULL, dloss = NULL){
  nv <- sqrt(sum(v^2))
  p <- X%*%v/nv
  n <- length(p)
  o <- order(p)
  po <- p[o]
  yo <- y[o]
  S1 <- kndksum(po, numeric(n)+1, po, h, n, n, length(betas)-1, betas)-rep(c(betas[1], 0), each = n)
  Sx <- kndksum(po, po, po, h, n, n, length(betas)-1, betas)-cbind(po*betas[1], 0)
  Sy <- kndksum(po, yo, po, h, n, n, length(betas)-1, betas)-cbind(yo*betas[1], 0)
  Sxx <- kndksum(po, po^2, po, h, n, n, length(betas)-1, betas)-cbind(po^2*betas[1], 0)
  Sxy <- kndksum(po, po*yo, po, h, n, n, length(betas)-1, betas)-cbind(po*yo*betas[1], 0)
  N <- ((Sxx[,1]*Sy[,1]-Sx[,1]*Sxy[,1])+(S1[,1]*Sxy[,1]-Sx[,1]*Sy[,1])*po)
  D <- (S1[,1]*Sxx[,1]-Sx[,1]^2)
  hy <- N/D
  if(is.null(loss)){
    dL <- 2*(hy-yo)
  }
  else dL <- dloss(yo, hy)
  dLD <- dL/D
  S_SydLD <- kndksum(po, Sy[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(Sy[,1]*dLD, 0)
  S_xSydLD <- kndksum(po, Sy[,1]*po*dLD, po, h, n, n, length(betas)-1, betas)-cbind(Sy[,1]*po*dLD, 0)
  S_SxxmxSxdLD <- kndksum(po, (Sxx[,1]-po*Sx[,1])*dLD, po, h, n, n, length(betas)-1, betas)-cbind((Sxx[,1]-po*Sx[,1])*dLD, 0)
  S_xS1mSxdLD <- kndksum(po, (S1[,1]*po-Sx[,1])*dLD, po, h, n, n, length(betas)-1, betas)-cbind((S1[,1]*po-Sx[,1])*dLD, 0)
  S_xSxydLD  <- kndksum(po, Sxy[,1]*po*dLD, po, h, n, n, length(betas)-1, betas)-cbind(Sxy[,1]*po*dLD, 0)
  S_SxydLD  <- kndksum(po, Sxy[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(Sxy[,1]*dLD, 0)
  S_S1dLD  <- kndksum(po, S1[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(S1[,1]*dLD, 0)
  S_SxxdLD  <- kndksum(po, Sxx[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(Sxx[,1]*dLD, 0)
  S_SxdLD  <- kndksum(po, Sx[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(Sx[,1]*dLD, 0)
  S_hyS1dLD <- kndksum(po, hy*S1[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(hy*S1[,1]*dLD, 0)
  S_hySxxdLD <- kndksum(po, hy*Sxx[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(hy*Sxx[,1]*dLD, 0)
  S_hySxdLD <- kndksum(po, hy*Sx[,1]*dLD, po, h, n, n, length(betas)-1, betas)-cbind(hy*Sx[,1]*dLD, 0)
  dp <- 2*po*S_SydLD[,1]-S_xSydLD[,1]+po/h*S_xSydLD[,2]-po^2/h*S_SydLD[,2]-yo/h*S_SxxmxSxdLD[,2]+yo*S_xS1mSxdLD[,1]
  dp <- dp-po*yo/h*S_xS1mSxdLD[,2]-1/h*S_xSxydLD[,2]+po/h*S_SxydLD[,2]-S_SxydLD[,1]-2*po*S_hyS1dLD[,1]+po^2/h*S_hyS1dLD[,2]
  dp <- dp+1/h*S_hySxxdLD[,2]+2*S_hySxdLD[,1]-2*po/h*S_hySxdLD[,2]
  dp <- dp+dLD*(1/h*Sy[,2]*(Sx[,1]*po-Sxx[,1])+Sy[,1]*(po/h*Sx[,2]-1/h*Sxx[,2]-Sx[,1])-1/h*Sxy[,2]*(po*S1[,1]-Sx[,1]))
  dp <- dp+dLD*(Sxy[,1]*(1/h*Sx[,2]-po/h*S1[,2]+S1[,1])+hy/h*(S1[,2]*Sxx[,1]+Sxx[,2]*S1[,1]-2*Sx[,2]*Sx[,1]))
  c(dp)%*%(X[o,]/nv-(po)%*%t(v)/nv^2)
}


fancy_PPR_initialisation <- function(X, y, hmult, betas){
  Xt <- apply(X, 2, function(x){
    h <- sd(x)/length(x)^.2*hmult
    ret <- kLLreg(x, y, h, betas)
    ret/sqrt(mean(ret-y)^2)
  })
  glmnet::glmnet(Xt, y, intercept = FALSE, lambda = 1e-5)$beta@x
}

kPPR <- function(X, y, nterms = 1, hmult = 1, betas = NULL, loss = NULL, dloss = NULL, initialisation = 'fancy'){
  n <- nrow(X)
  d <- ncol(X)
  mu <- mean(y)
  mu_X <- colMeans(X)
  r <- y-mu
  X <- sweep(X, 2, mu_X, '-')
  hs <- numeric(nterms)
  vs <- matrix(0,nterms,d)
  fitted <- matrix(0,nterms,n)
  if(is.null(betas)) betas = c(1,1)
  for(tm in 1:nterms){
    if(initialisation=='fancy') v0 <- fancy_PPR_initialisation(X, r, hmult, betas)
    else if(initialisation=='lm') v0 <- glmnet::glmnet(X, r, intercept = FALSE, lambda = 1e-5)$beta@x
    else if(initialisation=='random') v0 <- rnorm(d)
    else if(is.function(initialisation)) v0 <- initialisation(X, r)
    else stop('argument "initialisation" must be a function of X and y or one of "fancy", "lm" and "random".')
    v0 <- v0/sqrt(sum(v0^2))
    h <- sd(X%*%v0)/n^.2*hmult
    vs[tm,] <- optim(v0, f_ppr, df_ppr, X, r, h, betas, loss, dloss, method = 'BFGS')$par
    vs[tm,] <- vs[tm,]/sqrt(sum(vs[tm,]^2))
    fitted[tm,] <- kLLreg(X%*%vs[tm,], r, h, betas)
    r <- r - fitted[tm,]
    hs[tm] <- h
  }
  sol <- list(mu = mu, mu_X = mu_X, y = y, X = X, hs = hs, vs = vs, fitted = fitted, betas = betas)
  class(sol) <- "kPPR"
  sol
}


predict.kPPR <- function(model, Xtest){
  Xtest <- sweep(Xtest, 2, model$mu_X, '-')
  betas <- model$betas
  n <- nrow(model$X)
  ntest <- nrow(Xtest)
  yhat <- model$mu
  for(tm in 1:length(model$hs)){
    h <- model$hs[tm]
    p <- model$X%*%model$vs[tm,]
    o <- order(p)
    ptest <- Xtest%*%model$vs[tm,]
    otest <- order(ptest)
    sK <- ksum(p[o], numeric(n)+1, ptest[otest], h, n, ntest, length(betas)-1, betas)
    sKx <- ksum(p[o], p[o], ptest[otest], h, n, ntest, length(betas)-1, betas)
    sKy <- ksum(p[o], model$fitted[tm,o], ptest[otest], h, n, ntest, length(betas)-1, betas)
    sKx2 <- ksum(p[o], p[o]^2, ptest[otest], h, n, ntest, length(betas)-1, betas)
    sKxy <- ksum(p[o], p[o]*model$fitted[tm,o], ptest[otest], h, n, ntest, length(betas)-1, betas)
    yhat <- yhat + (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*ptest[otest])/(sK*sKx2-sKx^2))[rank(ptest)]
  }
  yhat
}
