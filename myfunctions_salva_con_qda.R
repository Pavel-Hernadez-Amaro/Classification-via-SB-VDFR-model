# packages ----------------------------------------------------------------

library(MASS)
library(fda)
library(mvtnorm)
library(fpca)
library(fdapace)
source("fpca-modifica.R")


### functions

# data generation ---------------------------------------------------------

data.gen.splines.uneven <- function(n, p, grid, knots, mean, err, coef){
  
  ########################################################################################
  ##                                                                                    ##
  #       genera un campione di curve da una base spline con nodi irregolari
  # 
  #       genera una muestra de curvas a partir de una base spline con nodos irregulares
  ##                                                                                    ##
  ########################################################################################
  
  #     n: number of subjects
  #     p: number of observation points
  #  grid: evaluation grid
  # knots: vector of knots placement
  #  mean: mean vector for the coefficients
  #   err: standard error of measurement error
  #  coef: variance vector for the coefficients (the covariance matrix is supposed to be diagonal)
  
  x     = matrix(nrow = n, ncol = p)
  x.sm  = matrix(nrow = n, ncol = p)    
  
  basis = create.bspline.basis(rangeval=grid[c(1,p)], breaks=knots)
  
  for(i in 1:n){
    coefs    = rmvnorm(1, mean, coef*diag(rep(1,length(knots)+2)))
    x.sm[i,] = as.vector(coefs%*%t(eval.basis(grid, basis))) # smooth data (B-splines evaluated)
    x[i,]    = as.vector(coefs%*%t(eval.basis(grid, basis))) + rnorm(p, sd = err) # noisy data (B-splines evaluated + an error)
  }
  
  out = list(x,x.sm)
  names(out) = c("data", "smooth.data")
  return(out)
}


# missing generation ------------------------------------------------------

data.agg.miss <- function(x){
  
  ########################################################################################
  ##                                                                                    ##
  #       genera starting points ed ending points da una distribuzione uniforme          # 
  ##      
  #       genera una matriz con datos faltantes (las curvas partially observed)
  ##      
  #       todas las curvas tienen el intervalo común timepoints
  ########################################################################################
  
  # x: data matrix, for example "rbind()" of the previous data generation
  
  p    = ncol(x)
  n    = nrow(x)
  init = round(runif(n-1, 1, round(p/3)))          # generiamo i punti di inizio e fine
  fin  = round(runif(n-1, round(2*p/3+1), p))
  ep   = rbind(c(1,p),cbind(init, fin))
  
  for(i in 2:n){                                   # generiamo gli NAs
    aa = ep[i,1]
    bb = ep[i,2]
    x[i, 1:aa] = NA
    x[i, bb:p] = NA
  }
  
  obs = matrix(nrow = n, ncol = 2)
  for(i in 1:n){
    obs[i,] = c(min(which(!is.na(x[i,]))), max(which(!is.na(x[i,]))))
  }
  min = max(obs[,1])
  max = min(obs[,2])
  t = min:max
  m = length(t)
  
  domain=rbind(cbind(1,p),cbind(ep[2:n,1]+1,ep[2:n,2]-1))
  
  out = list(x, t, m, domain)
  names(out) = c("x.miss", "timepoints", "interval.length", "domain")
  return(out)
}

data.agg.miss.betadecay <- function(x, rate){
  
  ########################################################################################
  ##                                                                                    ##
  #         genera starting points ed ending points da una distribuzione beta            # 
  ##                                                                                    ##
  ########################################################################################
  
  # x: data matrix, for example "rbind()" of the previous data generation
  
  p    = ncol(x)
  n    = nrow(x)
  init = round((1-rbeta(n-1, 1, rate))*(p/3))          # generiamo i punti di inizio e fine
  fin  = round(2*p/3 + rbeta(n-1, 1, rate)*(p/3)) 
  ep   = rbind(c(1,p),cbind(init, fin))
  
  for(i in 2:n){                                       # generiamo gli NAs
    aa = ep[i,1]
    bb = ep[i,2]
    x[i, 1:aa] = NA
    x[i, bb:p] = NA
  }
  
  obs = matrix(nrow = n, ncol = 2)
  for(i in 1:n){
    obs[i,] = c(min(which(!is.na(x[i,]))), max(which(!is.na(x[i,]))))
  }
  min = max(obs[,1])
  max = min(obs[,2])
  t = min:max
  m = length(t)
  
  out = list(x, t, m)
  names(out) = c("x.miss", "timepoints", "interval.length")
  return(out)
}


# Functional PCA for partially observed data ------------------------------

pca.nomiss.class2 <- function(x, y, maxi, t.par, lambda.sm){      
  
  ########################################################################################
  ##                                                                                    ##
  #    performa functional PCA + QDA sugli scores, dopo aver effettuato lo smoothing     #
  #    dei dati, in un dato intervallo, per un numero di componenti massimo a piacere    #
  #                    *** NB: non funziona con dati mancanti ***                        #
  #
  #    realiza un PCA + QDA funcional en los scores, después de suavizar los datos,
  #    en un intervalo dado, para un número máximo de componentes a voluntad *** NB: no funciona con datos faltantes ***                                                                                  ##
  ########################################################################################
  
  #         x: data matrix 
  #         y: grouping variable
  #      maxi: maximum number of components
  #     t.par: the working inteval
  # lambda.sm: smoothing parameter
  
  x.cut    = x[,t.par]
  m        = length(t.par)
  
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m) # smoothing
  curv.Lfd = int2Lfd(2)
  fdob     = fdPar(bbasis,curv.Lfd,lambda.sm) 
  x.sm     = t(eval.fd(1:m, smooth.basis(1:m,t(x.cut),fdob)$fd))
  
  pca      = eigen(var(scale(x.sm, scale=F)))                     # principal component analysis
  scores.t = t(t(pca$vectors[,1:maxi]) %*% t(scale(x.sm, scale=F))) # dim should be dim(x.cut)[1] x maxi
  
  CV.err = c()                                                    # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.t[,1:i])
    mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best    = which(CV.err == min(CV.err, na.rm=T))[1]              # storing
  errcv.t = CV.err[best] 
  
  
  out = list(errcv.t,best,scores.t,pca$vectors[,1:maxi])
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
}

pca.miss.james <- function(x, y, maxi, t.par, k){      
  
  ########################################################################################
  ##                                                                                    ##
  #    performa "James et al.(2000)" + QDA sugli scores, per un dato numbero di basi,    #
  #         in un dato intervallo, per un numero di componenti massimo a piacere         #
  ##                                                                                    ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  #     k: number of basis for fpca.mle
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  data.y         = list()                    # preparing               
  data.curve     = list()
  data.timeindex = list()
  for(i in 1:n){
    data.timeindex[[i]] = which(apply(t(matrix(x.cut[i,])), 2, is.na)==F)
    data.y[[i]]         = x.cut[i,data.timeindex[[i]]]
    data.curve[[i]]     = rep(i,length(data.timeindex[[i]]))
  }
  data.y         = c(data.y, recursive=T)
  data.curve     = c(data.curve, recursive=T)
  data.timeindex = c(data.timeindex, recursive=T) 
  dati.m         = cbind(data.curve, data.y, data.timeindex)
  
  # principal component analysis
  environment(myfpca.mle) <- asNamespace("fpca")
  pca.james = tryCatch(myfpca.mle(data.m = dati.m, M.set = k, r.set = maxi, grids = seq(0,1,length.out=m),
                                ini.method = "EM", max.step = 1), error=function(e) NA)
  # computing scores 
  scores.james = tryCatch(modifica.fpca(data.m = dati.m, grids.u = 1:m, muhat = pca.james$fitted_mean,
                                        eigenvals = pca.james$eigenvalues, eigenfuncs = pca.james$eigenfunctions,
                                        sig2hat = pca.james$error_var, K = maxi), error=function(e) NA)
  
  CV.err = NA                               # classification with cross validation
  best = NA
  err.cv.james = NA
  if(!all(is.na(scores.james))){
    for(i in 1:maxi){
      scores    = cbind(scores.james[,1:i])
      mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
      CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
    }
    
    best      = which(CV.err == min(CV.err, na.rm=T))[1]      # storing
    err.cv.james = CV.err[best] 
  }
  
  if(!all(is.na(pca.james))) out4 = t(pca.james$eigenfunctions[1:maxi,]) else out4 = NA
  out = list(err.cv.james, best, scores.james, out4)
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
  
}

pca.miss.yao2 <- function(x, y, maxi, t.par, bwmean, bwcov){     # pca WITH missing data, Yao method
  
  ########################################################################################
  ##                                                                                    ##
  #     performa "Yao et al.(2005)" + QDA sugli scores, per due bandwidth fissate,       #
  #          in un dato intervallo, per un numero di componenti massimo a piacere        #
  ##                                                                                    ##
  ########################################################################################
  
  #      x: data matrix
  #      y: grouping variable
  #   maxi: maximum number of components
  #  t.par: working interval
  # bwmean: bandwidth for mean smoothing
  #  bwcov: bandwidth for covariance smoothing
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  data.y         = list()                   # preparing                   
  data.timeindex = list()
  for(i in 1:n){
    data.timeindex[[i]] = which(apply(t(matrix(x.cut[i,])), 2, is.na)==F)
    data.y[[i]]         = x.cut[i,data.timeindex[[i]]]
  }
  
  # principal component analysis
  pca.yao = tryCatch(FPCA(Ly = data.y, Lt = data.timeindex, optns=list(methodXi="CE", maxK=maxi,
                          nRegGrid=length(t.par), methodMuCovEst = "smooth", methodSelectK=maxi,
                          userBwCov=bwcov, userBwMu=bwmean)),
                     error=function(e) NA)
  
  scores.yao = pca.yao$xiEst
  
  CV.err = NA                               # classification with cross validation
  best = NA
  err.cv.yao = NA
  if(!all(is.na(scores.yao))){
    for(i in 1:maxi){
      scores    = cbind(scores.yao[,1:i])
      mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
      CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
    }
    
    best      = which(CV.err == min(CV.err, na.rm=T))[1]      # storing
    err.cv.yao = CV.err[best] 
  }
  
  if(!all(is.na(pca.yao)) && dim(pca.yao$phi[,1:maxi])[1]==length(t.par) ) out4 = pca.yao$phi[,1:maxi] else out4 = matrix(NA,nrow=length(t.par), ncol=maxi)
  out = list(err.cv.yao, best, scores.yao, out4)
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
  
  
}

pca.miss.buja2 <- function(x, y, maxi, t.par, lambda.buja){     
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, per un dato livello di       #
  #           smoothing - *uguale per tutte le PC* -  in un dato intervallo,             #
  #                 per un numero di componenti massimo a piacere                        #
  ##                                                                                    ##
  ########################################################################################
  
  #           x: data matrix
  #           y: grouping variable
  #        maxi: maximum number of components
  #       t.par: working inteval
  # lambda.buja: smoothing parameter
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  lambda.buja = lambda.buja                                      # setting parameters
  tol = 1e-5
  tt  = 1
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  
  for(k in 1:maxi){
    
    vv1 = rep(1, m)
    v1  = vv1
    u1  = c()
    v1. = c()
    z1  = c()
    
    while(tt>tol){
      for(i in 1:n){                                # scores step
        obs   = !is.na(x.cut[i,])
        u1[i] = x.cut[i,obs] %*% v1[obs]
      }
      for(j in 1:m){                                # eigenvectors step
        obs   = !is.na(x.cut[,j])
        z1[j] = t(x.cut[obs,j]) %*% u1[obs]
      }
      
      fdob = fdPar(bbasis,curv.Lfd,lambda.buja)          # smoothing step
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
    }
    
    tt = 1
    x.cut = x.cut - u1 %*% t(v1)                    # matrix deflation
    eigenfun.buja[,k] <- v1
    scores.buja[,k] <- u1
  }
  
  CV.err = c()                                      # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.buja[,1:i])
    mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best      = which(CV.err == min(CV.err, na.rm=T))[1]    # storing
  err.cv.buja = CV.err[best]
  
  
  out = list(err.cv.buja, best, scores.buja, eigenfun.buja)
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
  
}


# main function -----------------------------------------------------------

all.together.splines.uneven <- function(knots, mean1, err1, mean2, err2, n1, n2, p, grid, maxi, maxl,
                                        k, lambda.sm, lambda.buja, coef, bwmean, bwcov, unif.miss=T,
                                        rate=NA){
  
  ########################################################################################
  ##                                                                                    ##
  #            genera i dati e performa la classificazione con tutti i metodi,           #
  #             per *ogni* estensione del dominio, con relativi parametri di             #
  #                 smoothing *fissi*, per un tipo di missing a piacere                  #
  #                   e per un numero di componenti massimo a piacere                    #
  ##                                                                                    ##
  ########################################################################################
  
  #       knots: vector of knots placement
  #       mean1: mean vector for the coefficients, group 1
  #        err1: standard error of measurement error, group 1
  #       mean2: mean vector for the coefficients, group 2
  #        err2: standard error of measurement error, group 2
  #        coef: variance vector for the coefficients (the covariance matrix is supposed to be diagonal)
  #          n1: number of subjects, group 1
  #          n2: number of subjects, group 2
  #           p: number of observation points
  #        grid: evaluation grid
  #        maxi: maximum number of components
  #        maxl: maximum number of extensions
  #           k: number of bases (for "James et al." method)
  #   lambda.sm: presmoothing parameter (for standard fPCA with complete data)
  # lambda.buja: smoothing parameter (for "Huang et al." method)
  #      bwmean: bandwidth for mean smoothing (for "Yao et al." method)
  #       bwcov: bandwidth for covariance smoothing (for "Yao et al." method)
  
  
  ### data generation
  
  x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  xbind    <- rbind(x1$data, x2$data)
  x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  y        <- c(rep(0,n1), rep(1,n2))  
  x        <- if(unif.miss==T) data.agg.miss(xbind) else data.agg.miss.betadecay(xbind) 
  
  err.cv.t     <- c()
  err.cv.james <- c()
  err.cv.yao   <- c()
  err.cv.buja  <- c()
  
  ### preparing
  
  eigenfun.t     <- array(dim=c(p, maxi, maxl+1))
  eigenfun.james <- array(dim=c(p, maxi, maxl+1))
  eigenfun.yao   <- array(dim=c(p, maxi, maxl+1))
  eigenfun.buja  <- array(dim=c(p, maxi, maxl+1))
  
  scores.t     <- array(dim=c(n, maxi, maxl+1))
  scores.james <- array(dim=c(n, maxi, maxl+1))
  scores.yao   <- array(dim=c(n, maxi, maxl+1))
  scores.buja  <- array(dim=c(n, maxi, maxl+1))
  
  best.t     <- c()
  best.james <- c()
  best.buja  <- c()
  best.yao   <- c()
  
  ### data analysis
  
  for(l in 0:maxl){
    
    t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
    
    pcanomiss <- pca.nomiss.class2(xbind, y, maxi, t.par, lambda.sm)
    pcajames  <- pca.miss.james(x$x.miss, y, maxi, t.par, k)
    pcayao    <- pca.miss.yao2(x$x.miss, y, maxi, t.par, bwmean, bwcov)
    pcabuja   <- pca.miss.buja2(x$x.miss, y, maxi, t.par, lambda.buja)
    
    
    ### storing
    
    err.cv.t[l+1]     = pcanomiss$err
    err.cv.james[l+1] = pcajames$err
    err.cv.yao[l+1]   = pcayao$err
    err.cv.buja[l+1]  = pcabuja$err
    
    eigenfun.t[t.par,,l+1]     = pcanomiss$eigenfunctions
    eigenfun.james[t.par,,l+1] = pcajames$eigenfunctions
    eigenfun.yao[t.par,,l+1]   = pcayao$eigenfunctions
    eigenfun.buja[t.par,,l+1]  = pcabuja$eigenfunctions
    
    scores.t[,,l+1]     = pcanomiss$scores
    scores.james[,,l+1] = pcajames$scores
    scores.yao[,,l+1]   = pcayao$scores
    scores.buja[,,l+1]  = pcabuja$scores
    
    best.t[l+1]     = pcanomiss$best
    best.james[l+1] = pcajames$best
    best.yao[l+1]   = pcayao$best
    best.buja[l+1]  = pcabuja$best
    
  }
  
  out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, eigenfun.t, eigenfun.james,
             eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja,
             best.t, best.james, best.yao, best.buja)
  names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "eigenfun.t",
                 "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james",
                 "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja")
  return(out)
}


