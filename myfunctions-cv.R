library(caret)
# source("C:/Users/Marco/Desktop/Dottorato - III anno/Partially Observed FDA -- the day after tomorrow/Functions/fpca-modifica.R")
# source("C:/Users/Marco/Desktop/Dottorato - III anno/Partially Observed FDA -- the day after tomorrow/Functions/pred.missfd2.R")

pca.miss.james.cv <- function(x, y, maxi, t.par, k){        
  
  ########################################################################################
  ##                                                                                    ##
  #        performa "James et al.(2000)" + QDA sugli scores, in un dato intervallo,      #
  #         per un numero di componenti massimo a piacere, scegliendo il numero di       #
  #          basi ***con la CV implementata nel pacchetto "fpca"***, cercando            #
  #                   su una griglia di valori definita dall'utente                      #
  ##                     
  ##      realizar "James et al.(2000)" + QDA sobre las puntuaciones, en un intervalo dado,
  #       para un número máximo de componentes a voluntad, eligiendo el número de bases
  #       ***con el CV implementado en el paquete "fpca"** *, buscando en una cuadrícula de valores definida por el usuario
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  #     k: number of basis for fpca.mle (MUST BE A VECTOR ---> CV)
  
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
  k.opt = pca.james$selected_model[1]
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
  out = list(err.cv.james, best, scores.james, out4, k.opt)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "opt.basis")
  return(out)
  
}

pca.miss.yao2.cv <- function(x, y, maxi, t.par){     
  
  ########################################################################################
  ##                                                                                    ##
  #         performa "Yao et al.(2005)" + QDA sugli scores, in un dato intervallo,       #
  #       per un numero di componenti massimo a piacere, scegliendo il valore delle      #
  #           bandwidths ***con la CV implementata nel pacchetto "fdapace"***            #
  ##            
  #       realizar "Yao et al.(2005)" + QDA sobre las puntuaciones, en un intervalo dado,
  #       para un número máximo de componentes a voluntad, eligiendo el valor de los anchos de banda
  #       ***con el CV implementado en el paquete "fdapace"* **
  ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working interval
  
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
                                                                       nRegGrid=length(t.par), methodMuCovEst = "smooth", 
                                                                       methodBwMu="GCV", methodBwCov="GCV",methodSelectK=maxi,
                                                                       kFoldMuCov=5)),
                     error=function(e) NA)
  
  scores.yao = pca.yao$xiEst
  bw.m = pca.yao$bwMu
  bw.c = pca.yao$bwCov
  
  
  CV.err_aux_pred=NA
  CV.err_aux=NA
  optimal_cut_yao=NA
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
    
    # for (i in 1:N) {
    #   
    #   scores_aux    = cbind(scores.yao[-i,1:5])
    #   mod_aux       = tryCatch(qda(y[-i] ~ scores_aux, CV=T), error=function(e) NA)
    #   
    #   for (o in 1:length(max.col(mod_aux$posterior))) {
    #     
    #     aux=max.col(mod_aux$posterior)[o]
    #     
    #     prob_post_sample=mod_aux$posterior[,aux]
    #   
    #   }
    #   
    #   prob_post_sample=as.matrix(prob_post_sample)
    #   
    #   aux=ROC(prob_post_sample,y[-i],plot = NULL)
    #   best_pos=which.max(aux$res[,1]+aux$res[,2])
    #   optimal_cut_yao[i]=aux$res[best_pos,5]
    #   
    #   qda(y[i] ~ t(cbind(scores.yao[i,1:5])), CV=T)
    #   
    # 
    #         
    # }
    
  }
  
  if(!all(is.na(pca.yao)) && dim(pca.yao$phi[,1:maxi])[1]==length(t.par) ) out4 = pca.yao$phi[,1:maxi] else out4 = matrix(NA,nrow=length(t.par), ncol=maxi)
  out = list(err.cv.yao, best, scores.yao, out4, bw.m, bw.c)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "opt.bwm", "opt.bwc")
  return(out)
  
  
}

pca.miss.buja3 <- function(x, y, maxi, t.par, lambda.buja){    
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, per un dato livello di       #
  #             smoothing - *diverso per ogni PC* -  in un dato intervallo,              #
  #                 per un numero di componenti massimo a piacere                        #
  # 
  # 
  ##      realiza "Huang et al. (2008)" + QDA en puntajes, para un nivel dado de suavizado - *diferente para cada PC* -
  #       en un intervalo dado, para un número máximo de componentes a voluntad                                                  
  ##
  ########################################################################################
  
  #           x: data matrix
  #           y: grouping variable
  #        maxi: maximum number of components
  #       t.par: working inteval
  # lambda.buja: smoothing parameter (MUST BE A VECTOR)
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  tol = 1e-5
  tt  = 1
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  flds <- createFolds(y, k=5, list = TRUE, returnTrain = FALSE)
  
  for(h in 1:maxi){
    
    lambda.buja.now = lambda.buja[h] 
    
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
      
      fdob = fdPar(bbasis,curv.Lfd,lambda.buja.now)          # smoothing step
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
    };tt = 1
    
    eigenfun.buja[,h] <- v1
    scores.buja[,h] <- u1
    
    x.cut = x.cut - u1 %*% t(v1)   
  }
  
  
  CV.err = c()                                      # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.buja[,1:i])
    mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best      = which(CV.err == min(CV.err, na.rm=T))[1]    # storing
  err.cv.buja = CV.err[best]
  
  
  out = list(err.cv.buja, best, scores.buja, eigenfun.buja, lambda.buja)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "lambdas used")
  return(out)
  
}

pca.miss.buja2.cv <- function(x, y, maxi, t.par, lambda.buja){  
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, in un dato intervallo,       #
  #       per un numero di componenti massimo a piacere, scegliendo il valore del        #
  #         parametro di smoothing per ogni componente ***con la CV implementata         #
  #                    minimizzando l'errore di ricostruzione***                         #
  ##      
  #       realizar "Huang et al.(2008)" + QDA sobre las puntuaciones, en un intervalo dado,
  #       para un número máximo de componentes a voluntad, eligiendo el valor del parámetro de suavizado
  #       para cada componente *** con el CV implementado minimizando el error de reconstrucción***
  # 
  ##
  ########################################################################################
  
  #           x: data matrix
  #           y: grouping variable
  #        maxi: maximum number of components
  #       t.par: working inteval
  # lambda.buja: smoothing parameter (MUST BE A VECTOR ---> CV)
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  tol = 1e-5
  tt  = 1
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  lambda.opt = c()
  
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  flds <- createFolds(y, k=5, list = TRUE, returnTrain = FALSE)
  
  for(h in 1:maxi){
    
    cvs <- c()
    
    for(a in 1:length(lambda.buja)){
      
      lambda.buja.now = lambda.buja[a] 
      
      cvscore <- c()
      
      for(w in 1:5){
        
        test <- x.cut[flds[[w]],]
        train  <- x.cut[-flds[[w]],]
        
        vv1 = rep(1, m)
        v1  = vv1
        u   = vector(length=n)
        u1  = c()
        u2  = c()
        v1. = c()
        z1  = c()
        
        while(tt>tol){
          for(i in 1:nrow(train)){                                # scores step
            obs   = !is.na(train[i,])
            u1[i] = train[i,obs] %*% v1[obs]
          }
          for(j in 1:m){                                # eigenvectors step
            obs   = !is.na(train[,j])
            z1[j] = t(train[obs,j]) %*% u1[obs]
          }
          
          fdob = fdPar(bbasis,curv.Lfd,lambda.buja.now)          # smoothing step
          v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
          v1.  = v1./sqrt(sum(v1.^2))
          tt   = sum((v1.-v1)^2)
          v1   = v1.
        };tt = 1
        
        for(i in 1:nrow(test)){                                # scores step
          obs   = !is.na(test[i,])
          u2[i] = test[i,obs] %*% v1[obs]
        }
        
        #u[-flds[[w]]] <- u1
        #u[flds[[w]]] <- u2
        
        obj <- matrix(nrow=length(flds[[w]]), ncol=m)
        
        for(i in 1:length(flds[[w]])){
          for(k in 1:m){
            obj[i,k] <- (x.cut[flds[[w]][i],k] - u2[i]*v1[k])^2
          }
        }
        cvscore[w] <- sum((apply(obj, 2, sum, na.rm=T))/(n*m))
      }
      
      cvs[a] <- mean(cvscore)
      
    }
    
    bestlambda <- lambda.buja[which.min(cvs)]
    
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
      
      fdob = fdPar(bbasis,curv.Lfd,bestlambda)          # smoothing step
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
    };tt = 1
    
    lambda.opt[h] = bestlambda
    
    eigenfun.buja[,h] <- v1
    scores.buja[,h] <- u1
    
    x.cut = x.cut - u1 %*% t(v1)   
  }
  
  
  CV.err = c()                                      # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.buja[,1:i])
    mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best      = which(CV.err == min(CV.err, na.rm=T))[1]    # storing
  err.cv.buja = CV.err[best]
  
  
  out = list(err.cv.buja, best, scores.buja, eigenfun.buja, lambda.opt)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "opt.lambda")
  return(out)
  
}

pca.miss.buja.gcv <- function(x, y, maxi, t.par, max.rs){     
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, con un livello di            #
  #     smoothing - *scelto con Generalized Cross Validation per ogni PC* -  in          #
  #        un dato intervallo, per un numero di componenti massimo a piacere             #
  ##                 
  #       realiza "Huang et al.(2008)" + QDA en las puntuaciones, con un nivel de suavizado
  #       - *elegido con Validación Cruzada Generalizada para cada PC* - en un intervalo dado,
  #       para un número máximo de componentes a voluntad
  ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  
  tol = 1e-5
  tt  = 1
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  best.lam = 1
  
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
      sol = NA
      rs  = 0
      while(is.na(sol) & rs<max.rs){
        sol = tryCatch(optim(runif(1,0.001, best.lam^2), gcvfun, m=m, z1=z1, bbasis=bbasis, 
                             curv.Lfd=curv.Lfd)$par, error=function(e) NA)
        rs = rs + 1
      }
      best.lam = sol
      if(!is.na(best.lam)) fdob = fdPar(bbasis,curv.Lfd,best.lam) else fdob = fdPar(bbasis,curv.Lfd,0.0001)
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
      print(best.lam)
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

gcvfun <- function(m, z1, bbasis, curv.Lfd, alp){ # THIS FUNCTION IS USED IN THE pca.miss.buja.gcv()
  fdob = fdPar(bbasis,curv.Lfd,alp) 
  gcv  = smooth.basis(1:m,z1,fdob)$gcv
  return(gcv)
}

pca.miss.kraus <- function(x, y, maxi, t.par, max.rs, sd.rs){      
  
  ########################################################################################
  ##                                                                                    ##
  #              performa "Kraus (2005)" + QDA sugli scores, in un dato                  #
  #             intervallo, per un numero di componenti massimo a piacere                #
  ##                                                                                    ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  mu = mean.missfd(x.cut)
  R  = tryCatch(var.missfd(x.cut), error=function(e) NA)
  
  if(!is.na(R)){
    
    eig.R = eigen.missfd(R)
    
    scores.kraus <- matrix(nrow=n, ncol=maxi)
    for(i in 1:n){
      for(j in 1:maxi){
        scores.kraus[i,j] = pred.score.missfd(x.cut[i,],phi=eig.R$vectors[,j],x=x.cut,
                                              max.rs=max.rs, sd.rs=sd.rs)[1]
        
      }
    }
    
    CV.err = NA                               # classification with cross validation
    best = NA
    err.cv.kraus = NA
    if(!all(is.na(scores.kraus))){
      for(i in 1:maxi){
        scores    = cbind(scores.kraus[,1:i])
        mod       = tryCatch(qda(y ~ scores, CV=T), error=function(e) NA)
        CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
      }
      
      best      = which(CV.err == min(CV.err, na.rm=T))[1]      # storing
      err.cv.kraus = CV.err[best] 
    }
    
    eigenfunctions.kraus = eig.R$vectors[,1:maxi]
  } else {
    err.cv.kraus = NA 
    best = NA
    scores.kraus = matrix(nrow=n, ncol=maxi)
    eigenfunctions.kraus = matrix(nrow=m,ncol=maxi)
  }
  
  out = list(err.cv.kraus, best, scores.kraus, eigenfunctions.kraus)
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
  
}

flda.miss.james.cv <- function(x, y, t.par, k){      
  
  ########################################################################################
  ##                                                                                    ##
  #            performa "James & Hastie (2001)", in un dato intervallo,                  #
  #        per un numbero di basi scelto con CV sul misclassification error              #
  ##                                                                                    ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  # t.par: working interval
  #     k: number of basis for fldafit (MUST BE A VECTOR ---> CV)
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  flds <- createFolds(y, k=5, list = TRUE, returnTrain = FALSE)
  cvs <- c()
  
  for(a in 1:length(k)){
    
    cvscore <- c()
    
    for(w in 1:5){
      
      test   <- x.cut[flds[[w]],]
      train  <- x.cut[-flds[[w]],]
      ytest  <- y[flds[[w]]]
      ytrain <- y[-flds[[w]]]
      
      train.y         = list()                    # preparing training set            
      train.curve     = list()
      train.timeindex = list()
      for(i in 1:nrow(train)){
        train.timeindex[[i]] = which(apply(t(matrix(train[i,])), 2, is.na)==F)
        train.y[[i]]         = train[i,train.timeindex[[i]]]
        train.curve[[i]]     = rep(i,length(train.timeindex[[i]]))
      }
      train.y         = c(train.y, recursive=T)
      train.curve     = c(train.curve, recursive=T)
      train.timeindex = c(train.timeindex, recursive=T) 
      train.m         = list(train.y, train.curve, train.timeindex, ytrain+1)
      names(train.m)  = c("y", "curve", "timeindex", "class")
      
      test.y         = list()                    # preparing testing set            
      test.curve     = list()
      test.timeindex = list()
      for(i in 1:nrow(test)){
        test.timeindex[[i]] = which(apply(t(matrix(test[i,])), 2, is.na)==F)
        test.y[[i]]         = test[i,test.timeindex[[i]]]
        test.curve[[i]]     = rep(i,length(test.timeindex[[i]]))
      }
      test.y         = c(test.y, recursive=T)
      test.curve     = c(test.curve, recursive=T)
      test.timeindex = c(test.timeindex, recursive=T) 
      test.m         = list(test.y, test.curve, test.timeindex, ytest+1)
      names(test.m)  = c("y", "curve", "timeindex", "class")
      
      nbasis = k[a]
      flda.james = fldafit(data=train.m, grid=t.par, q=nbasis)
      pred = fldapred(flda.james, data=test.m)
      err = sum(pred$class.pred!=ytest+1)/nrow(test)
      cvscore[w] <- err
      
    }
    
    cvs[a] <- mean(cvscore)
    
  }
  
  bestbasis <- k[which.min(cvs)]
  err.cv.flda <- cvs[which.min(cvs)]*n
  
  
  out = list(err.cv.flda, bestbasis)
  names(out) = c("err", "best basis")
  return(out)
  
}
###

pca.miss.james.cv_lda <- function(x, y, maxi, t.par, k){        
  
  ########################################################################################
  ##                                                                                    ##
  #        performa "James et al.(2000)" + QDA sugli scores, in un dato intervallo,      #
  #         per un numero di componenti massimo a piacere, scegliendo il numero di       #
  #          basi ***con la CV implementata nel pacchetto "fpca"***, cercando            #
  #                   su una griglia di valori definita dall'utente                      #
  ##                     
  ##      realizar "James et al.(2000)" + QDA sobre las puntuaciones, en un intervalo dado,
  #       para un número máximo de componentes a voluntad, eligiendo el número de bases
  #       ***con el CV implementado en el paquete "fpca"** *, buscando en una cuadrícula de valores definida por el usuario
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  #     k: number of basis for fpca.mle (MUST BE A VECTOR ---> CV)
  
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
  k.opt = pca.james$selected_model[1]
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
      mod       = tryCatch(lda(y ~ scores, CV=T), error=function(e) NA)
      CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
    }
    
    best      = which(CV.err == min(CV.err, na.rm=T))[1]      # storing
    err.cv.james = CV.err[best] 
  }
  
  if(!all(is.na(pca.james))) out4 = t(pca.james$eigenfunctions[1:maxi,]) else out4 = NA
  out = list(err.cv.james, best, scores.james, out4, k.opt)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "opt.basis")
  return(out)
  
}

pca.miss.yao2.cv_lda <- function(x, y, maxi, t.par){     
  
  ########################################################################################
  ##                                                                                    ##
  #         performa "Yao et al.(2005)" + QDA sugli scores, in un dato intervallo,       #
  #       per un numero di componenti massimo a piacere, scegliendo il valore delle      #
  #           bandwidths ***con la CV implementata nel pacchetto "fdapace"***            #
  ##            
  #       realizar "Yao et al.(2005)" + QDA sobre las puntuaciones, en un intervalo dado,
  #       para un número máximo de componentes a voluntad, eligiendo el valor de los anchos de banda
  #       ***con el CV implementado en el paquete "fdapace"* **
  ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working interval
  
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
                                                                       nRegGrid=length(t.par), methodMuCovEst = "smooth", 
                                                                       methodBwMu="GCV", methodBwCov="GCV",methodSelectK=maxi,
                                                                       kFoldMuCov=5)),
                     error=function(e) NA)
  
  scores.yao = pca.yao$xiEst
  bw.m = pca.yao$bwMu
  bw.c = pca.yao$bwCov
  
  
  
  CV.err = NA                               # classification with cross validation
  best = NA
  err.cv.yao = NA
  if(!all(is.na(scores.yao))){
    for(i in 1:maxi){
      scores    = cbind(scores.yao[,1:i])
      mod       = tryCatch(lda(y ~ scores, CV=T), error=function(e) NA)
      CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
    }
    
    best      = which(CV.err == min(CV.err, na.rm=T))[1]      # storing
    err.cv.yao = CV.err[best] 
  }
  
  if(!all(is.na(pca.yao)) && dim(pca.yao$phi[,1:maxi])[1]==length(t.par) ) out4 = pca.yao$phi[,1:maxi] else out4 = matrix(NA,nrow=length(t.par), ncol=maxi)
  out = list(err.cv.yao, best, scores.yao, out4, bw.m, bw.c)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "opt.bwm", "opt.bwc")
  return(out)
  
  
}

pca.miss.buja3_lda <- function(x, y, maxi, t.par, lambda.buja){    
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, per un dato livello di       #
  #             smoothing - *diverso per ogni PC* -  in un dato intervallo,              #
  #                 per un numero di componenti massimo a piacere                        #
  # 
  # 
  ##      realiza "Huang et al. (2008)" + QDA en puntajes, para un nivel dado de suavizado - *diferente para cada PC* -
  #       en un intervalo dado, para un número máximo de componentes a voluntad                                                  
  ##
  ########################################################################################
  
  #           x: data matrix
  #           y: grouping variable
  #        maxi: maximum number of components
  #       t.par: working inteval
  # lambda.buja: smoothing parameter (MUST BE A VECTOR)
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  tol = 1e-5
  tt  = 1
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  flds <- createFolds(y, k=5, list = TRUE, returnTrain = FALSE)
  
  for(h in 1:maxi){
    
    lambda.buja.now = lambda.buja[h] 
    
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
      
      fdob = fdPar(bbasis,curv.Lfd,lambda.buja.now)          # smoothing step
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
    };tt = 1
    
    eigenfun.buja[,h] <- v1
    scores.buja[,h] <- u1
    
    x.cut = x.cut - u1 %*% t(v1)   
  }
  
  
  CV.err = c()                                      # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.buja[,1:i])
    mod       = tryCatch(lda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best      = which(CV.err == min(CV.err, na.rm=T))[1]    # storing
  err.cv.buja = CV.err[best]
  
  
  out = list(err.cv.buja, best, scores.buja, eigenfun.buja, lambda.buja)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "lambdas used")
  return(out)
  
}

pca.miss.buja2.cv_lda <- function(x, y, maxi, t.par, lambda.buja){  
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, in un dato intervallo,       #
  #       per un numero di componenti massimo a piacere, scegliendo il valore del        #
  #         parametro di smoothing per ogni componente ***con la CV implementata         #
  #                    minimizzando l'errore di ricostruzione***                         #
  ##      
  #       realizar "Huang et al.(2008)" + QDA sobre las puntuaciones, en un intervalo dado,
  #       para un número máximo de componentes a voluntad, eligiendo el valor del parámetro de suavizado
  #       para cada componente *** con el CV implementado minimizando el error de reconstrucción***
  # 
  ##
  ########################################################################################
  
  #           x: data matrix
  #           y: grouping variable
  #        maxi: maximum number of components
  #       t.par: working inteval
  # lambda.buja: smoothing parameter (MUST BE A VECTOR ---> CV)
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  tol = 1e-5
  tt  = 1
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  lambda.opt = c()
  
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  flds <- createFolds(y, k=5, list = TRUE, returnTrain = FALSE)
  
  for(h in 1:maxi){
    
    cvs <- c()
    
    for(a in 1:length(lambda.buja)){
      
      lambda.buja.now = lambda.buja[a] 
      
      cvscore <- c()
      
      for(w in 1:5){
        
        test <- x.cut[flds[[w]],]
        train  <- x.cut[-flds[[w]],]
        
        vv1 = rep(1, m)
        v1  = vv1
        u   = vector(length=n)
        u1  = c()
        u2  = c()
        v1. = c()
        z1  = c()
        
        while(tt>tol){
          for(i in 1:nrow(train)){                                # scores step
            obs   = !is.na(train[i,])
            u1[i] = train[i,obs] %*% v1[obs]
          }
          for(j in 1:m){                                # eigenvectors step
            obs   = !is.na(train[,j])
            z1[j] = t(train[obs,j]) %*% u1[obs]
          }
          
          fdob = fdPar(bbasis,curv.Lfd,lambda.buja.now)          # smoothing step
          v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
          v1.  = v1./sqrt(sum(v1.^2))
          tt   = sum((v1.-v1)^2)
          v1   = v1.
        };tt = 1
        
        for(i in 1:nrow(test)){                                # scores step
          obs   = !is.na(test[i,])
          u2[i] = test[i,obs] %*% v1[obs]
        }
        
        #u[-flds[[w]]] <- u1
        #u[flds[[w]]] <- u2
        
        obj <- matrix(nrow=length(flds[[w]]), ncol=m)
        
        for(i in 1:length(flds[[w]])){
          for(k in 1:m){
            obj[i,k] <- (x.cut[flds[[w]][i],k] - u2[i]*v1[k])^2
          }
        }
        cvscore[w] <- sum((apply(obj, 2, sum, na.rm=T))/(n*m))
      }
      
      cvs[a] <- mean(cvscore)
      
    }
    
    bestlambda <- lambda.buja[which.min(cvs)]
    
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
      
      fdob = fdPar(bbasis,curv.Lfd,bestlambda)          # smoothing step
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
    };tt = 1
    
    lambda.opt[h] = bestlambda
    
    eigenfun.buja[,h] <- v1
    scores.buja[,h] <- u1
    
    x.cut = x.cut - u1 %*% t(v1)   
  }
  
  
  CV.err = c()                                      # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.buja[,1:i])
    mod       = tryCatch(lda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best      = which(CV.err == min(CV.err, na.rm=T))[1]    # storing
  err.cv.buja = CV.err[best]
  
  
  out = list(err.cv.buja, best, scores.buja, eigenfun.buja, lambda.opt)
  names(out) = c("err", "best.comp", "scores", "eigenfunctions", "opt.lambda")
  return(out)
  
}

pca.miss.buja.gcv_lda <- function(x, y, maxi, t.par, max.rs){     
  
  ########################################################################################
  ##                                                                                    ##
  #       performa "Huang et al.(2008)" + QDA sugli scores, con un livello di            #
  #     smoothing - *scelto con Generalized Cross Validation per ogni PC* -  in          #
  #        un dato intervallo, per un numero di componenti massimo a piacere             #
  ##                 
  #       realiza "Huang et al.(2008)" + QDA en las puntuaciones, con un nivel de suavizado
  #       - *elegido con Validación Cruzada Generalizada para cada PC* - en un intervalo dado,
  #       para un número máximo de componentes a voluntad
  ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  eigenfun.buja = matrix(nrow=m, ncol=maxi)
  scores.buja   = matrix(nrow=n, ncol=maxi)
  
  tol = 1e-5
  tt  = 1
  bbasis   = create.bspline.basis(rangeval=m,norder=4,breaks=1:m)
  curv.Lfd = int2Lfd(2)
  best.lam = 1
  
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
      sol = NA
      rs  = 0
      while(is.na(sol) & rs<max.rs){
        sol = tryCatch(optim(runif(1,0.001, best.lam^2), gcvfun, m=m, z1=z1, bbasis=bbasis, 
                             curv.Lfd=curv.Lfd)$par, error=function(e) NA)
        rs = rs + 1
      }
      best.lam = sol
      if(!is.na(best.lam)) fdob = fdPar(bbasis,curv.Lfd,best.lam) else fdob = fdPar(bbasis,curv.Lfd,0.0001)
      v1.  = eval.fd(1:m, smooth.basis(1:m,z1,fdob)$fd)
      v1.  = v1./sqrt(sum(v1.^2))
      tt   = sum((v1.-v1)^2)
      v1   = v1.
      print(best.lam)
    }
    
    tt = 1
    x.cut = x.cut - u1 %*% t(v1)                    # matrix deflation
    eigenfun.buja[,k] <- v1
    scores.buja[,k] <- u1
  }
  
  CV.err = c()                                      # classification with cross validation
  for(i in 1:maxi){
    scores    = cbind(scores.buja[,1:i])
    mod       = tryCatch(lda(y ~ scores, CV=T), error=function(e) NA)
    CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
  }
  
  best      = which(CV.err == min(CV.err, na.rm=T))[1]    # storing
  err.cv.buja = CV.err[best]
  
  
  out = list(err.cv.buja, best, scores.buja, eigenfun.buja)
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
  
}

pca.miss.kraus_lda <- function(x, y, maxi, t.par, max.rs, sd.rs){      
  
  ########################################################################################
  ##                                                                                    ##
  #              performa "Kraus (2005)" + QDA sugli scores, in un dato                  #
  #             intervallo, per un numero di componenti massimo a piacere                #
  ##                                                                                    ##
  ########################################################################################
  
  #     x: data matrix
  #     y: grouping variable
  #  maxi: maximum number of components
  # t.par: working inteval
  
  m      = length(t.par)
  n      = nrow(x)
  x.cut  = scale(x[,t.par], scale=F)
  
  mu = mean.missfd(x.cut)
  R  = tryCatch(var.missfd(x.cut), error=function(e) NA)
  
  if(!is.na(R)){
    
    eig.R = eigen.missfd(R)
    
    scores.kraus <- matrix(nrow=n, ncol=maxi)
    for(i in 1:n){
      for(j in 1:maxi){
        scores.kraus[i,j] = pred.score.missfd(x.cut[i,],phi=eig.R$vectors[,j],x=x.cut,
                                              max.rs=max.rs, sd.rs=sd.rs)[1]
        
      }
    }
    
    CV.err = NA                               # classification with cross validation
    best = NA
    err.cv.kraus = NA
    if(!all(is.na(scores.kraus))){
      for(i in 1:maxi){
        scores    = cbind(scores.kraus[,1:i])
        mod       = tryCatch(lda(y ~ scores, CV=T), error=function(e) NA)
        CV.err[i] = if(class(mod)=="list") sum(mod$class!=y) else NA
      }
      
      best      = which(CV.err == min(CV.err, na.rm=T))[1]      # storing
      err.cv.kraus = CV.err[best] 
    }
    
    eigenfunctions.kraus = eig.R$vectors[,1:maxi]
  } else {
    err.cv.kraus = NA 
    best = NA
    scores.kraus = matrix(nrow=n, ncol=maxi)
    eigenfunctions.kraus = matrix(nrow=m,ncol=maxi)
  }
  
  out = list(err.cv.kraus, best, scores.kraus, eigenfunctions.kraus)
  names(out) = c("err", "best", "scores", "eigenfunctions")
  return(out)
  
}

####
SB_class = function(x, y, domain, grid, c1, c2){ # numeric approximation of the best cut (same cut as the ROC in EPI but much slower)

#x: data matrix
#y: grouping variable
#grid: grid where evaluate the B-splines # same as the one passed to data.gen.splines.uneven()
#domain: data domain 
#k: number of basis for SB_VDFR model (MUST BE A 2D VECTOR)

# err=NULL
# tag=1
# optimal_cut_grid=seq(0,1,10^-tag)
# err_old=100
# err_new=0


# k=10
# folds=createFolds(y, k=k, list = TRUE, returnTrain = FALSE)

#err=NULL
# probs=array(dim=c(length(folds),dim(x)[1]*(k-1)/k))
# probs_test=array(dim=c(length(folds),dim(x)[1]/k))
# optimal_cut_aux=array(dim=c(length(folds)))


# for (i in 1:k) {
# 
#   # print(i)
#   
#   test <- x[folds[[i]],]
#   train  <- x[-folds[[i]],]
#   
#   SB_VDFR_Class=Data2B_simpson_Classification(train, domain[-folds[[i]],], grid, nbasis=c(c1,c2),sub = 25)
#   E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
#   res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y[-folds[[i]]], family = binomial())
#   
#   probs[i,]=res_SB_VDFR$fit$fitted.values
#   aux=ROC(probs[i,],y[-folds[[i]]],plot = NULL)
#   best_pos=which.max(aux$res[,1]+aux$res[,2])
#   optimal_cut_aux[i]=aux$res[best_pos,5]
#   
#   res_clas=rep(1,length(y[folds[[i]]]))
# 
#   SB_VDFR_Class_test=Data2B_simpson_Classification(test, domain[folds[[i]],], grid, nbasis=c(c1,c2),sub = 25)
#   E_SB_VDFR_test=B2XZG_1d(SB_VDFR_Class_test$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
#   res_SB_VDFR_test=XZG2theta(X = E_SB_VDFR_test$X, Z = E_SB_VDFR_test$Z, G = E_SB_VDFR_test$G, T = E_SB_VDFR_test$T, y = y[folds[[i]]], family = binomial())
#   probs_test[i,]=res_SB_VDFR_test$fit$fitted.values
#   
#   res_clas[which(res_SB_VDFR_test$fit$fitted.values<optimal_cut_aux[i])]=0
#   
#   err[i]  = sum(res_clas!=y[folds[[i]]])
# 
#   }

# where_opt=which.min(err)
# optimal_cut=optimal_cut_aux[where_opt]

SB_VDFR_Class=Data2B_simpson_Classification(x, domain, grid, nbasis=c(c1,c2),sub = 25)
E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())

aux=ROC(res_SB_VDFR$fit$fitted.values,y,plot = NULL)
best_pos=which.max(aux$res[,1]+aux$res[,2])
optimal_cut=aux$res[best_pos,5]

res_clas=rep(1,length(y))

res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut)]=0
err_opt  = sum(res_clas!=y)
probs=res_SB_VDFR$fit$fitted.values
# probs_opt=probs[where_opt,]

# res_clas_save=rep(1,length(y))
# res_clas_save[which(res_SB_VDFR$fit$fitted.values<optimal_cut_save)]=0
# err_save  = sum(res_clas!=y)
# probs_opt=probs[where_opt,]

# SB_VDFR_Class=Data2B_simpson_Classification(x, domain, grid, nbasis=c(c1,c2),sub = 25)
# E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
# res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())
# 
# probs[i,]=res_SB_VDFR$fit$fitted.values
# aux=ROC(probs[i,],y,plot = NULL)
# best_pos=which.max(aux$res[,1]+aux$res[,2])
# optimal_cut_aux[i]=aux$res[best_pos,5]

# probs[i,,j]=res_SB_VDFR$fit$fitted.values
# aux=ROC(probs[i,,j],y,plot = NULL)
# best_pos=which.max(aux$res[,1]+aux$res[,2])
# optimal_cut_aux[i,j]=aux$res[best_pos,5]

# res_clas=rep(1,length(y))
# res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut_aux[i])]=0

# err[i]  = sum(res_clas!=y)
# err[i,j]  = sum(res_clas!=y)




#### numeric approximation of the best cutoff (deprecated because results are the same as ROC in EPI but much slower)

# while (err_old>err_new) {
# 
# 
#   for (i in 1:length(optimal_cut_grid)) {
# 
#     res_clas=rep(1,length(y))
#     res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut_grid[i])]=0
#     err[i]  = sum(res_clas!=y)
# 
#   }
# 
#   aux=which.min(err)
# 
#   optimal_cut_old=optimal_cut_grid[aux]
# 
#   err_old=err[aux]
# 
#   if (err[aux-1]<err[aux+1]) {
# 
#     optimal_cut_grid=seq(optimal_cut_grid[aux-1],optimal_cut_grid[aux],10^-(tag+1))
# 
#   }else{
#     optimal_cut_grid=seq(optimal_cut_grid[aux],optimal_cut_grid[aux+1],10^-(tag+1))
#   }
# 
#   for (i in 1:length(optimal_cut_grid)) {
# 
#     res_clas=rep(1,length(y))
#     res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut_grid[i])]=0
#     err[i]  = sum(res_clas!=y)
# 
#   }
# 
#   aux_new=which.min(err)
# 
#   optimal_cut_new=optimal_cut_grid[aux_new]
# 
#   err_new=err[aux_new]
# 
#   tag=tag+1
# 
# 
# }

out = list(probs, err_opt, optimal_cut)
names(out) = c("probs", "err_SB_VDFR","optimal_cut_SB" )
return(out)

}

SB_class_CV = function(x, y, domain, grid, c1, c2, sub=500){
  
  #x: data matrix
  #y: grouping variable
  #grid: grid where evaluate the B-splines # same as the one passed to data.gen.splines.uneven()
  #domain: data domain 
  
  N=length(y)
  
  AUC=res_clas_12=res_clas_34=prob_pred=optimal_cut_12=optimal_cut_34=err_opt_12=err_opt_34=err_opt_sample=array(dim = N)
  prob=res_clas_sample=array(dim = c(N,N-1))
  
  for (i in 1:N) {
    
    # print(c("i = ",i))
    
    SB_VDFR_Class=Data2B_simpson_Classification(x[-i,], domain[-i,], grid, nbasis=c(c1,c2),sub = sub)
    E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
    res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y[-i], family = binomial())
    
    aux=ROC(res_SB_VDFR$fit$fitted.values,y[-i],plot = NULL)

    best_pos_34=which.max(aux$res[,3]+aux$res[,4]) # which.max(aux$res[,1]+aux$res[,2])
    optimal_cut_34[i]=aux$res[best_pos_34,5]
    
    best_pos_12=which.max(aux$res[,1]+aux$res[,2]) # which.max(aux$res[,1]+aux$res[,2])
    optimal_cut_12[i]=aux$res[best_pos_12,5]
    
    
    AUC[i]=1-aux$AUC
    
    #####
    res_clas_sample[i,]=rep(1,length(y[-i]))
    res_clas_sample[i,which(res_SB_VDFR$fit$fitted.values<optimal_cut_34[i])]=0
    err_opt_sample[i]=sum(res_clas_sample[i,]!=y[-i])
    prob[i,]=res_SB_VDFR$fit$fitted.values
    #####
    
    #
    SB_VDFR_Class_CV=Data2B_simpson_Classification(x[i,], domain[i,], grid, nbasis=c(c1,c2),sub = sub) #t(as.matrix(domain[i,]))
    nu=SB_VDFR_Class_CV$B %*% res_SB_VDFR$theta
    
    prob_pred[i]= exp(nu)/(1+ exp(nu))
    
    if (prob_pred[i]<optimal_cut_34[i]) {
      
      res_clas_34[i]=0
    }else{
      res_clas_34[i]=1
    }
    
    err_opt_34[i]  = as.double(res_clas_34[i]!=y[i])
    
    if (prob_pred[i]<optimal_cut_12[i]) {
      
      res_clas_12[i]=0
    }else{
      res_clas_12[i]=1
    }
    
    err_opt_12[i]  = as.double(res_clas_12[i]!=y[i])
    
    
  }

  aux=c(sum(err_opt_12),sum(err_opt_34))
  aux_2=cbind(res_clas_12,res_clas_34)
  
  min_track=which.min(aux)
  
  aux_3=cbind(optimal_cut_12,optimal_cut_34)
  
  res_clas_final=aux_2[,min_track]
  miss_classification_error=aux[min_track]
  optimal_cut=aux_3[,min_track]
  
  out = list(prob, err_opt_sample, optimal_cut, miss_classification_error, prob_pred, res_clas_final, min_track, AUC)

  names(out) = c("prob_sample", "err_SB_VDFR_Sample", "optimal_cut", "err_SB_VDFR", "prob_pred", "clas_pred", "min_track", "AUC" )

  return(out)
  
}

all.together.splines.uneven.cv <- function(knots, mean1, err1, mean2, err2, n1, n2, p, grid, maxi, maxl,
                                           k, lambda.sm, lambda.buja, coef, unif.miss=T,
                                           rate=NA){
  
  ########################################################################################
  ##                                                                                    ##
  #     genera i dati e performa la classificazione con tutti i metodi, per *ogni*       #
  #            estensione del dominio, con parametri di smoothing ***scelti              #
  #            con relativi metodi di CV ad ogni estensione***, per un tipo              #
  #       di missing a piacere e per un numero di componenti massimo a piacere           #
  #
  #                                                                               
  #     genera los datos y realiza la clasificación con todos los métodos, para *cada* extensión del dominio,
  #     con parámetros de suavizado ***elegidos con métodos de CV relativos para cada extensión***,
  #     para un tipo de faltante a voluntad y para un máximo número de componentes a voluntad
  ##
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
  #           k: number of bases (for "James et al." method)                   (MUST BE A VECTOR ---> CV)
  #   lambda.sm: presmoothing parameter (for standard fPCA with complete data) (MUST BE A VECTOR ---> CV)
  # lambda.buja: smoothing parameter (for "Huang et al." method)               (MUST BE A VECTOR ---> CV)

   
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
    pcajames  <- pca.miss.james.cv(x$x.miss, y, maxi, t.par, k)
    pcayao    <- pca.miss.yao2.cv(x$x.miss, y, maxi, t.par)
    pcabuja   <- pca.miss.buja2.cv(x$x.miss, y, maxi, t.par, lambda.buja)
    
    
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


all.together.splines.last_kneip <- function(knots, mean1, err1, mean2, err2, n1, n2, p, grid, maxi, maxl, c1, c2, sub, c_kneip, k, lambda.sm, lambda.buja, coef, unif.miss=T, rate=NA){
  
  ########################################################################################
  ##                                                                                    ##
  #       genera i dati e performa la classificazione con tutti i metodi, per *ogni*     #
  #       estensione del dominio, con parametri di smoothing ***scelti con relativi      #
  #     metodi di CV in un dato dominio e poi tenuti fissi ad ogni estensione***, per    #
  #     un tipo di missing a piacere e per un numero di componenti massimo a piacere     #
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
  #           k: number of bases (for "James et al." method)                   (MUST BE A VECTOR ---> CV)
  #   lambda.sm: presmoothing parameter (for standard fPCA with complete data) (MUST BE A VECTOR ---> CV)
  # lambda.buja: smoothing parameter (for "Huang et al." method)               (MUST BE A VECTOR ---> CV)
  
  ### data generation
  
  x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  xbind    <- rbind(x1$data, x2$data)
  x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  y        <- c(rep(0,n1), rep(1,n2))  
  # nu=abs(rowSums(x$x.miss,na.rm=TRUE)/100) # POSSIBLE NU THAT IS DATA RELATED (BAD RESULTS IN 5 ITERATIONS)
  # y=rbinom(n1+n2,1,nu) 
  if(unif.miss==TRUE) {
    
    x=data.agg.miss(xbind)
    
  }else {
    print("Beta")
    x=data.agg.miss.betadecay(xbind,1) 
  }
  
  ############## Reconstruction of curves
  # 
  # X_list=t_list=list()
  # 
  # AUC=res_clas_12=res_clas_34=prob_pred=optimal_cut_12=optimal_cut_34=err_opt_12=err_opt_34=array(dim = N)
  # 
  # for (i in 1:dim(x$x.miss)[1]) {
  #   
  #   X_list[[i]]=x$x.miss[i,x$domain[i,1]:x$domain[i,2]]
  #   t_list[[i]]=grid[x$domain[i,1]:x$domain[i,2]]
  # }
  # 
  # reconst_result <- reconstructKneipLiebl(Ly = X_list, Lu = t_list, method = 'Error>0_AlignYES')
  # 
  # X_reconst_mat  <- t(matrix(unlist(reconst_result[['Y_reconst_list']]), ncol=100))
  # 
  # for (i in 1:length(y)) {
  #   
  #   X_aux=X_reconst_mat[-i,]
  #   y_aux=y[-i]
  #   
  #   fit_Kneip <- pfr( y_aux~ lf(X_aux, argvals = grid, bs="ps",k=c_kneip),family=binomial())
  #   
  #   aux=ROC(fit_Kneip$fitted.values,y_aux,plot = NULL)
  #   
  #   best_pos_34=which.max(aux$res[,3]+aux$res[,4]) # which.max(aux$res[,1]+aux$res[,2])
  #   optimal_cut_34[i]=aux$res[best_pos_34,5]
  #   
  #   best_pos_12=which.max(aux$res[,1]+aux$res[,2]) # which.max(aux$res[,1]+aux$res[,2])
  #   optimal_cut_12[i]=aux$res[best_pos_12,5]
  #   
  #   AUC[i]=1-aux$AUC
  #   
  #   prob_pred[i]=predict(fit_Kneip, newdata = list(X_aux=t(as.matrix(X_reconst_mat[i,]))),type='response')
  #   
  #   if (prob_pred[i]<optimal_cut_34[i]) {
  #     
  #     res_clas_34[i]=0
  #   }else{
  #     res_clas_34[i]=1
  #   }
  #   
  #   err_opt_34[i]  = as.double(res_clas_34[i]!=y[i])
  #   
  #   if (prob_pred[i]<optimal_cut_12[i]) {
  #     
  #     res_clas_12[i]=0
  #   }else{
  #     res_clas_12[i]=1
  #   }
  #   
  #   err_opt_12[i]  = as.double(res_clas_12[i]!=y[i])
  #   
  #   
  # }
  # 
  # aux=c(sum(err_opt_12),sum(err_opt_34))
  # aux_2=cbind(res_clas_12,res_clas_34)
  # 
  # min_track=which.min(aux)
  # 
  # aux_3=cbind(optimal_cut_12,optimal_cut_34)
  # 
  # res_clas_final=aux_2[,min_track]
  # miss_classification_error=aux[min_track]
  # optimal_cut=aux_3[,min_track]
  # 
  
  #######################
  
  
  t.par  = 51:100
  # 
  # pcanomiss_lda <- pca.nomiss.class2_lda(xbind, y, maxi, t.par, lambda.sm)
  # pcajames_lda  <- pca.miss.james.cv_lda(x$x.miss, y, maxi, t.par, k)
  # pcayao_lda    <- pca.miss.yao2.cv_lda(x$x.miss, y, maxi, t.par)
  # pcabuja_lda   <- pca.miss.buja2.cv_lda(x$x.miss, y, maxi, t.par, lambda.buja)
  # 
  # opt.basis_lda <- pcajames_lda$opt.basis
  # opt.bwm_lda   <- pcayao_lda$opt.bwm
  # opt.bwc_lda   <- pcayao_lda$opt.bwc
  # opt.lam_lda   <- pcabuja_lda$opt.lambda
  # 
  # err.cv.t_lda     <- c()
  # err.cv.james_lda <- c()
  # err.cv.yao_lda   <- c()
  # err.cv.buja_lda  <- c()
  # 
  # eigenfun.t_lda     <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.james_lda <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.yao_lda   <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.buja_lda  <- array(dim=c(p, maxi, maxl+1))
  # 
  # scores.t_lda     <- array(dim=c(n, maxi, maxl+1))
  # scores.james_lda <- array(dim=c(n, maxi, maxl+1))
  # scores.yao_lda   <- array(dim=c(n, maxi, maxl+1))
  # scores.buja_lda  <- array(dim=c(n, maxi, maxl+1))
  # 
  # best.t_lda     <- c()
  # best.james_lda <- c()
  # best.buja_lda  <- c()
  # best.yao_lda   <- c()
  # 
  ####
  
  pcanomiss <- pca.nomiss.class2(xbind, y, maxi, t.par, lambda.sm)
  pcajames  <- pca.miss.james.cv(x$x.miss, y, maxi, t.par, k)
  pcayao    <- pca.miss.yao2.cv(x$x.miss, y, maxi, t.par)
  pcabuja   <- pca.miss.buja2.cv(x$x.miss, y, maxi, t.par, lambda.buja)
  
  opt.basis <- pcajames$opt.basis
  opt.bwm   <- pcayao$opt.bwm
  opt.bwc   <- pcayao$opt.bwc
  opt.lam   <- pcabuja$opt.lambda
  
  
  err.cv.t     <- c()
  err.cv.james <- c()
  err.cv.yao   <- c()
  err.cv.buja  <- c()
  err_SB_VDFR = c()
  
  # prob_SB_VDFR = array(dim=c(length(y), maxl+1))
  clas_pred=optimal_cut_SB_VDFR=AUC_SB_VDFR = array(dim=c(length(y),maxl+1))
  prob_SB_VDFR_sample = array(dim=c(length(y),length(y)-1, maxl+1))
  prob_SB_VDFR_pred = array(dim=c(length(y),maxl+1))
  min_track=array(dim=c(maxl+1))
  
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
  
  domain=matrix(nrow = length(y),ncol=2)
  
  for(l in 0:maxl){
    
    print(c("l = ", l))
    
    domain_original=x$domain
    
    t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
    
    pcanomiss <- pca.nomiss.class2(x.smbind, y, maxi, t.par, lambda.sm)
    pcajames  <- pca.miss.james(x$x.miss, y, maxi, t.par, opt.basis)
    pcayao    <- pca.miss.yao2(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
    pcabuja   <- pca.miss.buja3(x$x.miss, y, maxi, t.par, opt.lam)
    # 
    # pcanomiss_lda <- pca.nomiss.class2_lda(x.smbind, y, maxi, t.par, lambda.sm)
    # pcajames_lda  <- pca.miss.james_lda(x$x.miss, y, maxi, t.par, opt.basis)
    # pcayao_lda    <- pca.miss.yao2_lda(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
    # pcabuja_lda   <- pca.miss.buja3_lda(x$x.miss, y, maxi, t.par, opt.lam)
    
    for (ind in 1:length(y)) {
      
      if (domain_original[ind,1]>t.par[1]) {
        domain[ind,1]=domain_original[ind,1]-t.par[1]+1
      }else{domain[ind,1]=1}
      
      if (domain_original[ind,2]<t.par[length(t.par)]) {
        domain[ind,2]=domain_original[ind,2]-t.par[length(t.par)]+length(t.par)
      }else{domain[ind,2]=length(t.par)}
      
    }
    
    res=SB_class_CV(x$x.miss[,t.par], y, domain, grid[t.par], c1, c2, sub=sub) # new expriment with grid[t.par], before:grid
    
    clas_pred[,l+1]=res$clas_pred
    err_SB_VDFR[l+1]  = res$err_SB_VDFR
    min_track[l+1]=res$min_track
    optimal_cut_SB_VDFR[,l+1] = res$optimal_cut
    AUC_SB_VDFR[,l+1] = res$AUC
    prob_SB_VDFR_sample[,,l+1] =res$prob_sample
    prob_SB_VDFR_pred[,l+1] =res$prob_pred
    
    
    # SB_VDFR_Class=Data2B_simpson_Classification(x$x.miss[,t.par], domain, grid, nbasis=c(c1,c2),sub = 1000)
    # E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
    # res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())
    # 
    # probs=res_SB_VDFR$fit$fitted.values
    # 
    # aux=ROC(probs,y,plot = NULL)
    # 
    # best_pos=which.max(aux$res[,1]+aux$res[,2])
    # 
    # optimal_cut=aux$res[best_pos,5]
    
    # pred=prediction(probs,y)
    # cost=performance(aux_pred,"err")
    # optimal_cut=pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
    
    # res_clas=rep(1,length(y))
    # res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut)]=0
    
    ### storing
    
    # err_SB_VDFR[l+1]  = sum(res_clas!=y)
    # 
    # optimal_cut_SB_VDFR[l+1] = optimal_cut
    # prob_SB_VDFR[,l+1] =probs
    
    err.cv.t[l+1]     = pcanomiss$err
    err.cv.james[l+1] = pcajames$err
    err.cv.yao[l+1]   = pcayao$err
    err.cv.buja[l+1]  = pcabuja$err
    
    # err.cv.t_lda[l+1]     = pcanomiss_lda$err
    # err.cv.james_lda[l+1] = pcajames_lda$err
    # err.cv.yao_lda[l+1]   = pcayao_lda$err
    # err.cv.buja_lda[l+1]  = pcabuja_lda$err
    
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
    
    #####
    
    # eigenfun.t_lda[t.par,,l+1]     = pcanomiss_lda$eigenfunctions
    # eigenfun.james_lda[t.par,,l+1] = pcajames_lda$eigenfunctions
    # eigenfun.yao_lda[t.par,,l+1]   = pcayao_lda$eigenfunctions
    # eigenfun.buja_lda[t.par,,l+1]  = pcabuja_lda$eigenfunctions
    # 
    # scores.t_lda[,,l+1]     = pcanomiss_lda$scores
    # scores.james_lda[,,l+1] = pcajames_lda$scores
    # scores.yao_lda[,,l+1]   = pcayao_lda$scores
    # scores.buja_lda[,,l+1]  = pcabuja_lda$scores
    # 
    # best.t_lda[l+1]     = pcanomiss_lda$best
    # best.james_lda[l+1] = pcajames_lda$best
    # best.yao_lda[l+1]   = pcayao_lda$best
    # best.buja_lda[l+1]  = pcabuja_lda$best
    
  }
  
  
  # out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR, optimal_cut_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja,err.cv.t_lda, err.cv.james_lda, err.cv.yao_lda, err.cv.buja_lda, eigenfun.t_lda, eigenfun.james_lda, eigenfun.yao_lda, eigenfun.buja_lda, scores.t_lda, scores.james_lda, scores.yao_lda, scores.buja_lda, best.t_lda, best.james_lda, best.yao_lda, best.buja_lda)
  # names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR", "optimal_cut_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja", "err.cv.t_lda", "err.cv.james_lda", "err.cv.yao_lda", "err.cv.buja_lda", "eigenfun.t_lda","eigenfun.james_lda", "eigenfun.yao_lda", "eigenfun.buja_lda", "scores.t_lda", "scores.james_lda","scores.yao_lda", "scores.buja_lda", "best.t_lda", "best.james_lda", "best.yao_lda", "best.buja_lda")
  
  # out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR_sample, prob_SB_VDFR_pred, optimal_cut_SB_VDFR, AUC_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja,err.cv.t_lda, err.cv.james_lda, err.cv.yao_lda, err.cv.buja_lda, eigenfun.t_lda, eigenfun.james_lda, eigenfun.yao_lda, eigenfun.buja_lda, scores.t_lda, scores.james_lda, scores.yao_lda, scores.buja_lda, best.t_lda, best.james_lda, best.yao_lda, best.buja_lda)
  # names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR_sample", "prob_SB_VDFR_pred", "optimal_cut_SB_VDFR", "AUC_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja", "err.cv.t_lda", "err.cv.james_lda", "err.cv.yao_lda", "err.cv.buja_lda", "eigenfun.t_lda","eigenfun.james_lda", "eigenfun.yao_lda", "eigenfun.buja_lda", "scores.t_lda", "scores.james_lda","scores.yao_lda", "scores.buja_lda", "best.t_lda", "best.james_lda", "best.yao_lda", "best.buja_lda")
  
  out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR_sample, prob_SB_VDFR_pred, optimal_cut_SB_VDFR, min_track, clas_pred, AUC_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja)
  names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR_sample", "prob_SB_VDFR_pred", "optimal_cut_SB_VDFR", "min_track", "clas_pred", "AUC_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja")
  
  return(out)
}

all.together.splines.last <- function(knots, mean1, err1, mean2, err2, n1, n2, p, grid, maxi, maxl, c1, c2, sub, k, lambda.sm, lambda.buja, coef, unif.miss=T, rate=NA){
  
  ########################################################################################
  ##                                                                                    ##
  #       genera i dati e performa la classificazione con tutti i metodi, per *ogni*     #
  #       estensione del dominio, con parametri di smoothing ***scelti con relativi      #
  #     metodi di CV in un dato dominio e poi tenuti fissi ad ogni estensione***, per    #
  #     un tipo di missing a piacere e per un numero di componenti massimo a piacere     #
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
  #           k: number of bases (for "James et al." method)                   (MUST BE A VECTOR ---> CV)
  #   lambda.sm: presmoothing parameter (for standard fPCA with complete data) (MUST BE A VECTOR ---> CV)
  # lambda.buja: smoothing parameter (for "Huang et al." method)               (MUST BE A VECTOR ---> CV)
  
  ### data generation
  
  x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  xbind    <- rbind(x1$data, x2$data)
  x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  y        <- c(rep(0,n1), rep(1,n2))  
  # nu=abs(rowSums(x$x.miss,na.rm=TRUE)/100) # POSSIBLE NU THAT IS DATA RELATED (BAD RESULTS IN 5 ITERATIONS)
  # y=rbinom(n1+n2,1,nu) 
  if(unif.miss==TRUE) {
    
    x=data.agg.miss(xbind)
    
    }else {
    print("Beta")
    x=data.agg.miss.betadecay(xbind,1) 
    }
  
  #######################
  
  
  t.par  = 51:100
  # 
  # pcanomiss_lda <- pca.nomiss.class2_lda(xbind, y, maxi, t.par, lambda.sm)
  # pcajames_lda  <- pca.miss.james.cv_lda(x$x.miss, y, maxi, t.par, k)
  # pcayao_lda    <- pca.miss.yao2.cv_lda(x$x.miss, y, maxi, t.par)
  # pcabuja_lda   <- pca.miss.buja2.cv_lda(x$x.miss, y, maxi, t.par, lambda.buja)
  # 
  # opt.basis_lda <- pcajames_lda$opt.basis
  # opt.bwm_lda   <- pcayao_lda$opt.bwm
  # opt.bwc_lda   <- pcayao_lda$opt.bwc
  # opt.lam_lda   <- pcabuja_lda$opt.lambda
  # 
  # err.cv.t_lda     <- c()
  # err.cv.james_lda <- c()
  # err.cv.yao_lda   <- c()
  # err.cv.buja_lda  <- c()
  # 
  # eigenfun.t_lda     <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.james_lda <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.yao_lda   <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.buja_lda  <- array(dim=c(p, maxi, maxl+1))
  # 
  # scores.t_lda     <- array(dim=c(n, maxi, maxl+1))
  # scores.james_lda <- array(dim=c(n, maxi, maxl+1))
  # scores.yao_lda   <- array(dim=c(n, maxi, maxl+1))
  # scores.buja_lda  <- array(dim=c(n, maxi, maxl+1))
  # 
  # best.t_lda     <- c()
  # best.james_lda <- c()
  # best.buja_lda  <- c()
  # best.yao_lda   <- c()
  # 
  ####
  
  pcanomiss <- pca.nomiss.class2(xbind, y, maxi, t.par, lambda.sm)
  pcajames  <- pca.miss.james.cv(x$x.miss, y, maxi, t.par, k)
  pcayao    <- pca.miss.yao2.cv(x$x.miss, y, maxi, t.par)
  pcabuja   <- pca.miss.buja2.cv(x$x.miss, y, maxi, t.par, lambda.buja)
  
  opt.basis <- pcajames$opt.basis
  opt.bwm   <- pcayao$opt.bwm
  opt.bwc   <- pcayao$opt.bwc
  opt.lam   <- pcabuja$opt.lambda
  
  
  err.cv.t     <- c()
  err.cv.james <- c()
  err.cv.yao   <- c()
  err.cv.buja  <- c()
  err_SB_VDFR = c()
  
  # prob_SB_VDFR = array(dim=c(length(y), maxl+1))
  clas_pred=optimal_cut_SB_VDFR=AUC_SB_VDFR = array(dim=c(length(y),maxl+1))
  prob_SB_VDFR_sample = array(dim=c(length(y),length(y)-1, maxl+1))
  prob_SB_VDFR_pred = array(dim=c(length(y),maxl+1))
  min_track=array(dim=c(maxl+1))
  
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
  
  domain=matrix(nrow = length(y),ncol=2)
  
  for(l in 0:maxl){
    
    print(c("l = ", l))
    
    domain_original=x$domain
    
    t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
    
    pcanomiss <- pca.nomiss.class2(x.smbind, y, maxi, t.par, lambda.sm)
    pcajames  <- pca.miss.james(x$x.miss, y, maxi, t.par, opt.basis)
    pcayao    <- pca.miss.yao2(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
    pcabuja   <- pca.miss.buja3(x$x.miss, y, maxi, t.par, opt.lam)
    # 
    # pcanomiss_lda <- pca.nomiss.class2_lda(x.smbind, y, maxi, t.par, lambda.sm)
    # pcajames_lda  <- pca.miss.james_lda(x$x.miss, y, maxi, t.par, opt.basis)
    # pcayao_lda    <- pca.miss.yao2_lda(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
    # pcabuja_lda   <- pca.miss.buja3_lda(x$x.miss, y, maxi, t.par, opt.lam)
    
    for (ind in 1:length(y)) {

      if (domain_original[ind,1]>t.par[1]) {
        domain[ind,1]=domain_original[ind,1]-t.par[1]+1
      }else{domain[ind,1]=1}

      if (domain_original[ind,2]<t.par[length(t.par)]) {
        domain[ind,2]=domain_original[ind,2]-t.par[length(t.par)]+length(t.par)
      }else{domain[ind,2]=length(t.par)}
      
    }

    res=SB_class_CV(x$x.miss[,t.par], y, domain, grid[t.par], c1, c2, sub=sub) # new expriment with grid[t.par], before:grid
    
    clas_pred[,l+1]=res$clas_pred
    err_SB_VDFR[l+1]  = res$err_SB_VDFR
    min_track[l+1]=res$min_track
    optimal_cut_SB_VDFR[,l+1] = res$optimal_cut
    AUC_SB_VDFR[,l+1] = res$AUC
    prob_SB_VDFR_sample[,,l+1] =res$prob_sample
    prob_SB_VDFR_pred[,l+1] =res$prob_pred
    
    
    # SB_VDFR_Class=Data2B_simpson_Classification(x$x.miss[,t.par], domain, grid, nbasis=c(c1,c2),sub = 1000)
    # E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
    # res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())
    # 
    # probs=res_SB_VDFR$fit$fitted.values
    # 
    # aux=ROC(probs,y,plot = NULL)
    # 
    # best_pos=which.max(aux$res[,1]+aux$res[,2])
    # 
    # optimal_cut=aux$res[best_pos,5]
    
    # pred=prediction(probs,y)
    # cost=performance(aux_pred,"err")
    # optimal_cut=pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
    
    # res_clas=rep(1,length(y))
    # res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut)]=0
    
    ### storing
    
    # err_SB_VDFR[l+1]  = sum(res_clas!=y)
    # 
    # optimal_cut_SB_VDFR[l+1] = optimal_cut
    # prob_SB_VDFR[,l+1] =probs
    
    err.cv.t[l+1]     = pcanomiss$err
    err.cv.james[l+1] = pcajames$err
    err.cv.yao[l+1]   = pcayao$err
    err.cv.buja[l+1]  = pcabuja$err

    # err.cv.t_lda[l+1]     = pcanomiss_lda$err
    # err.cv.james_lda[l+1] = pcajames_lda$err
    # err.cv.yao_lda[l+1]   = pcayao_lda$err
    # err.cv.buja_lda[l+1]  = pcabuja_lda$err

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

    #####

    # eigenfun.t_lda[t.par,,l+1]     = pcanomiss_lda$eigenfunctions
    # eigenfun.james_lda[t.par,,l+1] = pcajames_lda$eigenfunctions
    # eigenfun.yao_lda[t.par,,l+1]   = pcayao_lda$eigenfunctions
    # eigenfun.buja_lda[t.par,,l+1]  = pcabuja_lda$eigenfunctions
    # 
    # scores.t_lda[,,l+1]     = pcanomiss_lda$scores
    # scores.james_lda[,,l+1] = pcajames_lda$scores
    # scores.yao_lda[,,l+1]   = pcayao_lda$scores
    # scores.buja_lda[,,l+1]  = pcabuja_lda$scores
    # 
    # best.t_lda[l+1]     = pcanomiss_lda$best
    # best.james_lda[l+1] = pcajames_lda$best
    # best.yao_lda[l+1]   = pcayao_lda$best
    # best.buja_lda[l+1]  = pcabuja_lda$best
    
  }
  
  
  # out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR, optimal_cut_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja,err.cv.t_lda, err.cv.james_lda, err.cv.yao_lda, err.cv.buja_lda, eigenfun.t_lda, eigenfun.james_lda, eigenfun.yao_lda, eigenfun.buja_lda, scores.t_lda, scores.james_lda, scores.yao_lda, scores.buja_lda, best.t_lda, best.james_lda, best.yao_lda, best.buja_lda)
  # names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR", "optimal_cut_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja", "err.cv.t_lda", "err.cv.james_lda", "err.cv.yao_lda", "err.cv.buja_lda", "eigenfun.t_lda","eigenfun.james_lda", "eigenfun.yao_lda", "eigenfun.buja_lda", "scores.t_lda", "scores.james_lda","scores.yao_lda", "scores.buja_lda", "best.t_lda", "best.james_lda", "best.yao_lda", "best.buja_lda")

  # out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR_sample, prob_SB_VDFR_pred, optimal_cut_SB_VDFR, AUC_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja,err.cv.t_lda, err.cv.james_lda, err.cv.yao_lda, err.cv.buja_lda, eigenfun.t_lda, eigenfun.james_lda, eigenfun.yao_lda, eigenfun.buja_lda, scores.t_lda, scores.james_lda, scores.yao_lda, scores.buja_lda, best.t_lda, best.james_lda, best.yao_lda, best.buja_lda)
  # names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR_sample", "prob_SB_VDFR_pred", "optimal_cut_SB_VDFR", "AUC_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja", "err.cv.t_lda", "err.cv.james_lda", "err.cv.yao_lda", "err.cv.buja_lda", "eigenfun.t_lda","eigenfun.james_lda", "eigenfun.yao_lda", "eigenfun.buja_lda", "scores.t_lda", "scores.james_lda","scores.yao_lda", "scores.buja_lda", "best.t_lda", "best.james_lda", "best.yao_lda", "best.buja_lda")

  out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR_sample, prob_SB_VDFR_pred, optimal_cut_SB_VDFR, min_track, clas_pred, AUC_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja)
  names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR_sample", "prob_SB_VDFR_pred", "optimal_cut_SB_VDFR", "min_track", "clas_pred", "AUC_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja")
  
  return(out)
}

all.together.splines.last_solo <- function(knots, mean1, err1, mean2, err2, n1, n2, p, grid, maxi, maxl, c1, c2, sub, k, lambda.sm, lambda.buja, coef, unif.miss=T, rate=NA){
  
  ########################################################################################
  ##                                                                                    ##
  #       genera i dati e performa la classificazione con tutti i metodi, per *ogni*     #
  #       estensione del dominio, con parametri di smoothing ***scelti con relativi      #
  #     metodi di CV in un dato dominio e poi tenuti fissi ad ogni estensione***, per    #
  #     un tipo di missing a piacere e per un numero di componenti massimo a piacere     #
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
  #           k: number of bases (for "James et al." method)                   (MUST BE A VECTOR ---> CV)
  #   lambda.sm: presmoothing parameter (for standard fPCA with complete data) (MUST BE A VECTOR ---> CV)
  # lambda.buja: smoothing parameter (for "Huang et al." method)               (MUST BE A VECTOR ---> CV)
  
  ### data generation
  
  x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  xbind    <- rbind(x1$data, x2$data)
  x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  y        <- c(rep(0,n1), rep(1,n2))  
  # nu=abs(rowSums(x$x.miss,na.rm=TRUE)/100) # POSSIBLE NU THAT IS DATA RELATED (BAD RESULTS IN 5 ITERATIONS)
  # y=rbinom(n1+n2,1,nu) 
  if(unif.miss==TRUE) {
    
    x=data.agg.miss(xbind)
    
  }else {
    print("Beta")
    x=data.agg.miss.betadecay(xbind,1) 
  }
  
  
  t.par  = 51:100
  # 
  # pcanomiss_lda <- pca.nomiss.class2_lda(xbind, y, maxi, t.par, lambda.sm)
  # pcajames_lda  <- pca.miss.james.cv_lda(x$x.miss, y, maxi, t.par, k)
  # pcayao_lda    <- pca.miss.yao2.cv_lda(x$x.miss, y, maxi, t.par)
  # pcabuja_lda   <- pca.miss.buja2.cv_lda(x$x.miss, y, maxi, t.par, lambda.buja)
  # 
  # opt.basis_lda <- pcajames_lda$opt.basis
  # opt.bwm_lda   <- pcayao_lda$opt.bwm
  # opt.bwc_lda   <- pcayao_lda$opt.bwc
  # opt.lam_lda   <- pcabuja_lda$opt.lambda
  # 
  # err.cv.t_lda     <- c()
  # err.cv.james_lda <- c()
  # err.cv.yao_lda   <- c()
  # err.cv.buja_lda  <- c()
  # 
  # eigenfun.t_lda     <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.james_lda <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.yao_lda   <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.buja_lda  <- array(dim=c(p, maxi, maxl+1))
  # 
  # scores.t_lda     <- array(dim=c(n, maxi, maxl+1))
  # scores.james_lda <- array(dim=c(n, maxi, maxl+1))
  # scores.yao_lda   <- array(dim=c(n, maxi, maxl+1))
  # scores.buja_lda  <- array(dim=c(n, maxi, maxl+1))
  # 
  # best.t_lda     <- c()
  # best.james_lda <- c()
  # best.buja_lda  <- c()
  # best.yao_lda   <- c()
  # 
  ####
  
  # pcanomiss <- pca.nomiss.class2(xbind, y, maxi, t.par, lambda.sm)
  # pcajames  <- pca.miss.james.cv(x$x.miss, y, maxi, t.par, k)
  # pcayao    <- pca.miss.yao2.cv(x$x.miss, y, maxi, t.par)
  # pcabuja   <- pca.miss.buja2.cv(x$x.miss, y, maxi, t.par, lambda.buja)
  # 
  # opt.basis <- pcajames$opt.basis
  # opt.bwm   <- pcayao$opt.bwm
  # opt.bwc   <- pcayao$opt.bwc
  # opt.lam   <- pcabuja$opt.lambda
  # 
  # 
  # err.cv.t     <- c()
  # err.cv.james <- c()
  # err.cv.yao   <- c()
  # err.cv.buja  <- c()
  err_SB_VDFR = c()
  
  # prob_SB_VDFR = array(dim=c(length(y), maxl+1))
  clas_pred=optimal_cut_SB_VDFR=AUC_SB_VDFR = array(dim=c(length(y),maxl+1))
  prob_SB_VDFR_sample = array(dim=c(length(y),length(y)-1, maxl+1))
  prob_SB_VDFR_pred = array(dim=c(length(y),maxl+1))
  min_track=array(dim=c(maxl+1))
  
  # eigenfun.t     <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.james <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.yao   <- array(dim=c(p, maxi, maxl+1))
  # eigenfun.buja  <- array(dim=c(p, maxi, maxl+1))
  # 
  # scores.t     <- array(dim=c(n, maxi, maxl+1))
  # scores.james <- array(dim=c(n, maxi, maxl+1))
  # scores.yao   <- array(dim=c(n, maxi, maxl+1))
  # scores.buja  <- array(dim=c(n, maxi, maxl+1))
  # 
  # best.t     <- c()
  # best.james <- c()
  # best.buja  <- c()
  # best.yao   <- c()
  
  ### data analysis
  
  domain=matrix(nrow = length(y),ncol=2)
  
  for(l in 0:maxl){
    
    print(c("l = ", l))
    
    domain_original=x$domain
    
    t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
    
    # pcanomiss <- pca.nomiss.class2(x.smbind, y, maxi, t.par, lambda.sm)
    # pcajames  <- pca.miss.james(x$x.miss, y, maxi, t.par, opt.basis)
    # pcayao    <- pca.miss.yao2(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
    # pcabuja   <- pca.miss.buja3(x$x.miss, y, maxi, t.par, opt.lam)
    # 
    # pcanomiss_lda <- pca.nomiss.class2_lda(x.smbind, y, maxi, t.par, lambda.sm)
    # pcajames_lda  <- pca.miss.james_lda(x$x.miss, y, maxi, t.par, opt.basis)
    # pcayao_lda    <- pca.miss.yao2_lda(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
    # pcabuja_lda   <- pca.miss.buja3_lda(x$x.miss, y, maxi, t.par, opt.lam)
    
    for (ind in 1:length(y)) {
      
      if (domain_original[ind,1]>t.par[1]) {
        domain[ind,1]=domain_original[ind,1]-t.par[1]+1
      }else{domain[ind,1]=1}
      
      if (domain_original[ind,2]<t.par[length(t.par)]) {
        domain[ind,2]=domain_original[ind,2]-t.par[length(t.par)]+length(t.par)
      }else{domain[ind,2]=length(t.par)}
      
    }
    
    res=SB_class_CV(x$x.miss[,t.par], y, domain, grid[t.par], c1, c2, sub=sub) # new expriment with grid[t.par], before:grid
    
    clas_pred[,l+1]=res$clas_pred
    err_SB_VDFR[l+1]  = res$err_SB_VDFR
    min_track[l+1]=res$min_track
    optimal_cut_SB_VDFR[,l+1] = res$optimal_cut
    AUC_SB_VDFR[,l+1] = res$AUC
    prob_SB_VDFR_sample[,,l+1] =res$prob_sample
    prob_SB_VDFR_pred[,l+1] =res$prob_pred
    
    # SB_VDFR_Class=Data2B_simpson_Classification(x$x.miss[,t.par], domain, grid, nbasis=c(c1,c2),sub = 1000)
    # E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
    # res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())
    # 
    # probs=res_SB_VDFR$fit$fitted.values
    # 
    # aux=ROC(probs,y,plot = NULL)
    # 
    # best_pos=which.max(aux$res[,1]+aux$res[,2])
    # 
    # optimal_cut=aux$res[best_pos,5]
    
    # pred=prediction(probs,y)
    # cost=performance(aux_pred,"err")
    # optimal_cut=pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
    
    # res_clas=rep(1,length(y))
    # res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut)]=0
    
    ### storing
    
    # err.cv.t[l+1]     = pcanomiss$err
    # err.cv.james[l+1] = pcajames$err
    # err.cv.yao[l+1]   = pcayao$err
    # err.cv.buja[l+1]  = pcabuja$err
    
    # err.cv.t_lda[l+1]     = pcanomiss_lda$err
    # err.cv.james_lda[l+1] = pcajames_lda$err
    # err.cv.yao_lda[l+1]   = pcayao_lda$err
    # err.cv.buja_lda[l+1]  = pcabuja_lda$err
    
    # eigenfun.t[t.par,,l+1]     = pcanomiss$eigenfunctions
    # eigenfun.james[t.par,,l+1] = pcajames$eigenfunctions
    # eigenfun.yao[t.par,,l+1]   = pcayao$eigenfunctions
    # eigenfun.buja[t.par,,l+1]  = pcabuja$eigenfunctions
    # 
    # scores.t[,,l+1]     = pcanomiss$scores
    # scores.james[,,l+1] = pcajames$scores
    # scores.yao[,,l+1]   = pcayao$scores
    # scores.buja[,,l+1]  = pcabuja$scores
    # 
    # best.t[l+1]     = pcanomiss$best
    # best.james[l+1] = pcajames$best
    # best.yao[l+1]   = pcayao$best
    # best.buja[l+1]  = pcabuja$best
    
    #####
    
    # eigenfun.t_lda[t.par,,l+1]     = pcanomiss_lda$eigenfunctions
    # eigenfun.james_lda[t.par,,l+1] = pcajames_lda$eigenfunctions
    # eigenfun.yao_lda[t.par,,l+1]   = pcayao_lda$eigenfunctions
    # eigenfun.buja_lda[t.par,,l+1]  = pcabuja_lda$eigenfunctions
    # 
    # scores.t_lda[,,l+1]     = pcanomiss_lda$scores
    # scores.james_lda[,,l+1] = pcajames_lda$scores
    # scores.yao_lda[,,l+1]   = pcayao_lda$scores
    # scores.buja_lda[,,l+1]  = pcabuja_lda$scores
    # 
    # best.t_lda[l+1]     = pcanomiss_lda$best
    # best.james_lda[l+1] = pcajames_lda$best
    # best.yao_lda[l+1]   = pcayao_lda$best
    # best.buja_lda[l+1]  = pcabuja_lda$best
    
  }
  
 
  out = list(err_SB_VDFR, prob_SB_VDFR_sample, prob_SB_VDFR_pred, optimal_cut_SB_VDFR, min_track, clas_pred, AUC_SB_VDFR)
  names(out) = c("err_SB_VDFR", "prob_SB_VDFR_sample", "prob_SB_VDFR_pred", "optimal_cut_SB_VDFR", "min_track", "clas_pred", "AUC_SB_VDFR")

  # out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR_sample, prob_SB_VDFR_pred, optimal_cut_SB_VDFR, min_track, clas_pred, AUC_SB_VDFR, eigenfun.t, eigenfun.james,eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja, best.t, best.james, best.yao, best.buja)
  # names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR_sample", "prob_SB_VDFR_pred", "optimal_cut_SB_VDFR", "min_track", "clas_pred", "AUC_SB_VDFR", "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja")
  
  return(out)
}

# all.together.splines.last <- function(knots, mean1, err1, mean2, err2, n1, n2, p, grid, c1, c2, maxi,# maxl,
#                                       k, lambda.sm, lambda.buja, coef, unif.miss=T,
#                                       rate=NA){
#   
#   ########################################################################################
#   ##                                                                                    ##
#   #       genera i dati e performa la classificazione con tutti i metodi, per *ogni*     #
#   #       estensione del dominio, con parametri di smoothing ***scelti con relativi      #
#   #     metodi di CV in un dato dominio e poi tenuti fissi ad ogni estensione***, per    #
#   #     un tipo di missing a piacere e per un numero di componenti massimo a piacere
#   #
#   #
#   #     genera los datos y realiza la clasificación con todos los métodos, para *cada*
#   #     extensión del dominio, con parámetros de suavizado ***elegidos con métodos CV relativos en un dominio dado
#   #     y luego mantenidos fijos en cada extensión***, para un tipo de falta a voluntad y para un número máximo de componentes
#   #     a voluntad
#   #
#   ##                                                                                    ##
#   ########################################################################################
#   
#   #       knots: vector of knots placement
#   #       mean1: mean vector for the coefficients, group 1
#   #        err1: standard error of measurement error, group 1
#   #       mean2: mean vector for the coefficients, group 2
#   #        err2: standard error of measurement error, group 2
#   #        coef: variance vector for the coefficients (the covariance matrix is supposed to be diagonal)
#   #          n1: number of subjects, group 1
#   #          n2: number of subjects, group 2
#   #           p: number of observation points
#   #        grid: evaluation grid
#   #        maxi: maximum number of components
#   #        maxl: maximum number of extensions
#   #           k: number of bases (for "James et al." method)                   (MUST BE A VECTOR ---> CV)
#   #   lambda.sm: presmoothing parameter (for standard fPCA with complete data) (MUST BE A VECTOR ---> CV)
#   # lambda.buja: smoothing parameter (for "Huang et al." method)               (MUST BE A VECTOR ---> CV)
#   
#   ### data generation
#   
#   x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
#   x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
#   xbind    <- rbind(x1$data, x2$data)
#   x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
#   y        <- c(rep(0,n1), rep(1,n2))
#   x        <- if(unif.miss==T) data.agg.miss(xbind) else data.agg.miss.betadecay(xbind)
#   
#   t.par  = 51:100 # MAYBE range(x$timepoints)
#   
#   pcanomiss <- pca.nomiss.class2(xbind, y, maxi, t.par, lambda.sm)
#   pcajames  <- pca.miss.james.cv(x$x.miss, y, maxi, t.par, k)
#   pcayao    <- pca.miss.yao2.cv(x$x.miss, y, maxi, t.par)
#   pcabuja   <- pca.miss.buja2.cv(x$x.miss, y, maxi, t.par, lambda.buja)
#   
#   # res_SB=SB_class(x$x.miss, y, x$domain, grid, c=c(c1,c2))
#     
#   
#   opt.basis <- pcajames$opt.basis
#   opt.bwm   <- pcayao$opt.bwm
#   opt.bwc   <- pcayao$opt.bwc
#   opt.lam   <- pcabuja$opt.lambda
#   
#   
#   # err.cv.t     <- c()
#   # err.cv.james <- c()
#   # err.cv.yao   <- c()
#   # err.cv.buja  <- c()
#   # 
#   # eigenfun.t     <- array(dim=c(p, maxi))
#   # eigenfun.james <- array(dim=c(p, maxi))
#   # eigenfun.yao   <- array(dim=c(p, maxi))
#   # eigenfun.buja  <- array(dim=c(p, maxi))
#   # 
#   # scores.t     <- array(dim=c(n, maxi))
#   # scores.james <- array(dim=c(n, maxi))
#   # scores.yao   <- array(dim=c(n, maxi))
#   # scores.buja  <- array(dim=c(n, maxi))
#   # 
#   # best.t     <- c()
#   # best.james <- c()
#   # best.buja  <- c()
#   # best.yao   <- c()
#   
#   ### data analysis
#     
#     t.par = 1:p  #round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
#     
#     pcanomiss <- pca.nomiss.class2(x.smbind, y, maxi, t.par, lambda.sm) # this does not work with missing data hence recive x.smbind (why smoothed data?)
#     pcajames  <- pca.miss.james(x$x.miss, y, maxi, t.par, opt.basis)
#     pcayao    <- pca.miss.yao2(x$x.miss, y, maxi, t.par, opt.bwm, opt.bwc)
#     pcabuja   <- pca.miss.buja3(x$x.miss, y, maxi, t.par, opt.lam)
#     
#     SB_VDFR_Class=Data2B_simpson_Classification(x$x.miss, x$domain, grid, nbasis=c(c1,c2),sub = 25)
#     E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
#     res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())
# 
#     probs=res_SB_VDFR$fit$fitted.values
# 
#     aux=ROC(probs,y,plot = NULL)
#     
#     best_pos=which.max(aux$res[,1]+aux$res[,2])
#     
#     optimal_cut=aux$res[best_pos,5]
#     
#     # pred=prediction(probs,y)
#     # cost=performance(aux_pred,"err")
#     # optimal_cut=pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
# 
#     res_clas=rep(1,length(y))
#     res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut)]=0
#     
# 
#     # aux=confusionMatrix(as.factor(res_clas), as.factor(y), positive = "1")
#     
#     ### storing
#     
#     err.cv.t     = pcanomiss$err
#     err.cv.james = pcajames$err
#     err.cv.yao   = pcayao$err
#     err.cv.buja  = pcabuja$err
#     err_SB_VDFR  = sum(res_clas!=y)
#     # err_SB_VDFR  = res_SB$err_SB_VDFR
#     
#     eigenfun.t     = pcanomiss$eigenfunctions
#     eigenfun.james = pcajames$eigenfunctions
#     eigenfun.yao   = pcayao$eigenfunctions
#     eigenfun.buja  = pcabuja$eigenfunctions
#     
#     scores.t     = pcanomiss$scores
#     scores.james = pcajames$scores
#     scores.yao   = pcayao$scores
#     scores.buja  = pcabuja$scores
#     
#     best.t     = pcanomiss$best
#     best.james = pcajames$best
#     best.yao   = pcayao$best
#     best.buja  = pcabuja$best
#     
#     optimal_cut_SB_VDFR = optimal_cut
#     prob_SB_VDFR =probs
# 
#     # optimal_cut_SB_VDFR = res_SB$optimal_cut_SB
#     # prob_SB_VDFR =res_SB$probs
#   
#   out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err_SB_VDFR, prob_SB_VDFR, optimal_cut_SB_VDFR, eigenfun.t, eigenfun.james,
#              eigenfun.yao, eigenfun.buja, scores.t, scores.james, scores.yao, scores.buja,
#              best.t, best.james, best.yao, best.buja)
#   names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err_SB_VDFR", "prob_SB_VDFR", "optimal_cut_SB_VDFR", "eigenfun.t",
#                  "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "scores.t", "scores.james",
#                  "scores.yao", "scores.buja", "best.t", "best.james", "best.yao", "best.buja")
#   return(out)
# }




all.together.splines.last.totalcv <- function(knots, mean1, err1, mean2, err2, n1, n2, p,
                                              grid, maxi, maxl, k, lambda.buja, coef, 
                                              unif.miss=T, rate=NA){
  
  ########################################################################################
  ##                                                                                    ##
  #         genera i dati e performa la classificazione con tutti i metodi,              #
  #         per *ogni* estensione del dominio, con parametri di smoothing                #
  #   ***scelti con relativi metodi di CV in in ogni estensione del dominio***, per      #
  #    un tipo di missing a piacere e per un numero di componenti massimo a piacere      #
  ##
  #     genera los datos y realiza la clasificación con todos los métodos, para *cada* extensión del dominio,
  #     con parámetros de suavizado ***elegidos con métodos de CV relativos en cada extensión del dominio***,
  #     para un tipo de missing a voluntad y para un número máximo de componentes a voluntad
  #
  ########################################################################################
  
  ### data generation
  
  x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  xbind    <- rbind(x1$data, x2$data)
  x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  y        <- c(rep(0,n1), rep(1,n2))  
  x        <- if(unif.miss==T) data.agg.miss(xbind) else data.agg.miss.betadecay(xbind, rate) 
  
  
  ### empty arrays for results
  
  
  err.cv.t     <- c()
  err.cv.james <- c()
  err.cv.yao   <- c()
  err.cv.buja  <- c()
  err.cv.buja.gcv  <- c()
  err.cv.kraus <- c()
  
  eigenfun.t     <- array(dim=c(p, maxi, maxl+1))
  eigenfun.james <- array(dim=c(p, maxi, maxl+1))
  eigenfun.yao   <- array(dim=c(p, maxi, maxl+1))
  eigenfun.buja  <- array(dim=c(p, maxi, maxl+1))
  eigenfun.buja.gcv  <- array(dim=c(p, maxi, maxl+1))
  eigenfun.kraus <- array(dim=c(p, maxi, maxl+1))
  
  scores.t     <- array(dim=c(n, maxi, maxl+1))
  scores.james <- array(dim=c(n, maxi, maxl+1))
  scores.yao   <- array(dim=c(n, maxi, maxl+1))
  scores.buja  <- array(dim=c(n, maxi, maxl+1))
  scores.buja.gcv  <- array(dim=c(n, maxi, maxl+1))
  scores.kraus <- array(dim=c(n, maxi, maxl+1))
  
  best.t     <- c()
  best.james <- c()
  best.yao   <- c()
  best.buja  <- c()
  best.buja.gcv  <- c()
  best.kraus <- c()
  
  param.james <- c()
  param.yao.m <- c()
  param.yao.c <- c()  
  param.buja  <- matrix(ncol=maxi, nrow=maxl+1)

  
  ### data analysis
  
  for(l in 0:maxl){
    
    t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
    
    pcanomiss <- pca.nomiss.class2(x.smbind, y, maxi, t.par, lambda.sm=1e-8)
    pcajames  <- pca.miss.james.cv(x$x.miss, y, maxi, t.par, k)
    pcayao    <- pca.miss.yao2.cv(x$x.miss, y, maxi, t.par)
    pcabuja   <- pca.miss.buja2.cv(x$x.miss, y, maxi, t.par, lambda.buja)
    pcabujagcv <- pca.miss.buja.gcv(x$x.miss, y, maxi, t.par, max.rs)
    pcakraus  <- pca.miss.kraus(x$x.miss, y, maxi, t.par, max.rs, sd.rs)
    ### storing
    
    err.cv.t[l+1]     = pcanomiss$err
    err.cv.james[l+1] = pcajames$err
    err.cv.yao[l+1]   = pcayao$err
    err.cv.buja[l+1]  = pcabuja$err
    err.cv.buja.gcv[l+1]  = pcabujagcv$err
    err.cv.kraus[l+1] = pcakraus$err
    
    eigenfun.t[t.par,,l+1]     = pcanomiss$eigenfunctions
    eigenfun.james[t.par,,l+1] = pcajames$eigenfunctions
    eigenfun.yao[t.par,,l+1]   = pcayao$eigenfunctions
    eigenfun.buja[t.par,,l+1]  = pcabuja$eigenfunctions
    eigenfun.buja.gcv[t.par,,l+1]  = pcabujagcv$eigenfunctions
    eigenfun.kraus[t.par,,l+1] = pcakraus$eigenfunctions
    
    scores.t[,,l+1]     = pcanomiss$scores
    scores.james[,,l+1] = pcajames$scores
    scores.yao[,,l+1]   = pcayao$scores
    scores.buja[,,l+1]  = pcabuja$scores
    scores.buja.gcv[,,l+1]  = pcabujagcv$scores
    scores.kraus[,,l+1] = pcakraus$scores
    
    best.t[l+1]     = pcanomiss$best
    best.james[l+1] = pcajames$best
    best.yao[l+1]   = pcayao$best
    best.buja[l+1]  = pcabuja$best
    best.buja.gcv[l+1]  = pcabujagcv$best
    best.kraus[l+1] = pcakraus$best
    
    param.james[l+1] = pcajames$opt.basis
    param.yao.m[l+1] = pcayao$opt.bwm
    param.yao.c[l+1] = pcayao$opt.bwc
    param.buja[l+1,] = pcabuja$opt.lambda
    
    }
  
  out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, err.cv.buja.gcv, err.cv.kraus, eigenfun.t, 
             eigenfun.james, eigenfun.yao, eigenfun.buja, eigenfun.buja.gcv,eigenfun.kraus, scores.t, 
             scores.james, scores.yao, scores.buja, scores.buja.gcv, scores.kraus, best.t, best.james, 
             best.yao, best.buja, best.buja.gcv, best.kraus, param.james, param.yao.m, param.yao.c, param.buja)
  names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja", "err.cv.buja.gcv","err.cv.kraus",
                 "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja", "eigenfun.buja.gcv", "eigenfun.kraus",
                 "scores.t", "scores.james", "scores.yao", "scores.buja", "scores.buja.gcv", "scores.kraus", "best.t",
                 "best.james", "best.yao", "best.buja", "best.buja.gcv", "best.kraus", "best.basis", "best.bwmean",
                 "best.bwcov", "best.lambda")
  return(out)
}


