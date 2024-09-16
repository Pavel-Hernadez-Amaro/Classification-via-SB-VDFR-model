# FUNCTION FOR USING SIMIULATION_2

SB_class_cv_2 = function(x, y, domain){

  #x: data matrix
  #y: grouping variable
  #domain: data domain 
  
  # IN THE SIMULATION STUDY 2 WE CREATE THE CURVES SO WE DO NOT NEED THE GRID HENCE WE WILL USE Data2B_simpson_no_vc instead of the classification version

  SB_VDFR_Class=Data2B_simpson_Classification_2(x, domain, nbasis=c(30,30),sub = 25)
  E_SB_VDFR=B2XZG_1d(SB_VDFR_Class$B,c=c(30)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
  res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = binomial())
  
  prod_SB=SB_VDFR_Class$Phi$B
  Beta=as.vector(prod_SB %*% res_SB_VDFR$theta) 
  
  aux=ROC(res_SB_VDFR$fit$fitted.values,y,plot = NULL)
  best_pos=which.max(aux$res[,1]+aux$res[,2])
  optimal_cut=aux$res[best_pos,5]
  
  res_clas=rep(1,length(y))
  
  res_clas[which(res_SB_VDFR$fit$fitted.values<optimal_cut)]=0
  err_opt  = sum(res_clas!=y)
  probs=res_SB_VDFR$fit$fitted.values

  out = list(probs, err_opt, optimal_cut, Beta)
  names(out) = c("probs", "err_SB_VDFR","optimal_cut_SB", "Beta" )
  
  return(out)
}


all.together.splines.last.totalcv_sim_2 <- function(N=100, p=150, case, maxi, maxl, k, lambda.buja,unif.miss=T, rate=NA){
  
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
  
  # x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  # x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  # xbind    <- rbind(x1$data, x2$data)
  # x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  # y        <- c(rep(0,n1), rep(1,n2))  
  # x        <- if(unif.miss==T) data.agg.miss(xbind) else data.agg.miss.betadecay(xbind, rate) 
  
  
  
  M_a =   round(runif(N-1, 1, round(p/3))) #round(runif(N,1,p-10),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
  M_b =   round(runif(N-1, round(2*p/3+1), p)) #round(runif(N,M_a+10,p),digits = 0)

  M_a=c(1,M_a)
  M_b=c(p,M_b)
  
  # c(max(M_a),min(M_b)) # THIS IS THE COMMON DOMAIN
    
  # M_a = rnegbin(N,24,10) # Para 1000 poner (240,2)
  # 
  # M_b = round(runif(N,M_a+10,J),digits = 0)
  
  if (sum(M_a<M_b)!=length(M_b)) {
    stop("M_a is bigger than M_b")
  }
  
  
  if (max(M_b)>p) {
    # print("si")
    M_b[which(M_b>p)]=p
  }
  
  
  if (min(M_a)<=1) {
    # print("si")
    # M[which(M<=10)]=round(runif(length(which(M<=10)),31,max(M)))
    M_a[which(M_a<=1)]=1
  }
  
  T=max(M_b)
  
  t=1:T
  
  o=order(M_a)
  
  M_a=M_a[o]
  
  M_b=M_b[o]
  
  M=cbind(M_a,M_b)
  
  X_se=X_se_all=matrix(NA,N,p) # NOISY
  
  X_s=X_s_all=matrix(NA,N,p) # NOT NOISY
  
  for (i in 1:N) {
    
    u=rnorm(1)
    
    temp=matrix(NA,10,p)
    
    for (k in 1:10) {
      
      v_i1=rnorm(1,0,4/k^2)
      v_i2=rnorm(1,0,4/k^2)
      
      temp[k,]=v_i1*sin(2*pi*k*(1:T)/100)+v_i2*cos(2*pi*k*(1:p)/100)
    }
    
    B=apply(temp,2,sum)
    
    B=B+u
    
    X_s_all[i,]=B
    X_se_all[i,]=B+rnorm(p,0,1) # WE ADD NOISE
    
    X_s[i,M_a[i]:M_b[i]]=X_s_all[i,M_a[i]:M_b[i]]
    X_se[i,M_a[i]:M_b[i]]=X_se_all[i,M_a[i]:M_b[i]]
    
  }
  
  ###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE
  
  Beta=array(dim = c(2,T))  
  nu_1=nu_2=nu_1_noisy=nu_2_noisy=y=rep(0,N)
  
  # TRUE FUNCTIONAL COEFFICIENTS
  
  Beta[1,1:p]=((10*t/p)-5)/100
  Beta[2,1:p]=-1*(5-40*((t/p)-0.5)^2)/100
  # Beta[3,1:T]=(5-10*((T-t)/T))/10
  #
  
  dim_beta=dim(Beta)[1]
  
  for (i in 1:N) {
    
    # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)  
    
    nu_1[i]=sum(X_s_all[i,]*Beta[1,],na.rm = 1)/(M_b[i]-M_a[i]+1) # NOISY
    nu_2[i]=sum(X_s_all[i,]*Beta[2,],na.rm = 1)/(M_b[i]-M_a[i]+1) # NOISY
    
    nu_1_noisy[i]=sum(X_se_all[i,]*Beta[1,],na.rm = 1)/(M_b[i]-M_a[i]+1) # NOISY
    nu_2_noisy[i]=sum(X_se_all[i,]*Beta[2,],na.rm = 1)/(M_b[i]-M_a[i]+1) # NOISY
    
  }
  
  
  # y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE GAUSSIAN MODEL 
  
  y_1=rbinom(N,1,exp(nu_1)/(1+exp(nu_1)))
  y_2=rbinom(N,1,exp(nu_2)/(1+exp(nu_2)))
  y_1_noisy=rbinom(N,1,exp(nu_1_noisy)/(1+exp(nu_1_noisy)))
  y_2_noisy=rbinom(N,1,exp(nu_2_noisy)/(1+exp(nu_2_noisy)))
  
  if (case==1) {
    y=y_1
    X=X_se
    X_no_miss=X_se_all
  }
  if (case==2) {
    y=y_2
    X=X_se
    X_no_miss=X_se_all
  }
  if (case==3) {
    y=y_1_noisy
    X=X_se
    X_no_miss=X_se_all
  }
  if (case==4) {
    y=y_2_noisy
    X=X_se
    X_no_miss=X_se_all
  }
  
  
  ### empty arrays for results
  
  
  err.cv.t     <- c()
  err.cv.james <- c()
  err.cv.yao   <- c()
  err.cv.buja  <- c()
  err_SB_VDFR = c()
  
  optimal_cut_SB_VDFR = c()
  prob_SB_VDFR = array(dim=c(length(y), maxl+1))
  Beta_estimated = array(dim=c(p, maxl+1))
  
  eigenfun.t     <- array(dim=c(p, maxi, maxl+1))
  eigenfun.james <- array(dim=c(p, maxi, maxl+1))
  eigenfun.yao   <- array(dim=c(p, maxi, maxl+1))
  eigenfun.buja  <- array(dim=c(p, maxi, maxl+1))

  scores.t     <- array(dim=c(N, maxi, maxl+1))
  scores.james <- array(dim=c(N, maxi, maxl+1))
  scores.yao   <- array(dim=c(N, maxi, maxl+1))
  scores.buja  <- array(dim=c(N, maxi, maxl+1))

  best.t     <- c()
  best.james <- c()
  best.yao   <- c()
  best.buja  <- c()

  param.james <- c()
  param.yao.m <- c()
  param.yao.c <- c()  
  param.buja  <- matrix(ncol=maxi, nrow=maxl+1)
  
  err.cv.t_lda     <- c()
  err.cv.james_lda <- c()
  err.cv.yao_lda   <- c()
  err.cv.buja_lda  <- c()
  
  eigenfun.t_lda     <- array(dim=c(p, maxi, maxl+1))
  eigenfun.james_lda <- array(dim=c(p, maxi, maxl+1))
  eigenfun.yao_lda   <- array(dim=c(p, maxi, maxl+1))
  eigenfun.buja_lda  <- array(dim=c(p, maxi, maxl+1))
  
  scores.t_lda     <- array(dim=c(N, maxi, maxl+1))
  scores.james_lda <- array(dim=c(N, maxi, maxl+1))
  scores.yao_lda   <- array(dim=c(N, maxi, maxl+1))
  scores.buja_lda  <- array(dim=c(N, maxi, maxl+1))
  
  best.t_lda     <- c()
  best.james_lda <- c()
  best.buja_lda  <- c()
  best.yao_lda   <- c()
  
  param.james_lda <- c()
  param.yao.m_lda <- c()
  param.yao.c_lda <- c()  
  param.buja_lda  <- matrix(ncol=maxi, nrow=maxl+1)
  
  ### data analysis
  
  for(l in 0:maxl){
    
    print(l)
    
    domain_original=M
    
    t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
    
    pcanomiss_lda <- pca.nomiss.class2_lda(X_no_miss, y, maxi, t.par, lambda.sm=1e-8)
    pcajames_lda  <- pca.miss.james.cv_lda(X, y, maxi, t.par, k)
    pcayao_lda    <- pca.miss.yao2.cv_lda(X, y, maxi, t.par)
    pcabuja_lda   <- pca.miss.buja2.cv_lda(X, y, maxi, t.par, lambda.buja)
    
    pcanomiss <- pca.nomiss.class2(X_no_miss, y, maxi, t.par, lambda.sm=1e-8)
    pcajames  <- pca.miss.james.cv(X, y, maxi, t.par, k)
    pcayao    <- pca.miss.yao2.cv(X, y, maxi, t.par)
    pcabuja   <- pca.miss.buja2.cv(X, y, maxi, t.par, lambda.buja)

    for (ind in 1:length(y)) {
      
      if (domain_original[ind,1]>t.par[1]) {
        domain[ind,1]=domain_original[ind,1]-t.par[1]+1
      }else{domain[ind,1]=1}
      
      if (domain_original[ind,2]<t.par[length(t.par)]) {
        domain[ind,2]=domain_original[ind,2]-t.par[length(t.par)]+length(t.par)
      }else{domain[ind,2]=length(t.par)}
      
    }

  # c(max(domain_original[,1]),min(domain_original[,2]))

    res=SB_class_cv_2(X[,t.par], y, domain)

### storing
    
    err_SB_VDFR[l+1]  = res$err_SB_VDFR
    
    optimal_cut_SB_VDFR[l+1] = res$optimal_cut_SB
    prob_SB_VDFR[,l+1] =res$probs
    
    Beta_estimated[t.par,l+1]=res$Beta
    
    err.cv.t_lda[l+1]     = pcanomiss_lda$err
    err.cv.james_lda[l+1] = pcajames_lda$err
    err.cv.yao_lda[l+1]   = pcayao_lda$err
    err.cv.buja_lda[l+1]  = pcabuja_lda$err
    
    eigenfun.t_lda[t.par,,l+1]     = pcanomiss_lda$eigenfunctions
    eigenfun.james_lda[t.par,,l+1] = pcajames_lda$eigenfunctions
    eigenfun.yao_lda[t.par,,l+1]   = pcayao_lda$eigenfunctions
    eigenfun.buja_lda[t.par,,l+1]  = pcabuja_lda$eigenfunctions
    
    scores.t_lda[,,l+1]     = pcanomiss_lda$scores
    scores.james_lda[,,l+1] = pcajames_lda$scores
    scores.yao_lda[,,l+1]   = pcayao_lda$scores
    scores.buja_lda[,,l+1]  = pcabuja_lda$scores
    
    best.t_lda[l+1]     = pcanomiss_lda$best
    best.james_lda[l+1] = pcajames_lda$best
    best.yao_lda[l+1]   = pcayao_lda$best
    best.buja_lda[l+1]  = pcabuja_lda$best
    
    param.james_lda[l+1] = pcajames_lda$opt.basis
    param.yao.m_lda[l+1] = pcayao_lda$opt.bwm
    param.yao.c_lda[l+1] = pcayao_lda$opt.bwc
    param.buja_lda[l+1,] = pcabuja_lda$opt.lambda
#    
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

    param.james[l+1] = pcajames$opt.basis
    param.yao.m[l+1] = pcayao$opt.bwm
    param.yao.c[l+1] = pcayao$opt.bwc
    param.buja[l+1,] = pcabuja$opt.lambda
    
  }
  
  
  out = list(err.cv.t, err.cv.james, err.cv.yao, err.cv.buja, eigenfun.t, 
             eigenfun.james, eigenfun.yao, eigenfun.buja, scores.t, 
             scores.james, scores.yao, scores.buja, best.t, best.james, 
             best.yao, best.buja, param.james, param.yao.m, param.yao.c, param.buja,
             err.cv.t_lda, err.cv.james_lda, err.cv.yao_lda, err.cv.buja_lda, eigenfun.t_lda, 
             eigenfun.james_lda, eigenfun.yao_lda, eigenfun.buja_lda, scores.t_lda, 
             scores.james_lda, scores.yao_lda, scores.buja_lda, best.t_lda, best.james_lda, 
             best.yao_lda, best.buja_lda, param.james_lda, param.yao.m_lda, param.yao.c_lda, param.buja_lda,
             err_SB_VDFR, optimal_cut_SB_VDFR, prob_SB_VDFR, Beta_estimated)
  names(out) = c("err.cv.t", "err.cv.james", "err.cv.yao", "err.cv.buja",
                 "eigenfun.t", "eigenfun.james", "eigenfun.yao", "eigenfun.buja",
                 "scores.t", "scores.james", "scores.yao", "scores.buja", "best.t",
                 "best.james", "best.yao", "best.buja", "best.basis", "best.bwmean",
                 "best.bwcov", "best.lambda",
                 "err.cv.t_lda", "err.cv.james_lda", "err.cv.yao_lda", "err.cv.buja_lda",
                 "eigenfun.t_lda", "eigenfun.james_lda", "eigenfun.yao_lda", "eigenfun.buja_lda",
                 "scores.t_lda", "scores.james_lda", "scores.yao_lda", "scores.buja_lda", "best.t_lda",
                 "best.james_lda", "best.yao_lda", "best.buja_lda", "best.basis_lda", "best.bwmean_lda",
                 "best.bwcov_lda", "best.lambda_lda", "err_SB_VDFR", "optimal_cut_SB_VDFR", "prob_SB_VDFR", "Beta_estimated")
  return(out)
}

