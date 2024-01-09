# NECESSARY FUNCTIONS IN OUR NEW METHODOLOGY

####### FUNCTION TO TRANSFORM ALL CURVES TO BE FITTED BY PFR GOLDSMITH


############## FUNCTION TO CREATE THE CORRECT B-SPLINE BASIS

bspline <-function(X., XL., XR., NDX., BDEG.){
  dx <- (XR. - XL.)/NDX.
  knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
  B <- spline.des(knots, X., BDEG.+1, 0*X.)$design
  res <- list(B = B, knots = knots)
  res
}

############## THE FUNCTION OF THE PACKAGE SOP AND SOPExamples WITH SOME NECESSARY ADJUSTMENS

sop.fit=function (y, X, Z, weights = NULL, G = NULL, vcstart = NULL, 
                  etastart = NULL, mustart = NULL, offset = NULL, family = gaussian(), 
                  control = sop.control()) 
{
  deviance <- function(C, G, w, sigma2, ssr, edf) {
    log_det_C <- determinant(C)$modulus
    log_det_G <- determinant(G)$modulus
    deviance <- log_det_C + log_det_G + sum(log(sigma2 * 
                                                  1/w)) + ssr/sigma2 + edf
    deviance
  }
  control <- do.call("sop.control", control)
  if (missing(X)) 
    stop("Missing argument: 'X' must be provided")
  if (missing(y)) 
    stop("Missing argument: 'y' must be provided")
  if (missing(Z)) 
    stop("Missing argument: 'Z' must be provided")
  if (!is.null(vcstart)) {
    if (length(vcstart) != (length(G) + 1)) {
      stop("The length of 'vcstart' should be equal to the length of 'G' + 1)")
    }
  }
  trace <- control$trace
  ynames <- if (is.matrix(y)) {
    rownames(y)
  }
  else {
    names(y)
  }
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- NCOL(X)
  EMPTY <- nvars == 0
  nc <- ncol(X)
  ncomp <- length(G)
  nq <- ncol(Z)
  if (!is.null(vcstart)) {
    la <- vcstart
  }
  else {
    la <- rep(1, len = ncomp + 1)
  }
  devold <- 1e+10
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  prior.weights <- weights
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  if (is.null(G)) {
    stop("Missing argument: 'G' must be provided")
  }
  Xnames <- colnames(X)
  na.ind <- is.na(y)
  y.tmp <- y
  y[na.ind] <- 1
  weights <- weights * (!na.ind)
  X <- as.matrix(X)
  start <- NULL
  variance <- family$variance
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) {
    stop("illegal `family' argument")
  }
  valideta <- family$valideta
  if (is.null(valideta)) {
    valideta <- function(eta) TRUE
  }
  validmu <- family$validmu
  if (is.null(validmu)) {
    validmu <- function(mu) TRUE
  }
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (NCOL(y) > 1) {
    stop("y must be univariate")
  }
  eta <- if (!is.null(etastart)) {
    etastart
  }
  else if (!is.null(start)) {
    if (length(start) != nvars) {
      stop(gettextf("Length of start should equal %d and correspond to initial coefs.", 
                    nvars))
    }
    else {
      coefold <- start
      offset + as.vector(if (nvars == 1) 
        c(X, Z) * start
        else c(X, Z) %*% start)
    }
  }
  else {
    family$linkfun(mustart)
  }
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) {
    stop("Can't find valid starting values: please specify some")
  }
  for (i in 1:control$maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offset) + (y - mu)/deriv
    w <- as.vector(deriv^2/family$variance(mu))
    w <- w * weights
    z[!weights] <- 0
    mat <- construct.matrices(X, Z, z, w,GLAM = FALSE) #############################
    if (trace) 
      start1 <- proc.time()[3]
    for (it in 1:control$maxit) {
      Ginv <- 0
      for (ii in 1:ncomp) {
        Ginv <- Ginv + 1/la[ii + 1] * G[[ii]]
      }
      GG <- 1/Ginv
      V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., 
                           mat$ZtZ.)
      D <- diag(c(rep(0, nc), Ginv))
      H <- (1/la[1]) * V + D
      Hinv <- try(solve(H))
      if (inherits(Hinv, "try-error")) {
        Hinv <- ginv(H)
      }
      b <- (1/la[1]) * Hinv %*% mat$u
      b.fixed <- b[1:nc]
      b.random <- b[-(1:nc)]
      aux <- GG - diag(Hinv[-(1:nc), -(1:nc)])
      ed.sum <- 0
      ied <- taus <- NULL
      for (i.fun in 1:ncomp) {
        d1.d <- (1/la[i.fun + 1]) * G[[i.fun]]
        ed1 <- sum(aux * d1.d)
        ed1 <- ifelse(ed1 <= 1e-10, 1e-10, ed1)
        tau1 <- sum(b.random^2 * G[[i.fun]])/ed1
        tau1 <- ifelse(tau1 <= 1e-10, 1e-10, tau1)
        taus <- c(taus, tau1)
        ied <- c(ied, ed1)
      }
      ssr <- mat$yty. - t(c(b.fixed, b.random)) %*% (2 * 
                                                       mat$u - V %*% b)
      dev <- deviance(H, diag(GG), w[w != 0], la[1], ssr, 
                      sum(b.random^2 * Ginv))[1]
      if (family$family == "gaussian" | family$family == 
          "Gamma" | family$family == "quasipoisson") {
        sig2 <- as.numeric((ssr/(length(y[weights != 
                                            0]) - sum(ied) - nc)))
      }
      else {
        sig2 <- 1
      }
      lanew <- c(sig2, taus)
      dla <- abs(devold - dev)
      if (trace) {
        cat(sprintf("%1$3d %2$10.6f", it, dev))
        cat(sprintf("%8.3f", ied), "\n")
      }
      if (dla < control$epsilon) 
        break
      la <- lanew
      if (la[1]<1e-6) {
        la[1]=1e-6
      }
      devold <- dev
    }
    if (trace) {
      end1 <- proc.time()[3]
      cat("Timings:\nSOP", (end1 - start1), "seconds\n")
    }
    eta.old <- eta
    eta <- X %*% b.fixed + Z %*% b.random + offset
    mu <- linkinv(eta)
    tol <- sum((eta - eta.old)^2)/sum(eta^2)
    if (tol < control$epsilon | (family$family == "gaussian" & 
                                 family$link == "identity")) 
      break
  }
  end <- proc.time()[3]
  mu.eta <- family$mu.eta
  mu.eta.val <- mu.eta(eta)
  linear.predictor <- eta
  mu <- linkinv(eta)
  names(mu) <- ynames
  names(linear.predictor) <- ynames
  residuals <- family$dev.resids(y, mu, weights)
  s <- attr(residuals, "sign")
  if (is.null(s)) {
    s <- sign(y - mu)
  }
  residuals <- sqrt(pmax(residuals, 0)) * s
  names(residuals) <- ynames
  residuals[na.ind] <- NA
  names(b.fixed) <- Xnames
  names(b.random) <- colnames(Z)
  names(la) <- c("ssr", names(G))
  names(ied) <- c(names(G))
  dev.residuals <- family$dev.resids(y, mu, weights)
  dev.residuals[na.ind] <- NA
  deviance <- sum(dev.residuals, na.rm = TRUE)
  null.deviance <- glm(y ~ offset(offset), family = family, 
                       weights = prior.weights)$deviance
  out <- list(tol.ol = tol, it.ol = i, tol.il = dla, it.in = it, 
              vc = la, edf = ied)
  fit <- list()
  fit$b.fixed <- b.fixed
  fit$b.random <- b.random
  fit$fitted.values <- mu
  fit$linear.predictor <- linear.predictor
  fit$residuals <- residuals
  fit$X <- X
  fit$Z <- Z
  fit$G <- G
  fit$y <- y.tmp
  fit$weights <- weights
  fit$family <- family
  fit$out <- out
  fit$deviance <- deviance
  fit$null.deviance <- null.deviance
  fit$Vp <- Hinv
  class(fit) <- "sop"
  invisible(fit)
}


############### INTERNAL FUNCTIONS OF THE SOP PACKAGE

H= function (X, A) 
{
  d <- dim(A)
  M <- matrix(A, nrow = d[1])
  XM <- X %*% M
  array(XM, c(nrow(XM), d[-1]))
}
#################

fit.SOP = function (y, X, Z, Lambda, diagonal = FALSE, family = gaussian(), 
                    offset = NULL, weights = NULL, maxit = 200, thr = 0.001, 
                    trace = TRUE) 
{
  if (trace) 
    start <- proc.time()[3]
  if (is.null(offset)) 
    offset <- rep(0, length(y))
  if (is.null(weights)) 
    weights <- rep(1, length(y))
  if (is.null(X)) {
    stop("The design matrix associated with the fixed part of the model can not be NULL")
  }
  if (diagonal) {
    g <- construct.capital.lambda(Lambda)
  }
  else {
    g <- construct.capital.lambda.matrices(Lambda)
  }
  la = c(1, rep(1, length = length(g)))
  devold = 1e+10
  np <- c(ncol(X), ncol(Z))
  mustart <- etastart <- NULL
  nobs <- length(y)
  eval(family$initialize)
  mu <- mustart
  eta <- family$linkfun(mustart)
  if (trace) {
    cat("Effective dimensions\n")
    cat("-------------------------\n")
    cat(sprintf("%1$3s %2$12s", "It.", "Deviance"), sep = "")
    cat(sprintf("%12s", names(g)), sep = "")
    cat("\n")
  }
  for (it in 1:maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offset) + (y - mu)/deriv
    w <- as.vector(deriv^2/family$variance(mu))
    w <- w * weights
    mat <- construct.matrices(X, Z, z, w, FALSE)
    for (it in 1:maxit) {
      if (diagonal) {
        Ginv <- vector(length = length(g[[1]]))
        for (i in 1:length(g)) {
          Ginv <- Ginv + (1/la[i + 1]) * g[[i]]
        }
        G <- 1/Ginv
        V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., 
                             mat$ZtZ.)
        D <- Matrix::bdiag(diag(rep(0, np[1]), ncol = np[1]), 
                           diag(Ginv))
      }
      else {
        Ginv = 0
        for (i in 1:length(g)) {
          Ginv <- Ginv + (1/la[i + 1]) * g[[i]]
        }
        Ginv <- Matrix(Ginv)
        G <- try(solve(Ginv))
        if (class(G) == "try-error") {
          G <- ginv(Ginv)
        }
        V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., 
                             mat$ZtZ.)
        D <- Matrix::bdiag(diag(rep(0, np[1]), nrow = np[1], 
                                ncol = np[1]), Ginv)
      }
      H <- (1/la[1]) * V + D
      Hinv <- try(solve(H))
      if (class(Hinv) == "try-error") {
        Hinv <- ginv(as.matrix(H))
      }
      b <- (1/la[1]) * Hinv %*% mat$u
      b.fixed <- b[1:np[1]]
      b.random <- b[-(1:np[1])]
      if (diagonal) {
        aux <- G - diag(Hinv[-(1:np[1]), -(1:np[1])])
        ssv <- ed <- tau <- vector(mode = "list", length = length(g))
        for (i in 1:length(g)) {
          g.inv.d <- (1/la[i + 1]) * g[[i]]
          ed[[i]] <- sum(aux * g.inv.d)
          ed[[i]] <- ifelse(ed[[i]] <= 1e-10, 1e-10, 
                            ed[[i]])
          ssv[[i]] <- sum(b.random^2 * g[[i]])
          tau[[i]] <- ssv[[i]]/ed[[i]]
          tau[[i]] <- ifelse(tau[[i]] <= 1e-10, 1e-10, 
                             tau[[i]])
        }
        ssr = mat$yty. - t(c(b.fixed, b.random)) %*% 
          (2 * mat$u - V %*% b)
        dev <- deviance(H, diag(G), w[w != 0], la[1], 
                        ssr, sum(b.random^2 * Ginv))[1]
      }
      else {
        aux <- G - Hinv[-(1:np[1]), -(1:np[1])]
        updates <- lapply(1:length(g), function(i, g, 
                                                la, b.random, aux) {
          g.inv.d <- (1/la[i + 1]) * g[[i]]
          ed <- sum(colSums(t(aux) * g.inv.d))
          tau <- as.vector((t(b.random) %*% g[[i]] %*% 
                              b.random)/ed)
          tau <- ifelse(tau <= 1e-10, 1e-10, tau)
          res <- list(ed = ed, tau = tau)
          res
        }, g = g, la = la, b.random = b.random, aux = aux)
        ed <- as.list(unlist(updates)[seq(1, 2 * length(g), 
                                          by = 2)])
        tau <- as.list(unlist(updates)[seq(2, 2 * length(g), 
                                           by = 2)])
        ssr = mat$yty. - t(c(b.fixed, b.random)) %*% 
          (2 * mat$u - V %*% b)
        dev <- deviance(H, G, w[w != 0], la[1], ssr, 
                        (b.random) %*% Ginv %*% b.random)[1]
      }
      if (family$family == "gaussian" | family$family == 
          "quasipoisson") {
        phi <- as.numeric((ssr/(length(y[w != 0]) - sum(unlist(ed)) - 
                                  np[1])))
      }
      else {
        phi <- 1
      }
      lanew = c(phi, unlist(tau))
      dla = abs(devold - dev)
      if (trace) {
        cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
        cat(sprintf("%12.3f", unlist(ed)), sep = "")
        cat("\n")
      }
      if (dla < thr) 
        break
      la <- lanew
      if (la[1]<1e-6) {
        la[1]=1e-6
      }
      devold <- dev
    }
    eta.old <- eta
    eta <- X %*% b.fixed + Z %*% b.random + offset
    mu <- family$linkinv(eta)
    tol <- sum((eta - eta.old)^2)/sum(eta^2)
    if (tol < 1e-06 | family$family == "gaussian") 
      break
  }
  if (trace) {
    end <- proc.time()[3]
    cat("All process", end - start, "seconds\n")
  }
  res <- list()
  res$dat <- list(y = y, X = X, Z = Z)
  res$family = family
  res$ed <- unlist(ed)
  res$tot_ed <- sum(unlist(ed))
  res$vc <- unlist(tau)
  res$phi <- phi
  res$coeff <- c(b.fixed, b.random)
  res$eta <- eta
  res$mu <- mu
  res$deviance <- dev
  res$convergence <- ifelse(it < maxit, TRUE, FALSE)
  res$Vp <- Hinv
  res
}

################################## INNER_PRODUCT NECESSARY FUNCTIONS:

fdchk <- function(fdobj) {
  
  #  check the class of FDOBJ and extract coefficient matrix
  
  if (inherits(fdobj, "fd")) {
    coef  <- fdobj$coefs
  } else {
    if (inherits(fdobj, "basisfd")) {
      coef  <- diag(rep(1,fdobj$nbasis - length(fdobj$dropind)))
      fdobj <- fd(coef, fdobj)
    } else { 
      stop("FDOBJ is not an FD object.")
    }
  }
  
  #  extract the number of replications and basis object
  
  coefd <- dim(as.matrix(coef))
  if (length(coefd) > 2) stop("Functional data object must be univariate")
  nrep     <- coefd[2]
  basisobj <- fdobj$basis
  
  return(list(nrep, fdobj))
  
}

knotmultchk <- function(basisobj, knotmult) {
  type <- basisobj$type
  if (type == "bspline") {
    # Look for knot multiplicities in first basis
    params  <- basisobj$params
    nparams <- length(params)
    norder  <- basisobj$nbasis - nparams
    if (norder == 1) {
      knotmult <- c(knotmult, params)
    } else {
      if (nparams > 1) {
        for (i in 2:nparams) 
          if (params[i] == params[i-1]) knotmult <- c(knotmult, params[i])
      }
    }
  }
  return(knotmult)
}

inprod=function (fdobj1, fdobj2 = NULL, Lfdobj1 = int2Lfd(0), Lfdobj2 = int2Lfd(0), rng = range1, wtfd = 0) {
  
  result1 <- fdchk(fdobj1)
  nrep1 <- result1[[1]]
  fdobj1 <- result1[[2]]
  coef1 <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1 <- basisobj1$type
  range1 <- basisobj1$rangeval
  if (is.null(fdobj2)) {
    tempfd <- fdobj1
    tempbasis <- tempfd$basis
    temptype <- tempbasis$type
    temprng <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- create.bspline.basis(temprng, 1, 1)
    }
    else {
      if (temptype == "fourier") 
        basis2 <- create.fourier.basis(temprng, 1)
      else basis2 <- create.constant.basis(temprng)
    }
    fdobj2 <- fd(1, basis2)
  }
  result2 <- fdchk(fdobj2)
  nrep2 <- result2[[1]]
  fdobj2 <- result2[[2]]
  coef2 <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2 <- basisobj2$type
  range2 <- basisobj2$rangeval
  if (rng[1] < range1[1] || rng[2] > range1[2]) 
    stop("Limits of integration are inadmissible.")
  # if (is.fd(fdobj1) && is.fd(fdobj2) && type1 == "bspline" && 
  #     type2 == "bspline" && is.eqbasis(basisobj1, basisobj2) && 
  #     is.integer(Lfdobj1) && is.integer(Lfdobj2) && length(basisobj1$dropind) == 
  #     0 && length(basisobj1$dropind) == 0 && wtfd == 0 && all(rng == 
  #                                                             range1)) {
  #   inprodmat <- inprod.bspline(fdobj1, fdobj2, Lfdobj1$nderiv, 
  #                               Lfdobj2$nderiv)
  #   return(inprodmat)
  # }
  Lfdobj1 <- int2Lfd(Lfdobj1)
  Lfdobj2 <- int2Lfd(Lfdobj2)
  iter <- 0
  rngvec <- rng
  knotmult <- numeric(0)
  if (type1 == "bspline") 
    knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline") 
    knotmult <- knotmultchk(basisobj2, knotmult)
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < 
                           rng[2]]
    rngvec <- c(rng[1], knotmult, rng[2])
  }
  if ((all(c(coef1) == 0) || all(c(coef2) == 0))) 
    return(matrix(0, nrep1, nrep2))
  JMAX <- 15
  JMIN <- 5
  EPS <- 1e-04
  inprodmat <- matrix(0, nrep1, nrep2)
  nrng <- length(rngvec)
  for (irng in 2:nrng) {
    rngi <- c(rngvec[irng - 1], rngvec[irng])
    if (irng > 2) 
      rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) 
      rngi[2] <- rngi[2] - 1e-10
    iter <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1, JMAXP)
    h[2] <- 0.25
    s <- array(0, c(JMAXP, nrep1, nrep2))
    sdim <- length(dim(s))
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2)
    if (!is.numeric(wtfd)) {
      wtd <- eval.fd(rngi, wtfd, 0)
      fx2 <- matrix(wtd, dim(wtd)[1], dim(fx2)[2]) * fx2
    }
    s[1, , ] <- width * matrix(crossprod(fx1, fx2), nrep1, 
                               nrep2)/2
    tnm <- 0.5
    for (iter in 2:JMAX) {
      tnm <- tnm * 2
      if (iter == 2) {
        x <- mean(rngi)
      }
      else {
        del <- width/tnm
        x <- seq(rngi[1] + del/2, rngi[2] - del/2, del)
      }
      fx1 <- eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- eval.fd(x, fdobj2, Lfdobj2)
      if (!is.numeric(wtfd)) {
        wtd <- eval.fd(wtfd, x, 0)
        fx2 <- matrix(wtd, dim(wtd)[1], dim(fx2)[2]) * 
          fx2
      }
      chs <- width * matrix(crossprod(fx1, fx2), nrep1, 
                            nrep2)/tnm
      s[iter, , ] <- (s[iter - 1, , ] + chs)/2
      if (iter >= 5) {
        ind <- (iter - 4):iter
        ya <- s[ind, , ]
        ya <- array(ya, c(5, nrep1, nrep2))
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y <- ya[ns, , ]
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5 - m)) {
            ho <- xa[i]
            hp <- xa[i + m]
            w <- (cs[i + 1, , ] - ds[i, , ])/(ho - hp)
            ds[i, , ] <- hp * w
            cs[i, , ] <- ho * w
          }
          if (2 * ns < 5 - m) {
            dy <- cs[ns + 1, , ]
          }
          else {
            dy <- ds[ns, , ]
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        }
        else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN) 
          break
      }
      s[iter + 1, , ] <- s[iter, , ]
      h[iter + 1] <- 0.25 * h[iter]
      if (iter == JMAX) 
        warning("Failure to converge.")
    }
    inprodmat <- inprodmat + ss
  }
  if (length(dim(inprodmat) == 2)) {
    return(as.matrix(inprodmat))
  }
  else {
    return(inprodmat)
  }
}

# SIMPSON NUMERIC INTEGRATION VERSION 1.0 SEE BURDEN AND FAIRES (2016)

Simpson <- function(fdobj1, fdobj2=NULL, fdobj3=NULL, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),rng = rng, sub=25, wtfd = 0){
  
  #  Check FDOBJ1 and get no. replications and basis object
  
  result1   <- fdchk(fdobj1)
  nrep1     <- result1[[1]]
  fdobj1    <- result1[[2]]
  coef1     <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1     <- basisobj1$type
  range1    <- basisobj1$rangeval
  
  #  Default FDOBJ2 to a constant function, using a basis that matches
  #  that of FDOBJ1 if possible.
  
  if (is.null(fdobj2)) {
    tempfd    <- fdobj1
    tempbasis <- tempfd$basis
    temptype  <- tempbasis$type
    temprng   <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- create.bspline.basis(temprng, 1, 1)
    } else {
      if (temptype == "fourier") basis2 <- create.fourier.basis(temprng, 1)
      else                       basis2 <- create.constant.basis(temprng)
    }
    fdobj2 <- fd(1,basis2)
  }
  
  #  Check FDOBJ2 and get no. replications and basis object
  
  result2   <- fdchk(fdobj2)
  nrep2     <- result2[[1]]
  fdobj2    <- result2[[2]]
  coef2     <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2     <- basisobj2$type
  range2    <- basisobj2$rangeval
  
  
  if (is.null(fdobj3)) {
    return(inprod(fdobj1, fdobj2, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                  rng = rng, wtfd = 0))
  }
  
  # check ranges
  
  if (rng[1] < range1[1] || rng[2] > range1[2]) stop(
    "Limits of integration are inadmissible.")

  #  check LFDOBJ1 and LFDOBJ2
  
  Lfdobj1 <- int2Lfd(Lfdobj1)
  Lfdobj2 <- int2Lfd(Lfdobj2)
  # Lfdobj3 <- int2Lfd(Lfdobj3)
  
  #  set iter
  
  iter <- 0
  
  # The default case, no multiplicities.
  
  rngvec <- rng
  
  #  check for any knot multiplicities in either argument
  
  knotmult <- numeric(0)

  if (type1 == "bspline") knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline") knotmult <- knotmultchk(basisobj2, knotmult)

  #  Modify RNGVEC defining subinvervals if there are any
  #  knot multiplicities.
  
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
  }
  
  #  check for either coefficient array being zero
  
  if ((all(c(coef1) == 0) || all(c(coef2) == 0)))
    return(matrix(0,nrep1,nrep2*length(fdobj3)))
  
   n=2*sub # THIS IS THE INTERVALS FOR THE SIMPSON METHOD

  nrng <- length(rngvec)
  for (irng  in  2:nrng) {  
    rngi <- c(rngvec[irng-1],rngvec[irng])
    
    #  change range so as to avoid being exactly on multiple knot values
    
    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10

    width = (rngi[2] - rngi[1])/n
    
    # the first iteration uses just the endpoints
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2)
    fx3 <- matrix(fdobj3,nrow = dim(fx2)[1],ncol = length(fdobj3),byrow = TRUE)
    fx2 <- Rten2(fx2,fx3)
    
    XI0 = matrix(crossprod(fx1,fx2),nrep1,nrep2*length(fdobj3))
    
    XI1=0
    XI2=0
    
    for (i in 1:(n-1)) {

### THESE LINES OF CODE ARE THE CORE OF THE INNER PRODUCT. NOTICE HOW WE PERFORM THE INTEGRATION ONLY IN THE t VARIABLE BUT THEN WE CREATE THE BIDEMENSIONAL BASIS IN EVERY ITERATION
      
      x = rngi[1] + i*width
      fx1 <- eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- eval.fd(x, fdobj2, Lfdobj2)
      fx3 <- matrix(fdobj3,nrow = dim(fx2)[1],ncol = length(fdobj3),byrow = TRUE)
      fx2 <- Rten2(fx2,fx3)
      
      Fx=matrix(crossprod(fx1,fx2),nrep1,nrep2*length(fdobj3))
      
      if (i%%2==0) {
        
        XI2= XI2 + Fx
      }else{
        XI1= XI1 + Fx
      }
      
    }
    
    XI= width*(XI0 + 2*XI2 + 4*XI1)/3
  }
  
  if(length(dim(XI) == 2)) {
    #  coerce inprodmat to be nonsparse
    return(as.matrix(XI))
  } else {
    #  allow inprodmat to be sparse if it already is
    return(XI)
  }
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL:

Data2B_simpson=function(X, M, nbasis=c(30,30,30), bdeg=c(3,3,3),sub=25, lim=NULL){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  c3=nbasis[3]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }

    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(1, M[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))
    
     # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,1:M[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,1:M[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)
    
    c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_Phi <- paste("Phi", i, sep = "_")
    L_Phi[[i]]=assign(nam_Phi, bspline(1:M[i], XL, XR, c_t, bdeg[2]))
    
  }
  
  ####### HERE WE CREATE THE MARGINAL BASIS FOR THE T VARIABLE In B(t,T)
  
  xlim_T <- c(min(M), max(M))
  XL_T=xlim_T[1]-0.001
  XR_T=xlim_T[2]+0.001
  
  c_T=c3-bdeg[3] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  if (is.null(lim)) {
  
    B_T=bspline(M, XL_T, XR_T, c_T, bdeg[3])
      
  }else{
    
    B_T=bspline(M, lim[1]-0.001, lim[2]+0.001, c_T, bdeg[3])
  }
  
  
  # matplot(B_T$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT  
                  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS
  for (i in 1:N) {
    
    # DATA BASIS
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
    
    # t MARGINAL BASIS
    
    breaks=L_Phi[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_Phi_aux <- paste("Phi_aux", i, sep = "_")
    L_Phi_aux[[i]]=assign(nam_Phi_aux, create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2))))
    
  }
  
  # PERFORMING THE INNER PRODUCT
  
  
  for (i in 1:N) {
    
    PROD=Simpson(L_X_aux[[i]],L_Phi_aux[[i]],B_T$B[i,],rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
    }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, B_T=B_T, B_Phi=L_Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT FOR ALL SUBJECTS:

Data2B_simpson_no_vc_old=function(X, M, nbasis=c(30,31), bdeg=c(3,3),sub=25){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }
    
    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(1, M[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,1:M[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,1:M[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)
    
  }
  
  rng_t= c(1, max(M))
  
  XL_t=rng_t[1]-0.001
  XR_t=rng_t[2]+0.001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(1:max(M), XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT  
  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS
  for (i in 1:N) {
    
    # DATA BASIS
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
    
  }
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n=length(breaks)
  
  
  Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
  
  for (i in 1:N) {
    

    PROD=Simpson(L_X_aux[[i]],Phi_aux,rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT FOR ALL SUBJECTS:

Data2B_simpson_no_vc=function(X, M, nbasis=c(30,31), bdeg=c(3,3),sub=25){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }
    
    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(1, M[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,1:M[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,1:M[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)
    
  }
  
  rng_t= c(1, max(M))
  
  XL_t=rng_t[1]-0.001
  XR_t=rng_t[2]+0.001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(1:max(M), XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT

  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT  
  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS
  for (i in 1:N) {
    
    # DATA BASIS
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
  
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
  
  }
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n_knots=length(breaks)
  
  
  # Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n_knots-2):(n_knots+2)))
  
  for (i in 1:N) {
  

    # if (i!=N) {

      ind_knots=max(which(Phi$knots<=(M[i]+0.001)))

    # }
    #else{
    #   ind_knots=length(Phi$knots)
    # }

   # IDEA 1: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE)

    # breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_1=c(breaks[1:5],rep(0,n_knots-n),breaks[6:(n-5)],breaks[(n-4):n])
    # 
    # n_breaks_no_vc_1=length(breaks_no_vc_1)
    # 
    # Phi_aux_1=create.bspline.basis(breaks=breaks_no_vc_1,norder=bdeg[2]+1,dropin=c(1:5,(n_breaks_no_vc_1-2):(n_breaks_no_vc_1+2)))
      
   # IDEA 2: TOMAR LOS NODOS DEL INTERVALO

    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_2=c(breaks[1:2],rep(0,n_knots-n),breaks[3:(n-2)],breaks[(n-1):n])
    # 
    # n_breaks_no_vc_2=length(breaks_no_vc_2)
    # 
    # Phi_aux_2=create.bspline.basis(breaks=breaks_no_vc_2,norder=bdeg[2]+1,dropin=c(1:7,(n_breaks_no_vc_2):(n_breaks_no_vc_2+2)))
      
    # IDEA 3: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE) PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    tag=3
    Phi_aux_3=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
      
      
    # IDEA 4: TOMAR LOS NODOS DEL INTERVALO PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
      
      # breaks=Phi$knots[4:(ind_knots)]
      # dife=diff(breaks)[1]
      # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
      # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
      # n=length(breaks)
      # 
      # tag=4
      # 
      # Phi_aux_4=create.bspline.basis(breaks=breaks, norder=bdeg[2]+1,dropin=c(1:2,(n):(n+2)))
      
    
    PROD=Simpson(L_X_aux[[i]],Phi_aux_3,rng = c(1,M[i]),sub = sub)/M[i]
    
    if (tag==3 || tag==4) {
      PROD=cbind(PROD,matrix(0, nrow = dim(PROD)[1], ncol=nbasis[2]-dim(PROD)[2]))
    }
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}

# matplot(eval.basis(1:max(M),Phi_aux_2),type="l")

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT AND PARTIALLY OBSERVED CURVES:

Data2B_simpson_Classification=function(X, M, grid, nbasis=c(80,80), bdeg=c(3,3),sub=1000){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M_b)
  # M is the matrix with columns M_a and M_b of numbers of observations dim(M) == N x 2 
  
  if (is.null(dim(X))) {
    # print("ENTRÓ")
    X=t(as.matrix(X))
    M=t(as.matrix(M))
  }    
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  M_a=M[,1]
  M_b=M[,2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M_a)!=N | length(M_b)!=N) {
    stop("length of 'M_a' and 'M_b' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M_b[i]-M_a[i]+1) {
      stop("Incorrect numbers of NAs in column ",i, " of 'X'")
    }
    
    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(grid[M_a[i]], grid[M_b[i]]) #c(range(grid)[1],range(grid)[2])
    
    XL=rng[i,1]-0.0001
    XR=rng[i,2]+0.0001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(grid[M_a[i]:M_b[i]], XL, XR, c, bdeg[1]))
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,M_a[i]:M_b[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,M_a[i]:M_b[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
  }
  
  ############### HERE WE CREATE THE BASIS FOR THE t VARIABLE In B(t)
  
  # rng_t= c(grid[min(M_a)], grid[max(M_b)]) # THE KNOTS FOR THE FUNCTIONAL COEFFICIENTS HAVE TO BE MAXIMUM 

  rng_t= range(grid) # THE KNOTS FOR THE FUNCTIONAL COEFFICIENTS HAVE TO BE MAXIMUM 
  
  XL_t=rng_t[1]-0.0001
  XR_t=rng_t[2]+0.0001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  # Phi=bspline(grid[min(M_a):max(M_b)], XL_t, XR_t, c_t, bdeg[2])
  
  Phi=bspline(grid, XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  

  # PERFORMING THE INNER PRODUCT

  res=matrix(nrow=N,ncol=c2)
  
  for (i in 1:N) {
    
    # Phi_short=Phi$B[(M_a[i]-min(M_a)+1):(M_b[i]-min(M_a)+1),]

    Phi_short=Phi$B[M_a[i]:M_b[i],]
    
    aux_del=which(colSums(Phi_short)==0)
    
    if (length(aux_del)!=0) {
      
      Phi_short=as.matrix(Phi_short[,-aux_del])
      
    }
    
    # range(grid[i,],na.rm=1)
    # 
    # range(Phi$knots)
    
    min_knot=max(which(Phi$knots<=range(grid[M_a[i]:M_b[i]],na.rm=1)[1]))-3
    
    max_knot=min(which(Phi$knots>=range(grid[M_a[i]:M_b[i]],na.rm=1)[2]))+3
    
    if (min_knot<1) {
      min_knot=1
    }

    if (max_knot>length(Phi$knots)) {
      max_knot=length(Phi$knots)
    }
    
    new_knots=Phi$knots[min_knot:max_knot]
    
    # B_short <- spline.des(new_knots, grid[M_a[i]:M_b[i]], bdeg[2]+1, 0*grid[M_a[i]:M_b[i]])$design
    # 
    # all.equal(B_short,Phi_short)
    # 
    # breaks=new_knots
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # basis_Phi = create.bspline.basis(breaks = breaks, norder=bdeg[1]+1, nbasis = dim(B_short)[2]+10, dropin=c(1:5,(n-2):(n+2)))
    # 
    # Phi_short_aux=eval.basis(grid[M_a[i]:M_b[i]], basis_Phi)#[aux_sup:aux_inf,-aux_del]
    # 
    # max(abs(Phi_short_aux-Phi_short))==0 # THESE CHECKS ARE TO DETERMINE IF THE TRANSFORMATION INTO THE FDA B-splines ARE CORRECT

    # Phi_X=eval.basis(grid[i,1:M[i]], L_X_aux[[i]])

    B_X_a=spline.des(L_X[[i]]$knots, rng[i,1], bdeg[1]+1, 0*rng[i,1])$design
    B_Phi_a=spline.des(new_knots, rng[i,1], bdeg[2]+1, 0*rng[i,1])$design
    
    B_X_b=spline.des(L_X[[i]]$knots, rng[i,2], bdeg[1]+1, 0*rng[i,2])$design
    B_Phi_b=spline.des(new_knots, rng[i,2], bdeg[2]+1, 0*rng[i,2])$design
    
    # max(abs(B_X-L_X[[i]]$B[1,]))==0
    
    n=2*sub # THIS IS THE INTERVALS FOR THE SIMPSON METHOD
    
    width = (rng[i,2] - rng[i,1])/n
    
    XIa = t(B_X_a)%*% B_Phi_a
    XIb = t(B_X_b)%*% B_Phi_b
    XI1=0
    XI2=0
    
    
    for (j in 1:(n-1)) {
      
      ### THESE LINES OF CODE ARE THE CORE OF THE INNER PRODUCT.
      
      x = rng[i,1] + j*width
      
      fx1=spline.des(L_X[[i]]$knots, x, bdeg[1]+1, 0*x)$design
      
      fx2=spline.des(new_knots, x, bdeg[2]+1, 0*x)$design
      
      Fx=t(fx1) %*% fx2
      
      
      
      if (j%%2==0) {
        
        XI2= XI2 + Fx
      }else{
        XI1= XI1 + Fx
      }
      
      
      
    }
    
    XI = matrix(0,nrow = c1, ncol = c2)
    
    
    if (length(aux_del)!=0) {
      
      XI[,-aux_del]= width*(XIa + XIb + 2*XI2 + 4*XI1)/3
      
    }else{
      
      XI= width*(XIa + XIb + 2*XI2 + 4*XI1)/3
      
    }
    
    
    
    
    res[i,]=t(L_theta[[i]]) %*% XI # ESTO ES LO MISMO QUE RELLENAR CON CEROS UNA MATRIZ DE 
                                   # PRODUCTOS INTERIORES DE TAMAÑO COMPLETO 
                                   # Y LUEGO MULTIPLICARLA POR LA MATRIZ A
  }
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT AND PARTIALLY OBSERVED CURVES:

Data2B_simpson_Classification_knots=function(X, M, grid, nbasis=c(20,20), bdeg=c(3,3),sub=25){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M_b)
  # M is the matrix with columns M_a and M_b of numbers of observations dim(M) == N x 2 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  M_a=M[,1]
  M_b=M[,2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M_a)!=N | length(M_b)!=N) {
    stop("length of 'M_a' and 'M_b' must be equal to 'N'")
  }
  
  for (i in 1:N) {

    # print(i)
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M_b[i]-M_a[i]+1) {
      stop("Incorrect numbers of NAs in column ",i, " of 'X'")
    }

    ############### HERE WE CREATE THE BASIS FOR THE DATA, ALL CURVES WILL HAVE THE SAME KNOTS

    rng_t= c(grid[min(M_a)], grid[max(M_b)]) # THE KNOTS FOR THE FUNCTIONAL DATA ARE MAXIMUM 
    
    XL_t=rng_t[1]-0.0001
    XR_t=rng_t[2]+0.0001
    
    c_t=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    Phi_X=bspline(grid[min(M_a):max(M_b)], XL_t, XR_t, c_t, bdeg[1])
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT

    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)

    B=Phi_X$B[M_a[i]:M_b[i],]
    
    ind_col=which(colSums(B)!=0)
    
    B=B[,ind_col]
    
    # basis = create.bspline.basis(rangeval=grid[c(1,p)], breaks=Phi_X$knots[ind_knots_X])
    # 
    # B=eval.basis(grid, basis)
    
    
    aux= B #Phi_X$B[aux_knots_x,aux_knots]#(aux_knots[1]-bdeg[1]):(aux_knots[length(aux_knots)]+bdeg[1])]

    aux_2=B2XZG_1d(aux,2,dim(B)[2])

    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,M_a[i]:M_b[i]])

    A[i,((c1*(i-1))+1):(i*c1)]=c(aux_3$theta,rep(0,c1-length(aux_3$theta)))

    L_theta[[i]]=aux_3$theta

    nam_y <- paste("y_h", i, sep = "_")

    L_y[[i]]=assign(nam_y, aux%*%L_theta[[i]])

    error[i]=mean(abs((X[i,M_a[i]:M_b[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA

  }

  
  ############### HERE WE CREATE THE BASIS FOR THE t VARIABLE In B(t)
  
  rng_t= c(grid[min(M_a)], grid[max(M_b)]) # THE KNOTS FOR THE FUNCTIONAL COEFFICIENTS HAVE TO BE MAXIMUM 
  
  XL_t=rng_t[1]-0.0001
  XR_t=rng_t[2]+0.0001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(grid[min(M_a):max(M_b)], XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n_knots=length(breaks)
  
  
  # Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n_knots-2):(n_knots+2)))
  
  for (i in 1:N) {
    
    # print(i)
    
    # if (i!=N) {
    
    ind_knots=max(which(Phi$knots<=(grid[M_b[i]]+0.0001) & Phi$knots>=(grid[M_a[i]]-0.0001)))
    
    ind_knots_X=max(which(Phi_X$knots<=(grid[M_b[i]]+0.0001) & Phi_X$knots>=(grid[M_a[i]]-0.0001)))
    
    ind_knots_X=which(Phi_X$knots<=(grid[M_b[i]]+0.0001) & Phi_X$knots>=(grid[M_a[i]]-0.0001))
    
    aux_knots=ind_knots_X[which(ind_knots_X<c1 & ind_knots_X>bdeg[1])]
    
    aux_knots_x=which(grid>=Phi_X$knots[aux_knots[1]] & grid<=Phi_X$knots[aux_knots[length(aux_knots)]]) 
    
    # }
    #else{
    #   ind_knots=length(Phi$knots)
    # }
    
    # IDEA 1: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE)
    
    # breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_1=c(breaks[1:5],rep(0,n_knots-n),breaks[6:(n-5)],breaks[(n-4):n])
    # 
    # n_breaks_no_vc_1=length(breaks_no_vc_1)
    # 
    # Phi_aux_1=create.bspline.basis(breaks=breaks_no_vc_1,norder=bdeg[2]+1,dropin=c(1:5,(n_breaks_no_vc_1-2):(n_breaks_no_vc_1+2)))
    
    # IDEA 2: TOMAR LOS NODOS DEL INTERVALO
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_2=c(breaks[1:2],rep(0,n_knots-n),breaks[3:(n-2)],breaks[(n-1):n])
    # 
    # n_breaks_no_vc_2=length(breaks_no_vc_2)
    # 
    # Phi_aux_2=create.bspline.basis(breaks=breaks_no_vc_2,norder=bdeg[2]+1,dropin=c(1:7,(n_breaks_no_vc_2):(n_breaks_no_vc_2+2)))
    
    # IDEA 3: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE) PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    breaks=Phi$knots[(ind_knots_X[1]-bdeg[2]):(ind_knots_X[length(ind_knots_X)]+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    tag=3
    Phi_aux_3=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    breaks=Phi_X$knots[(ind_knots_X[1]-bdeg[2]):(ind_knots_X[length(ind_knots_X)]+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    Phi_aux_X=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    # IDEA 4: TOMAR LOS NODOS DEL INTERVALO PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # tag=4
    # 
    # Phi_aux_4=create.bspline.basis(breaks=breaks, norder=bdeg[2]+1,dropin=c(1:2,(n):(n+2)))
    
    
    PROD=Simpson(Phi_aux_X,Phi_aux_3,rng = c(grid[aux_knots_x[1]],grid[aux_knots_x[length(aux_knots_x)]]),sub = sub)/length(aux_knots_x)
    
    if (tag==3 || tag==4) {
      if (dim(PROD)[2]!=nbasis[2]) {
        
        PROD=cbind(PROD,matrix(0, nrow = dim(PROD)[1], ncol=nbasis[2]-dim(PROD)[2]))
      }
      
      if (dim(PROD)[1]!=nbasis[1]) {
        
        PROD=rbind(PROD,matrix(0, nrow = nbasis[1]-dim(PROD)[1], ncol=dim(PROD)[2]))
      }
    }
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT AND PARTIALLY OBSERVED CURVES:

Data2B_simpson_Classification_knots_2=function(X, M, grid, knots, nbasis=c(20,20), bdeg=c(3,3),sub=25){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M_b)
  # M is the matrix with columns M_a and M_b of numbers of observations dim(M) == N x 2 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  M_a=M[,1]
  M_b=M[,2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M_a)!=N | length(M_b)!=N) {
    stop("length of 'M_a' and 'M_b' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    # print(i)
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M_b[i]-M_a[i]+1) {
      stop("Incorrect numbers of NAs in column ",i, " of 'X'")
    }
    
    ############### HERE WE CREATE THE BASIS FOR THE DATA, ALL CURVES WILL HAVE THE SAME KNOTS
    
    
    rng_t= c(grid[min(M_a)], grid[max(M_b)]) # THE KNOTS FOR THE FUNCTIONAL DATA ARE MAXIMUM 
    
    XL_t=rng_t[1]-0.0001
    XR_t=rng_t[2]+0.0001
    
    c_t=c1-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    # breaks=knots
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    
    # basis = create.bspline.basis(rangeval=grid[c(1,p)], breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    basis = create.bspline.basis(rangeval=grid[c(1,p)], breaks=knots)
    
    B=eval.basis(grid, basis)
    
    

    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    ind_knots_X=which(knots<=(grid[M_b[i]]+0.0001) & knots>=(grid[M_a[i]]-0.0001))
    
    aux_knots=ind_knots_X[which(ind_knots_X<c1 & ind_knots_X>bdeg[1])]
    
    aux_knots_x=which(grid>=knots[aux_knots[1]] & grid<=knots[aux_knots[length(aux_knots)]]) 
    
    aux=B[aux_knots_x,aux_knots]#(aux_knots[1]-bdeg[1]):(aux_knots[length(aux_knots)]+bdeg[1])]
    
    L_X[[i]]=aux
    
    aux_2=B2XZG_1d(aux,2,length(aux_knots))
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,aux_knots_x])
    
    A[i,((c1*(i-1))+1):(i*c1)]=c(aux_3$theta,rep(0,c1-length(aux_3$theta)))
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, aux%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,aux_knots_x])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
  }
  
  
  ############### HERE WE CREATE THE BASIS FOR THE t VARIABLE In B(t)
  
  rng_t= c(grid[min(M_a)], grid[max(M_b)]) # THE KNOTS FOR THE FUNCTIONAL COEFFICIENTS HAVE TO BE MAXIMUM 
  
  XL_t=rng_t[1]-0.0001
  XR_t=rng_t[2]+0.0001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(grid[min(M_a):max(M_b)], XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n_knots=length(breaks)
  
  
  # Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n_knots-2):(n_knots+2)))
  
  for (i in 1:N) {
    
    # print(i)
    
    # if (i!=N) {
    
    ind_knots=max(which(Phi$knots<=(grid[M_b[i]]+0.0001) & Phi$knots>=(grid[M_a[i]]-0.0001)))
    
    ind_knots_Phi=which(Phi$knots<=(grid[M_b[i]]+0.0001) & Phi$knots>=(grid[M_a[i]]-0.0001))
    
    ind_knots_X=which(knots<=(grid[M_b[i]]+0.0001) & knots>=(grid[M_a[i]]-0.0001))
    
    aux_knots=ind_knots_X[which(ind_knots_X<c1 & ind_knots_X>bdeg[1])]
    
    aux_knots_x=which(grid>=knots[aux_knots[1]] & grid<=knots[aux_knots[length(aux_knots)]]) 
    
    # }
    #else{
    #   ind_knots=length(Phi$knots)
    # }
    
    # IDEA 1: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE)
    
    # breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_1=c(breaks[1:5],rep(0,n_knots-n),breaks[6:(n-5)],breaks[(n-4):n])
    # 
    # n_breaks_no_vc_1=length(breaks_no_vc_1)
    # 
    # Phi_aux_1=create.bspline.basis(breaks=breaks_no_vc_1,norder=bdeg[2]+1,dropin=c(1:5,(n_breaks_no_vc_1-2):(n_breaks_no_vc_1+2)))
    
    # IDEA 2: TOMAR LOS NODOS DEL INTERVALO
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_2=c(breaks[1:2],rep(0,n_knots-n),breaks[3:(n-2)],breaks[(n-1):n])
    # 
    # n_breaks_no_vc_2=length(breaks_no_vc_2)
    # 
    # Phi_aux_2=create.bspline.basis(breaks=breaks_no_vc_2,norder=bdeg[2]+1,dropin=c(1:7,(n_breaks_no_vc_2):(n_breaks_no_vc_2+2)))
    
    # IDEA 3: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE) PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    breaks=Phi$knots[(ind_knots_Phi[1]-bdeg[2]):(ind_knots_Phi[length(ind_knots_Phi)]+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    tag=3
    Phi_aux_3=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    breaks=knots[(ind_knots_X[1]):(ind_knots_X[length(ind_knots_X)])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    Phi_aux_X=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    # IDEA 4: TOMAR LOS NODOS DEL INTERVALO PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # tag=4
    # 
    # Phi_aux_4=create.bspline.basis(breaks=breaks, norder=bdeg[2]+1,dropin=c(1:2,(n):(n+2)))
    
    
    PROD=Simpson(Phi_aux_X,Phi_aux_3,rng = c(grid[aux_knots_x[1]],grid[aux_knots_x[length(aux_knots_x)]]),sub = sub)/length(aux_knots_x)
    
    if (tag==3 || tag==4) {
      if (dim(PROD)[2]!=nbasis[2]) {
        
        PROD=cbind(PROD,matrix(0, nrow = dim(PROD)[1], ncol=nbasis[2]-dim(PROD)[2]))
      }
      
      if (dim(PROD)[1]!=nbasis[1]) {
        
        PROD=rbind(PROD,matrix(0, nrow = nbasis[1]-dim(PROD)[1], ncol=dim(PROD)[2]))
      }
    }
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}


####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT AND PARTIALLY OBSERVED CURVES:

Data2B_simpson_Classification_2=function(X, M, nbasis=c(30,31), bdeg=c(3,3),sub=25){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M_b)
  # M is the matrix with columns M_a and M_b of numbers of observations dim(M) == N x 2 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  M_a=M[,1]
  M_b=M[,2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M_a)!=N | length(M_b)!=N) {
    stop("length of 'M_a' and 'M_b' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M_b[i]-M_a[i]+1) {
      stop("Incorrect numbers of NAs in column ",i, " of 'X'")
    }
    
    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(M_a[i], M_b[i]) #c(range(grid)[1],range(grid)[2])
    
    XL=rng[i,1]-0.0001
    XR=rng[i,2]+0.0001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(M_a[i]:M_b[i], XL, XR, c, bdeg[1]))
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,M_a[i]:M_b[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,M_a[i]:M_b[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
  }
  
  ############### HERE WE CREATE THE BASIS FOR THE t VARIABLE In B(t)
  
  rng_t= c(min(M_a), max(M_b)) # THE KNOTS FOR THE FUNCTIONAL COEFFICIENTS HAVE TO BE MAXIMUM 
  
  XL_t=rng_t[1]-0.0001
  XR_t=rng_t[2]+0.0001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(min(M_a):max(M_b), XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT  
  for (i in 1:N) {
    
    # DATA BASIS
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
    
  }
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n_knots=length(breaks)
  
  
  # Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n_knots-2):(n_knots+2)))
  
  for (i in 1:N) {
    
    
    # if (i!=N) {
    
    ind_knots=max(which(Phi$knots<=(M_b[i]+0.0001)))
    
    # }
    #else{
    #   ind_knots=length(Phi$knots)
    # }
    
    # IDEA 1: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE)
    
    # breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_1=c(breaks[1:5],rep(0,n_knots-n),breaks[6:(n-5)],breaks[(n-4):n])
    # 
    # n_breaks_no_vc_1=length(breaks_no_vc_1)
    # 
    # Phi_aux_1=create.bspline.basis(breaks=breaks_no_vc_1,norder=bdeg[2]+1,dropin=c(1:5,(n_breaks_no_vc_1-2):(n_breaks_no_vc_1+2)))
    
    # IDEA 2: TOMAR LOS NODOS DEL INTERVALO
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_2=c(breaks[1:2],rep(0,n_knots-n),breaks[3:(n-2)],breaks[(n-1):n])
    # 
    # n_breaks_no_vc_2=length(breaks_no_vc_2)
    # 
    # Phi_aux_2=create.bspline.basis(breaks=breaks_no_vc_2,norder=bdeg[2]+1,dropin=c(1:7,(n_breaks_no_vc_2):(n_breaks_no_vc_2+2)))
    
    # IDEA 3: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE) PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    tag=3
    Phi_aux_3=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    
    # IDEA 4: TOMAR LOS NODOS DEL INTERVALO PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # tag=4
    # 
    # Phi_aux_4=create.bspline.basis(breaks=breaks, norder=bdeg[2]+1,dropin=c(1:2,(n):(n+2)))
    
    
    PROD=Simpson(L_X_aux[[i]],Phi_aux_3,rng = c(M_a[i],M_b[i]),sub = sub)/(M_b[i]-M_a[i]+1)
    
    if (tag==3 || tag==4) {
      PROD=cbind(PROD,matrix(0, nrow = dim(PROD)[1], ncol=nbasis[2]-dim(PROD)[2]))
    }
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}
#
################### SAME AS THE PREVIOUS FUNCTION BUT USING THE ADAPTIVE APPORACH, SEE RODRÍGUEZ ET AL (2019)

Data2B_simpson_ad=function(X, M, nbasis=c(30,30,30), bdeg=c(3,3,3), sub=25, ndb=25){
  
  require(fda)
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  c3=nbasis[3]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }

    
    ###############
    
    rng[i,]= c(1, M[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1]
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))
    
    ######### 
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d_ad(aux,2,c1,ndb = ndb)
    
    aux_3=XZG2theta_1d_ad(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,1:M[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,1:M[i]])-L_y[[i]]))
    
    ###############
    
    c_t=c2-bdeg[2]
    
    nam_Phi <- paste("Phi", i, sep = "_")
    L_Phi[[i]]=assign(nam_Phi, bspline(1:M[i], XL, XR, c_t, bdeg[2]))
    
  }
  
  #######
  
  xlim_T <- c(min(M), max(M))
  XL_T=xlim_T[1]-0.001
  XR_T=xlim_T[2]+0.001
  
  c_T=c3-bdeg[3]
  
  B_T=bspline(M, XL_T, XR_T, c_T, bdeg[3])
  
  #################
  
  
  for (i in 1:N) {
    
    #
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
    
    #
    
    breaks=L_Phi[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_Phi_aux <- paste("Phi_aux", i, sep = "_")
    L_Phi_aux[[i]]=assign(nam_Phi_aux, create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2))))
    
  }
    
  for (i in 1:N) {
    
    PROD=Simpson(L_X_aux[[i]],L_Phi_aux[[i]],B_T$B[i,],rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, B_T=B_T, B_Phi=L_Phi)
  
}


################### THE NEXT FUNCTIONS GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL:

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 2D CASE, SEE LEE AND DURBÁN REGUERA (2010)

B2XZG <-function (B, pord=c(2,2),c=c(10,10)) {
    
    c1=c[1]
    c2=c[2]
    
    c1c2 = ncol(B)
    if (c1c2!=c1*c2) {
      stop("c1 * c2 must me equal to the number of colums of B")
    }
    
    D_1 = diff(diag(c1), differences=pord[1])
    D_2 = diff(diag(c2), differences=pord[2])
    
    P1.svd = svd(crossprod(D_1))
    P2.svd = svd(crossprod(D_2))
    
    U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
    U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
    d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
    
    U_2s = (P2.svd$u)[,1:(c2-pord[2])] # eigenvectors
    U_2n=((P2.svd$u)[,-(1:(c2-pord[2]))])
    d2 = (P2.svd$d)[1:(c2-pord[2])]  # eigenvalues
    
    
    T_n=kronecker(U_1n,U_2n)
    
    AUX_1=kronecker(U_1n,U_2s)
    AUX_2=kronecker(U_1s,U_2n)
    AUX_3=kronecker(U_1s,U_2s)
    
    T_s=cbind(AUX_1,AUX_2,AUX_3)
    
    
    Z = B%*%T_s
    X = B%*%T_n
    
    ####
    
    d_1s=diag(P1.svd$d)[1:(c1-pord[1]),1:(c1-pord[1])]
    d_2s=diag(P2.svd$d)[1:(c2-pord[2]),1:(c2-pord[2])]
    
    T_1=kronecker(diag(pord[1]),d_2s)
    T_2=matrix(0,nrow = pord[2]*(c1-pord[1]),ncol = pord[2]*(c1-pord[1]))
    T_3=kronecker(diag(c1-pord[1]),d_2s)
    
    T_21=cbind(T_1,matrix(0,nrow=dim(T_1)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(T_1)[2]))
    T_22=cbind(matrix(0,nrow=dim(T_2)[1],ncol=dim(T_1)[2]),T_2,matrix(0,nrow=dim(T_2)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(T_1)[2]-dim(T_2)[2]))
    T_23=cbind(matrix(0,nrow=((c2-pord[2])*(c1-pord[1])),ncol=(c1*c2-pord[1]*pord[2])-dim(T_3)[2]),T_3)
    
    H_1=matrix(0,nrow = pord[1]*(c2-pord[2]),ncol = pord[1]*(c2-pord[2]))
    H_2=kronecker(d_1s,diag(pord[2]))
    H_3=kronecker(d_1s,diag(c2-pord[2]))
    
    H_11=cbind(H_1,matrix(0,nrow=dim(H_1)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(H_1)[2]))
    H_12=cbind(matrix(0,nrow=dim(H_2)[1],ncol=dim(H_1)[2]),H_2,matrix(0,nrow=dim(H_2)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(H_1)[2]-dim(H_2)[2]))
    H_13=cbind(matrix(0,nrow=((c2-pord[2])*(c1-pord[1])),ncol=(c1*c2-pord[1]*pord[2])-dim(H_3)[2]),H_3)
    
    L_2=rbind(T_21,T_22,T_23)
    L_1=rbind(H_11,H_12,H_13)
    
    t_2=diag(L_1)
    t_1=diag(L_2)
    
    G=list(t_1,t_2)
    names(G)=c("t_1","t_2")
    
    T=cbind(T_n,T_s)
    
    ####
    
    list(X = X, Z = Z, G=G, T=T, d1 = d1, d2=d2, D_1 = D_1, D_2=D_2, U_1n=U_1n, U_1s=U_1s, U_2n=U_2n, U_2s=U_2s, T_n=T_n, T_s=T_s, t_1=t_1, t_2=t_2)
  }

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 1D CASE, SEE LEE AND DURBÁN REGUERA (2010)
B2XZG_1d <-function (B, pord=c(2),c=c(10)) {
    
    c1=c[1]
    
    
    
    
    D_1 = diff(diag(c1), differences=pord[1])
    
    P1.svd = svd(crossprod(D_1))
    
    U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
    U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
    d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
    
    T_n=(U_1n)
    
    T_s=U_1s
    
    
    Z = B%*%T_s
    X = B%*%T_n
    
    ####
    G=list(d1)

    
    # G=list(t_1,t_2)
    # names(G)=c("t_1","t_2")
    # 
    T=cbind(T_n,T_s)
    
    ####
    
    list(X = X, Z = Z, G=G, T=T, d1 = d1, D_1 = D_1, U_1n=U_1n, U_1s=U_1s, T_n=T_n, T_s=T_s)
  }

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 1D CASE WITH THE ADAPTIVE APPROACH, SEE LEE AND DURBÁN REGUERA (2010)
B2XZG_1d_ad <-function (B, pord=c(2),c=c(10), ndb=25) {
  
  c1=c[1]
  
  D_1 = diff(diag(c1), differences=pord[1])
  
  P1.svd = svd(crossprod(D_1))
  
  U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
  U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
  d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
  
  T_n=(U_1n)
  
  T_s=U_1s
  
  Z = B%*%T_s
  X = B%*%T_n
  
  bdeg=3
  
  v.comp_1 <- 1:ncol(U_1s)/ncol(U_1s)
  C_1 <- t(SOPExamples:::bbase(v.comp_1, min(v.comp_1), max(v.comp_1), ndx = ndb-bdeg, bdeg = bdeg)$B)
  # One random (smooth) component
  Lambda_1 <- vector("list", length = 1) 
  # Precision matrices associated with that component
  Lambda_1[[1]] <- vector("list", length = nrow(C_1))
  for (i in 1:nrow(C_1)) {
    Lambda_1[[1]][[i]] <-  t(P1.svd$u[,1:(c1-pord[1])]) %*% t(D_1) %*% diag(C_1[i,]) %*% D_1 %*% P1.svd$u[,1:(c1-pord[1])]
    }
  
  
  ####
  G=Lambda_1
  T=cbind(T_n,T_s)
  ####
  
  list(X = X, Z = Z, G=G, T=T, d1 = d1, D_1 = D_1, U_1n=U_1n, U_1s=U_1s, T_n=T_n, T_s=T_s)
}

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 1D CASE WHEN WE HAVE TWO VARIABLES, SEE LEE AND DURBÁN REGUERA (2010)
B2XZG_1d_2_var_intercept <-function (B, pord=c(2,2),c=c(20,20)) {
  
  # with intercept
  
  c1=c[1]
  
  D_1 = diff(diag(c1), differences=pord[1])
  
  P1.svd = svd(crossprod(D_1))
  
  U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
  U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
  d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
  
  # T_n1=as.matrix((U_1n)[,-1])
  T_n1=as.matrix(U_1n)
  
  T_s1=U_1s
  
  c2=c[2]
  
  D_2 = diff(diag(c2), differences=pord[2])
  
  P2.svd = svd(crossprod(D_2))
  
  U_2s = (P2.svd$u)[,1:(c2-pord[2])] # eigenvectors
  U_2n=((P2.svd$u)[,-(1:(c2-pord[2]))])
  d2 = (P2.svd$d)[1:(c2-pord[2])]  # eigenvalues
  
  # T_n2=as.matrix((U_2n)[,-1])
  T_n2=as.matrix(U_2n)
  
  T_s2=U_2s
  
  # P=matrix(0,nrow = (1+dim(crossprod(D_1))[1]+dim(crossprod(D_2))[1]),ncol=(1+dim(crossprod(D_1))[1]+dim(crossprod(D_2))[1]))
  # P[2:(c1+1),2:(c1+1)]=crossprod(D_1)
  # P[(c1+2):(c1+c2+1),(c1+2):(c1+c2+1)]=crossprod(D_2)
  # 
  # aux=svd(P)
  # 
  # U_1=aux$u[2:(c1+1),2:(c1+1)]
  # U_1n=U_1[,-(1:(c1-pord[1]))]
  # U_1s=U_1[,1:(c1-pord[1])]
  # 
  # U_2=aux$u[(c1+2):(c1+c2+1),(c1+2):(c1+c2+1)]
  # U_2n=U_2[,-(1:(c2-pord[2]))]
  # U_2s=U_2[,1:(c2-pord[2])]
  # 
  # U_armada=matrix(0,nrow = dim(P)[1],ncol=dim(P)[2])
  # U_armada[2:(c1+1),2:(c1+1)]=P1.svd$u
  # U_armada[(c1+2):(c1+c2+1),(c1+2):(c1+c2+1)]=P2.svd$u
  # 
  # all.equal(aux$u,U_armada)
  # 
  # 
  # all.equal(aux$d,c(P1.svd$d,P2.svd$d))
  
  
  
  T_n=rbind(cbind(T_n1,matrix(0,nrow = dim(T_n1)[1],ncol=dim(T_n2)[2])),cbind(matrix(0,nrow = dim(T_n2)[1],ncol=dim(T_n1)[2]),T_n2))
  T_n=rbind(diag(dim(T_n)[2]+1)[1,],cbind(0,T_n))
  
  T_s=rbind(cbind(T_s1,matrix(0,nrow = dim(T_s1)[1],ncol=dim(T_s2)[2])),cbind(matrix(0,nrow = dim(T_s2)[1],ncol=dim(T_s1)[2]),T_s2))
  T_s=rbind(0,T_s)
  
  # T=cbind(T_n,T_s)
  
  T=cbind(T_n[,-c(pord[1],pord[1]+pord[2])],T_s)
  
  # all.equal(t(T)%*%T,diag(c1+c2-1)) # THIS IS TRUE
  # all.equal(T%*%t(T),diag(c1+c2+1)) # BUT THIS IS NOT
  
  Z = B%*%T_s
  X = B%*%T_n
  
  X_wrong=X
  T_wrong=cbind(T_n,T_s)
  
  X=X[,-c(pord[1],pord[1]+pord[2])]
  
  
  
  ####
  
  
  G=list(c(d1,rep(0,length(d2))),c(rep(0,length(d1)),d2))
  
  
  # G=list(t_1,t_2)
  # names(G)=c("t_1","t_2")
  # 
  
  #T=rbind(cbind(T_n1,T_s1,matrix(0,nrow = dim(T_n1)[1],ncol=dim(T_n2)[2]+dim(T_s2)[2])),cbind(matrix(0,nrow = dim(T_n2)[1],ncol=dim(T_n1)[2]+dim(T_s1)[2]),T_n2,T_s2))
  
  
  ####
  
  list(X = X, Z = Z, G=G, T=T, T_wrong=T_wrong, X_wrong=X_wrong)
}


##################### FUNCTIONS FOR ESTIMATING THE COEFFICIENTS OF THE MODEL

# THIS FUNCTION FIT THE MODEL USING OUR APPROACH 
# AND RECOVER THE ORIGINAL THETA COEFFICIENTS THAT ARE NECESSARY FOR CALCULATE THE ESTIMATED FUNCTIONAL COEFFICIENT BETA FOR THE 2D CASE.
XZG2theta=function(X, Z, G, T, y, family=gaussian(), offset = NULL){
  require(SOP)
  
  if (dim(X)[1]!=dim(Z)[1]) {
    stop("'X' and 'Z'must have same numbers of rows")
  }
  
  # if ( dim(Z)[2]!=length(G[[1]]) || dim(Z)[2]!=length(G[[2]]) ) {
  #   stop("The number of columns of 'Z' must be equal to the length of 'G'")
  # }
  
  w = as.vector(rep(1,dim(X)[1]))
  
  fit=sop.fit(X = X, Z = Z, G = G, 
              y = y, family = family,
              control = list(trace = FALSE),offset = offset)
  
  if (dim(fit$Vp)[1]-dim(T)[1]==0) {
  
    theta_aux=c(fit$b.fixed,fit$b.random) 
    
  }else{
  
    aux=dim(fit$Vp)[1]-dim(T)[1]
    
    theta_aux=c(fit$b.fixed[-(1:aux)],fit$b.random) # ESTE [1:4] FUE AGREGADO PARA QUE AL AGREGAR OTRAS VARIABLES LINEALES SOLO COGIERA LA PARTE FUNCIONAL 
    
  }
  
  theta=T %*% theta_aux
  
  
  
  if (dim(fit$Vp)[1]==dim(T)[1]) {
    
    covar_theta=T %*% fit$Vp %*% t(T)
    std_error_theta=sqrt(diag(T %*% fit$Vp %*% t(T)))
    std_error_non_functional=NULL
    p_values=NULL
  }else{
    
    covar_theta=T %*% fit$Vp[-(1:aux),-(1:aux)] %*% t(T)
    std_error_theta=sqrt(diag(T %*% fit$Vp[-(1:aux),-(1:aux)] %*% t(T)))
    std_error_non_functional=sqrt(diag(fit$Vp[1:aux,1:aux]))
    WALD= (fit$b.fixed[(1:aux)]/std_error_non_functional)  
    p_values=2*pnorm(abs(WALD),lower.tail = FALSE)
    }
  
  
  
  list (fit=fit, theta=theta, std_error_theta=std_error_theta, std_error_non_functional=std_error_non_functional, covar_theta=covar_theta, p_values=p_values)
}

# THIS FUNCTION FIT THE MODEL USING OUR APPROACH 
# AND RECOVER THE ORIGINAL THETA COEFFICIENTS THAT ARE NECESSARY FOR CALCULATE THE ESTIMATED FUNCTIONAL COEFFICIENT BETA FOR THE 1D CASE.
XZG2theta_1d=function(X, Z, G, T, y, family=gaussian()){
  require(SOP)
  
  fit=sop.fit(X = X, Z = Z, G = G, 
                y = y, family = family, weights = rep(1,dim(X)[1]), 
                control = list(trace = FALSE))
  
  theta_aux=c(fit$b.fixed,fit$b.random)
  
  theta=T %*% theta_aux
  
  list (fit=fit,theta=theta)
}

# THIS FUNCTION FIT THE MODEL USING OUR APPROACH 
# AND RECOVER THE ORIGINAL THETA COEFFICIENTS THAT ARE NECESSARY FOR CALCULATE THE ESTIMATED FUNCTIONAL COEFFICIENT BETA FOR THE 1D CASE 
# WITH THE ADAPTIVE APPROACH.
XZG2theta_1d_ad=function(X, Z, G, T, y, family=gaussian()){
  require(SOP)
  require(SOPExamples)
  
  fit=fit.SOP(X = X, Z = Z, Lambda = G, 
              y = y, family = family ,trace = 0, diagonal = 0)
  
  theta_aux=fit$coeff
  
  theta=T %*% theta_aux
  
  list (fit=fit,theta=theta)
}
