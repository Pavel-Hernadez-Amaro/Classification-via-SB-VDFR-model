modifica.fpca <-function (data.m, grids.u, muhat, eigenvals, eigenfuncs, sig2hat, 
          K) 
{
  temp <- table(data.m[, 1])
  n <- length(temp)
  m.l <- as.vector(temp)
  result <- matrix(0, n, K)
  N <- length(grids.u)
  evalmat <- if(K>1) diag(eigenvals[1:K]) else eigenvals[1:K]
  current <- 0
  eigenfuncs.u <- t(eigenfuncs)
  data.u <- matrix(as.numeric(as.vector(data.m[, -1])), nrow = nrow(data.m[, 
                                                                           -1]), ncol = ncol(data.m[, -1]))
  for (i in 1:n) {
    Y <- as.vector(data.u[(current + 1):(current + m.l[i]), 
                          1])
    meastime <- data.u[(current + 1):(current + m.l[i]), 
                       2]
    gridtime <- ceiling(meastime)
    muy <- muhat[gridtime]
    Phiy <- matrix(eigenfuncs.u[gridtime, 1:K], ncol = K)
    Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
    temp.y <- matrix(Y - muy)
    result[i, ] <- evalmat %*% t(Phiy) %*% solve(Sigy, temp.y)
    current <- current + m.l[i]
  }
  return(result)
}

###

myfpca.mle <- function (data.m, M.set, r.set, ini.method = "EM", basis.method = "bs", 
          sl.v = rep(0.5, 10), max.step = 50, grid.l = seq(0, 1, 0.01), 
          grids = seq(0, 1, 0.002)) 
{
  tol <- 0.001
  cond.tol <- 1e+10
  if (length(sl.v) < max.step) {
    sl.v <- c(sl.v, rep(1, max.step - length(sl.v)))
  }
  else {
    sl.v <- sl.v[1:max.step]
  }
  if (basis.method == "bs") {
    basis.method <- "poly"
  }
  M.self <- 20
  basis.EM <- "poly"
  if (ini.method == "EM") {
    no.EM <- FALSE
  }
  else {
    no.EM <- TRUE
  }
  iter.EM.num <- 50
  sig.EM <- 1
  nmax <- max(table(data.m[, 1]))
  L1 <- min(as.numeric(data.m[, 3])) # MIN OBSERVATION POINT
  L2 <- max(as.numeric(data.m[, 3])) # MAX OBSERVATION POINT
  data.list <- fpca.format(data.m)
  n <- length(data.list)
  if (n == 0) {
    print("error: no subject has more than one measurements!")
    return(0)
  }
  data.list.new <- data.list
  for (i in 1:n) {
    cur <- data.list[[i]][[1]]
    temp <- (cur[, 2] - L1)/(L2 - L1)
    temp[temp < 1e-05] <- 1e-05
    temp[temp > 0.99999] <- 0.99999
    cur[, 2] <- temp
    data.list.new[[i]][[1]] <- cur
  }
  temp <- myLocLin.mean(data.list.new, n, nmax, grids)    ###   <----------------
  data.list.new <- temp[[1]]
  fitmu <- temp[[2]]
  result <- NULL
  for (k in 1:length(r.set)) {
    r.c <- r.set[k]
    #print(paste("r=", r.c))
    result.c <- fpca.fit(M.set, r.c, data.list.new, n, nmax, 
                         grid.l, grids, iter.EM.num, sig.EM, ini.method, basis.method, 
                         sl.v, max.step, tol, cond.tol, M.self, no.EM, basis.EM)
    if (is.vector(result.c[[1]])) {
      print("warning: see warning code")
    }
    grids.new <- grids * (L2 - L1) + L1
    temp <- result.c$eigenvec
    M.set.u <- M.set[M.set >= r.c]
    if (length(M.set.u) == 0) {
      temp <- NULL
    }
    else {
      for (j in 1:length(M.set.u)) {
        temp[[j]] <- temp[[j]]/sqrt(L2 - L1)
      }
    }
    result.c$eigenvec <- temp
    result.c$eigenval <- result.c$eigenval * (L2 - L1)
    result.c$grids <- grids.new
    result[[k]] <- result.c
  }
  mod.sele <- fpca.cv(result, M.set, r.set, tol = tol)
  temp <- mod.sele[[3]]
  index.r <- temp[1]
  index.M <- temp[2]
  cv.result <- mod.sele[[1]]
  con.result <- mod.sele[[2]]
  result.sele <- result[[index.r]]
  eigenf <- result.sele$eigenvec
  eigenf.sele <- eigenf[[index.M]][, , 1]
  eigenv <- result.sele$eigenval
  eigenv.sele <- eigenv[1, , index.M]
  sig <- result.sele$sig
  sig.sele <- sig[1, index.M]
  temp.model <- c(M.set[index.M], r.set[index.r])
  names(temp.model) <- c("M", "r")
  rownames(eigenf.sele) <- paste("eigenfunction", 1:r.set[index.r])
  names(eigenv.sele) <- paste("eigenvalue", 1:r.set[index.r])
  temp <- list(selected_model = temp.model, eigenfunctions = eigenf.sele, 
               eigenvalues = eigenv.sele, error_var = sig.sele^2, fitted_mean = fitmu, 
               grid = grids.new, cv_scores = cv.result, converge = con.result)
  return(temp)
}

environment(myfpca.mle) <- asNamespace('fpca')
  
myLocLin.mean <- function (data.list, n, nmax, grids = seq(0, 1, 0.002)) 
  {
    data.f <- Format.data(data.list, n, Nmax = nmax)
    Obs <- data.f[[1]]
    T <- data.f[[2]]
    N <- data.f[[3]]
    data <- TranMtoV(Obs, T, N)
    y <- data[, 2]
    t <- data[, 3]
    hmu.cv <- h.select(t, y, method = "cv")
    fitmu <- mysm.regression(t, y, h = hmu.cv, poly.index = 1,    ###      <-----------
                           eval.points = t)$estimate
    fitmu.s <- mysm.regression(t, y, h = hmu.cv, poly.index = 1,  ###      <-----------
                             eval.points = grids)$estimate
    data.s <- data
    data.s[, 2] <- data[, 2] - fitmu
    data.result <- data.list
    for (i in 1:n) {
      data.c <- matrix(data.s[data.s[, 1] == i, ], ncol = 3)
      data.result[[i]][[1]] <- data.c[, 2:3]
    }
    return(list(data.result, fitmu.s))
}

environment(myLocLin.mean) <- asNamespace('fpca')

mysm.regression <- function (x, y, h, design.mat = NA, model = "none", weights = NA, 
                             group = NA, ...) 
{
  if (!all(is.na(group))) 
    return(sm.ancova(x, y, group, h, model, weights = weights, 
                     ...))
  x.name <- deparse(substitute(x))
  if (isMatrix(x)) 
    x.names <- dimnames(x)[[2]]
  y.name <- deparse(substitute(y))
  opt <- sm.options(list(...))
  data <- sm.check.data(x = x, y = y, weights = weights, group = group, 
                        ...)
  x <- data$x
  y <- data$y
  weights <- data$weights
  group <- data$group
  nobs <- data$nobs
  ndim <- data$ndim
  opt <- data$options
  replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs)/ndim))
  rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = nobs, 
                  ndim = ndim)
  if (!((model %in% "none") | (model %in% "no effect") | (model %in% 
                                                          "no.effect") | (model %in% "linear"))) 
    stop("invalid setting for model argument.", call. = FALSE)
  if (model != "none") 
    replace.na(opt, test, TRUE)
  if (missing(h)) 
    h <- h.select(x = x, y = y, weights = weights, ...)
  else {
    if (length(h) != ndim) 
      stop("length(h) does not match size of x")
  }
  if (opt$nbins > 0) {
    if (!all(weights == 1) & opt$verbose > 0) 
      cat("Warning: weights overwritten by binning\n")
    if (!all(is.na(opt$h.weights))) 
      stop("use of h.weights is incompatible with binning - set nbins=0")
    bins <- binning(x, y, nbins = opt$nbins)
    x <- bins$x
    y <- bins$means
    weights <- bins$x.freq
    rawdata$devs <- bins$devs
    nx <- length(y)
  }
  else nx <- nobs
  replace.na(opt, h.weights, rep(1, nx))
  if (opt$panel && !require(rpanel)) {
    opt$panel <- FALSE
    cat("The rpanel package is not available.\n")
  }
  if (ndim == 1) {
    replace.na(opt, xlab, x.name)
    replace.na(opt, ylab, y.name)
    replace.na(opt, ngrid, 50)
    opt$period <- opt$period[1]
    if (opt$pch == ".") 
      replace.na(opt, cex, 1)
    else replace.na(opt, cex, 2/log(rawdata$nobs))
    if (!opt$panel) 
      est <- mysm.regression.1d(x, y, h, design.mat, model,     ###     <---------
                              weights, rawdata, options = opt)
    else {
      rp.smooth1(x, y, h, design.mat, model, weights, 
                 rawdata, opt)
    }
  }
  else {
    replace.na(opt, ngrid, 20)
    dimn <- x.names
    name.comp <- if (!is.null(dimn) & !all(dimn == "")) 
      dimn
    else {
      if (!is.null(attributes(x)$names)) 
        attributes(x)$names
      else outer(x.name, c("[1]", "[2]"), paste, sep = "")
    }
    replace.na(opt, xlab, name.comp[1])
    replace.na(opt, ylab, name.comp[2])
    replace.na(opt, zlab, y.name)
    if (all(is.na(opt$period))) 
      opt$period <- rep(NA, 2)
    if (!(length(opt$period) == 2)) 
      stop("the length of period should match the number of covariates.")
    if (opt$panel) 
      rp.smooth2(x, y, h, model, weights, rawdata, opt)
    else est <- sm.regression.2d(x, y, h, model, weights, 
                                 rawdata, options = opt)
  }
  if (opt$panel) 
    invisible()
  else {
    est$data <- list(x = x, y = y, opt$nbins, freq = weights)
    est$call <- match.call()
    invisible(est)
  }
}

environment(mysm.regression) <- asNamespace('sm')

mysm.regression.1d <- function (x, y, h, design.mat = NA, model = "none", weights = rep(1, 
                                                                                        length(x)), rawdata, options = list()) 
{
  opt <- sm.options(options)
  replace.na(opt, ngrid, 50)
  replace.na(opt, xlim, range(rawdata$x))
  replace.na(opt, ylim, range(rawdata$y))
  replace.na(opt, display, "none")
  replace.na(opt, col, "black")
  replace.na(opt, col.band, "cyan")
  replace.na(opt, col.points, "black")
  replace.na(opt, se, FALSE)
  hmult <- opt$hmult
  if (model == "none") {
    opt$band <- FALSE
    opt$test <- FALSE
  }
  else replace.na(opt, band, TRUE)
  band <- opt$band
  if (opt$add | opt$display %in% "none") 
    opt$panel <- FALSE
  r <- list(x = NA, y = NA, model.y = NA, se = NA, sigma = NA, 
            h = h * hmult, hweights = opt$h.weights, weights = weights)
  if (!opt$add & !(opt$display %in% "none")) 
    plot(rawdata$x, rawdata$y, xlab = opt$xlab, ylab = opt$ylab, 
         xlim = opt$xlim, ylim = opt$ylim, type = "n")
  if (!(opt$display %in% "none")) {
    opt1 <- opt
    opt1$test <- FALSE
    r <- smplot.regression(x, y, design.mat, h, r, model, 
                           weights, rawdata, options = opt1)
  }
  if (opt$test) 
    rtest <- sm.regression.test(x, y, design.mat, h, model, 
                                weights, rawdata, options = opt)
  if (!(any(is.na(opt$eval.points)))) 
    r <- sm.regression.eval.1d(x, y, design.mat, h, model, 
                               weights, rawdata, options = opt)
  else if ((opt$display %in% "none") & (model == "none")) {
    opt$eval.points <- seq(min(x), max(x), length = opt$ngrid)
    r <- sm.regression.eval.1d(x, y, design.mat, h, model, 
                               weights, rawdata, options = opt)
  }
  if (opt$test) 
    r <- list(eval.points = r$eval.points, estimate = r$estimate, 
              model.y = r$model.y, se = r$se, sigma = r$sigma, 
              h = r$h, hweights = r$hweights, weights = weights, 
              model = rtest$model, p = rtest$p, q.squared = rtest$q.squared)
  r
}
environment(mysm.regression.1d) <- asNamespace('sm')