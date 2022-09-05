
# main code -------------------------------------------------------------

dg_fcn <- function(x, mu, d, A) {
  # gradient of generator g_{\theta}(z) over theta
  p <- length(x)
  p2 <- 2*p +p*(p-1)/2
  
  x0 <- c(A%*%(x-mu)) /d
  
  dmat <- matrix(0, p, p2)
  dmat[,1:p] <- -A/d          #mu   
  dmat[,p+1:p] <- -mdiag(x0)  #log d
  
  if (p>1) {
    for (j in 2:p) {
      dmat[j, 2*p+(j-2)*(j-1)/2 +1:(j-1)] <- (x-mu)[1:(j-1)] /d[j]  #A
    }
  }
  
  list(dmat=dmat)
}

trunc_basis <- function(x, knots, degr=2, grad=T) {
  # x is scalar
  k <- length(knots)
  if (grad) {
    if (degr==2) {
      vec <- rep(0, 2*k)
      dvec <- rep(0, 2*k)
      for (j in 1:k) {
        if (x>knots[j]) {
          vec[j] <- x-knots[j]
          dvec[j] <- 1
          vec[j+k] <- (x-knots[j])^2 
          dvec[j+k] <- 2*(x-knots[j])
        }
      }
      vec <- c(x, vec, x^2)  # remove intercept
      dvec <- c(1, dvec, 2*x)
    } else { # degr=1
      vec <- rep(0, k)
      dvec <- rep(0, k)
      for (j in 1:k) {
        if (x>knots[j]) {
          vec[j] <- x-knots[j]
          dvec[j] <- 1
        }
      }
      vec <- c(x, vec) # remove intercept
      dvec <- c(1, dvec)
    }
  } else {
    if (degr==2) {
      vec <- rep(0, 2*k)
      for (j in 1:k) {
        if (x>knots[j]) {
          vec[j] <- x-knots[j]
          vec[j+k] <- (x-knots[j])^2 
          
        }
      }
      vec <- c(x, vec, x^2)  # remove intercept
    } else { # degr=1
      vec <- rep(0, k)
      for (j in 1:k) {
        if (x>knots[j]) {
          vec[j] <- x-knots[j]
        }
      }
      vec <- c(x, vec) #remove intercept
    }
    dvec <- NULL
  }
  
  list(vec=vec, dvec=dvec)  
}

linear_basis <- function(x, knots, sub_gamma, mat, dmat) {
  # linear spline basis (p-dim input)
  # x is vector of length p
  p <- length(x)
  k <- dim(knots)[1]
  k2 <- 2*k+2
  
  if (is.null(sub_gamma)) {
    for (i in 1:p) {
      out_basisse <- trunc_basis(x[i], knots[,i], degr=2, grad=F)
      mat[(i-1)*k2 + (1:k2)] <- out_basisse$vec
    }
    dmat <- NULL
  } else {
    for (i in 1:p) {
      out_basisse <- trunc_basis(x[i], knots[,i], degr=2)
      mat[(i-1)*k2 + (1:k2)] <- out_basisse$vec
      dmat[1, i] <- out_basisse$dvec %*% 
        sub_gamma[(i-1)*k2 + (1:k2)]
    }
  }
  
  list(mat=mat, dmat=dmat)  
}

cross_basis <- function(x, knots, sub_gamma, mat, dmat) {
  # interaction spline basis (bivariate input)
  # x is vector of length 2
  k <- dim(knots)[1]
  p <- length(x)
  
  # mat, dmat are large empty matrices
  # we do not generate those matrices here to save running time
  pos <- 0
  if (is.null(sub_gamma)) {
    for (i in 1:(p-1)) {
      out_basisse1 <- trunc_basis(x[i], knots[,i], degr=1, grad=F)
      for (j in (i+1):p) {
        out_basisse2 <- trunc_basis(x[j], knots[,j], degr=1, grad=F)
        
        vec1 <- out_basisse1$vec
        vec2 <- out_basisse2$vec
        
        mat[pos*(1+k)^2 + (1:(1+k)^2)] <-  # each iter produce (1+k)^2 new entries
          c(vec1*vec2[1], vec1[1]*vec2[-1], tprod(vec1[-1], vec2[-1]))
        pos <- pos + 1
      }
    }
    dmat <- NULL
  } else { # sub_gamma is not NULL, grad will be computed
    for (i in 1:(p-1)) {
      out_basisse1 <- trunc_basis(x[i],knots[,i], degr=1)
      for (j in (i+1):p) {
        out_basisse2 <- trunc_basis(x[j],knots[,j], degr=1)
        
        vec1 <- out_basisse1$vec
        vec2 <- out_basisse2$vec
        
        dvec1 <- out_basisse1$dvec
        dvec2 <- out_basisse2$dvec
        
        #
        mat[pos*(1+k)^2 + (1:(1+k)^2)] <-  # each iter produce (1+k)^2 new entries
          c(vec1*vec2[1], vec1[1]*vec2[-1], tprod(vec1[-1], vec2[-1]))
        # tprod equivalent to c(t(outer(vec1, vec2)))
        dmat[, i] <-  dmat[, i] + 
          c(dvec1*vec2[1], dvec1[1]*vec2[-1], tprod(dvec1[-1], vec2[-1])) %*% 
          sub_gamma[pos*(1+k)^2 + (1:(1+k)^2)]
        dmat[, j] <-  dmat[, j] + 
          c(vec1*dvec2[1], vec1[1]*dvec2[-1], tprod(vec1[-1], dvec2[-1])) %*% 
          sub_gamma[pos*(1+k)^2 + (1:(1+k)^2)]
        
        pos <- pos + 1
      }
    }
  } 
  
  
  list(mat=mat, dmat=dmat)
}

G_update <- function(theta, x, z, gamma, knots=K, Loss="rKL", requires_grad=F) {
  # Pass raw data through the spline "network" and generate
  # ...spline basis for real and fake data using current estimates. 
  # Gradient for the generator is computed in the process (chain rule).
  
  #-----------------------------------------------------------------------------
  # *theta* is the generator parameter (in choleskey form)
  # *x* is the real Gaussian+comtamination data 
  # *z* is standard Gaussian noise
  # *gamma* is the discriminator parameter \gamma obtained in previous run
  # *knots* is the spline knots. Default is object "K" which is a global object.
  # *Loss* could be "rKL", "JS", and "Hinge". General $f$-divergences can be
  # ...added and customized conveniently.
  #-----------------------------------------------------------------------------
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  m <- dim(z)[1]
  
  # deprecated
  # if (n != m) {
  #   stop("fake data size should match real data size.")
  #}
  
  # Obtained current estimated location and cov from choleskey factors
  mu2 <- theta[1:p]
  d2 <- exp(theta[p+1:p])
  A2 <- vec.lmat(theta[-(1:(2*p))], p)
  
  k <- dim(knots)[1]
  
  # Transform/standardize real data using current estimates
  p_comb <- choose(p,2) # numebr of combinations of pairs
  xz <- t( A2%*%(t(x)-mu2) /d2 )
  k2 <- p*(2*k+2) + p_comb*(k+1)^2
  p2 <- 2*p + p*(p-1)/2
  
  # Spline basis (and gradient of theta) for standardized real data
  f1 <- matrix(0, n, 1+k2)
  grad0 <- matrix(0, n, p2)

  empty_bvmat <- rep(0, (1+k)^2*p_comb) # pre-allocation
  empty_bvdmat <- matrix(0, 1, p) # pre-allocation
  cross_gamma <- gamma[1 + p*(2*k+2) + (1:((1+k)^2*p_comb))]
  
  empty_uvmat <- rep(0, p*(2*k+2)) # pre-allocation
  empty_uvdmat <- matrix(0, 1, p) # pre-allocation
  if (requires_grad) {
    linear_gamma <- gamma[1 + (1:(p*(2*k+2)))]
  } else {
    linear_gamma <- NULL
  }
  
  for (i in 1:n) {
    # spline basis
    out.ubase <- linear_basis(xz[i,], knots, linear_gamma, 
                          empty_uvmat, empty_uvdmat)
    out.bbase <- cross_basis(xz[i,], knots, cross_gamma,
                          empty_bvmat, empty_bvdmat)
    f1[i,] <- c(1, out.ubase$mat, out.bbase$mat)
    
    if (requires_grad) {
      out.dg <- dg_fcn(x[i,], mu2, d2, A2)   #x[i], not xz[i]
      
      grad0[i,] <- 0 + out.ubase$dmat %*% out.dg$dmat + 
        out.bbase$dmat%*%out.dg$dmat
    }
  }
  
  # Spline basis for fake data
  f0 <- matrix(0, m, 1+k2)
  linear_gamma <- NULL # do not need gradient for fake data
  cross_gamma <- NULL
  
  for (i in 1:m) {
    # spline basis
    out.ubase <- linear_basis(z[i,], knots, linear_gamma, 
                          empty_uvmat, empty_uvdmat)
    out.bbase <- cross_basis(z[i,], knots, cross_gamma,
                          empty_bvmat, empty_bvdmat)
    f0[i,] <- c(1, out.ubase$mat, out.bbase$mat)
    #--------------------------------------------------------------------------
    # Note that no theta gradient is computed for the fake data batch.
    # This is because theta only appears in the back transformation of
    # ...the real data. In fact, the "generator" is trying to "generate" 
    # ...back standard Gaussian data using the real data and the estimated 
    # ...location and shape so the fake gaussian noise is indistinguishable 
    # ...from the standardized real data. 
    # For details see Section 5.1 in https://arxiv.org/2112.12919.pdf
    #--------------------------------------------------------------------------
  }
  
  if (requires_grad) {
    # Compute gradient for theta and evaluate loss
    mat <- rbind(f1, f0)
    tr <- c(rep(1,n), rep(0,m))
    eta <- c(mat%*%gamma)
    
    if (Loss == "rKL") {
      d <- exp(-eta[tr==1])
      val <- 1-mean(exp(-eta[tr==1])) - mean(eta[tr==0])
    } else if (Loss == "Hinge") {
      d <- (1-eta[tr==1] > 0)
      val <- 2-mean(min(eta[tr==1], 1)) - mean(min(-eta[tr==0], 1))
    } else if (Loss == "JS") {
      val <- 2*log(2)-
        mean(log(1+exp(-eta[tr==1]))) - 
        mean(log(1+exp(eta[tr==0])))
      d <- 1/(1+exp(eta[tr==1]))
    } 
    
    grad <- apply(d*grad0, 2, mean) # chain rule
  } else {
    grad <- NULL
    val <- NA
  }
  
  list(gradient=grad, f1=f1, f0=f0, value=val)
}

D_update <- function(x, z, pen, gamma, Loss="rKL", L1_pen=FALSE, solver="MOSEK") {
  # CVX framework 
  # Package "CVXR" and "Rmosek" required
  
  #------------------------------------------------------------------------
  # *x* is the spline basis of real data
  # *z* is the spline basis of generated data
  # *pen* is $\lambda$ the penalty tuning parameter
  # *gamma* is the $\gamma$ parameter in previous run. It is returned if
  # ...the optimizer failed for this run.
  # *Loss* could be "rKL", "JS", and "Hinge". General $f$-divergences can be
  # ...easily added. The smooth hinge is supported but theory is not studied.
  # *L1_pen* is a logical parameter, indicating l1 penalty or l2 penalty.
  #------------------------------------------------------------------------
  
  gam <- Variable(dim(x)[2]) # gam for \gamma, the discriminator parameter
  
  wx <- x%*%gam
  wz <- z%*%gam

  if (Loss == "rKL") {
    ewx <- exp(-wx)
    #val <- max(mean(ewx) + mean(wz) -1, -10)
    val <- max(mean(ewx), -10) + max(mean(wz) -1, -10) 
    # clipped at a (large) constant level to prevent unbounded optim problem in
    # early phase of training when fake and real data could be separable. 
  } else if (Loss == "Hinge") {
    val <- mean(max_elemwise(0, 1-wx)) + mean(max_elemwise(0, 1+wz)) 
  } else if (Loss == "JS") {
   val <- mean(CVXR::logistic(-wx)) + mean(CVXR::logistic(wz))
  } 
  
  if (L1_pen == FALSE) {
    # L2 penalty on gamma for linear and interactions (different lambda)
    penalty <- pen[1]*cvxr_norm(gam[2:(4*P + 1)], 2) + 
     pen[2]*cvxr_norm(gam[-(1:(4*P + 1))], 2)
  } else {
    # L1 penalty on gamma
    penalty <- pen*cvxr_norm(gam[-1], 1) #recommanded
  }
  
  if (! solver %in% CVXR::installed_solvers()) {
    solver <- "ECOS"
  }
  
  result <- solve(Problem(Minimize(val + penalty)), 
                  solver=solver) 
  
  if (result$status == "optimal") {
    coef <- result$getValue(gam)
  } else {
      coef <- gamma 
  }
  return(list(coef=coef, optimized=result, value=result$value))
}




fGAN <- function(x, pen=0, L1_pen=F,
                 kendall_samp=1e2, 
                 iter=100, init_mu=NULL, init_Sigma=NULL, 
                 decay_lr=FALSE, lr=.1,  a=.1, 
                 batch_x=NULL, batch_z=NULL, 
                 mu0=NULL, Sigma0=NULL, optim=optim,
                 Loss=Loss) {
  
  # Model settings 
  n <- dim(x)[1]
  p <- dim(x)[2]
  

  if (is.null(batch_x)) batch_x <- n
  if (is.null(batch_z)) batch_z <- batch_x
  
  
  # Initialization
  # START ---------------------------------------------------------------------
  
  if (is.null(init_mu)) {
    # Median init
    mu <- apply(x, 2, median) 
  } else {
    mu <- init_mu
  }
  
  if (is.null(init_Sigma)) {
  # Kendall's tau estimator for Cov init, using only a small size of data
  S_hat <- SS(x)  # diag(apply(x, 2, mad))
  K_hat <- cor(x[1:kendall_samp,], method = 'kendall', use='pairwise')
  ssig <- S_hat^.5 %*% K_hat %*% S_hat^.5
  
  # or use MCD, "rrcov" package requried
  # ssig <- rrcov::CovMcd(x)@cov
  
  } else {
    ssig <- init_Sigma
  }

  
  ## Choleskey decomp on initial covariance matrix "ssig"
  out.chol <- chol(ssig) 
  d <- diag(out.chol) # diagonal of the upper triangular matrix
  B <- t(out.chol /d) # "standardized" upper triangular matrix
  A <- solve(B) # get the inverse for back transformation.
  
  # ---------------------------------------------------------------------
  # A %*% ((x - mu)/d) gives back standard normal if
  # x ~ N(mu, solve(A) %*% diag(d^2) %*% t(solve(A))).
  # We transform read data back to standard scale and shape.
  # For details see Section 5.1 in https://arxiv.org/pdf/2112.12919.pdf 
  # ---------------------------------------------------------------------
  
  init <- c(mu, log(d), lmat.vec(A)) # initial estimates (choleskey parameters)
  theta <- init
  ## Get fixed knots for splines (while be transformed in training)
  knots <- matrix(rep(c(-1, 0 ,1), p), ncol=p)
  ## Or use "apply(z, 2, quantile, c(.2,.5,.8))" (make sure 'z' exists)

  k <- dim(knots)[1] # number of knots
  K <<- knots # Preserved global variable 
  
  # Discriminator does not require initialization 
  gamma <- rep(NA, 1+p*(2*k+2)+ choose(p,2)*(k+1)^2) 

  if (! optim %in% CVXR::installed_solvers()) {
    warning("Convex solver not installed, using default ECOS solver instead. Speed could be substantially slower. MOSEK solver is recommanded and a free academic license is required. See https://www.mosek.com/products/academic-licenses/")
    optim <- "ECOS"
  }
  
  
  # END Initialization --------------------------------------------------------
  
  
  # Training
  Theta <- matrix(NA, nrow=iter+1, ncol=length(theta))
  Theta[1, ] <- theta
  Lr <- lr # object *Lr* could change, object *lr* is protected.

  cat("\nTotal number of discriminator parameters are:", length(gamma))
  cat("\nFollowing numbers are estimation erros measured in L2 (Location), Operator, and Frobenius norms")
  
  for (i in 1:iter) {
    # Discriminator update 
    # START -------------------------------------------------------------------
    ## Generate noise, sample random batch from real data.
    z <- matrix(rnorm(batch_z*p), batch_z, p, byrow=T) 
    ind <- sample(n, batch_x)
    
    ## Pass data to generator and discriminator, get spline basis.
    out_basis <- G_update(theta, x[ind,], z, gamma, requires_grad = F)
    x_basis <- out_basis$f1
    z_basis <- out_basis$f0
    
   
    #------------------------------------------------------------------
    # Requires CVXR convex optimization framework (R package "CVXR"),
    # see https://cvxr.rbind.io 
    # Requires MOSEK solver, using the default ECOS solver will
    # ...substantially slow down the training.
    # For free academic license of MOSEK, 
    # see https://www.mosek.com/products/academic-licenses/
    #------------------------------------------------------------------
    D <- D_update(x_basis, z_basis, gamma=gamma, pen=pen,
                  Loss=Loss, L1_pen=L1_pen, solver=optim)
    gamma <- D$coef
    
    # END ---------------------------------------------------------------------
    
    
    # Generator update 
    # START -------------------------------------------------------------------
    ## Decaying learning rate
    if (decay_lr == T) {
      lr <- Lr*exp(-.05*i)
    }

    ## Generate random noise
    z <- matrix(rnorm(batch_z*p), batch_z, p, byrow=T)
    
    ## Pass real data and noise to generator, 
    ## update generator parameters with gradient descent.
    G <- G_update(theta, x[ind,], z, gamma, Loss=Loss, requires_grad=T) 
    grad <- G$gradient
    grad <- soft(grad, a) # gradient clip
    theta <- theta - lr*grad
  
    Theta[i+1, ] <- theta
    
    # END ---------------------------------------------------------------------
    
    # Error evaluation if true values are known
    if (is.null(mu0) | is.null(Sigma0)) {
      cat("\n[", as.character(Sys.time()),
          "Iter", i, "] ",
          "Gloss:", format(G$value, digits=5), 
          "Dloss:", format(D$value, digits=5)
      )
    } else {
      temp_theta <- trans.theta(theta, p) # choleskey components back to cov matrix
      cat("\n[", as.character(Sys.time()),
          "Iter", i, "] ",
          "l2:", format(sqrt(sum((temp_theta$mu - mu0)^2)), digits=5), 
          "F:", format(sqrt(sum((temp_theta$ssig - Sigma0)^2)), digits=5), 
          "Op:", format(svd(temp_theta$ssig - Sigma0)$d[1], digits=5),
          "Gloss:", format(G$value, digits=5), 
          "Dloss:", format(D$value, digits=5)
      )      
    }
    cat(" | cvx: time =", D$optimized$solve_time)
  }
  
  # Return training results 
  theta <- colMeans(tail(Theta, 3)) # take average on last 3 steps
  temp_theta <-  trans.theta(theta, p) 
  return(list(gamma=gamma, # discriminator linear coeff.
              thet=theta, # estimated generator parameters in choleskey form (a long vector)
              Thet=Theta, # generator parameter traj
              init_thet=init,# initial generator parameter
              real_spline_basis=x_basis, # spline basis for real data
              fake_spline_basis=z_basis, # spline basis for fake data
              mu_hat=temp_theta$mu, # estimated location
              Sigma_hat=temp_theta$ssig # # estimated covariance matrix
              )
         )
}



# Error bypass version for implementation on HPC cluster. 
fGAN_error_bypass <- 
  purrr::possibly(
    fGAN,
    otherwise=list(gamma = NA, thet = NA, Thet = NA, init_thet=NA)
  )

svd_error_bypass <- purrr::possibly(svd, otherwise=list(d=NA, u=NA, v=NA))
