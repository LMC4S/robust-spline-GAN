# resolve model settings and run experiments

train <- function(eps, p, n, Cov="Ar", sd=1, mu=NULL,  # model spec
                  Q="A",  # contamination spec
                  fold=1, # repeat times
                  pen_lvl=1, # penalty level 
                  L1_pen=F, # l1 penalty or l2 (logical)
                  kendall_size=500, # sample size for initial cov esti. 
                  lr=.1, # learning rate for generator
                  iter=100, # iteration times
                  decay_lr=F, # decaying lr or not(logical)
                  optim="MOSEK", # convex optimizer for discriminator
                  Loss=Loss, # loss. rKL, JS, or Hinge 
                  batch_x=n, batch_z=NULL, # training batch size 
                  seed=0) {
  
  n1 <- floor((1-eps)*n) # uncontaminated data size
  n2 <- n - n1 # outlier size
  
  # model and data preparation
  if (length(sd) == 1) {
    sd <- rep(sd, p)
  } else if (length(sd) != p) {
    stop("sd must be a scalar or a p-dim vector")
  }
  
  if (L1_pen) {
    # lambda for l1
    pen0 <- sqrt(log(p)/batch_x) 
  } else {
    # l2 lambdas are different for linear and interactions.
    pen0 <- c(sqrt(p/batch_x), sqrt(p^2/batch_x)) 
  } 
  
  if (Cov == "Ar") {
    # Autoregression structure
    Cor_ar <- diag(rep(0, p))
    for (i in 1:p) {
      for (j in 1:p) {
        Cor_ar[i, j] <- 1/2^abs(i-j)
      }
    }
    Sigma_ar <- diag(sd) %*% Cor_ar %*% diag(sd)
    Sigma0 <- Sigma_ar
  } else {
    Sigma0 <- diag(sd^2)
  } 
  
  if (is.null(mu)) {
    mu <- rep(0, p)
  } else if (length(mu) == 1) {
    mu <- rep(mu, p)
  } else if (length(mu) != p) {
    stop("mu must be a scalar or a p-dim vector")
  }
  

  if (Q == "A") {
    scale_Q <- sqrt(1/3) 
    sign_mu_Q <- ifelse(seq(1:p) %% 2 == 0, -1, 1) 
    range_mu_Q <- rep(2.25, p) + mu 
    mu_Q <- range_mu_Q * sign_mu_Q
  } else if (Q == "B") {
    scale_Q <- sqrt(5)
    mu_Q <- rep(5, p) + mu
  } else {
    stop("Q must be of type 'A' or 'B'")
  }

  
  L2 <- c()
  op <- c()
  fro <- c()
  gamma <- list()
  X <- list()
  Thet <- list()
  spline_basis <- list()
  i <- 1
  set.seed(seed)
  while (sum(!is.na(L2)) < fold) {
    cat("\nNow on repeat", sum(!is.na(L2))+1,
        "of", fold,  "-------------------------\n")

    x <- rbind(
      Samples(n1, mu, Sigma0), 
      Samples(n2, mu_Q, diag(rep(scale_Q, p)))
    )
    # for Cauchy contamination try the following 
    # Samples_mCauchy(n2, mu_Q, scale_Q)
    
    res <- fGAN(x=x, decay_lr=decay_lr, lr=lr, iter=iter,
                          pen=pen_lvl*pen0,
                          L1_pen = L1_pen,
                          batch_x=batch_x, batch_z=batch_z,
                          kendall_samp=kendall_size,
                          mu0=mu, Sigma0=Sigma0, 
                          optim=optim,
                          Loss=Loss)
    
    if (any(is.na(res$thet))) {
      cat("\n NA produced, possibly because penalty level is not large enough, skipping to next run...\n")

      thet <- NA
      init_thet <- NA
      gamma <- NA
      L2 <- c(L2)  
      op <- c(op)
      fro <- c(fro)
    } else {
      gamma[[i]] <- res$gamma
      thet <- trans.thet(res$thet, p)
      init_thet <- trans.thet(res$init_thet, p)
      L2 <- c(L2, sqrt( sum( (thet$mu - mu)^2 ) ) )  
      op <- c(op, max(svd_error_bypass(thet$ssig - Sigma0)$d))
      fro <- c(fro, sqrt(sum(svd_error_bypass(thet$ssig - Sigma0)$d^2)))
      Thet[[i]] <- thet
      X[[i]] <- x
      i <- i+1
    }


  }
  
  return(list( L2=L2, op=op, fro=fro,
               Sigma0=Sigma0, Mu=mu,
               X=X, Theta_hat=Thet, gamma=gamma
               )
         )
}
