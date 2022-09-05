
# utilities ----------------------------

lmat.vec <- function(mat) {
  # Lower triangular matrix to vector
  # Diagonal is removed
  k <- dim(mat)[1]
  
  vec <- NULL
  for (j in 2:k)
    vec <- c(vec, mat[j, 1:(j-1)])
  return(vec)
}

vec.lmat <- function(vec, k) {
  # Vector to lower triangular matrix
  # Diagonal is one
  
  mat <- diag(rep(1,k))
  
  cnt <- 0 
  for (j in 2:k) {
    mat[j, 1:(j-1)] <- vec[cnt+1:(j-1)]
    cnt <- cnt+j-1
  }
  return(mat)
}

tprod <- function(vec1, vec2) {
  # Vectorized outer prod
  # Exactly same as "c(t(outer(vec1, vec2)))"
  
  k1 <- length(vec1)
  k2 <- length(vec2)
  
  vec1[rep(1:k1, each=k2)] * vec2[rep(1:k2, times=k1)]
}


mdiag <- function(u) {
  # u is vector
  # avoid diag(.3) = 0 x 0 matrix
  if (length(u)>1) {
    diag(u)
  } else {
    matrix(u, 1, 1)
  }  
}


Samples <- function(n, mu, sigma) {
  # Generate p-dim multivariate Gaussian sample
  
  p <- length(mu) #p >=2
  z <- matrix(rnorm(p*n), n, p)
  L <- chol(sigma)
  d <- diag(L)
  
  return(t(mu + (t(L/d)) %*%  (t(z)*d)))
}


Samples_mCauchy <- function(n, mu, scale) {
  # Generate p-dim multivariate Cauchy with 0 correlation 
  # and same scale on each coordinate.
  
  p <- length(mu) # p >=2
  
  z <- matrix(rnorm(n*p, mean=0, sd=scale), n, p)
  chi2 <- rchisq(n, df=1)
  x <- z / sqrt(chi2)
    
  return(t(mu+t(x))) # n x p
}

soft <- function(u, a) {
  # Gradient Clipped at a
  # u is vector, a is scalar
  
  if (is.na(max(u))) {
    print(u)
    u <- rep(0, length(u))
  }
  
  if (max(u) > max(-u)) {
    if (max(u) > a) {
      u <- u /max(u) * a
    }
  } else { 
    if (max(-u) > a) {
      u <- u /max(-u) * a
    }
  } 
  return(u)
}


trans.theta <- function(thet, p) {
  # Transfom vectorized parameters to 
  #   corresponding location vector and
  #   covariance matrix.
  
  mu_temp <- thet[1:p]
  d_temp <- exp(thet[p+1:p])
  if (p>1) {
    A2 <- vec.lmat(thet[-(1:(2*p))], p)
  } else { #p=1, not supported 
    A2 <- matrix(1, 1, 1) 
  } 
  ssig_temp <- solve(A2) %*% diag(d_temp^2) %*% t(solve(A2))
  return(list(mu=mu_temp, ssig=ssig_temp))
}


trans.thet <- function(thet, p) {
  # Legacy name
  mu_temp <- thet[1:p]
  d_temp <- exp(thet[p+1:p])
  if (p>1) {
    A2 <- vec.lmat(thet[-(1:(2*p))], p)
  } else { #p=1, not supported 
    A2 <- matrix(1, 1, 1) 
  } 
  ssig_temp <- solve(A2) %*% diag(d_temp^2) %*% t(solve(A2))
  return(list(mu=mu_temp, ssig=ssig_temp))
}



kendall <- function(x) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  tao <- matrix(0, p, p)
  for (j in 1:p) {
    for (k in 1:p) {
      temp <- 0
      for (i in (1:(n-1))) {
        for (ii in ((i+1):n)) {
          temp <- temp + sign( (x[i,j] - x[ii, j])*(x[i,k]-x[ii,k]))
        }
      }
      tao[j,k] <- 2/(n*(n-1)) * temp
    }
  }
  return(sin(pi/2 *tao))
}

SS <- function(x) {
  # For Kendall's tau estimatior
  n <- dim(x)[1]
  p <- dim(x)[2]
  s <- rep(0, p)
  for (j in 1:p) {
    s[j] <- median(x[, j]^2)
  }
  return(diag(s))
}

cor2cov <- function(R, S) {
  # Correlation matrix and standard deviation vector to covariance matrix  
  sweep(sweep(R, 1, S, "*"), 2, S, "*")
}


expit <- function(x) {
  1/(1+exp(-x))
}


# Svd to matrix norms
svd_to_norm <- function(Sigma, true_Sigma) {
  svd_res <- svd_error_bypass(true_Sigma - Sigma)$d
  return(c(svd_res[1], sqrt(sum(svd_res^2)), max(abs(true_Sigma - Sigma))))
}

norm_to_mean_sd <- function(vec) {
  return(c(mean(vec, na.rm=T), sd(vec, na.rm=T)))
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
# https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r

norm_to_mean_sd_chr <- function(vec) {
  return(
    paste(specify_decimal(mean(vec, na.rm=T), 4),
          " (", specify_decimal(sd(vec, na.rm=T), 4), ")",
          sep="")
  )
}
