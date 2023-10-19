Nowcasting <- function(data){
  
  # Hyperparameters for Gamma prior of delta
  a.delta <- b.delta <- 10^(-5)
  
  # Prior for overdispersion parameter phi
  a.phi <- b.phi <- 10^(-5)
  
  # nu prior parameter for the penalty
  nu <- 3
  
  D <- max(data$d) # Maximum delay
  TT <- max(data$t) # Maximum time
  
  # Rows for not yet reported
  nyr <- which(data$Reported == "Not yet reported")
  
  # Time and delay
  t <- unique(data$t)
  d <- unique(data$d)
  
  # Model matrices
  
  Kt = 40 # Number of B-splines for time dimension
  Kd = 10 # Number of B-splines for delay dimension
  
  # B-spline basis matrix
  Bt <- cubicbs(t,lower = min(t),upper = max(t),K = Kt)$Bmatrix
  Bd <- cubicbs(d,lower = min(d),upper = max(d),K = Kd)$Bmatrix
  
  # Two dimensional B-spline matrix
  B <- kronecker(Bd, Bt)
  
  y <- data$Cases[-nyr] # Reported cases
  
  # Difference order of the penalty
  penorder <- 2
  
  # Penalty for column (delay dimension)
  Dd <- diag(Kt)
  for (k in 1:penorder) Dd <- diff(Dd)
  Pd <- t(Dd) %*% Dd
  Pd <- Pd + diag(1e-12, Kt)
  
  # Penalty for row (time dimension)
  Dt <- diag(Kd)
  for (k in 1:penorder) Dt <- diff(Dt)
  Pt <- t(Dt) %*% Dt
  Pt <- Pt + diag(1e-12, Kd)
  
  # Full design matrix
  X1 <- B
  
  X_nyr <- X1[nyr,] # Design matrix for not yet reported cases
  X <- X1[-nyr,] # Design matrix for reported cases
  
  # Precision matrix for B-spline parameters
  Pv <- function(v) exp(v[1])*(kronecker(diag(1,Kd),Pd)) + 
    exp(v[2])*(kronecker(Pt,diag(1,Kt)))
  
  # Precision matrix for parameter xi
  Qv <- Pv
  
  # Negative binomial GLM with log-link
  mu.nb <- function(xi) exp(as.numeric(X %*% xi))
  var.nb <- function(xi, v) mu.nb(xi) + (1/exp(v[3]))*(mu.nb(xi)^2)
  W.nb <- function(xi, v) diag(((exp(as.numeric(X %*% xi)))^2)*(1/var.nb(xi, v)))
  D.nb <- function(xi) diag(1/mu.nb(xi))
  M.nb <- function(xi) diag(y - mu.nb(xi))
  V.nb <- function(xi, v) diag(mu.nb(xi) * (1/var.nb(xi, v) - 
                                              (mu.nb(xi)/(var.nb(xi, v)^2)) * 
                                              (1 + 2*mu.nb(xi)*(1/exp(v[3])))))
  gamma.nb <- function(xi, v) exp(v[3]) * log(mu.nb(xi)/(mu.nb(xi) + exp(v[3])))
  bgamma.nb <- function(xi, v) - (exp(v[3])^2) * log(exp(v[3])/(exp(v[3]) + mu.nb(xi)))
  
  # Log conditional posterior of xi given v
  log_pxi <- function(xi, v) {
    value <- (1/exp(v[3])) * sum((y * gamma.nb(xi, v)) - bgamma.nb(xi, v)) - .5 * t(xi) %*% Qv(v[1:2]) %*% xi
    return(value)
  }
  
  # Gradient of parameter xi
  Grad.logpxi <- function(xi,v){
    value <- t(X)%*%W.nb(xi,v)%*%D.nb(xi)%*%(y - mu.nb(xi)) - Qv(v[1:2])%*%xi
    as.numeric(value)
  }
  
  # Hessian of parameter xi
  Hess.logpxi <- function(xi,v){
    value <- t(X)%*%M.nb(xi)%*%V.nb(xi, v)%*%X - 
      t(X)%*%W.nb(xi, v)%*%X - Qv(v[1:2])
    value
  }
  
  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm
  NR_xi <- function(xi0, v){
    
    epsilon <- 1e-05 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxi <- as.numeric((-1) * solve(Hess.logpxi(xi0, v),
                                     Grad.logpxi(xi0, v)))
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, v)
      while (log_pxi(xi.new, v) <= logpxi.current) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    xistar <- xi0 
    return(xistar)
  }
  
  # Initial values for log-penalty and log-overdispersion parameter
  v_init = c(1,1,1)
  
  # Initial estimate for xi
  xi_init <- NR_xi(xi0 = rep(0,dim(X)[2]), v = v_init)
  
  # Log-conditional (joint) posterior of overdispersion parameter and penalty on the time dimension
  log_pvcond <- function(vpar){
    
    v <- c(vpar[1], -3 ,vpar[2])
    
    e1 <- eigen(Pv(v[1:2]),only.values = T)$values
    e2 <- eigen(-t(X)%*%(M.nb(xi = xi_init)%*%V.nb(xi = xi_init, v = v)-W.nb(xi = xi_init, v = v))%*%X +
                  Qv(v[1:2]),only.values = T)$values
    
    value <- sum((1/exp(v[3]))*((y * gamma.nb(xi = xi_init, v)) - bgamma.nb(xi = xi_init, v)) +
                   lgamma(y + exp(v[3])) - lgamma(exp(v[3]))) +
      0.5*sum(sapply(e1[e1>0], log)) - 0.5 * sum((xi_init * Qv(v[1:2])) %*% xi_init) +
      a.phi*v[3] - b.phi*exp(v[3]) - 0.5*sum(sapply(e2[e2>0],log)) +
      0.5*nu*(v[1]+v[2]) - (0.5*nu + a.delta)*(log(b.delta + 0.5*nu*exp(v[1]))+log(b.delta + 0.5*nu*exp(v[2])))
    
    return(value)
  }
  
  vstar <- optim(par = c(1,1), fn = log_pvcond,method = "Nelder-Mead", control = list(fnscale = -1))$par
  
  # Conditional posterior mode of v
  v_mode <- c(vstar[1],-3, vstar[2])
  
  # Mode a posteriori estimate for xi
  xi_mode <- NR_xi(xi0 = xi_init, v = v_mode)
  
  ## Nowcast for not yet reported
  mu_nyr <- exp(X_nyr%*%xi_mode)
  nowcast <- data
  nowcast[nyr,"Cases"] <- mu_nyr
  
  # Nowcasted cases (reported + nowcast)
  cases.now <- aggregate(Cases ~ t, data = nowcast, FUN = function(x) ceiling(sum(x)))
  colnames(cases.now) <- c("t", "y")
  
  output <- list(data = data,
                 cases.now = cases.now)
  output
}
