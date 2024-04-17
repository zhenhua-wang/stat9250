load("H:/My Drive/Simulations class/Homeworks/HW3_shalini/GPD_samples.RData")
n <- length(samples)
min(samples)

d_log_l <- function(param){
  d_out <- array(NA, dim = 2)
  out <- ( 1 + param[2] *(samples - mu)/param[1] )
  d_out[1] <- -n/param[1] + ( 1/param[2] + 1) * sum( out^(-1) * param[2] * (samples - mu)/ param[1]^2 )
  d_out[2] <- (1/param[2]^2) * sum(log(out))  - ( 1/param[2] + 1) * sum( out^(-1)*(samples - mu)/ param[1] )
  return(d_out)
}

H_log_l <- function(param){
  H_mat <- matrix(NA, nrow = 2, ncol = 2)
  out <- ( 1 + param[2] *(samples - mu)/param[1] )
  H_mat[1,1] <- n/param[1]^2 + (1/param[2]+1) * param[2]^2/param[1]^4 * sum( out^(-2)*(samples-mu) ) - 
    2 * (1/param[2]+1) * param[2]/param[1]^3 * sum(out^(-1)*(samples-mu))
  H_mat[2,1] <- H_mat[1,2] <- sum( out^(-1) * (samples-mu)/param[1]^2) - (1/param[2]+1) * param[2]/param[1]^3 * sum(out^(-2)*(samples-mu)^2)
  H_mat[2,2] <- - 2/param[2]^3 * sum(log(out)) + 2/param[2]^2 * sum( out^(-1)*(samples - mu)/param[1] ) +
    (1/param[2] + 1) * sum(out^(-2) * (samples - mu)^2/param[1]^2)
  return(H_mat)
}

param_diff <- 1
i <- 1

mu <- 10.5
param <- c(1, 0.1) #sigma and xi respectively # we keep mu a constant 


while(param_diff>1e-10){
  d_l <- d_log_l(param)
  H <- H_log_l(param)
  param_new <- param - solve(H) %*% d_l
  #lik_new <- log_l(param_new)
  param_diff <- max(abs(param - param_new))
  param <- param_new
  if(i%%100 == 0) print(i)
  i <- i+1
  print(paste0("i = ",i))
  print(param_new)
}

est_params_BFGS <- param_new

I_mat <- -H_log_l(param_new)
std_errors <- sqrt(diag(solve(I_mat)))

alpha <- 0.05  # Confidence level
z_score <- qnorm(1 - alpha/2)
CI_NR <- cbind(param_new - z_score * std_errors,
                 param_new + z_score * std_errors)


# optim function BFGS -----------------------------------------------------

n_log_l <- function(param){
  n <- length(samples)
  mu <- 10.5
  if(all(param > 0.00001)){
    out <- ( 1 + param[2] *(samples - mu)/param[1] )
    out = -n*log( param[1] ) - (1/param[2] + 1) * sum(log(out))
    return(-out)
  }else
    return(+1e+6)
}

init_param <- c(1, 0.1) #sigma and xi respectively # we keep mu a constant 

fit <- optim(init_param, n_log_l, method = "BFGS", hessian = TRUE)
est_params_BFGS <- fit$par

I_mat <- fit$hessian
std_errors <- sqrt(diag(solve(I_mat)))

# Calculate confidence intervals
alpha <- 0.05  # Confidence level
z_score <- qnorm(1 - alpha/2)
CI_BFGS <- cbind(est_params_BFGS - z_score * std_errors,
                 est_params_BFGS + z_score * std_errors)

fit$counts
fit$convergence


# CG optim ----------------------------------------------------------------

fit <- optim(init_param, n_log_l, method = "CG", hessian = TRUE)
est_params_CG <- fit$par

I_mat <- fit$hessian
std_errors <- sqrt(diag(solve(I_mat)))

# Calculate confidence intervals
alpha <- 0.05  # Confidence level
z_score <- qnorm(1 - alpha/2)
CI_CG <- cbind(est_params_CG - z_score * std_errors,
               est_params_CG + z_score * std_errors)


fit$counts
