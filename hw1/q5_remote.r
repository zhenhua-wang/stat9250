# Install and load packages
library(doMC)
library(foreach)

## register workers
args = commandArgs(trailingOnly=TRUE)
num_cores <- as.numeric(args[1])
registerDoMC(num_cores)

# Define a function
myfunc <- function(v_s, i_v, iter)
{
  d_vals <- round(i_v %% 256)
  cos_vals <- -cos(2 * d_vals)
  v_mat <- v_s * cos_vals;
  return(v_mat/cos(iter))
}

N1 <- 1e3; N2 <- 2e3; N_tot <- 64
vi_v <- rep(NA, N1)
vd_s <- matrix(NA, N1, N2)

set.seed(123)
for (i in 1:N1){
  vi_v[i] = i + rnorm(1, sd = sqrt(i)*0.01);
  for (j in 1:N2)
    vd_s[i,j] = j + i
}

# Use foreach and print time
ptm <- proc.time()
Res <- foreach(i = 1:num_cores, .combine = c) %dopar% {
  Res_worker <- NULL
  for (iter in seq(i, N_tot, by = num_cores)){
    res_mat <- myfunc(vd_s, vi_v, iter)
    Res_worker = c(Res_worker, mean(res_mat))
  }
  Res_worker
}
print(proc.time() - ptm)

# Print the result
print(Res)
