
#Function for simulating process as listed in 1
sim_process <- function(n){
  z_vec <- rt(n, 10)
  X <- rep(NA, n)
  X[1] <- z_vec[1]
  for (i in 2:n) {
    X[i] <- 0.9 * X[i - 1] + z_vec[i]
  }
  return(X)
}



