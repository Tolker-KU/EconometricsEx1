
sim_process <- function(n){
  z_vec = rt(n, 10)
  return(z_vec + c(0, head(z_vec, -1) * 0.9))
}

test <- 1 + 1
