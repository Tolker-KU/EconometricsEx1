
sim_process <- function(n){
  z_vec = rt(n, 10)
  return(z_vec + c(0, head(z_vec, -1) * 0.9))
}

<<<<<<< HEAD
test <- 1 + 1
=======
>>>>>>> 485359772b4574af629fa62d663864e3f1afed2b