trueTheta <- pi/4; trueSigma2 <- 1;
n <- 5000; times <- seq(0, by = 1, length = n)
firstStep <- 0.5
nStart <- floor(n / 3)
base_pow <- 0.6
get_grad_steps <- function(gradient_power, cst = 8){
  firstStep <- 0.5
  nStart <- 300
  steps <- c(rep(firstStep, nStart), sapply(cst * firstStep * (1:(n - nStart) )^(-gradient_power),
                                            function(x) min(x, firstStep)))
  matrix(steps , ncol = 2, nrow = n)
}