rm(list=ls())

perceptron_1passo <- function(dados, pesos, eta) {
  x <- dados[1:2]
  y <- dados[3]
  x <- c(x, 1)
  yhat <- 1.0*(x %*% pesos)
  
  if(yhat >= 0) {
    yhat <- 1.0*(+1)
  } else {
    yhat <- 1.0*(-1)
  }
  
  erro <- y - yhat
  dw <- eta*erro*x[1:3]
  pesos <- pesos + dw
  
  return(pesos)
}

eta <- 0.1
dados <- as.matrix(c(0.9, 1.9, 1), nrow = 3, ncol = 1)
pesos <- as.matrix(c(-0.9, -1.7, -0.5))

pesos_1p <- perceptron_1passo(dados, pesos, eta)