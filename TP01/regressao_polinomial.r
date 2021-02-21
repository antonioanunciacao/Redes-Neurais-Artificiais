rm(list=ls())
library('corpcor')

# Geração dos dados sinteticos:
N <- 10 # Tamanho da Amostra
n <- 8 # Maior grau da regressão polinomial

x <- runif(n = N, min = -15, max = 10)
yr <- (0.5*x^2+3*x+10) + rnorm(n = length(x), mean = 0, sd = 4)

xgrid <- seq(-15, 15, 0.1)
ygrid <- 0.5*xgrid^2+3*xgrid+10

H <- c(1)
Hgrid <- c(1)

for(i in 1:n) {
  H <- cbind(x^i, H)
  w <- pseudoinverse(H) %*% yr
  yhat <- H %*% w
  Hgrid <- cbind(xgrid^i, Hgrid)
  yhatgrid <- Hgrid %*% w
  
  nome <- paste(as.character(i))
  img <- paste("regPol_p",nome,"_N",N,'.png',sep = "",collapse = "")
  png(file = img, width=720, height=500)
  
  plot(x, yr, col='red', xlim=c(-20,20), ylim=c(-20,100), xlab = 'x', ylab = 'y')
  par(new = T)
  plot(x, yhat, col='black', xlim=c(-20,20), ylim=c(-20,100), xlab = '', ylab = '', axes = F)
  par(new = T)
  plot(xgrid, ygrid, col='green', type = 'l', xlim=c(-20,20), ylim=c(-20,100), xlab = '', ylab = '', axes = F)
  par(new = T)
  plot(xgrid, yhatgrid, col='blue', type = 'l', xlim=c(-20,20), ylim=c(-20,100), xlab = '', ylab = '', axes = F)
  
  dev.off()
}
