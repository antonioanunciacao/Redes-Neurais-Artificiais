rm(list=ls())
dev.off()
library('plot3D')
library('roccv')

# Multi Layer Perceptron:
MLPerceptron <- function(xin, yd, eta, tol, maxepocas, neuronios, xtest, ytest) {
  dimxin <- dim(xin)
  N <- dimxin[1]
  n <- dimxin[2]

  wo <- matrix( runif( (n+1)*neuronios, -0.5, 0.5), nrow =neuronios, ncol=n+1 )
  wt <- matrix(runif(neuronios+1)-0.5, nrow = 1)
  xin <- cbind(1, xin)
  xtest <- cbind(1, xtest)
  
  nepocas <- 0
  eepoca <- tol + 1
  
  evec <- matrix(0, nrow = 1, ncol = maxepocas)
  eTestvec <- matrix(0, nrow = 1, ncol = maxepocas)
  while((nepocas < maxepocas) && (eepoca > tol)) {
    erro <- errotest <- 0
    xseq <- sample(N)
    
    for(i in 1:N) {
      irand <- xseq[i]
      
      z1 <- wo %*% xin[irand, ]
      a1 <- rbind(1, tanh(z1))
      
      z2 <- wt %*% a1
      #yhati <- tanh(z2)
      yhati <- z2
      
      e <- yd[irand]-yhati
      deltaE2 <- -1*e
      dwt <- eta*deltaE2 %*% t(a1)
      
      dwo <- matrix(0,dim(wo)[1], dim(wo)[2])
      for(i in 1:dim(wo)[1]) {
        dwo[i,] <- ( eta*deltaE2*wt[,i+1]*( 1/cosh(z1[i,])^2 ) ) %*% t(xin[irand, ])
      }

      wt <- wt - dwt
      wo <- wo - dwo
      erro <- erro + e*e
    }
    
    xtestseq <- sample(dim(xtest)[1])
    for(i in 1:dim(xtest)[1]) {
      irandtest <- xtestseq[i]
      Z1test <- wo %*% xtest[irandtest, ]
      A1test <- tanh(Z1test)
      Yhattest <- wt %*% rbind(1,A1test)
      Predict <- Yhattest
      etest <- ytest[irandtest] - Predict
      errotest <- errotest + etest*etest
    }
    
    nepocas <- nepocas + 1
    
    evec[nepocas] <- erro/N
    eTestvec[nepocas] <- errotest/dim(xtest)[1]
    
    eepoca <- evec[nepocas]

    cat("Erro[", nepocas, "]: ", evec[nepocas], "\n")
  }
  retlist <- list(wo, wt, evec[1:nepocas], eTestvec[1:nepocas])
  return(retlist)
}

MLPredict <- function(xin, model) {
  W1 <- model[[1]]
  W2 <- model[[2]]
  X <- cbind(1,xin)
  Predict <- matrix(0, dim(xin)[1])
  
  for(i in 1:dim(X)[1]) {
    Z <- W1 %*% X[i,]
    A <- tanh(Z)
    Yhat <- W2 %*% rbind(1,A)
    Predict[i] <- Yhat
  }
  return(Predict)
}

s1 <- s2 <- 0.4
nc <- 100
x11<-cbind(matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,2)), ncol=nc, nrow=2)), 1)
x22<-cbind(matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,4)), ncol=nc, nrow=2)), 1)
x12<-cbind(matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,4)), ncol=nc, nrow=2)), -1)
x21<-cbind(matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,2)), ncol=nc, nrow=2)), -1)
X1 <- rbind(x11, x22)
X2 <- rbind(x12, x21)

png("Distribuicao.png", 550, 380)
  plot(x11[,1],x11[,2], col = 'green', xlim = c(0,6), ylim = c(0,6), xlab = 'X1', ylab = 'X2')
  par(new=T)
  plot(x12[,1],x12[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
  par(new=T)
  plot(x21[,1],x21[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
  par(new=T)
  plot(x22[,1],x22[,2], col = 'green', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
dev.off()

X1_index <- sample (c(1:200), size=0.7*dim(X1)[1], replace =F)
X2_index <- sample (c(1:200), size=0.7*dim(X2)[1], replace =F)

XTrain <- rbind(X1[ X1_index,], X2[ X2_index,])
XTest  <- rbind(X1[-X1_index,], X2[-X2_index,])

# Hiperparametros
eta <- 0.01
tol <- 10^-3
maxepocas <- 10^3
neuronios <- 3

# Treinamento
model <- MLPerceptron(XTrain[,1:2], XTrain[,3], eta, tol, maxepocas, neuronios, XTest[,1:2], XTest[,3])

Ypredict <- MLPredict(XTest[,1:2], model)
result <- cbind(XTest[,3], sign(Ypredict))

YpredicTrain <- MLPredict(XTrain[,1:2], model)
resultTrain <- cbind(XTrain[,3], sign(YpredicTrain))

acuraciaTest <- 1 - abs(mean(XTest[,3] - sign(Ypredict)))
sdTest <- sd(XTest[,3] - sign(Ypredict))
acuraciaTrain <- 1 - abs(mean(XTrain[,3] - sign(YpredicTrain)))
sdTrain <- sd(XTrain[,3] - sign(YpredicTrain))
tabResultados <- rbind(cbind(acuraciaTrain, sdTrain), cbind(acuraciaTest, sdTest))
seqi<-seq(0, 6, 6/100)
seqj<-seq(0, 6, 6/100)

Class <-matrix(0,nrow=length(seqi),ncol=length(seqj)) 
ci<-0
for(i in seqi){
  ci<-ci+1
  cj<-0
  for(j in seqj){
    cj<-cj+1
    xi<- as.matrix(t(c(i,j)))
    Class[ci,cj]<- MLPredict(xi, model)
  }
}
png("superficieXOR.png", 550, 380)
  plot(x11[,1],x11[,2], col = 'green', xlim = c(0,6), ylim = c(0,6), xlab = 'X1', ylab = 'X2', main = "Superficie de Separação")
  par(new=T)
  plot(x12[,1],x12[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = NA, ylab = NA)
  par(new=T)
  plot(x21[,1],x21[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = NA, ylab = NA)
  par(new=T)
  plot(x22[,1],x22[,2], col = 'green', xlim = c(0,6), ylim = c(0,6), xlab = NA, ylab = NA)
  par(new=T)
  contour2D(Class, seqi, seqj, colkey = NULL, col = 'red', xlim = c(0,6), ylim = c(0,6), levels =0,  xlab = NA, ylab = NA, add = TRUE)
dev.off()

png("CurvaAprendizado.png", 550, 380)
  plot(seq(1:maxepocas), model[[3]], type ='l', col = 'green', xlim = c(0,maxepocas), ylim = c(0,1), xlab = 'Epoca', ylab = 'Erro', main = "Curva de Aprendizado")
  par(new=T)
  plot(seq(1:maxepocas), model[[4]], type = 'l', col = 'blue', xlim = c(0,maxepocas), ylim = c(0,1), xlab = NA, ylab = NA)
  legend(x=700, y=1, legend = c('Erro Treino','Erro Teste'), col = c('green','blue'), pch=c('-','-'))
dev.off()