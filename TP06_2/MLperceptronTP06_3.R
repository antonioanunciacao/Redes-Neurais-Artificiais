rm(list=ls())
dev.off()
library('plot3D')
library('roccv')
library('mlbench')

# Multi Layer Perceptron:
MLPerceptron <- function(xin, yd, eta, tol, maxepocas, neuronios, xtest, ytest, fold) {
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
    
    if(nepocas %% 100 == 0) cat("Erro[", fold, ",", nepocas,"]:", evec[nepocas], "\n")
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
    Yhat <- tanh(W2 %*% rbind(1,A))
    Predict[i] <- Yhat
  }
  return(Predict)
}

##################################################################
# Load Dataset
data(BreastCancer)
data <- data.matrix(BreastCancer)
data[is.na(data)] <- 0
dataV <- data[,2:10]
for(i in 1:dim(data)[1]){
  if(data[i,11] == 1) data[i,11] = 0
  else data[i,11] = 1
}

XBenign <- data[which( data[,11] == 0 ),]
XMalignant <- data[which( data[,11] == 1 ),]

kfolds <- 10
model_kfolds <- list()
resultTrain <- matrix(0, kfolds, 2)
resultTest <- matrix(0, kfolds, 2)
for(i in 1:kfolds){
  bach_size <- round(0.7*min(dim(XBenign)[1], dim(XMalignant)[1]))
  XBenign_index <- sample (c(1:dim(XBenign)[1]), size=bach_size, replace =F)
  XMalignant_index <- sample (c(1:dim(XMalignant)[1]), size=bach_size, replace =F)
  
  XBenignTest <- XBenign[-XBenign_index,]
  XBenign_index_test <- sample (c(1:dim(XBenignTest)[1]), size=dim(XMalignant[-XMalignant_index,]), replace =F)
  XTrain <- rbind(XBenign[XBenign_index,],XMalignant[XMalignant_index,])
  XTest <- rbind(XBenign[XBenign_index_test,],XMalignant[-XMalignant_index,])
  
  YTrain <- as.matrix(XTrain[,11])
  XTrain <- as.matrix(XTrain[,2:10])
  
  YTest <- XTest[,11]
  XTest <- XTest[,2:10]
  
  # Hiperparametros
  eta <- 0.01
  tol <- 10^-3
  maxepocas <- 10^3
  neuronios <- 3
  
  # Entradas da Rede
  model <- MLPerceptron(XTrain, YTrain, eta, tol, maxepocas, neuronios, XTest, YTest, i)
  
  YTrainpredict <- MLPredict(XTrain, model)
  YTestpredict <- MLPredict(XTest, model)
  
  resultTrain[i,1] <- mean(YTrain-YTrainpredict)
  resultTrain[i,2] <- sd(YTrain-YTrainpredict)
  resultTest[i,1] <- mean(YTest-YTestpredict)
  resultTest[i,2] <- sd(YTest-YTestpredict)
  
  model_kfolds[i] <- list(model) # [Wo, W1], [erroTreino], [erroTeste]
}
acuraciaTest <- t(as.matrix(colMeans(resultTest)))
acuraciaTest[1,1] <- 1-acuraciaTest[1,1]
acuraciaTrain <- t(as.matrix(colMeans(resultTrain)))
acuraciaTrain[1,1] <- 1-acuraciaTrain[1,1]
ResultadoFinal <- as.matrix(rbind(acuraciaTrain, acuraciaTest))

limit1 <- 0.6
png("erroBreastCEq.png", 550, 380)
  plot(seq(1:kfolds), resultTrain[,1], type='l', col = 'green', xlim = c(1,kfolds), ylim = c(0,limit1), xlab = 'fold', ylab = 'Erro')
  par(new=T)
  plot(seq(1:kfolds), resultTrain[,2], type='b', col = 'green', xlim = c(1,kfolds), ylim = c(0,limit1), xlab = '', ylab = '')
  par(new=T)
  plot(seq(1:kfolds), resultTest[,1], type='l', col = 'blue', xlim = c(1,kfolds), ylim = c(0,limit1), xlab = '', ylab = '')
  par(new=T)
  plot(seq(1:kfolds), resultTest[,2], type='b', col = 'blue', xlim = c(1,kfolds), ylim = c(0,limit1), xlab = '', ylab = '')
  legend(x=6.5, y=limit1, legend = c('Erro Treino','Desvio Padrao','Erro Teste', 'Desvio Padrao'), col = c('green', 'green','blue', 'blue'), pch=c('-','-.', '-', '-.'))
dev.off()

png("CurvaAprendizadoBreastCEq.png", 550, 380)
  plot(seq(1:maxepocas), model[[3]], type='l', col = 'green', xlim = c(0,maxepocas), ylim = c(0,0.1), xlab = 'Epoca', ylab = 'Erro')
  par(new=T)
  plot(seq(1:maxepocas), model[[4]], type='l', col = 'blue', xlim = c(0,maxepocas), ylim = c(0,0.1), xlab = '', ylab = '')
  legend(x=700, y=0.1, legend = c('Erro Treino','Erro Teste'), col = c('green','blue'), pch=c('-','-'))
dev.off()

png("CurvaAcuraciaBreastCEq.png", 550, 380)
  plot(seq(1:maxepocas), 1-model[[3]], type='l', col = 'green', xlim = c(0,maxepocas), ylim = c(0.90,1), xlab = 'Epoca', ylab = 'Acuracia')
  par(new=T)
  plot(seq(1:maxepocas), 1-model[[4]], type='l', col = 'blue', xlim = c(0,maxepocas), ylim = c(0.90,1), xlab = '', ylab = '')
  legend(x=600, y=0.94, legend = c('Acuracia Treino','Acuracia Teste'), col = c('green','blue'), pch=c('-','-'))
dev.off()