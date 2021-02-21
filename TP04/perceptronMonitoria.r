rm(list=ls())
library('plot3D')
library('roccv')

# Perceptron:
trainperceptron <- function(xin, yd, eta, tol, maxepocas, par) {
  dimxin <- dim(xin)
  N <- dimxin[1]
  n <- dimxin[2]
  
  if(par == 1) {
    wt <- as.matrix(runif(n+1)-0.5)
    xin <- cbind(1, xin)
  }
  else wt <- as.matrix(runif(n)-0.5)
  
  nepocas <- 0
  eepoca <- tol + 1
  
  evec <- matrix(nrow = 1, ncol = maxepocas)
  while((nepocas < maxepocas) && (eepoca > tol)) {
    ei2 <- 0
    xseq <- sample(N)
    for(i in 1:N) {
      irand <- xseq[i]
      yhati <- 1.0*((xin[irand, ] %*% wt) >= 0)
      ei <- yd[irand]-yhati
      dw <- eta*ei*xin[irand,]
      wt <- wt + dw
      ei2 <- ei2 + ei*ei
    }
    nepocas <- nepocas + 1
    evec[nepocas] <- ei2/N
    
    eepoca <- evec[nepocas]
  }
  retlist <- list(wt, evec[1:nepocas])
  return(retlist)
}

yperceptron  <- function(xin, model){
  ClassPredict <- cbind(1, xin) %*% model[[1]]
  return(ClassPredict)
}

# Definiçao dos dados em duas dimensões:
# Declarando as variancias:
s1 <- s2 <- 0.4
nc<-200

fold <- 10
eta <- 0.01
tol <- 10^-3
maxepocas <- 10^3
par <- 1

# Criando os dados:
xc1<-matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,2)), ncol=nc, nrow=2))
xc2<-matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,4)), ncol=nc, nrow=2))
dataClass1 <- matrix(nrow = dim(xc1)[1], ncol = 1,0)
dataClass2 <- matrix(nrow = dim(xc2)[1], ncol = 1,1)
X1 <- cbind(xc1, dataClass1)
X2 <- cbind(xc2, dataClass2)

xseqc1<-randomly_assign(dim(dataClass1)[1],fold)
xseqc2<-randomly_assign(dim(dataClass2)[1],fold)

ClassPredict <- matrix(0,2*dim(dataClass1)[1]*(1-1/fold),fold)
result <- matrix(0,2*dim(dataClass1)[1]*(1-1/fold),fold)

for (j in 1:fold){
  xc1train <- dataClass1[which(xseqc1 == j), ]
  xc2train <- dataClass2[which(xseqc2 == j), ]
  xc1test <- dataClass1[which(xseqc1 != j),]
  xc2test <- dataClass2[which(xseqc2 != j),]
  
  XTrain <- rbind(xc1train, xc2train)
  XTest <- rbind(xc1test, xc2test)
  
  dataTrain <- trainperceptron(XTrain[,1:(dim(XTrain)[2]-1)], XTrain[,dim(XTrain)[2]], eta, tol, maxepocas, par)
  ClassPredict[,j]<- yperceptron(XTest[,1:(dim(XTrain)[2]-1)], dataTrain)
  result[,j] <- ClassPredict[,j] - XTest[,dim(XTrain)[2]]
}
ValResult <- list(colMeans(result), ClassPredict)

