rm(list=ls())
dev.off()
library('plot3D')
library('roccv')
library("mlbench")

# Perceptron:
trainPerceptron <- function(xin, yd, eta, tol, maxepocas, par) {
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
      yhati <- binariza(1.0*((xin[irand, ] %*% wt)))
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

yPerceptron  <- function(xin, model){
  ClassPredict <- cbind(1, xin) %*% model[[1]]
  return(binariza(ClassPredict))
}

Validation <- function(dataClass1, dataClass2, fold, eta, tol, maxepocas, par){
  
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
    
    dataTrain <- trainPerceptron(XTrain[,1:(dim(XTrain)[2]-1)], XTrain[,dim(XTrain)[2]], eta, tol, maxepocas, par)
    ClassPredict[,j]<- yPerceptron(XTest[,1:(dim(XTrain)[2]-1)], dataTrain)
    result[,j] <- ClassPredict[,j] - XTest[,dim(XTrain)[2]]
  }
  ValResult <- list(colMeans(result), ClassPredict)
  
}

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

binariza <- function(y){
  for(i in 1:length(y)){
    if(y[i] >= 0) { y[i] <-1 }
    else { y[i] <- 0 }
  }
  return (y)
}

# Definiçao dos dados:
# Declarando as variancias:
#pega os dados da package mlbench
data(BreastCancer)
data <- data.matrix(BreastCancer[,2:10])
data[is.na(data)] <- 0
class <- as.numeric(BreastCancer$Class)
class[which(class == 2)] <- 0
X <- cbind(data, class)
X1 <- X[which(class == 1), ]
X2 <- X[which(class == 0), ]

k <- 10
eta <- 0.01
tol <- 10^-3
maxepocas <- 10^3
par <- 1

xseqc1<-randomly_assign(dim(X1)[1],k)
xseqc2<-randomly_assign(dim(X2)[1],k)

modelPerceptron <- list()
ClassPredict <- list()

result <- list()
for (j in 1:k){
  xc1train <- X1[which(xseqc1 == j), ]
  xc2train <- X2[which(xseqc2 == j), ]
  xc1test <- X1[which(xseqc1 != j),]
  xc2test <- X2[which(xseqc2 != j),]
  
  XTrain <- rbind(xc1train, xc2train)
  XTest <- rbind(xc1test, xc2test)
  
  modelPerceptron[[j]] <- trainPerceptron(XTrain[,1:(dim(XTrain)[2]-1)], XTrain[,dim(XTrain)[2]], eta, tol, maxepocas, par)
  ClassPredict[[j]] <- yPerceptron(XTest[,1:(dim(XTrain)[2]-1)], modelPerceptron[[j]])
  result[[j]] <- ClassPredict[[j]] - XTest[,dim(XTrain)[2]]
}
resultF <- array(0,k)
desvioF <- array(0,k)
for(i in 1:k){
  resultF[i] <- abs(mean(result[[i]]))
  desvioF[i] <- sd(result[[i]])
}

Acuracia <- array(0,k)
for(i in 1:k){
  Acuracia[i] <- 1 - resultF[i]
}
ResultadoFinal <- cbind(Acuracia, desvioF)

Matriz_pesos <- matrix(0,k, k+1)
for(i in 1:k){
  model <- modelPerceptron[[i]]
  model <- model[[1]]
  Matriz_pesos[i,] <- cbind(i, t(model))
}

write.table(ResultadoFinal, file = "ResultadoFinal.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

write.table(Matriz_pesos, file = "Matriz_pesos.txt", sep = "\t",
            row.names = TRUE, col.names = NA)