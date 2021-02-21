rm(list=ls())
library('plot3D')
library('roccv')

# Adaline:
trainAdaline <- function(xin, yd, eta, tol, maxepocas, par) {
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
      yhati <- 1.0*((xin[irand, ] %*% wt))
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
      yhati <- sign(1.0*((xin[irand, ] %*% wt)))
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
    
    dataTrain <- trainAdaline(XTrain[,1:(dim(XTrain)[2]-1)], XTrain[,dim(XTrain)[2]], eta, tol, maxepocas, par)
    ClassPredict[,j]<- yperceptron(XTest[,1:(dim(XTrain)[2]-1)], dataTrain)
    result[,j] <- ClassPredict[,j] - XTest[,dim(XTrain)[2]]
  }
  ValResult <- list(colMeans(result), ClassPredict)
  
}

# Definiçao dos dados em duas dimensões:
# Declarando as variancias:
s1 <- s2 <- 0.4
nc<-100

# Criando os dados:
xc1<-matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,2)), ncol=nc, nrow=2))
xc2<-matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,4)), ncol=nc, nrow=2))
xc1Class <- matrix(nrow = dim(xc1)[1], ncol = 1,-1)
xc2Class <- matrix(nrow = dim(xc2)[1], ncol = 1,1)
X1 <- cbind(xc1,xc1Class)
X2 <- cbind(xc2,xc2Class)
X <- rbind(X1, X2)

plot(xc1, col='red', xlim=c(0,6), ylim=c(0,6), xlab = 'XC1', ylab = 'XC2')
par(new = T)
plot(xc2, col='blue', xlim=c(0,6), ylim=c(0,6), xlab = '', ylab = '', axes = F)


k <- 10
eta <- 0.01
tol <- 10^-3
maxepocas <- 10^3
par <- 1

model_Adaline <- trainAdaline(X[,1:(ncol(X)-1)], X[,ncol(X)], eta, tol, maxepocas, par)
model_Perceptron <- trainPerceptron(X[,1:(ncol(X)-1)], X[,ncol(X)], eta, tol, maxepocas, par)
pesos_Adaline <- model_Adaline[[1]]

resultV <- Validation(X1, X2, k, eta, tol, maxepocas, par)

coeficientes <- -pesos_Adaline[1:2]/pesos_Adaline[3]
seqx = seq(0,6,6/100)
xi = as.matrix(seqx)
yi = array(0, length(xi))
for (i in 1:length(xi)){
  yi[i] = coeficientes[1]+coeficientes[2]*xi[i]
}
par(new = T)
plot(xi, yi, type='l', col='green', xlim=c(0,6), ylim=c(0,6), xlab = '', ylab = '', axes = F)

seqi<-seq(0, 6, 6/100)
seqj<-seq(0, 6, 6/100)
xy <- matrix(0,nrow=length(seqi),ncol=length(seqj))
ci<-0
for(i in seqi){
  ci<-ci+1
  cj<-0
  for(j in seqj){
    cj<-cj+1
    xy[ci,cj] <- sign(t(as.matrix(c(1,i,j))) %*% model_Adaline[[1]])
  }
}

# Plotando as densidades
persp3D(seqi, seqj, xy, counter=T, theta = 55, phi = 30, r = 40, d = 0.1, expand = 0.5,
        ltheta = 90, lphi = 180, shade = 0.4, ticktype = "detailed", nticks=5)

write.table(pesos_Adaline/pesos_Adaline[3], file = "pesos_adaline.txt", sep = "\t",
            row.names = TRUE, col.names = NA)