rm(list=ls())
library('plot3D')
library('roccv')
library('corpcor')

# Funções:

ELMmodel <- function(X, Y, nos) {
  Z <- replicate(nos, runif((dim(X)[2]+1),-0.5,0.5))
  Xaug <- cbind(replicate(dim(X)[1], 1), X)
  H <- tanh(Xaug %*% Z)
  W <- pseudoinverse(H) %*% Y
  parametros <- list(W, Z)
  return(parametros)
}

ELMpredict <- function(Xtest, Ytest, parametros){
  Xaug <- cbind(replicate(dim(Xtest)[1], 1), Xtest)
  H <- tanh(Xaug %*% parametros[[2]])
  Yhat <- sign(H %*% parametros[[1]])
  erro <- sum((Ytest-Yhat)^2)/4
  resultados <- list(Yhat, erro)
  return(resultados)
}

predict <- function(X, pesos, z){
  Xaug <- cbind(replicate(dim(X)[1], 1), X)
  H <- tanh(Xaug %*% z)
  Yhat <- sign(H %*% pesos)
  return(Yhat)
}

s1 <- s2 <- 0.4
nc <- 100
x11<-cbind(matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,2)), ncol=nc, nrow=2)), 1)
x22<-cbind(matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,4)), ncol=nc, nrow=2)), 1)
x12<-cbind(matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,4)), ncol=nc, nrow=2)), -1)
x21<-cbind(matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,2)), ncol=nc, nrow=2)), -1)
X1 <- rbind(x11, x22)
X2 <- rbind(x12, x21)

plot(x11[,1],x11[,2], col = 'green', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
par(new=T)
plot(x12[,1],x12[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
par(new=T)
plot(x21[,1],x21[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
par(new=T)
plot(x22[,1],x22[,2], col = 'green', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')


k <- 10
p <- 10
iter <- 0
W <- matrix(0,p,100)
z <- list()
erro <- matrix(0,100)
for(i in 1:10){
  xseq1<-randomly_assign(dim(X1)[1],k)
  xseq2<-randomly_assign(dim(X2)[1],k)
  
  for(j in 1:k){
    iter <- iter+1
    X1train <- X1[which(xseq1 == j),]
    X2train <- X2[which(xseq2 == j),]
    X1test <- X1[which(xseq1 != j),]
    X2test <- X2[which(xseq2 != j),]
    
    XTrain <- rbind(X1train, X2train)
    XTest <- rbind(X1test, X2test)
    
    model <- ELMmodel(XTrain[,1:2], XTrain[,3], p)
    W[,iter] <- model[[1]]
    z[[iter]] <- model[[2]]
    resultados <- ELMpredict(XTest[,1:2], XTest[,3], model)
    erro[iter] <- resultados[[2]]/dim(XTest)[1]
  }
}
index <- which(erro == min(erro))
index <- index[1]
peso_otimo <- W[, index]
z_otimo <- z[[index]]

seqi<-seq(0, 6, 6/100)
seqj<-seq(0, 6, 6/100)

Class <-matrix(0,nrow=length(seqi),ncol=length(seqj)) 
ci<-0
for(i in seqi){
  ci<-ci+1
  cj<-0
  for(j in seqj){
    cj<-cj+1
    xi<- t(as.matrix(c(i,j)))
    Class[ci,cj]<- predict(xi, peso_otimo, z_otimo)
  }
}

par(new=T)
contour2D(Class, seqi, seqj, colkey = NULL, xlim = c(0,6), ylim = c(0,6), xlab = 'X1', ylab = 'X2')
