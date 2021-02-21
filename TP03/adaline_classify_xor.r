rm(list=ls())
library('plot3D')

trainadaline <- function(xin, yd, eta, tol, maxepocas, par) {
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

# Criando os dados:
s1 <- s2 <- 0.4
nc <- 500
x11<-matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(-1,-1)), ncol=nc, nrow=2))
x12<-matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(1,1)), ncol=nc, nrow=2))
x21<-matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(-1,1)), ncol=nc, nrow=2))
x22<-matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(1,-1)), ncol=nc, nrow=2))
X <- rbind(x11, x12, x21, x22)

# Plotando os graficos
plot(x11[,1],x11[,2], col = 'green', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = 'X1', ylab = 'X2')
par(new=T)
plot(x12[,1],x12[,2], col = 'green', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')
par(new=T)
plot(x21[,1],x21[,2], col = 'blue', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')
par(new=T)
plot(x22[,1],x22[,2], col = 'blue', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')

# Separando os dados em um grid vertical e outra horizontal:
XVerticalC1 <- rbind(x11, x21) # Classe -1
XVerticalC1Class <- matrix(-1, nrow = dim(XVerticalC1)[1], ncol = 1)
XVerticalC2 <- rbind(x12, x22) # Classe +1
XVerticalC2Class <- matrix(1, nrow = dim(XVerticalC2)[1], ncol = 1)
XVertical <- rbind(XVerticalC1, XVerticalC2)
XVerticalClass <- rbind(XVerticalC1Class, XVerticalC2Class)

XHorizontalC1 <- rbind(x11, x22) # Classe -1
XHorizontalC1Class <- matrix(-1, nrow = dim(XHorizontalC1)[1], ncol = 1)
XHorizontalC2 <- rbind(x21, x12) # Classe +1
XHorizontalC2Class <- matrix(1, nrow = dim(XHorizontalC2)[1], ncol = 1)
XHorizontal <- rbind(XHorizontalC1, XHorizontalC2)
XHorizontalClass <- rbind(XHorizontalC1Class, XHorizontalC2Class)

# Separando dados em um grupo de teste e outro de treino:
MtrainX<- 0.9*dim(X)[1]
seqX <- sample(dim(X)[1])

XtrainVertical <- XVertical[seqX[(1:MtrainX)], ]
XTestVertical <- XVertical[seqX[((MtrainX+1):dim(X)[1])], ]
XClassTrainVertical <- XVerticalClass[seqX[(1:MtrainX)], ]
XClassTestVertical <- XVerticalClass[seqX[((MtrainX+1):dim(X)[1])], ]

XtrainHorizontal <- XHorizontal[seqX[(1:MtrainX)], ]
XTestHorizontal <- XHorizontal[seqX[((MtrainX+1):dim(X)[1])], ]
XClassTrainHorizontal <- XHorizontalClass[seqX[(1:MtrainX)], ]
XClassTestHorizontal <- XHorizontalClass[seqX[((MtrainX+1):dim(X)[1])], ]

# Parametros de controle do algoritmo:
eta <- 0.01
tol <- 0.001
maxepocas <- 100
par <- 1

vertical_classify <- trainadaline(XtrainVertical, XClassTrainVertical, eta, tol, maxepocas, par)
Horizontal_classify <- trainadaline(XtrainHorizontal, XClassTrainHorizontal, eta, tol, maxepocas, par)
pesos_verticais <- as.matrix(vertical_classify[[1]])
pesos_horizontais <- as.matrix(Horizontal_classify[[1]])

# Classificação final (Classificação Vertical X horizontal):

seqi<-seq(-3, 3, 6/100)
seqj<-seq(-3, 3, 6/100)

Class <-matrix(0,nrow=length(seqi),ncol=length(seqj)) 
ci<-0
for(i in seqi){
  ci<-ci+1
  cj<-0
  for(j in seqj){
    cj<-cj+1
    Class[ci,cj]<- (c(1,i,j)%*%pesos_verticais)*(c(1,i,j)%*%pesos_horizontais)
  }
}

plot(x11[,1],x11[,2], col = 'green', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')
par(new=T)
plot(x12[,1],x12[,2], col = 'green', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')
par(new=T)
plot(x21[,1],x21[,2], col = 'blue', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')
par(new=T)
plot(x22[,1],x22[,2], col = 'blue', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = '', ylab = '')
par(new=T)
contour2D(Class, seqi, seqj, colkey = NULL, xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), xlab = 'X1', ylab = 'X2', levels = -0.1:0.1)


plot(vertical_classify[[2]][0:5], type='l', xlab = 'epocas', ylab = 'ERRO')
plot(Horizontal_classify[[2]][0:5], type='l', xlab = 'epocas', ylab = 'ERRO')

XTest <- rbind(XTestVertical, XTestHorizontal)
XTest <- cbind(1, XTest)
XTestClass <- rbind(as.matrix(XClassTestVertical), as.matrix(XClassTestHorizontal))
erroT <- cbind((XTest%*%pesos_verticais)*(XTest%*%pesos_horizontais),XTestClass)
erroMQ <- mean((erroT[2]-erroT[1])^2)