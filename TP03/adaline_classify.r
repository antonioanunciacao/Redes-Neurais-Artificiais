rm(list=ls())
library('plot3D')

# Definiçao dos dados em duas dimensões:
# Declarando as variancias:
s1 <- s2 <- 0.4
nc<-200

# Criando os dados:
xc1<-matrix(rnorm(nc*2), ncol=2)*s1 + t(matrix((c(2,2)), ncol=nc, nrow=2))
xc2<-matrix(rnorm(nc*2), ncol=2)*s2 + t(matrix((c(4,4)), ncol=nc, nrow=2))

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

# Definição dos dados de treinamento e de teste:
Mtrain<- 0.9*dim(xc1)[1]
seqc1 <- sample(dim(xc1)[1])
seqc2 <- sample(dim(xc2)[1])

xc1train <- xc1[seqc1[(1:Mtrain)], ]
xc1tst <- xc1[seqc1[((Mtrain+1):dim(xc1)[1])], ]
xc1Class <- matrix(nrow = dim(xc1train)[1], ncol = 1,-1)
xc1Classtst <- matrix(nrow = dim(xc1tst)[1], ncol = 1,-1)

xc2train <- xc2[seqc2[(1:Mtrain)], ]
xc2tst <- xc2[seqc2[((Mtrain+1):dim(xc2)[1])], ]
xc2Class <- matrix(nrow = dim(xc2train)[1], ncol = 1,1)
xc2Classtst <- matrix(nrow = dim(xc2tst)[1], ncol = 1,1)

xcTrain <- rbind(xc1train, xc2train)
xcClass <- rbind(xc1Class, xc2Class)

xcTst <- rbind(xc1tst, xc2tst)
xcClassTst <- rbind(xc1Classtst, xc2Classtst)

eta <- 0.01
tol <- 10^-3
maxepocas <- 10^3
par <- 1

dataTrain <- trainadaline(xcTrain, xcClass, eta, tol, maxepocas, par)

seqi<-seq(0, 6, 6/100)
seqj<-seq(-0, 6, 6/100)

Class <-matrix(0,nrow=length(seqi),ncol=length(seqj)) 
ci<-0
for(i in seqi){
  ci<-ci+1
  cj<-0
  for(j in seqj){
    cj<-cj+1
    Class[ci,cj]<- c(1,i,j)%*%dataTrain[[1]]
  }
}

# Plotando os graficos
plot(xc1[,1],xc1[,2], col = 'red', xlim = c(0,6), ylim = c(0,6), xlab = 'X1', ylab = 'X2')
par(new=T)
plot(xc2[,1],xc2[,2], col = 'blue', xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '')
par(new=T)
contour2D(Class, seqi, seqj, colkey = NULL, xlim = c(0,6), ylim = c(0,6), xlab = '', ylab = '', levels= 0, axis = F)

pesos <- dataTrain[[1]]
validation <- cbind(cbind(1,xcTst)%*%pesos,xcClassTst)
erro <- (validation[,2]-validation[,1])^2
erro_mQ <- mean(erro)

plot(dataTrain[[2]][1:50], type = 'l', xlab = 'epoca', ylab = 'ERRO')