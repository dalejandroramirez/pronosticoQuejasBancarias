
  library(wavelets)
  library(waveslim)
  library(wavethresh)
###########################################
###### ONDALETA HAAR ######
layout(matrix(c(1:2), nrow=1, byrow=FALSE))

layout.show(2)
#par(mfrow=c(2,3))

set.seed(512)
e<-rnorm(512,0,1)
TT<-512; t<-1:512
fl12<-function(s){-1*cos(pi*(2*s/TT))}
at<-fl12(t) 
plot.ts(at)
xt<-rep(0,512)    
for(i in 1:511){
  xt[i+1]<-at[i+1]*xt[i]+e[i+1]   
}
length(xt)
plot.ts(xt,xlab="Tiempo",main="tvAR(1)")

acf(xt,lag.max = 35,xlab="Retraso",ylab=" ",main="ACF")
pacf(xt,lag.max = 35,xlab="Retraso",ylab=" ",main="PACF")
##########################################################
tamanho.T<-512 # filas
valor.J<-5 #2^{J-1}<sqrt{T}<2^{J}
tamanho.J<-2^(valor.J) # 32 columnas
##########################################################
#matrix.coef sera una matriz que
#-calcula una matrix donde se encuentran los coeficientes 
# de wavelets evaluados en un tiempo k
##########################################################
matrix.coef<-function(k){vetor<-rep(0,tamanho.T)
vetor[k]<-1
y=dwt(vetor,wf="d8",n.levels=5)
wav.dwt<-c(y$s5, y$d5,y$d4, y$d3, y$d2, y$d1) 
vetor.Psi<-wav.dwt[1:tamanho.J]
vetor.Psi} 
# coloca los primeros tamanho.J wavelet en el tiempo k.
######################################################################
#Matrix que contenga los valores de psi phi
#####################################################################
matrix.wav<-matrix(0, tamanho.T, tamanho.J) #512x32
dim(matrix.wav)
indiceU<-seq(1,tamanho.T,1)
matrix.wav<-t(apply(as.array(indiceU),1, matrix.coef))
dim(matrix.wav)
#############Construir la matriz psi
ThetaX<-matrix(0,tamanho.T-1,tamanho.J)
dim(ThetaX)
for(i in 1:511){
  ThetaX[i,]<-matrix.wav[i+1,]*xt[i]
}
dim(ThetaX)


Beta<-(solve(t(ThetaX)%*%ThetaX))%*%(t(ThetaX))%*%xt[2:512]
length(Beta)
hat<-matrix.wav%*%Beta
plot.ts(hat,xlab="Tiempo",ylab="a(t)",main="D8")
lines(at,col=4)








##############################################################
##optimizacion
###############################################################
aux<-e
xt=matrix(xt)
e1<-xt[2:512,]-hat[2:512,]*xt[1:511,]

#R <- matrix(rep(1:511, times =511), nrow = 511)
#C <- matrix(rep(1:511, each = 511), nrow = 511)
#pos <- abs(R - C) + 1
#AA=acf(e1)$acf
#A <- matrix(AA[pos], nrow = 511)
A=e1%*%t(e1)  #aqui esta el prroblemaa#
dim(A)
print(A)

Beta<-(solve(t(ThetaX)%*%solve(A)%*%ThetaX))%*%(t(ThetaX))%*%solve(A)%*%xt[2:512]
length(Beta)
hat1<-matrix.wav%*%Beta
plot.ts(hat1,xlab="Tiempo",ylab="a(t)",main="D8")
lines(at)
