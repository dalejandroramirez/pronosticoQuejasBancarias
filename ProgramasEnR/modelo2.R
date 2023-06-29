library(wavelets)
library(waveslim)
library(wavethresh)
##############################################################################
###### ONDALETA HAAR ######
set.seed(512)
e<-rnorm(512,0,1)
TT<-512; t<-1:512
fl12<-function(s){-1*cos(pi*(2*s/TT))}
fl11<-function(s){-1*sin(pi*(2*s/TT))}
at<-fl12(t)
bt<-fl11(t)
plot.ts(at)
plot.ts(bt)
xt<-rep(0,TT)    
for(i in 1:511){
  xt[i+1]<-at[i+1]*xt[i]+e[i+1]   
}
#for(i in 1:510){
#  xt[i+2]<-at[i+2]*xt[i+1]+bt[i+2]*xt[i]+e[i+2]   
#}

length(xt)
plot.ts(xt)

##################################################3
tamanho.T<-length(xt) # filas
p=1
v=tamanho.T-p
valor.J<-5 #2^{J-1}<sqrt{T}<2^{J}
tamanho.J<-2^(valor.J) # 32 columnas
acf(xt, lag.max = 100)
pacf(xt, lag.max = 100)   ### Porque con este modelo no deja??

#####################################################
#matrix.coef será una matriz que
#- calcula una matrix donde se encuentran los coeficient de wavelets evaluados
#- calcula la matrix Psi
######################################################
matrix.coef<-function(k){vetor<-rep(0,tamanho.T)

vetor[k]<-1
y=dwt(vetor,wf="haar",n.levels=5)
wav.dwt<-c(y$s5, y$d5,y$d4, y$d3, y$d2, y$d1) 
vetor.Psi<-wav.dwt[1:tamanho.J]
length(vetor.Psi)
vetor.Psi} # coloca los primeros tamanho.J wavelet en el tiempo k.


######################################################################
#Matrix que contenga los valores de psi phi
#####################################################################

matrix.wav<-matrix(0, tamanho.T, tamanho.J) #512 x32
dim(matrix.wav)
indiceU<-seq(1,tamanho.T,1)
matrix.wav<-t(apply(as.array(indiceU),1, matrix.coef))
dim(matrix.wav)                                  ##aqui tengo una duda
#############Construir la matriz psi
s=p*tamanho.J
ThetaY=matrix(0,v,s)
a=1
b=32
for(i in 1:p){
  ThetaX=matrix(0,v,tamanho.J)
for(j in 1:v){
  dim(matrix.wav)
  J=j+p
  P=p+j-i
  ThetaX[j,]<-matrix.wav[J,]*xt[P]   ###aca tengo una duda
}

  ThetaY[,a:b]=ThetaX
  a=b+1
  b=b+32
}
dim(ThetaY)
ThetaY<-ThetaY[,1:32]
t(ThetaY)%*%ThetaY
psu=p+1
Beta<-(solve(t(ThetaY)%*%ThetaY))%*%(t(ThetaY))%*%xt[psu:tamanho.T]
length(Beta)
dim(ThetaY)
#matrix.wav<-cbind(matrix.wav,matrix.wav,matrix.wav)
dim(matrix.wav)
hat<-matrix.wav[1:tamanho.T,]%*%Beta[1:tamanho.J]
plot.ts(hat)
lines(at)
###############################################################




