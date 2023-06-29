
library(wavelets)
library(waveslim)
library(wavethresh)
##############################################################################
###### ONDALETA HAAR ######
layout(matrix(c(1:2), nrow=1, byrow=FALSE))

layout.show(2)


TT<-512; t<-1:2157

fl12<-function(s){X[s]}
at<-rep(0,2176)
at[1:2157]<-fl12(t)
at=ts(at,start = c(2014,1),frequency = 365)
plot.ts(at)
##################################################3
tamanho.T<-2176 # filas
valor.J<-6 #2^{J-1}<sqrt{T}<2^{J}
tamanho.J<-2^(valor.J) # 32 columnas
#####################################################
#matrix.coef será una matriz que
#- calcula una matrix donde se encuentran los coeficient de wavelets evaluados
#- calcula la matrix Psi
######################################################
matrix.coef<-function(k){vetor<-rep(0,tamanho.T)
vetor[k]<-1
y=dwt(vetor,wf="haar",n.levels=6)     ##d8 es el otro 
wav.dwt<-c(y$s6,y$d6,y$d5,y$d4, y$d3, y$d2,y$d1) 
vetor.Psi<-wav.dwt[1:tamanho.J]
vetor.Psi} # coloca los primeros tamanho.J wavelet en el tiempo k.


######################################################################
#Matrix que contenga los valores de psi phi
#####################################################################
matrix.wav<-matrix(0, tamanho.T, tamanho.J) #512x32
dim(matrix.wav)
indiceU<-seq(1,tamanho.T,1)
matrix.wav<-t(apply(as.array(indiceU),1, matrix.coef))
dim(matrix.wav)


Beta<-(solve(t(matrix.wav)%*%matrix.wav))%*%(t(matrix.wav))%*%at
length(Beta)
hat<-matrix.wav%*%Beta
hat=ts(hat,start =c(2014,1),frequency = 365 )
plot.ts(hat,xlab="J=6",ylab="Haar")


#lines(at)
###############################################################################

##############################################################################
###### ONDALETA D8 ######
 t<-1:2176
#fl12<-function(s){-1*cos(pi*(2*s/TT))}
#at<-fl12(t) 
#plot.ts(at)
xt<-rep(0,2176)    
for(i in 1:2157){
  xt[i]<-X[i]   
}
length(xt)
xt=ts(xt,start = c(2014,1),frequency = 365)
plot.ts(xt)
##################################################3
tamanho.T<-2176 # filas
valor.J<-6 #2^{J-1}<sqrt{T}<2^{J}
tamanho.J<-2^(valor.J) # 32 columnas
#####################################################
#matrix.coef será una matriz que
#- calcula una matrix donde se encuentran los coeficient de wavelets evaluados
#- calcula la matrix Psi
######################################################
matrix.coef<-function(k){vetor<-rep(0,tamanho.T)
vetor[k]<-1
y=dwt(vetor,wf="d8",n.levels=6)
wav.dwt<-c(y$s6, y$d6,y$d5, y$d4, y$d3, y$d2) 
vetor.Psi<-wav.dwt[1:tamanho.J]
vetor.Psi} # coloca los primeros tamanho.J wavelet en el tiempo k.


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
for(i in 1:2175){
  ThetaX[i,]<-matrix.wav[i+1,]*xt[i]
}
dim(ThetaX)


Beta<-(solve(t(ThetaX)%*%ThetaX))%*%(t(ThetaX))%*%xt[2:2176]
length(Beta)
hat<-matrix.wav%*%Beta
hat=ts(hat,start = c(2014,1),frequency = 365)
plot.ts(hat)




