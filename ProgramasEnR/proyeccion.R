
library(wavelets)
library(waveslim)
library(wavethresh)
##############################################################################
###### ONDALETA HAAR ######
#layout(matrix(c(1:4), nrow=2, byrow=FALSE))

#layout.show(4)


TT<-512; t<-1:512
fl12<-function(s){-1*cos(pi*(2*s/TT))}
at<-fl12(t) 
plot.ts(at)
##################################################3
tamanho.T<-512 # filas
valor.J<-5 #2^{J-1}<sqrt{T}<2^{J}
tamanho.J<-2^(valor.J) # 32 columnas
#####################################################
#matrix.coef será una matriz que
#- calcula una matrix donde se encuentran los coeficient de wavelets evaluados
#- calcula la matrix Psi
######################################################
matrix.coef<-function(k){vetor<-rep(0,tamanho.T)
vetor[k]<-1
y=dwt(vetor,wf="haar",n.levels=6)
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


Beta<-(solve(t(matrix.wav)%*%matrix.wav))%*%(t(matrix.wav))%*%fl12(1:512)
length(Beta)
hat<-matrix.wav%*%Beta
plot.ts(hat,xlab="J=6",ylab="Haar")
lines(at)


