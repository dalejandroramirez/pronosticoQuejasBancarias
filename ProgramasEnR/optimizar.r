
library(wavelets)
library(waveslim)
library(wavethresh)
###########################################
###### ONDALETA HAAR ######
#par(mfrow=c(2,3))
L=1

TT<-512; t<-1:512
#set.seed(512)
e<-mvrnorm(TT,rep(0,L),diag(L))

fl12<-function(s){-1*cos(pi*(2*s/TT))}
at<-fl12(t) 
#plot.ts(at)
xt<-rep(0,TT)
xt<-matrix(0,TT,L)
for(j in 1:L){
  TTT=TT-1
for(i in 1:TTT){
  xt[i+1,j]<-at[i+1]*xt[i,j]+e[i+1,j]   
}}
dim(xt)
plot.ts(xt,xlab="Tiempo",main="tvAR(1)")

##########################################################
tamanho.T<-512 # filas
valor.J<-3 #2^{J-1}<sqrt{T}<2^{J}
tamanho.J<-2^(valor.J) # 32 columnas
##########################################################
#matrix.coef sera una matriz que
#-calcula una matrix donde se encuentran los coeficientes 
# de wavelets evaluados en un tiempo k
##########################################################
matrix.coef<-function(k){vetor<-rep(0,tamanho.T)
vetor[k]<-1
y=dwt(vetor,wf="d8",n.levels=6)
wav.dwt<-c(y$s6, y$d6,y$d5, y$d4, y$d3, y$d2,y$d1) 
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
M<-matrix(0,tamanho.J,L)
for (j in 1:L){
  ThetaX<-matrix(0,tamanho.T-1,tamanho.J)
  dim(ThetaX)
  for(i in 1:511){
  ThetaX[i,]<-matrix.wav[i+1,]*xt[i,j]
                  }  
dim(ThetaX)
Beta<-(solve(t(ThetaX)%*%ThetaX))%*%(t(ThetaX))%*%xt[2:tamanho.T,j]
length(xt[2:tamanho.T,1])
length(Beta)

hat<-matrix.wav%*%Beta
#plot.ts(hat,xlab="Tiempo",ylab="a(t)",main="D8")
M[,j]=Beta 

}


Beta=rep(0,tamanho.J)             ## usando el promedio
for(i in 1:tamanho.J){
  Beta[i]=mean(M[i,])
}
hat<-matrix.wav%*%Beta
print(max(abs(hat-at)))

for(i in 1:tamanho.J){           ##Usando la mediana
  Beta[i]=median(M[i,])  
}

hat<-matrix.wav%*%Beta
max(abs(hat-at))

MM=t(M)

#layout(matrix(c(1:2), nrow=1, byrow=FALSE))

#layout.show(2)


plot.ts(hat,xlab="Tiempo",ylab="a(t)",main="D8",sub="L=50",)
lines(at,col=2)
lines(hat,col=4)
legend("topleft", legend=c("Real", "Estimada"),
       col=c("red", "blue"), lty=1,cex=0.8,horiz = FALSE, bty="n")


boxplot(MM,main="Distribucion Coeficientes de D8")




