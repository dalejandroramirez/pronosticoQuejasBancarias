##Funcion de prueba de tamaño 512
## un modelo tvar(1)

set.seed(512)  #fija el aleatorio siguiente#
e<- rnorm(512,0,1)
TT<-512 ; t <- 1:512

f<- function(k){-1*cos(pi*(2*k/TT))}
at<-f(t)
plot.ts(at)
xt<-rep(0,512)    
for(i in 1:511){
  xt[i+1]<-at[i+1]*xt[i]+e[i+1]   
}
length(xt)
plot.ts(xt)
#####################################################################
T <- 512   ## tamaño de la muestra
J<- 5      ##nivel de resolucion optimo
dimEV<- 32  ## pues dimension de Vj= 2^5 
#####################################################################
###    realizar una funcion que me calcule las funcion de haar    ###
#####################################################################
M<-function(k){v<-rep(0,T)
v[k]<-1         ###aqui indico el tiempo en que quiero ver la ondaleta
wav.dwt<-c(dwt(v,wf="haar",n.levels=6)$s6, dwt(v,wf="haar",n.levels=6)$d6,   #¿No deberia ser el nivel 9?   
           dwt(v,wf="haar",n.levels=6)$d5, dwt(v,wf="haar",n.levels=6)$d4,     
           dwt(v,wf="haar",n.levels=6)$d3, dwt(v,wf="haar",n.levels=6)$d2)#,
v.Psi<-wav.dwt #[1:dimEV]                                                       #¿Porque no toma los otros?
v.Psi                                         ##esta normalizado ? como se cual es el coefciente correcto??
}
##############################################################################
### Matriz que contiene a psi y phi                                       ####      
##############################################################################
matrix.wav<-matrix(0, T, dimEV)
indice<-seq(1,T,1)
matrix.wav<-t(apply(as.array(indice),1, M))   ###  llena por columnas (el trasport es porque la funcion appy devulve una fila)
dim(matrix.wav)                               ### ese 1 porque quiero columna
#################################################################################
ThetaX<-matrix(0,T-1,dimEV)   ##columna llena de ceros ordern T-1
dim(ThetaX)
for(i in 1:511){
  ThetaX[i,]<-matrix.wav[i+1,]*xt[i]
}
dim(ThetaX)
print(ThetaX)
Beta<-(solve(t(ThetaX)%*%ThetaX))%*%(t(ThetaX))%*%xt[2:512]
length(Beta)
hat<-matrix.wav%*%Beta
print(hat)
plot.ts(hat)
lines(at)




dim(Beta)
dim(Cosa)
dim(ThetaX)
length(xt)


A= xt-hat*xt
print(A)
dim(A)
A=A%*%t(A)
A/T
print(A)
dim(A)
Beta<-(solve(t(ThetaX)%*%solve(A)%*%ThetaX))%*%(t(ThetaX))%*%solve(A)%*%xt[2:512]
length(Beta)
hat<-matrix.wav%*%Beta
print(hat)
length(xt)
plot.ts(hat)




        