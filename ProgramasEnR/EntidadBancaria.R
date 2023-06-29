
library(dplyr)
library(ggplot2)
library(PASWR)
library(urca)
library(MASS) # Boxcox function
library(Metrics)
#library(fpp3)  # encuenntro el lambda con BoxCox.numeric(X,interval = c(-1,1))
library(moments)  #prueba de normalidad de jarque.test( )
library(astsa)
#library(sarima)

###########################################################################3
Y <- Serie_Reclamos$ N ;   
Y <- rev(Y)

Fe<- Serie_Reclamos$FECHA   

plot(Fe,Y)

X=ts(Y,start=c(2014,1),frequency = 7)##Grafica de Y como serie de tiempo
                                            ##frecuencia me dice como estan organizados
                                            ##si lo dejo 365 esta ordenado por anho 1 ,2,3,4,5
boxplot(X ~ cycle(X))                       ## Si queremos comparar la distribución de las quejas para cada semana,

#X=ts(Y,start=c(2014,1))
##########################################################################
### boxplox anual por cada semana 
#################################################################


MIE<-X14[seq(1,363,7)]
JUE<-X14[seq(2,363,7)]
VIE<-X14[seq(3,363,7)]
SAB<-X14[seq(4,363,7)]
DOM<-X14[seq(5,363,7)]
LUN<-X14[seq(6,363,7)]
MAR<-X14[seq(7,364,7)]
BX14<-cbind(MIE,JUE,VIE,SAB,DOM,LUN,MAR)
#boxplot(BX14)

par(mfrow=c(2,3))

X14<-ts(X[1:363],start=c(2014,1),frequency = 7) ## primer dia es miercoles
NOMBRES=c("MIE","JUE","VIE","SAB","DOM","LUN","MAR")
boxplot(X14 ~ cycle(X14),main="2014",xlab="",names=NOMBRES,ylab = "")

X15=ts(X[364:727],start=c(2014,1),frequency = 7) ##  primer dia es jueves
NOMBRES=c("JUE","VIE","SAB","DOM","LUN","MAR","MIE")
boxplot(X15~ cycle(X15),main="2015",xlab="",ylab = "",names=NOMBRES)

X16=ts(X[728:1094],start=c(2014,1),frequency = 7) ##  primer dia es viernes
NOMBRES=c("VIE","SAB","DOM","LUN","MAR","MIE","JUE")
boxplot(X16 ~ cycle(X16),main="2016",xlab="",ylab = "",names=NOMBRES)

X17=ts(X[1095:1460],start=c(2014,1),frequency = 7) ## primer dia es Sabado
NOMBRES=c("SAB","DOM","LUN","MAR","MIE","JUE","VIE")
boxplot(X17 ~ cycle(X17),main="2017",xlab="",ylab = "",names=NOMBRES)

X18=ts(X[1461:1826],start=c(2014,1),frequency = 7) ## primer dia es lunes
NOMBRES=c("LUN","MAR","MIE","JUE","VIE","SAB","DOM")
boxplot(X18 ~ cycle(X18),main="2018",xlab="",ylab = "",names=NOMBRES)

X19=ts(X[1827:2157],start=c(2014,1),frequency = 7) ## primer dia es martes
NOMBRES=c("MAR","MIE","JUE","VIE","SAB","DOM","LUN")
boxplot(X19 ~ cycle(X19),main="2019",xlab="",ylab = "",names=NOMBRES)

##########################################################################

X=ts(Y,start=c(2014,1),frequency = 1)

##2014 tuvo 363 dias, Marzo 1 (festivo) 29 septiembre
##2015 tuvo 363 dias, Octubre 11 , 21 noviembre
##2016 366 bisiesto, 2017 y 2018 normal, 2019 hasta diciembre 2 (334) dias
##########################################################################
#             visualizacion de la serie 
##########################################################################
plot(X,xlab="Tiempo",ylab="Quejas")
EDA(X)   ## me arroja todos los datos estadisticos necesarios 
acf(X)
pacf(X)
hist(X)
length(X)
plot(density(X),col=4, main = "Densidad de probabilidad")  ##Genera la funcion de densidad de distribucion

########################################################################
   ### Calcular el periodograma en cada anho y toda la seria
########################################################################

par(mfrow=c(2,3))

X.per = spec.pgram(X[1:363], taper=0, log="no",main=2014,sub="",xlab="Frecuencias")   ##2014
X.per = spec.pgram(X[364:727], taper=0, log="no",main=2015,sub="",xlab="Frecuencias")  ##2015
X.per = spec.pgram(X[728:1094], taper=0, log="no",main=2016,sub="",xlab="Frecuencias")  #2016
X.per = spec.pgram(X[1095:1460], taper=0, log="no",main=2017,sub="",xlab="Frecuencias")  #2017
X.per = spec.pgram(X[1461:1826], taper=0, log="no",main=2018,sub="",xlab="Frecuencias")  ##2018
X.per = spec.pgram(X[1827:2157], taper=0, log="no",main=2019,sub="",xlab="Frecuencias")  ##2019
X.per = spec.pgram(X, taper=0, log="no",main="Periodograma")
#X.per = spec.pgram(X, taper=0, log="no",kernel("modified.daniell",c(5,7)))  #14-19
#############################################x#################

#X.per$spec
#par(new=T)
X.per = spec.pgram(X, taper=0, log="no")
#X.per = spec.pgram(X, taper=0, log="no",kernel("modified.daniell",c(5,7)))
#plot(X.per$spec,xlab = " ",ylab = " ")
##
##puntos de maximos 
#(52.21528,180.9488)  plot(X.per$freq[309],X.per$spec[309],xlab = " ",ylab = " ")
#(104.4306,44.55463)  plot(X.per$freq[618],X.per$spec[618],xlab = " ",ylab = " ")  
#(156.4769,6.30687)    X.per$freq[926]
################################################################
#                 w = (365/2160)k k=0,1,2,...,1080
################################################################
#  transformacion para obtener una serie estacionaria
###############################################################

X1 = decompose(X)
plot(X1, xlab='Año')

dqfX=ur.df(X,type="none",lags=7,selectlags = c("AIC"))  ##criterio de dickey fuller

summary(dqfX)
##Como tau 1 es menor que el valor del test entonces no es estacionario

#BoxCox.numeric(X,interval = c(-1,1))  ## encontrar lambda de box cox
lamda=0.675
X=ts(Y,start=c(2014,1),frequency = 365)

diffX=diff(sqrt(X))



dqfX=ur.df(diffX,type="none",lags=7,selectlags = c("AIC"))

summary(dqfX)
 plot(diffX)
 par(mfrow=c(1,2))
acf2(diffX,48)      #Calcula tambien el PACF

#pacf(diffX,lag.max = 45)

##Como tau1 es mayor que el test estadistico entonces 
##cumple el test de dicket fuller (raices unitarias)
##########################################################

##¿Que pasa si le quito el periodo?

diff7.diffX=diff(diffX,lag=7)


plot(diff7.diffX)

dqf7X=ur.df(diff7.diffX,type="none",lags=7,selectlags = c("AIC"))

summary(dqf7X)
#par(mfrow=c(1,2))
acf2(diff7.diffX,48)
#pacf(diff7.diffX,lag.max = 60)

###tambian pasa el test de estacionariedad 

####################################################################
# Modelo SARIMA
####################################################################


acf2(X, 48)
acf2(diff(X), 48)
acf2(diff(X,7))
acf2(diff(diff(sqrt(X), 7)),35)
#########################################################
#Modelos tentativos
##############################################################
mod1=sarima(X^(1/7),1,1,1,1,1,1,7)# fit model (ii)

hist(mod1$fit$residuals)


resud=mod1$fit$residuals
truehist(resud,col="red",xlab="",ylab="",main="Distribución Residuales")
lines(density(resud),col ="blue", lwd = 1, add = FALSE)

curve(dnorm(x, mean(mod1$fit$residuals), sd(mod1$fit$residuals)),col = "blue", lwd = 2, add = TRUE)
plot(mod1$fit$residuals)

###############################################################
#Pronostico
###############################################################
X1=X[2003:2157]
X2=X[1:2002]
X.test=X1^(1/7)
X.train=X2^(1/7)
forest=sarima.for(X.train, 155, 1, 1, 1, 1, 1, 1, 7) # forecast
X.forest=forest$pred

MAPE<-function(test,estimado){
  M1<-abs(estimado-test)
  M1<-M1/abs(test)
  M<-mean(M1)
  return(mean(M))
}
ECM<-function(test,estimado){
  M1<-abs(estimado-test)
  M<-mean(M1)
  return(mean(M))
}
MAPE(X.test,X.forest)
ECM(X.test,X.forest)

mape(X.test,X.forest)


## Δy(t)=Gamma∗y(t−1)+e(t)
##si pongo "None" :   la Hipotesis nula es Gamma=0 
#Si el estadistico de prueba esta dentro de las 3 regiones 
#tau1 : 1%  5%  10%  

#https://stats.stackexchange.com/questions/24072/interpreting-rs-ur-df-dickey-fuller-unit-root-test-results

#https://rpubs.com/joser/SeriesTemporalesBasicas
