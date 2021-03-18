#***************************************************************************************
#***************************** ITAM Diplomado Estadistica Aplicada *********************
#***************************** Modulo 2. Muestreo  *************************************
#***************************** Proyecto Practico  **************************************
#***************************************************************************************

#********* 0. Preparar Datos ********* 
#Clear console with cat("\014")
# Clear environment variables rm(list = ls())
df  <- read.csv(file = '007_ResElect_NACIONAL_seccs_PRES2012.csv',fileEncoding="latin1")
#fileEncoding='latin1'
head(df)
View(df)
nrow(df)
df<-df[complete.cases(df), ] #clear NA values
nrow(df) 
colnames(df) 
df$JVM<-df$PAN
df$EPN<-df$PRI+df$COALICIÓN.PRI.PVEM+df$PVEM
df$AMLO<-df$PRD+df$PT+df$MOVIMIENTO.CIUDADANO+df$COALICIÓN.PRD.PT.MC+df$COALICIÓN.PRD.PT+df$COALICIÓN.PRD.MC+df$COALICIÓN.PT.MC
df$QCT<-df$NUEVA.ALIANZA
df$NoReq=df$NO.REGISTRADOS
df$Nulos=df$NULOS
colnames(df)
TotalEPN=sum(df$EPN)/sum(df$TOTAL)
TotalJVM=sum(df$JVM)/sum(df$TOTAL)
TotalAMLO=sum(df$AMLO)/sum(df$TOTAL)
TotalQCT=sum(df$QCT)/sum(df$TOTAL)
TotalNoReq=sum(df$NoReq)/sum(df$TOTAL)
TotalNulos=sum(df$Nulos)/sum(df$TOTAL)
Candidatos=c('EPN','JVM','AMLO','QCT','NoReq','Nulos')
Totales=c(TotalEPN,TotalJVM,TotalAMLO,TotalQCT,TotalNoReq,TotalNulos)
names(Totales)<-Candidatos
Totales*100

colnames(df)
#***************************************************************************************
#***************************** 1. Muestro Aleatorio Simple *****************************
#***************************************************************************************
library(sampling)
require(samplingVarEst)

n=500
N=nrow(df)
#1.1 Seleccion aleatoria
s.SI1.U<-srswor1(n,N) 

#1.2 Var interes

VecY.s.SI.EPN         <- df$EPN[s.SI1.U==1] #Mascara True para X Candidato
VecY.s.SI.AMLO         <- df$AMLO[s.SI1.U==1]
VecY.s.SI.JVM         <- df$JVM[s.SI1.U==1]
VecY.s.SI.QCT         <- df$QCT[s.SI1.U==1]
VecY.s.SI.NoReq         <- df$NoReq[s.SI1.U==1]
VecY.s.SI.Nulos         <- df$Nulos[s.SI1.U==1]

#1.3 PKs

VecPk.s <- rep(n/N, times=n)

#1.4 Est Narain, Horvitz-Thompson

EstThetaEPN <- Est.Total.NHT(VecY.s.SI.EPN, VecPk.s)
EstThetaAMLO <- Est.Total.NHT(VecY.s.SI.AMLO, VecPk.s)
EstThetaJVM <- Est.Total.NHT(VecY.s.SI.JVM, VecPk.s)
EstThetaQCT <- Est.Total.NHT(VecY.s.SI.QCT, VecPk.s)
EstThetaNoReq <- Est.Total.NHT(VecY.s.SI.NoReq, VecPk.s)
EstThetaNulos <- Est.Total.NHT(VecY.s.SI.Nulos, VecPk.s)

# Ver estimaciones
EstThetaEPN
EstThetaAMLO
EstThetaJVM 
EstThetaQCT
EstThetaNoReq 
EstThetaNulos 

#1.5 Variacion estimada de cada Variacion
#1.5.1 Matrix de probabilidades de Inclusión
MatPkl.s           <- Pkl.Hajek.s(VecPk.s)
dim(MatPkl.s )
length(VecPk.s)

#1.5.2 Varianzas
EstVarEstThetaEPN    <- VE.HT.Total.NHT(VecY.s.SI.EPN, VecPk.s, MatPkl.s)
EstVarEstThetaAMLO    <- VE.HT.Total.NHT(VecY.s.SI.AMLO, VecPk.s, MatPkl.s)
EstVarEstThetaJVM    <- VE.HT.Total.NHT(VecY.s.SI.AMLO, VecPk.s, MatPkl.s)
EstVarEstThetaQCT    <- VE.HT.Total.NHT(VecY.s.SI.QCT, VecPk.s, MatPkl.s)
EstVarEstThetaNoReq    <- VE.HT.Total.NHT(VecY.s.SI.NoReq, VecPk.s, MatPkl.s)
EstVarEstThetaNulos    <- VE.HT.Total.NHT(VecY.s.SI.Nulos, VecPk.s, MatPkl.s)

#1.5.3 Error Estandar
StdErrEstThetaEPN    <- sqrt(EstVarEstThetaEPN)
StdErrEstThetaAMLO    <- sqrt(EstVarEstThetaAMLO)
StdErrEstThetaJVM    <- sqrt(EstVarEstThetaJVM)
StdErrEstThetaQCT    <- sqrt(EstVarEstThetaQCT)
StdErrEstThetaNoReq    <- sqrt(EstVarEstThetaNoReq)
StdErrEstThetaNulos    <- sqrt(EstVarEstThetaNulos)

#1.6 Error Absoluto o precision al 95%
alpha=0.05
AbsErrEstThetaEPN    <- StdErrEstThetaEPN*qnorm(1-alpha/2)
AbsErrEstThetaAMLO    <- StdErrEstThetaAMLO*qnorm(1-alpha/2)
AbsErrEstThetaJVM    <- StdErrEstThetaJVM*qnorm(1-alpha/2)
AbsErrEstThetaQCT    <- StdErrEstThetaQCT*qnorm(1-alpha/2)
AbsErrEstThetaNoReq    <- StdErrEstThetaNoReq*qnorm(1-alpha/2)
AbsErrEstThetaNulos    <- StdErrEstThetaNulos*qnorm(1-alpha/2)

#1.7 Intervalos de Confianza
LimInfICEstThetaEPN  <- EstThetaEPN - AbsErrEstThetaEPN
LimInfICEstThetaAMLO  <- EstThetaAMLO - AbsErrEstThetaAMLO
LimInfICEstThetaJVM  <- EstThetaJVM - AbsErrEstThetaJVM
LimInfICEstThetaQCT  <- EstThetaQCT - AbsErrEstThetaQCT
LimInfICEstThetaNoReq  <- EstThetaNoReq - AbsErrEstThetaNoReq
LimInfICEstThetaNulos  <- EstThetaNulos - AbsErrEstThetaNulos

LimSupICEstThetaEPN  <- EstThetaEPN + AbsErrEstThetaEPN
LimSupICEstThetaAMLO  <- EstThetaAMLO + AbsErrEstThetaAMLO
LimSupICEstThetaJVM  <- EstThetaJVM + AbsErrEstThetaJVM
LimSupICEstThetaQCT  <- EstThetaQCT + AbsErrEstThetaQCT
LimSupICEstThetaNoReq  <- EstThetaNoReq + AbsErrEstThetaNoReq
LimSupICEstThetaNulos  <- EstThetaNulos + AbsErrEstThetaNulos

#1.8 Coeficiente de Variacion Estimada (CVE)
CVEEstThetaEPN       <- StdErrEstThetaEPN/EstThetaEPN
CVEEstThetaAMLO       <- StdErrEstThetaAMLO/EstThetaAMLO
CVEEstThetaJVM       <- StdErrEstThetaJVM/EstThetaJVM
CVEEstThetaQCT       <- StdErrEstThetaQCT/EstThetaQCT
CVEEstThetaNoReq       <- StdErrEstThetaNoReq/EstThetaNoReq
CVEEstThetaNulos       <- StdErrEstThetaNulos/EstThetaNulos

#1.9 Output 1
OUTPUT1            <- c(EstThetaEPN, EstThetaAMLO, EstThetaJVM, EstThetaQCT,EstThetaNoReq, EstThetaNulos)
OUTPUT1            <- cbind(EstTheta = OUTPUT1, StdErr = c(StdErrEstThetaEPN, StdErrEstThetaAMLO, StdErrEstThetaJVM, StdErrEstThetaQCT, StdErrEstThetaNoReq, StdErrEstThetaNulos))
OUTPUT1            <- cbind(OUTPUT1, LInfCI95 = c(LimInfICEstThetaEPN, LimInfICEstThetaAMLO, LimInfICEstThetaJVM, LimInfICEstThetaQCT, LimInfICEstThetaNoReq, LimInfICEstThetaNulos))
OUTPUT1            <- cbind(OUTPUT1, LSupCI95 = c(LimSupICEstThetaEPN, LimSupICEstThetaAMLO, LimSupICEstThetaJVM, LimSupICEstThetaQCT,LimSupICEstThetaNoReq, LimSupICEstThetaNulos))
OUTPUT1            <- cbind(OUTPUT1, CVE = c(CVEEstThetaEPN, CVEEstThetaAMLO, CVEEstThetaJVM, CVEEstThetaQCT,CVEEstThetaNoReq, CVEEstThetaNulos))
OUTPUT1

#***************************************************************************************
#***************************** 2. Muestreo Proporcional *****************************
#***************************************************************************************


#2.1 Construir probabilidades de Muestreo Proporcional
colnames(df)
#LISTA.NOMINAL
#2.1.1
VecPk.U            <- Pk.PropNorm.U(n, df$LISTA.NOMINAL)
length(VecPk.U)
#2.2 MatPk
#MatPkl.U           <- Pkl.Hajek.U(VecPk.U)   #-> OVERFLOW ERROR

#2.2.1' Muestreo de Brewer

s.BrX.U            <- UPbrewer(VecPk.U) #DUDA!
length(s.BrX.U)

VecPk.sX           <- VecPk.U[s.BrX.U==1]
length(VecPk.sX)

#2.2 Matrix PKl
MatPkl.U           <- Pkl.Hajek.s(VecPk.sX) # Cambiando al VecPk.sX de la muestra
# DUDA DUDA DUDA
MatPkl.sX <-MatPkl.U
MatPkl.sX[1:5,1:5]
dim(MatPkl.sX)
#MatPkl.sX          <- MatPkl.U[s.BrX.U==1,s.BrX.U==1]

#2.3 Datos Muestrales
VecY.s.BrEPN         <- df$EPN[s.BrX.U==1]
VecY.s.BrAMLO         <- df$AMLO[s.BrX.U==1]
VecY.s.BrJVM         <- df$JVM[s.BrX.U==1]
VecY.s.BrQCT         <- df$QCT[s.BrX.U==1]
VecY.s.BrNoReq         <- df$NoReq[s.BrX.U==1]
VecY.s.BrNulos         <- df$Nulos[s.BrX.U==1]

#2.4 Narain, Horvitz-Thompson
EstThetaEPN         <- Est.Total.NHT(VecY.s.BrEPN, VecPk.sX)
EstThetaAMLO         <- Est.Total.NHT(VecY.s.BrAMLO, VecPk.sX)
EstThetaJVM          <- Est.Total.NHT(VecY.s.BrJVM, VecPk.sX)
EstThetaQCT          <- Est.Total.NHT(VecY.s.BrQCT, VecPk.sX)
EstThetaNoReq         <- Est.Total.NHT(VecY.s.BrNoReq, VecPk.sX)
EstThetaNulos         <- Est.Total.NHT(VecY.s.BrNulos, VecPk.sX)

#2.5 Calcular Varianzas con SYG
EstVarEstThetaEPN    <- VE.SYG.Total.NHT(VecY.s.BrEPN, VecPk.sX, MatPkl.sX)
EstVarEstThetaAMLO    <- VE.SYG.Total.NHT(VecY.s.BrAMLO, VecPk.sX, MatPkl.sX)
EstVarEstThetaJVM    <- VE.SYG.Total.NHT(VecY.s.BrJVM, VecPk.sX, MatPkl.sX)
EstVarEstThetaQCT    <- VE.SYG.Total.NHT(VecY.s.BrQCT, VecPk.sX, MatPkl.sX)
EstVarEstThetaNoReq    <- VE.SYG.Total.NHT(VecY.s.BrNoReq, VecPk.sX, MatPkl.sX)
EstVarEstThetaNulos    <- VE.SYG.Total.NHT(VecY.s.BrNulos, VecPk.sX, MatPkl.sX)

#2.6 Error Estandar
StdErrEstThetaEPN    <- sqrt(EstVarEstThetaEPN)
StdErrEstThetaAMLO    <- sqrt(EstVarEstThetaAMLO)
StdErrEstThetaJVM    <- sqrt(EstVarEstThetaJVM)
StdErrEstThetaQCT    <- sqrt(EstVarEstThetaQCT)
StdErrEstThetaNoReq    <- sqrt(EstVarEstThetaNoReq)
StdErrEstThetaNulos    <- sqrt(EstVarEstThetaNulos)

#2.7 Error Absoluto al 95%
alpha              <- 0.05
AbsErrEstThetaEPN    <- StdErrEstThetaEPN*qnorm(1-alpha/2)
AbsErrEstThetaAMLO    <- StdErrEstThetaAMLO*qnorm(1-alpha/2)
AbsErrEstThetaJVM    <- StdErrEstThetaJVM*qnorm(1-alpha/2)
AbsErrEstThetaQCT    <- StdErrEstThetaQCT*qnorm(1-alpha/2)
AbsErrEstThetaNoReq    <- StdErrEstThetaNoReq*qnorm(1-alpha/2)
AbsErrEstThetaNulos    <- StdErrEstThetaNulos*qnorm(1-alpha/2)

#2.8 Intervalos de Confianza
LimInfICEstThetaEPN  <- EstThetaEPN - AbsErrEstThetaEPN
LimInfICEstThetaAMLO  <- EstThetaAMLO - AbsErrEstThetaAMLO
LimInfICEstThetaJVM  <- EstThetaJVM - AbsErrEstThetaJVM
LimInfICEstThetaQCT  <- EstThetaQCT - AbsErrEstThetaQCT
LimInfICEstThetaNoReq  <- EstThetaNoReq - AbsErrEstThetaNoReq
LimInfICEstThetaNulos  <- EstThetaNulos - AbsErrEstThetaNulos

LimSupICEstThetaEPN  <- EstThetaEPN + AbsErrEstThetaEPN
LimSupICEstThetaAMLO  <- EstThetaAMLO + AbsErrEstThetaAMLO
LimSupICEstThetaJVM  <- EstThetaJVM + AbsErrEstThetaJVM
LimSupICEstThetaQCT  <- EstThetaQCT + AbsErrEstThetaQCT
LimSupICEstThetaNoReq  <- EstThetaNoReq + AbsErrEstThetaNoReq
LimSupICEstThetaNulos  <- EstThetaNulos + AbsErrEstThetaNulos

#2.9 Coeficiente de Variacion Estimado

CVEEstThetaEPN       <- StdErrEstThetaEPN/EstThetaEPN
CVEEstThetaAMLO       <- StdErrEstThetaAMLO/EstThetaAMLO
CVEEstThetaJVM       <- StdErrEstThetaJVM/EstThetaJVM
CVEEstThetaQCT       <- StdErrEstThetaQCT/EstThetaQCT
CVEEstThetaNoReq      <- StdErrEstThetaNoReq/EstThetaNoReq
CVEEstThetaNulos       <- StdErrEstThetaNulos/EstThetaNulos

#2.10 DEFF

VE.SYG.Total.NHT.EPN=VE.SYG.Total.NHT(VecY.s.BrEPN, VecPk.s, Pkl.Hajek.s(VecPk.s))
VE.SYG.Total.NHT.AMLO=VE.SYG.Total.NHT(VecY.s.BrAMLO, VecPk.s, Pkl.Hajek.s(VecPk.s))
VE.SYG.Total.NHT.JVM=VE.SYG.Total.NHT(VecY.s.BrJVM, VecPk.s, Pkl.Hajek.s(VecPk.s))
VE.SYG.Total.NHT.QCT=VE.SYG.Total.NHT(VecY.s.BrQCT, VecPk.s, Pkl.Hajek.s(VecPk.s))
VE.SYG.Total.NHT.NoReq=VE.SYG.Total.NHT(VecY.s.BrNoReq, VecPk.s, Pkl.Hajek.s(VecPk.s))
VE.SYG.Total.NHT.Nulos=VE.SYG.Total.NHT(VecY.s.BrNulos, VecPk.s, Pkl.Hajek.s(VecPk.s))

deffEstThetaEPN      <- EstVarEstThetaEPN/VE.SYG.Total.NHT.EPN
deffEstThetaAMLO      <- EstVarEstThetaAMLO/VE.SYG.Total.NHT.AMLO
deffEstThetaJVM      <- EstVarEstThetaJVM/VE.SYG.Total.NHT.JVM
deffEstThetaQCT      <- EstVarEstThetaQCT/VE.SYG.Total.NHT.QCT
deffEstThetaNoReq      <- EstVarEstThetaNoReq/VE.SYG.Total.NHT.NoReq
deffEstThetaNulos      <- EstVarEstThetaNulos/VE.SYG.Total.NHT.Nulos

#2.11 OUTPUT

OUTPUT2            <- c(EstThetaEPN, EstThetaAMLO, EstThetaJVM, EstThetaQCT,EstThetaNoReq, EstThetaNulos)
OUTPUT2            <- cbind(EstTheta = OUTPUT2, StdErr = c(StdErrEstThetaEPN, StdErrEstThetaAMLO,
                                                           StdErrEstThetaJVM, StdErrEstThetaQCT,
                                                           StdErrEstThetaNoReq, StdErrEstThetaNulos))
OUTPUT2            <- cbind(OUTPUT2, LInf = c(LimInfICEstThetaEPN, LimInfICEstThetaAMLO,
                                                  LimInfICEstThetaJVM, LimInfICEstThetaQCT,
                                                  LimInfICEstThetaNoReq, LimInfICEstThetaNulos))
OUTPUT2            <- cbind(OUTPUT2, LSup = c(LimSupICEstThetaEPN, LimSupICEstThetaAMLO,
                                                  LimSupICEstThetaJVM, LimSupICEstThetaQCT,
                                                  LimSupICEstThetaNoReq, LimSupICEstThetaNulos))
OUTPUT2            <- cbind(OUTPUT2, CVE = c(CVEEstThetaEPN, CVEEstThetaAMLO, CVEEstThetaJVM, CVEEstThetaQCT,
                                             CVEEstThetaNoReq, CVEEstThetaNulos))
OUTPUT2            <- cbind(OUTPUT2, deff = c(deffEstThetaEPN, deffEstThetaAMLO, deffEstThetaJVM, deffEstThetaQCT,
                                              deffEstThetaNoReq, deffEstThetaNulos))
OUTPUT2

#3 Comparacion
summary(1/VecPk.s)
plot(sort(1/VecPk.s))
summary(1/VecPk.sX)
plot(sort(1/VecPk.sX))


#Valor Real
ThetaTrue             <- c(sum(df$EPN),sum(df$AMLO),sum(df$JVM),sum(df$QCT),
                           sum(df$NoReq),sum(df$Nulos))
ThetaTrue
OUTPUT2
