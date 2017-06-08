
rm(list=ls())
library(lme4)
library(doBy)
library(MuMIn)
library(lsmeans)
library(influence.ME)
library(multcomp)

#logit <- function(x){log(x/(1-x))} #link function logit.

#############################################################
#################1. Data preparattion########################
#############################################################

##############################
######## 1.1 Raw data ########
setwd("C:/Users/carlos/Documents/SCIENCE/PROYECTOS/ADAPTA/geneflow_experiment")
                            
field=read.table("geneflow_field_experiment.txt", header=T, sep="\t")
str(field)

field<-field[,1:8]
key=read.table("geneflow_field_experiment_KEY.txt", header=T, sep="\t")
str(key)
#################################################
######## 1.2 Variable code (block & nail) ########

code=paste(key$blo,key$nail, sep = "_")
key<-cbind(key,code)
key<-key[!(duplicated(key$code)| duplicated(key$code, fromLast=TRUE)) ,] #26 casos duplicados error toma datos clavo. Cuando se resuelva borrar esta l�nea de codigo 
           
code=paste(field$blo,field$nail, sep = "_")
field=data.frame(code,field)

##########################################
######### 1.3 Variable treatment ########

###Najarra

naj=subset(key, pop=="naj")

#F1
naj.F <- naj[grep("naj", naj$father), ]
vec=rep("F1", nrow(naj.F))
naj.F1=data.frame(naj.F,treat=vec)

#F2
naj.F <- naj[grep("mor", naj$father), ]
vec=rep("F2", nrow(naj.F))
naj.F2=data.frame(naj.F,treat=vec)

#F3
naj.F <- naj[grep("pe�", naj$father), ]
vec=rep("F3", nrow(naj.F))
naj.F3=data.frame(naj.F,treat=vec)

#F4
naj.F <- rbind ( 
    (naj[grep("rui", naj$father), ]) , (naj[grep("agi", naj$father), ]) 
               )
vec=rep("F4", nrow(naj.F))
naj.F4=data.frame(naj.F,treat=vec)

#F5
naj.F <- naj[grep("neg", naj$father), ]
               
vec=rep("F5", nrow(naj.F))
naj.F5=data.frame(naj.F,treat=vec)

naj2=rbind(naj.F1,naj.F2,naj.F3,naj.F4,naj.F5)

###Morrena

mor=subset(key, pop=="mor")

#F1
mor.F <- mor[grep("mor", mor$father), ]
vec=rep("F1", nrow(mor.F))
mor.F1=data.frame(mor.F,treat=vec)

#F2
mor.F <- mor[grep("naj", mor$father), ]
vec=rep("F2", nrow(mor.F))
mor.F2=data.frame(mor.F,treat=vec)

#F3
mor.F <- mor[grep("pe�", mor$father), ]
vec=rep("F3", nrow(mor.F))
mor.F3=data.frame(mor.F,treat=vec)

#F4
mor.F <- rbind ( 
    (mor[grep("rui", mor$father), ]) , (mor[grep("agi", mor$father), ]) 
               )
vec=rep("F4", nrow(mor.F))
mor.F4=data.frame(mor.F,treat=vec)

#F5
mor.F <-    rbind ( 
    (mor[grep("neg", mor$father), ]) , (mor[grep("zon", mor$father), ]) 
              )
vec=rep("F5", nrow(mor.F))
mor.F5=data.frame(mor.F,treat=vec)

mor2=rbind(mor.F1,mor.F2,mor.F3,mor.F4,mor.F5)

###ruinas

rui=subset(key, pop=="rui")

#F1
rui.F <- rui[grep("rui", rui$father), ]
vec=rep("F1", nrow(rui.F))
rui.F1=data.frame(rui.F,treat=vec)

#F2
rui.F <- rui[grep("agi", rui$father), ]
vec=rep("F2", nrow(rui.F))
rui.F2=data.frame(rui.F,treat=vec)

#F3
rui.F <- rui[grep("neg", rui$father), ]
vec=rep("F3", nrow(rui.F))
rui.F3=data.frame(rui.F,treat=vec)

#F4
rui.F <- rbind ( 
    (rui[grep("naj", rui$father), ]) , (rui[grep("mor", rui$father), ]) 
               )
vec=rep("F4", nrow(rui.F))
rui.F4=data.frame(rui.F,treat=vec)

#F5
rui.F <- rui[grep("pe�", rui$father), ] 
               
vec=rep("F5", nrow(rui.F))
rui.F5=data.frame(rui.F,treat=vec)

rui2=rbind(rui.F1,rui.F2,rui.F3,rui.F4,rui.F5)

###Aguilas

agi=subset(key, pop=="agi")

#F1
agi.F <- agi[grep("agi", agi$father), ]
vec=rep("F1", nrow(agi.F))
agi.F1=data.frame(agi.F,treat=vec)

#F2
agi.F <- agi[grep("rui", agi$father), ]
vec=rep("F2", nrow(agi.F))
agi.F2=data.frame(agi.F,treat=vec)

#F3
agi.F <- agi[grep("neg", agi$father), ]
vec=rep("F3", nrow(agi.F))
agi.F3=data.frame(agi.F,treat=vec)

#F4
agi.F <- rbind ( 
    (agi[grep("naj", agi$father), ]) , (agi[grep("mor", agi$father), ]) 
               )
vec=rep("F4", nrow(agi.F))
agi.F4=data.frame(agi.F,treat=vec)

#F5
agi.F <- rbind ( 
    (agi[grep("pe�", agi$father), ]), (agi[grep("zon", agi$father), ]) 
                )
               
vec=rep("F5", nrow(agi.F))
agi.F5=data.frame(agi.F,treat=vec)

agi2=rbind(agi.F1,agi.F2,agi.F3,agi.F4,agi.F5)

###sestil

ses=subset(key, pop=="ses")

#F1
ses.F <- ses[grep("ses", ses$father), ]
vec=rep("F1", nrow(ses.F))
ses.F1=data.frame(ses.F,treat=vec)

#F2
ses.F <- ses[grep("cam", ses$father), ]
vec=rep("F2", nrow(ses.F))
ses.F2=data.frame(ses.F,treat=vec)

#F3
ses.F <- ses[grep("zon", ses$father), ]
vec=rep("F3", nrow(ses.F))
ses.F3=data.frame(ses.F,treat=vec)

#F4
ses.F <- rbind ( 
    (ses[grep("naj", ses$father), ]) , (ses[grep("mor", ses$father), ]) 
               )
vec=rep("F4", nrow(ses.F))
ses.F4=data.frame(ses.F,treat=vec)

#F5
ses.F <- ses[grep("pe�", ses$father), ] 
               
vec=rep("F5", nrow(ses.F))
ses.F5=data.frame(ses.F,treat=vec)

ses2=rbind(ses.F1,ses.F2,ses.F3,ses.F4,ses.F5)

###Campanarios

cam=subset(key, pop=="cam")

#F1
cam.F <- cam[grep("cam", cam$father), ]
vec=rep("F1", nrow(cam.F))
cam.F1=data.frame(cam.F,treat=vec)

#F2
cam.F <- cam[grep("ses", cam$father), ]
vec=rep("F2", nrow(cam.F))
cam.F2=data.frame(cam.F,treat=vec)

#F3
cam.F <- cam[grep("zon", cam$father), ]
vec=rep("F3", nrow(cam.F))
cam.F3=data.frame(cam.F,treat=vec)

#F4
cam.F <- rbind ( 
    (cam[grep("naj", cam$father), ]) , (cam[grep("mor", cam$father), ]) 
               )
vec=rep("F4", nrow(cam.F))
cam.F4=data.frame(cam.F,treat=vec)

#F5
cam.F <- cam[grep("pe�", cam$father), ] 
               
vec=rep("F5", nrow(cam.F))
cam.F5=data.frame(cam.F,treat=vec)

cam2=rbind(cam.F1,cam.F2,cam.F3,cam.F4,cam.F5)


key=rbind(naj2,mor2,rui2,agi2, ses2, cam2)
moun= data.frame(
    ifelse(key$pop=="naj" | key$pop=="mor", "gua", 
              ifelse(key$pop=="ses" | key$pop=="cam", "gre", "bej")))
colnames(moun)="moun"

key=data.frame(key, moun)

#############################################################
#################   2. Germination   ########################
#############################################################

#####################################
########  2.1 Preparing data  #######
  
field2=merge(field, data.frame(code=key$code, mother=key$mother, father=key$father, seeds=key$seeds, treat= key$treat), by= "code", )

#key2=subset(ger,select=c(code,moun,pop,mother,father,seeds,treat) )
#key2=key2[!duplicated(key2), ]

###Nails whitout seedlings (no germination)

#seedlings=data.frame(field[,1])
#seedlings=seedlings[!duplicated(seedlings), ]
#nail=data.frame(code=key[,1])

#nail.zero=apply(nail, 1, function (x) x %in% seedlings) 
#length(which(nail.zero==FALSE)) #FALSE== nails without seedlings
#nail.zero=data.frame( code=nail, condition=nail.zero)
#nail.zero=subset(nail.zero, condition=="FALSE",code)
#zero=rep(0,nrow(nail.zero))

#nail.zero=data.frame(code=nail.zero, n.ger=zero)
#no.ger= merge(nail.zero, key, by= "code", sort=TRUE)

###Germination at each time interval
code2=paste(field2$code,field2$pos,sep="_") #bloque_clavo_posicion
code3=paste(field2$code,field2$pos,field2$time,sep="_") #Bloque_clavo_posicion_tiempo

field3=data.frame(code2,code3,field2)
field3<-field3[!duplicated(field3$code3),]

t1=subset(field3, time==1)
t2=subset(field3, time==2)
t3=subset(field3, time==3)
t4=subset(field3, time==4)

t1.conditional=t1[,1] %in% t2[,1] #TRUE == germinado/vivo; FALSE == germinado/muerto
t1=data.frame(t1, conditional=t1.conditional)

t2.conditional= t2[,1] %in% t1[,1] #TRUE == 0/0 germinado/vivo; FALSE ==no germiniado/germinado
t2=data.frame(t2, conditional=t2.conditional)

t3.1.conditional=t3[,1] %in% t1[,1]
t3.2.conditional=t3[,1] %in% t2[,1]
t3.conditional<-ifelse(t3.1.conditional=="FALSE" & t3.2.conditional=="FALSE", "FALSE","TRUE") #FALSE = Germination in T3 (no ger/no ger/germinado), TRUE== germinado en T1 o T2 y vivo T3 )
t3=data.frame(t3, conditional=t3.conditional)

t4.1.conditional=t4[,1] %in% t1[,1]
t4.2.conditional=t4[,1] %in% t2[,1]
t4.3.conditional=t4[,1] %in% t3[,1]
t4.conditional<-ifelse(t4.1.conditional=="FALSE" & t4.2.conditional=="FALSE" & t4.3.conditional=="FALSE","FALSE","TRUE") #FALSE = Germination in T4 (no ger/no ger/no ger/germinado), TRUE==Superv T4 (germinado en T1 o T2 o T3 y vivo T4)
t4=data.frame(t4, conditional=t4.conditional)

ger<-rbind(t1,subset(t2, conditional=="FALSE"),subset(t3, conditional=="FALSE"),subset(t4, conditional=="FALSE") )  #Dataframe con todas las semillas germinadas y el censo (Time)en el que lo hicieron.

ger.sum=aggregate(class~code, data=ger, FUN=sum)
colnames(ger.sum)=c("code","n.ger")
ger2= merge(key,ger.sum, by="code", sort = TRUE)
nrow(ger2) #967 clavos con germinaci�n
length(which(!duplicated(ger$code)=="TRUE"))  #Compruebe que coincide con base de datos ger (code no duplicados es 967)
no.ger<-ger2$seeds-ger2$n.ger
p.ger=ger2$n.ger/ger2$seeds
ger2=data.frame(ger2,no.ger,p.ger)

#ger2=subset(ger2, c(blo!="B3" | blo!="B8"))    #Eliminar bloque 3 y bloque 8

#ger3<-subset(ger2, subset= c(treat=="F1" | treat== "F2" | treat=="F3"))   #Seleccionar tratamientos F1, F2, F3

ger3<-subset(ger2, subset= c(treat!="F4" & treat!= "F5" )) # Extraer tratamientos F4 y F5
ger3<-droplevels(ger3) #Limpiar de la memoria de R los niveles que has eliminado con subset


#####Distribution of data within populations

library(doBy)
sample.size=summaryBy(n.ger ~ blo, data = ger3,
  FUN = function(x) { c(n = length(which(x>0)))} )
sample.size2=summaryBy(n.ger ~ pop*treat, data = ger3,
  FUN = function(x) { c(n = length(which(x>0)))} )
sum(sample.size$p.ger.n)

sample.mean=summaryBy(n.ger ~ pop*treat, data = ger3,
  FUN = mean )


####################################
#####  2.2 Plot Germination ########

#####Boxplot

par(mfrow=c(1,2))
boxplot(n.ger~treat*moun, data=ger3, main="Germination", ylab= " Germination (n seedlings) ", xlab="Treatment") 
boxplot(p.ger~treat, data=ger3, main="Germination", ylab= " Proportion of germination", xlab="Treatment") 

##### barplot Germination treatment

table.plot=summaryBy(p.ger ~ treat, data = ger3,
   FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )

par(mfrow=c(1,1))
mp <- barplot(table.plot[1:3,2], axes=FALSE, axisnames=FALSE,ylim=c(0,0.4),  col=gray.colors(3, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Germination", xlab="Treatment", ylab="Proportion of germination")
axis(1, labels=c("F1", "F2","F3"), at = mp)   #Eje horizontal
axis(2, at=seq(0 , 0.4, by=0.1))      #Eje vertical

segments(mp, table.plot[1:3,2], mp, table.plot[1:3,2] + table.plot[1:3,3], lwd=2) #Add vertical line
segments(mp - 0.1, table.plot[1:3,2] + table.plot[1:3,3], mp + 0.1, table.plot[1:3,2] + table.plot[1:3,3], lwd=2)    #Add  horizontal line

###### barplot treatment vs mountain 

table.plot=summaryBy(p.ger ~ moun*treat, data = ger3,
  FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )

par(mfrow=c(1,3))
mp <- barplot(table.plot[1:3,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 0.5), col=gray.colors(3, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Bejar", xlab="Treatment", ylab="Proportion of germination")
axis(1,labels=c("F1", "F2","F3"),at = mp)
axis(2, at=seq(0 , 0.5, by=0.1))

segments(mp, table.plot[1:3,3], mp, table.plot[1:3,3] + table.plot[1:3,4], lwd=1)
segments(mp - 0.1, table.plot[1:3,3] + table.plot[1:3,4], mp + 0.1, table.plot[1:3,3] + table.plot[1:3,4], lwd=1)

mp <- barplot(table.plot[4:6,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 0.5), col=gray.colors(3, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Gredos", xlab="Treatment", ylab="Proportion of germination")
axis(1,labels=c("F1", "F2","F3"),at = mp)
axis(2, at=seq(0 , 0.5, by=0.1))

segments(mp, table.plot[4:6,3], mp, table.plot[4:6,3] + table.plot[4:6,4], lwd=1)
segments(mp - 0.1, table.plot[4:6,3] + table.plot[4:6,4], mp + 0.1, table.plot[4:6,3] + table.plot[4:6,4], lwd=1)

mp <- barplot(table.plot[7:9,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 0.5), col=gray.colors(3, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Guadarrama", xlab="Treatment", ylab="Proportion of germination")
axis(1,labels=c("F1", "F2","F3"),at = mp)
axis(2, at=seq(0 , 0.5, by=0.1))

segments(mp, table.plot[7:9,3], mp, table.plot[7:9,3] + table.plot[7:9,4], lwd=1)
segments(mp - 0.1, table.plot[7:9,3] + table.plot[7:9,4], mp + 0.1, table.plot[7:9,3] + table.plot[7:9,4], lwd=1)

#####################################
#####  2.3 GLMM - germination #######

###GLMMs

matrix = cbind(ger3$n.ger, ger3$no.ger)   #Matriz exitos y fracasos para usarla de variable dependiente
library(lme4)

ger.mod1= glmer(matrix ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
summary(ger.mod1)
anova(ger.mod1)
car::Anova(ger.mod1)     #Obtener P-valores

ger.mod2= glmer(matrix ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
ger.mod3= glmer(matrix ~ treat +  (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
ger.mod4= glmer(matrix ~ moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
ger.mod0= glmer(matrix ~ + (1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)

anova(ger.mod0, ger.mod1, ger.mod2,ger.mod3, ger.mod4)

##Check the prediction of the best model

a=predict(ger.mod2,type="link") # Predices el valor para cada caso despues de transformarlo con la funci�n de vinculo (logit para binomial)
a=predict(ger.mod1,type="response")      # Predices el valor para cada caso sin transfomar

par(mfrow=c(1,2))
boxplot(a~ger3$treat, main="Predict")   #Respuesta predicha
boxplot(ger3$p.ger~ger3$treat,main="Response")         #Respuesta bruta

###Residuals

par(mfrow=c(2,2), mar=c(5,4,1,2))
Res=residuals(ger.mod1)
Fit= fitted(ger.mod1)

plot(Res~Fit,xlim=c(0,0.8),xlab="Fitted values", ylab="residuals", main="Res vs Fit")
abline(h=0)
qqnorm(Res)
qqline(Res)
hist(Res, xlab="residuals", main="ger.mod1")
plot(ger3$treat, Res, xlab="density", ylab="residuals")
abline(h=0)


###R2

#MEjor usar la funcion   sem.model.fits del paquete "piecewiseSEM"

piecewiseSEM::sem.model.fits (ger.mod2)
piecewiseSEM::sem.model.fits (ger.mod4)


###Post hoc lsmeans

#Post hoc comparisons of least square means. 
#post hoc Z-test of least squares means
#Post hoc least square means tests used to test for intergroup differences

library(lsmeans)
means.germination.n=lsmeans(ger.mod1, pairwise~treat, adjust="none")       #no ajusta pvalor
means.germination.t=lsmeans(ger.mod1, pairwise~treat, adjust="bonferroni")  #muy restrictivo
means.germination.t=lsmeans(ger.mod1, pairwise~treat, adjust="fdr")       #poco restrictivo


means.germination.t.int=lsmeans(ger.mod1, pairwise~treat | moun, adjust="fdr")
means.germination.n.int=lsmeans(ger.mod1, pairwise~treat | moun, adjust="none")


lsmip.germination.t.int=lsmip(ger.mod1, moun~treat,type="response", pch = c(0,1,2),
lty = 1:3, col = 1)

plot( (lsmeans(ger.mod1, ~ treat | moun)  )  , type = "response", level = .95,
xlab = "Predicted probability of survival")

plot( (lsmeans(ger.mod1, ~ treat )  )  , type = "response", level = .95,
xlab = "Predicted probability of survival")


##Outliers
library(influence.ME)
infl <- influence(ger.mod1,group="pop")
#Calculate Cook's distance:
cooks.distance(infl)
#Plot Cook's distance:
plot(infl, which = "cook")



####MuMIn

library(MuMIn)
options(na.action = "na.fail") #MuMIn

###Aic selection

ms<-dredge(ger.mod1,rank="AICc",extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))
# Te hace todos los modelos posibles y adem�s te a�ade el estadistico R2
plot(ms, xlab=c("dns","plt", "dns:plt")) # Gr�fico �til para ver que variables salen significativas en cada modelo
delta2=subset(ms, subset = delta < 2.2) #Te devuelve los modelos con delta < 2
avg<-model.avg(ms, subset = delta < 2.2)
summary(avg) #Model averaging based on an information criterion
avg<-model.avg(ger.mod2,ger.mod4)
predict.avg<-MuMIn:::predict.averaging(avg, type="response") #type="link" si quieres los valores transformados

boxplot(predict.avg~ger3$trea)

############################################################
################# 3. Proportion of survival   ##############
############################################################

###################################
######## 3.1 Data managing ########

###Estimar n�mero de plantulas vivas y muertas en T2 respecto a T1 por cada clavo

vivo=ifelse(t1$conditional=="TRUE",1,0)
t1=data.frame(t1,vivo)
ger.0= aggregate(class~code,t1,sum)
vivo=aggregate(vivo~code, t1, sum)
death=ger.0$class-vivo$vivo
survival.prop=vivo$vivo/ger.0$class
head(survival.prop)

survival=data.frame(ger.0, vivo=vivo[,2], death,survival.prop)

survival=merge(survival,key, by="code", all=FALSE)
survival=subset(survival, c(blo!="B3" & blo!="B8" ))

###################################
######## 3.2 Plot Survival ########

boxplot(survival$survival.prop~survival$treat, main="Survival", ylab= " Survival proportion", xlab="Treatment")

table.plot=summaryBy(survival.prop ~ treat, data = survival,
  FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )

par(mfrow=c(1,1))
mp <- barplot(table.plot[1:5,2], axes=FALSE, axisnames=FALSE, ylim=c(0, 0.5), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Survival", xlab="Treatment", ylab="Survival proportion")
axis(1, labels=c("F1", "F2","F3", "F4","F5"), at = mp)
axis(2, at=seq(0 , 0.5, by=0.1))

segments(mp, table.plot[1:5,2], mp, table.plot[1:5,2] + table.plot[1:5,3], lwd=2)
segments(mp - 0.1, table.plot[1:5,2] + table.plot[1:5,3], mp + 0.1, table.plot[1:5,2] + table.plot[1:5,3], lwd=2)

###Plot treatment vs pop
table.plot=summaryBy(survival.prop ~ moun*treat, data = survival,
  FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )

par(mfrow=c(1,3))
mp <- barplot(table.plot[1:5,3],axes=FALSE, axisnames=FALSE, ylim=c(0, 0.6), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Bejar", xlab="Treatment", ylab="Survival proportion")
axis(1,labels=c("F1", "F2","F3", "F4","F5"),at = mp)
axis(2, at=seq(0 , 0.6, by=0.2))

segments(mp, table.plot[1:5,3], mp, table.plot[1:5,3] + table.plot[1:5,4], lwd=1)
segments(mp - 0.1, table.plot[1:5,3] + table.plot[1:5,4], mp + 0.1, table.plot[1:5,3] + table.plot[1:5,4], lwd=1)

mp <- barplot(table.plot[6:10,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 0.6), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Gredos", xlab="Treatment", ylab="Survival proportion")
axis(1,labels=c("F1", "F2","F3", "F4","F5"),at = mp)
axis(2, at=seq(0 , 0.6, by=0.2))

segments(mp, table.plot[6:10,3], mp, table.plot[6:10,3] + table.plot[6:10,4], lwd=1)
segments(mp - 0.1, table.plot[6:10,3] + table.plot[6:10,4], mp + 0.1, table.plot[6:10,3] + table.plot[6:10,4], lwd=1)

mp <- barplot(table.plot[11:15,3], axes=FALSE, axisnames=FALSE, ylim=c(0,0.6), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Guadarrama", xlab="Treatment", ylab="Proportion of seedlings death")
axis(1,labels=c("F1", "F2","F3", "F4","F5"),at = mp)
axis(2, at=seq(0 , 0.6, by=0.2))

segments(mp, table.plot[11:15,3], mp, table.plot[11:15,3] + table.plot[11:15,4], lwd=1)
segments(mp - 0.1, table.plot[11:15,3] + table.plot[11:15,4], mp + 0.1, table.plot[11:15,3] + table.plot[11:15,4], lwd=1)

######################################
######## 3.3 GLMM -  Survival ########

##Models
matrix.surv = cbind(survival$vivo, survival$death)

surv.mod1= glmer(matrix.surv ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = survival)
surv.mod2= glmer(matrix.surv ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = survival)
surv.mod3= glmer(matrix.surv ~ treat +  (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = survival)
surv.mod4= glmer(matrix.surv ~ moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = survival)
surv.mod0= glmer(matrix.surv ~ + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = survival)

#AIC
anova(surv.mod0, surv.mod1, surv.mod2,surv.mod3, surv.mod4)

##R2
r.surv.mod0=r.squaredGLMM(surv.mod0) #library(MuMIn)
r.surv.mod1=r.squaredGLMM(surv.mod1)  
r.surv.mod2=r.squaredGLMM(surv.mod2)
r.surv.mod3=r.squaredGLMM(surv.mod3)
r.surv.mod4=r.squaredGLMM(surv.mod4)

##Post hoc least square means tests

means.survival.t=lsmeans(surv.mod2, pairwise~treat, adjust="tukey")
means.survival.n=lsmeans(surv.mod2, pairwise~treat, adjust="none")

a=predict(surv.mod3,type="response")
boxplot(a~survival$treat)


lsmip.survival.t.int=lsmip(surv.mod2, moun~treat,type="response", pch = c(0,1,2),
lty = 1:3, col = 1)

plot( (lsmeans(surv.mod2, ~ treat | moun)  )  , type = "response", level = .95,
xlab = "Predicted probability of survival")

plot( (lsmeans(surv.mod2, ~ treat )  )  , type = "response", level = .95,
xlab = "Predicted probability of survival")

##Post hoc Dunnet

set.seed(20140123)
dunnet.surv <- glht(surv.mod2, linfct=mcp(treat="Dunnett"))
summary(dunnet.surv)

##################################################################
###################### 4. SIZE ###################################
##################################################################

###################################
######## 4.1 Data managing ########

size=aggregate(size~code*pos, ger, max)
code2=paste(size$code,size$pos,sep="_")
str(code2)
size=data.frame(code2,size)
head(size)
size=merge(size,key, by="code", all=FALSE)
size=subset(size, c(blo!="B3" & blo!="B8"))

boxplot(size$size~size$treat, ylab="size (mm)")

###################################
######## 4.2 Plotting Size ########

table.plot=summaryBy(size ~ treat, data = size,
  FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )

par(mfrow=c(1,1))
mp <- barplot(table.plot[1:5,2], axes=FALSE, axisnames=FALSE, ylim=c(0, 7), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Size (mm)", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1, labels=c("F1", "F2","F3", "F4","F5"), at = mp)
axis(2, at=seq(0 , 7, by=1))

segments(mp, table.plot[1:5,2], mp, table.plot[1:5,2] + table.plot[1:5,3], lwd=2)
segments(mp - 0.1, table.plot[1:5,2] + table.plot[1:5,3], mp + 0.1, table.plot[1:5,2] + table.plot[1:5,3], lwd=2)

###Plot SIZE  treatment vs pop
table.plot=summaryBy(size ~ moun*treat, data = size,
  FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )
par(mfrow=c(1,3))
mp <- barplot(table.plot[1:5,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 9), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Bejar", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1,labels=c("F1", "F2","F3", "F4","F5"),at = mp)
axis(2, at=seq(0 , 9, by=1))

segments(mp, table.plot[1:5,3], mp, table.plot[1:5,3] + table.plot[1:5,4], lwd=1)
segments(mp - 0.1, table.plot[1:5,3] + table.plot[1:5,4], mp + 0.1, table.plot[1:5,3] + table.plot[1:5,4], lwd=1)


mp <- barplot(table.plot[6:10,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 9), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Gredos", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1,labels=c("F1", "F2","F3", "F4","F5"),at = mp)
axis(2, at=seq(0 , 9, by=1))

segments(mp, table.plot[6:10,3], mp, table.plot[6:10,3] + table.plot[6:10,4], lwd=1)
segments(mp - 0.1, table.plot[6:10,3] + table.plot[6:10,4], mp + 0.1, table.plot[6:10,3] + table.plot[6:10,4], lwd=1)

mp <- barplot(table.plot[11:15,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 9), col=gray.colors(5, start = 0.1, end = 1, gamma = 2.2, alpha = NULL), main="Guadarrama", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1,labels=c("F1", "F2","F3", "F4","F5"),at = mp)
axis(2, at=seq(0 , 9, by=1))

segments(mp, table.plot[11:15,3], mp, table.plot[11:15,3] + table.plot[11:15,4], lwd=1)
segments(mp - 0.1, table.plot[11:15,3] + table.plot[11:15,4], mp + 0.1, table.plot[11:15,3] + table.plot[11:15,4], lwd=1)

#################################
######## 4.3 GLMM - Size ########

size.mod1= lmer(log(size) ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size)
size.mod2= lmer(log(size) ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size)
size.mod3= lmer(log(size) ~ treat +  (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size)
size.mod4= lmer(log(size) ~ moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size)
size.mod0=lmer(log(size) ~ + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size)


###AIC

anova(size.mod0, size.mod1, size.mod2,size.mod3, size.mod4)

###R2
r.size.mod0=r.squaredGLMM(size.mod0)
r.size.mod1=r.squaredGLMM(size.mod1)  
r.size.mod2=r.squaredGLMM(size.mod2)
r.size.mod3=r.squaredGLMM(size.mod3)
r.size.mod4=r.squaredGLMM(size.mod4)

###Posthoc lsmeans

means.size.t=lsmeans(size.mod2, pairwise~treat, adjust="fdr")
means.size.n=lsmeans(size.mod2, pairwise~treat, adjust="none")

####Posthoc dunnet


set.seed(20140123)
dunnet.size <- glht(size.mod2, linfct=mcp(treat="Dunnett"))
summary(dunnet.size)