
rm(list=ls())
library(lme4)
library(doBy)
library(MuMIn)
options(na.action = "na.fail")
library(lsmeans)
library(influence.ME)
library(multcomp)

#logit <- function(x){log(x/(1-x))} #link function logit.

#############################################################
#################1. Data preparattion########################
#############################################################

##############################
######## 1.1 Raw data ########
path<-"D:/Proyecto_AdAptA/Gene Flow Silene (Mallorca)/Silene-ciliata-gene-flow-master/Data" 
#Archivo original en https://github.com/javimorente/Silene-ciliata-gene-flow/tree/master/Data
setwd(path)

field=read.table("geneflow_field_experiment.txt", header=T, sep="\t")
str(field)

field<-field[,1:8]
key=read.table("geneflow_field_experiment_KEY.txt", header=T, sep="\t")
str(key)
str(field)
#fila with NA´s, I found it looking at the levels of each factor ("")
field <- field[-4758, ]  

#################################################
######## 1.2 Variable code (block & nail) ########

code=paste(key$blo,key$nail, sep = "_")
key<-cbind(key,code)

duplicated(key$code)
key[duplicated(key$code), ]
duplicates_KEY <- as.data.frame(key[duplicated(key$code), ])
#write.table(duplicates_KEY, file="duplicates_KEY.txt")
#key<-key[!(duplicated(key$code)| duplicated(key$code, fromLast=TRUE)) ,] #26 casos duplicados error toma datos clavo. Cuando se resuelva borrar esta línea de codigo 


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
naj.F <- naj[grep("pen", naj$father), ]
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
mor.F <- mor[grep("pen", mor$father), ]
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
rui.F <- rui[grep("pen", rui$father), ] 
               
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
    (agi[grep("pen", agi$father), ]), (agi[grep("zon", agi$father), ]) 
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
ses.F <- ses[grep("pen", ses$father), ] 
               
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
cam.F <- cam[grep("pen", cam$father), ] 
               
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


code=paste(field$blo,field$nail, sep = "_")
field=data.frame(field, code)
#write.table (field1, file = "field1.txt")

#field2=merge(field, key, by="code")
field2=merge(field, data.frame(code=key$code, mother=key$mother, father=key$father, seeds=key$seeds, treat= key$treat), by= "code", )
field2<-droplevels(field2)
#se pierde un caso, un NA

#########################################SE PIERDEN 5 CASOSO!!!!
#field.conditional=field$code %in% field2$code
#field1<-data.frame(field, conditional=field.conditional)
#field1False<-subset(field1, conditional <= FALSE)
#field1False
#key<-as.data.frame(key)
#key1False<-subset(key, code == "B16_5")
#no está en key.
#Los he encontrado en los estadillos. No se debieron de pasar esos códigos.
#modifico el original y lo subo también al GitHub
##########################################


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

ger<-rbind(t1,subset(t2, conditional=="FALSE"),subset(t3, conditional=="FALSE"),subset(t4, conditional=="FALSE") )  
#####Dataframe con todas las semillas germinadas y el censo (Time) en el que lo hicieron.

ger.sum=aggregate(class~code, data=ger, FUN=sum) #suma germinadas por clavo
colnames(ger.sum)=c("code","n.ger")
ger2= merge(key,ger.sum, by="code", sort = TRUE)
nrow(ger2) #967 clavos con germinación
length(which(!duplicated(ger$code)=="TRUE"))  #Compruebe que coincide con base de datos ger (code no duplicados es 967)
no.ger<-ger2$seeds-ger2$n.ger #calculo semillas no germinadas
p.ger=ger2$n.ger/ger2$seeds #calculo proporcion gernminadas
ger2=data.frame(ger2,no.ger,p.ger)
ger2<-droplevels(ger2)
nrow(ger2)

#ger2=subset(ger2, c(blo!="B3" | blo!="B8"))    #Eliminar bloque 3 y bloque 8
#ger3<-subset(ger2, subset= c(treat=="F1" | treat== "F2" | treat=="F3"))   #Seleccionar tratamientos F1, F2, F3

ger3<-subset(ger2, subset= c(treat!="F4" & treat!= "F5" )) # Extraer tratamientos F4 y F5
ger3<-droplevels(ger3) #Limpiar de la memoria de R los niveles que has eliminado con subset

###########Crear y  Unir unacolumna con el tamaño de la madre.

size<-read.table('Matriz _general_fenologia_2014.txt', header=T)
tail(size)
head(size)

size$poblacion
levels(size$poblacion)[1] = "agi"
levels(size$poblacion)[2] = "cam"
levels(size$poblacion)[3] = "mor"
levels(size$poblacion)[4] = "naj"
levels(size$poblacion)[5] = "neg"
levels(size$poblacion)[6] = "pen"
levels(size$poblacion)[7] = "rui"
levels(size$poblacion)[8] = "ses"
levels(size$poblacion)[9] = "zon"

size$code = paste(size$poblacion, size$madre, sep=" ")
head(size)

size_madre<-subset(size,select=c(code, poblacion, madre, sierra, size))
head(size_madre) 
str(size_madre)#720 datos
#quitamos NA's
size_madre_noNA <- na.omit(size_madre) #encontramos 20, se quedan en 700 datos
str(size_madre_noNA)
#quitamos duplicados haciendo la media entre ellos.
size_madre_cleen=aggregate(size~code, size_madre_noNA, FUN=mean) #quitamos 375, se queda en 325 datos
str(size_madre_cleen)
head(size_madre_cleen)
#guardamos la tabla.
write.table(size_madre_cleen, file='size_madre_cleened.txt', col.names = TRUE)
size_madre <- read.table("size_madre_cleened.txt", header = T, dec=",")

names(size_madre)[1]<-paste("mother")#renombramos igual a la columna para que se entiendan.
head(size_madre)
head(ger3)

ger4= merge(ger3,size_madre, by="mother", all.x = T)
head(ger4) 
head(ger3)#comparar con ger3

#######################################################################################################################
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

#par(mfrow=c(1,1))
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

ger.mod1= glmer(matrix ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
summary(ger.mod1)
anova(ger.mod1)
#install.packages('car')
car::Anova(ger.mod1) 

##Evaluar modelo complementary-log-log para "zero-inflated"
ger.mod1.log= glmer(matrix ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial(link=cloglog),data = ger3)
summary(ger.mod1.log)
anova(ger.mod1.log)
car::Anova(ger.mod1.log)     #Obtener P-valores
#modelos raros, todo significativo y Std. Error iguales.

##Modelos sin "zero inflated"
ger.mod1= glmer(matrix ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3) #waring message: no convergen
#hecho arriba
ger.mod2= glmer(matrix ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3) #waring message: no convergen
ger.mod3= glmer(matrix ~ treat +  (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
ger.mod4= glmer(matrix ~ moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)
ger.mod0= glmer(matrix ~ + (1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3)

anova(ger.mod0, ger.mod1, ger.mod2,ger.mod3, ger.mod4) 
#AIC, select models 2 and 4
#ojo, mas abajo se quitan outlayers
#Intentar primero optimizar aleatorios, decidir estructura de aleatorios y luego optimizar efectos fijos

##Check the prediction of the best model

a=predict(ger.mod2,type="link") # Predices el valor para cada caso despues de transformarlo con la función de vinculo (logit para binomial)
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

#Mejor usar la funcion   sem.model.fits del paquete "piecewiseSEM"
install.packages("parwiceSEM")
library(parwiceSEM)
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


####MuMIn

library(MuMIn)
options(na.action = "na.fail") #MuMIn

###Aic selection

ms<-dredge(ger.mod1,rank="AICc",extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))
# Generate a set of models with combinations (subsets) of terms in the global model
# extra= c("R^2") calculate a coefficient of determination based on the likelihood-ratio test (R_LR²).
plot(ms, xlab=c("dns","plt", "dns:plt")) # Gráfico útil para ver que variables salen significativas en cada modelo
delta2=subset(ms, subset = delta < 7) #Te devuelve los modelos con delta < 2
avg<-model.avg(delta2, subset = delta < 7)
summary(avg) #Model averaging based on an information criterion
avg<-model.avg(ger.mod2,ger.mod4)#2 modelos seleccionados
predict.avg<-MuMIn:::predict.averaging(avg, type="response") #type="link" si quieres los valores transformados
#juntamos al data set una columna con las predicciones del modelo (predict.avg)
ger3=cbind(ger3,predict.avg)

boxplot(ger3$predict.avg~ger3$treat)
boxplot(ger3$predict.avg~ger3$treat*ger3$moun)
#intentar hacer un barplot.... Barras de herror?


######Modelos eliminando outliers#####################
######################################################

#outliers are those observations that lie outside 1.5*IQR,

par(mfrow=c(1,2))
bp<-boxplot(n.ger~treat*moun, data=ger3, main="Germination", ylab= "Germination (n seedlings)", xlab="Treatment") 

 ger3out<-subset(ger3, moun!="bej" | n.ger!=8)
 ger3out<-subset(ger3out, moun!="gua" | n.ger!=10)
 
 matrix = cbind(ger3out$n.ger, ger3out$no.ger)   #Matriz exitos y fracasos para usarla de variable dependiente
 
 ger.mod1.out= glmer(matrix ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3out)
 
 
 ms<-dredge(ger.mod1.out,rank="AICc",extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]])) #select models mod2 y mod4 (mont & mont+treat)
 delta2=subset(ms, subset = delta < 2) #Te devuelve los modelos con delta < 2
summary(model.avg(ms, subset = delta < 2)) 

bp<-boxplot(n.ger~treat, data=ger3out, main="Germination", ylab= "Germination (n seedlings)", xlab="Treatment") 

###R2

 ger.mod2.out= glmer(matrix ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3out)
 ger.mod4.out= glmer(matrix ~moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother), family=binomial,data = ger3out)
  
piecewiseSEM::sem.model.fits (ger.mod2.out)
piecewiseSEM::sem.model.fits (ger.mod4.out)

################################3hacer modelos optimizando efectos aleatorios y despues optimizado los fijos.

##################################################################
###################### 3. SIZE ###################################
##################################################################

###################################
######## 3.1 Data managing ########

#quitamos NA's
field3<-subset(field3, code2!="B1_75_NA")
size=aggregate(size~code2*pos, field3, FUN=max) #Data with NA
size=aggregate(size~code3, field3, FUN=max) #hace lo mismo, ya que code3 ya contiene la posicion
#tamaño de cada plantula en base de datos "field3", limpia de duplicados (4662 plantulas en todos los tiempos).
#queremos obtener el tamaño maximo de cada plantula buscando en todos los tiempos (code 2 repetidos).
head(size) #mismo numero de obs que ger, tiene sentido, ger = numero de germinaciones total en todos los tiempos.
size=merge(size,ger[,c(-11,-16)], by="code3")
#unimos a ger para que tenga toda la información de cada caso, quitamos tamaño semillas en tiempo de su germinación
#size=subset(size, c(blo!="B3" & blo!="B8"))
size<-subset(size, subset= c(treat!="F4" & treat!= "F5" ))   #quitamos F4 y F5
size<-droplevels(size)
boxplot(size$size~size$treat, ylab="size (mm)")

which(size$size==35) #Caso 1805 es un outlier.  Una planta muy grande que será un adulto o una planta rebrotada

size<-size[-1806,] #eliminar u
boxplot(size$size~size$treat, ylab="size (mm)")


###################################
######## 3.2 Plotting Size ########

table.plot=summaryBy(size ~ treat, data = size,
                     FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )

mp <- barplot(table.plot[1:3,2], axes=FALSE, axisnames=FALSE, ylim=c(0, 6), col=gray.colors(3, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Size (mm)", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1, labels=c("F1", "F2","F3"), at = mp)
axis(2, at=seq(0 , 7, by=1))

segments(mp, table.plot[1:3,2], mp, table.plot[1:3,2] + table.plot[1:3,3], lwd=2)
segments(mp - 0.1, table.plot[1:3,2] + table.plot[1:3,3], mp + 0.1, table.plot[1:3,2] + table.plot[1:3,3], lwd=2)

###Plot SIZE  treatment vs pop
table.plot=summaryBy(size ~ moun*treat, data = size,
                     FUN = function(x) { c(m = mean(x), se= sd(x)/sqrt(length(x)),s = sd(x)) } )
par(mfrow=c(1,3))
mp <- barplot(table.plot[1:3,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 9), col=gray.colors(3, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Bejar", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1,labels=c("F1", "F2","F3"),at = mp)
axis(2, at=seq(0 , 9, by=1))

segments(mp, table.plot[1:3,3], mp, table.plot[1:3,3] + table.plot[1:3,4], lwd=1)
segments(mp - 0.1, table.plot[1:3,3] + table.plot[1:3,4], mp + 0.1, table.plot[1:3,3] + table.plot[1:3,4], lwd=1)

mp <- barplot(table.plot[4:6,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 9), col=gray.colors(5, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL), main="Gredos", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1,labels=c("F1", "F2","F3"),at = mp)
axis(2, at=seq(0 , 9, by=1))

segments(mp, table.plot[4:6,3], mp, table.plot[4:6,3] + table.plot[4:6,4], lwd=1)
segments(mp - 0.1, table.plot[4:6,3] + table.plot[4:6,4], mp + 0.1, table.plot[4:6,3] + table.plot[4:6,4], lwd=1)

mp <- barplot(table.plot[7:9,3], axes=FALSE, axisnames=FALSE, ylim=c(0, 9), col=gray.colors(5, start = 0.1, end =0.9, gamma = 2.2, alpha = NULL), main="Guadarrama", xlab="Treatment", ylab="Rosette Diameter (mm)")
axis(1,labels=c("F1", "F2","F3"),at = mp)
axis(2, at=seq(0 , 9, by=1))

segments(mp, table.plot[7:9,3], mp, table.plot[7:9,3] + table.plot[7:9,4], lwd=1)
segments(mp - 0.1, table.plot[7:9,3] + table.plot[7:9,4], mp + 0.1, table.plot[7:9,3] + table.plot[7:9,4], lwd=1)

#################################
######## 3.3 LMM - Size ########

#Mixed linear models (Gaussian Error)with log as link function

#Optimizamos los términos fijos con REML = F (Ajustamos por MAximum likelihood)
size.mod0= lmer(log(size) ~ 1 + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.mod1= lmer(log(size) ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.mod2= lmer(log(size) ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.mod3= lmer(log(size) ~ moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.mod4= lmer(log(size) ~ treat + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)

AICc.size<-AICc(size.mod0, size.mod1, size.mod2, size.mod3,size.mod4)
AICc.size[order(AICc.size$AICc),]

#Ajustamos los modelos finales mediante REML.
size.mod2= lmer(log(size) ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=T)
size.mod4= lmer(log(size) ~ treat + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=T)

r.size.mod2=r.squaredGLMM(size.mod2)
r.size.mod4=r.squaredGLMM(size.mod4)
par(mfrow=c(1,2))
qqnorm(resid(size.mod2))  #Un poco de heterocedasticidad
plot(size.mod2)           #Asunción normalidad ligeramente violada

##Vamos a probar que tal se comportan los modelos con link function = idendity (Sin transformar)

size.iden.mod0= lmer(size ~ 1 + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.iden.mod1= lmer(size ~ treat*moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.iden.mod2= lmer(size ~ treat+moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.iden.mod3= lmer(size ~ moun + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)
size.iden.mod4= lmer(size ~ treat + (1|moun:pop)+(1|moun:pop:blo)+ (1|moun:pop:blo:mother),data = size, REML=F)

AICc.size.iden<-AICc(size.iden.mod0, size.iden.mod1, size.iden.mod2, size.iden.mod3,size.iden.mod4)
AICc.size.iden[order(AICc.size.iden$AICc),]

par(mfrow=c(1,2))
qqnorm(resid(size.iden.mod2))  #Los residuos son peores nos quedamos con los modelos anteiores.
plot(size.iden.mod2)


############################################################
################# 4. Proportion of survival   ##############
############################################################

###################################
######## 4.1 Data managing ########

library(survival)
library(survminer)
library(KMsurv)
library(help=KMsurv)

###Germination at each time interval
#data.frame ger created by S_geneflow_fieldexperiment
#ger contains germinations (i.e. class==1)
#class -> 1=alive, 2=dead

head(ger)
 
ger.surv<- cbind(ger, binomial = rep(1, nrow(ger))     )
 
t1=subset(field3, time==1)    #Plants censored at T1
t2=subset(field3, time==2)     #Plants censored at T2
t3=subset(field3, time==3)     #Plants censored at T3
t4=subset(field3, time==4)     #Plants censored at T4
 
#mortalidad
 
ger.t1<-subset(ger.surv, time==1)
mort.t1.t2<- ger.t1[,1] %in% t2[,1] #True <- alive ; False <- dead at t2
mort.t1.t3<- ger.t1[,1] %in% t3[,1] #True <- alive ; False <- dead at t3 
mort.t1.t4<- ger.t1[,1] %in% t4[,1] #True <- alive ; False <- dead at t4
ger.t1<-cbind(ger.t1, mort.t1.t2,mort.t1.t3, mort.t1.t4) #2657  plants  germinated at T1 (alive)
mort.t1.2<-subset(ger.t1, mort.t1.t2==FALSE)          #1774 dead plants
mort.t1.3<-subset(ger.t1, mort.t1.t2==TRUE & mort.t1.t3==FALSE)    #490 dead plants
mort.t1.4<-subset(ger.t1, mort.t1.t2==TRUE & mort.t1.t3==TRUE & mort.t1.t4==FALSE)  #217 dead plants
mort<-rbind(mort.t1.2,mort.t1.3,mort.t1.4)
mort<-cbind(mort, binomial=rep(2,nrow(mort)))
surviv<- rbind(ger.t1, mort[,-17])
surviv.t1<-cbind(surviv,time2= rep(c(1, 15, 45, 60), c(nrow(ger.t1),nrow(mort.t1.2),nrow(mort.t1.3),nrow(mort.t1.4))) )

ger.t2<-subset(ger.surv, time==2)    #285  plant germinated at T2
 mort.t2.t3<- ger.t2[,1] %in% t3[,1] #True <- alive ; False <- dead at t3 
mort.t2.t4<- ger.t2[,1] %in% t4[,1] #True <- alive ; False <- dead at t4 
ger.t2<-cbind(ger.t2, mort.t2.t3, mort.t2.t4)     
mort.t2.3<-subset(ger.t2, mort.t2.t3==FALSE)          #204 dead plants   at T3
mort.t2.4<-subset(ger.t2, mort.t2.t3==TRUE & mort.t2.t4==FALSE)  #48 dead plants    at T4
mort<-rbind(mort.t2.3,mort.t2.4)
mort<-cbind(mort, binomial=rep(2,nrow(mort)))
surviv<- rbind(ger.t2, mort[,-17])
surviv.t2<-cbind(surviv,time2= rep(c(15, 45, 60), c(nrow(ger.t2),nrow(mort.t2.3),nrow(mort.t2.4))) )

ger.t3<-subset(ger.surv, time==3)
mort.t3.t4<- ger.t3[,1] %in% t4[,1] #True <- alive ; False <- dead at t4 
ger.t3<-cbind(ger.t3, mort.t3.t4)     #45 plants germinated at T3
mort.t2.4<-subset(ger.t3, mort.t3.t4==FALSE)  #33 dead plants at T4
mort<-cbind(mort.t2.4, binomial=rep(2,nrow(mort.t2.4)))
surviv<- rbind(ger.t3, mort[,-17])
surviv.t3<-cbind(surviv,time2= rep(c( 45, 60), c(nrow(ger.t3),nrow(mort.t2.4))) )

surviv<-rbind(surviv.t1[,c(-20,-19,-18)],surviv.t2[,c(-19,-18)],surviv.t3[,c(-18)])    
#5753 datos --> 2657 + 285+45 + 1774 + 490+ 217+ 204+48+33
#binomial--> 1 alive, 2 dead

surviv2<-subset(surviv, subset= c(treat!="F4" & treat!= "F5" )) # Extraer tratamientos F4 y F5
surviv2<-droplevels(surviv2) #Limpiar de la memoria de R los niveles que has eliminado con subset

###############################################################
####4.2 The Kaplan-Meier estimate of the survival function####

##### Interaccion montaña y tratamiento en un modelo #######

######### survfit() #########

#Use a Surv object as the response on the left of the ~ operator and, if desired, terms separated by + operators on the right. One of the terms may be a strata object. For a single survival curve the right hand side should be ~ 1.

fit.surv <- survfit(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat+surviv2$moun)     #Survfit crea las curvas de supervivencia

plot(fit.surv , main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival function")
summary(fit.surv )
 summary(fit.surv)$surv
#The median survival time is defined #to be the time t0.5 such that S(t0.5) = 0.5. Given an estimate of the survival function using Kaplan-Meier, this may be obtained graphically by drawing a horizontal line at 0.5.

print(survfit(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat+surviv2$moun), show.rmean=TRUE)

#########Tests for two or more samples###### 

  #Given two or more samples, is there a difference between the survival times? Setting up hypotheses for this problem, Then the test statistic is given by X2 = Z0ˆ-1Z,
  #which, under the null hypothesis, is distributed as a 2 distribution with n degrees of freedom
	#• H0 : h1(t) = h2(t) = · · · = hn(t) for all t.
	#• HA : hi(t0) 6= hj(t0) for at least one pair i, j and time t0.


  #The first argument is a survival object against a categorical covariate variable that is typically a variable designating which groups correspond to which survival times.
  #The second argument shown, rho, designates the weights. To give greater weight to the first part of the survival curves, use rho larger than 0. To give weight to the later part of the survival curves, use rho smaller than 0. 
  #The output of survdiff is relatively self-explanatory. A X2 statistic is computed along with a p-value.   #The test statistic is given by X2 = Z0ˆ-1Z, which, under the null hypothesis, is distributed as a X2 distribution with n degrees of freedom.

survdiff(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat+surviv2$moun, rho=0)
  
  #The test statistic is given by X2 = Z0ˆ-1Z, which, under the null hypothesis, is distributed as a X2 distribution with n degrees of freedom.

##########An accelerated failure-time (AFT) model##########

# The function survreg() is used for AFT modeling. 
# The argument dist has several options (’weibull’, ’exponential’, ’gaussian’, ’logistic’, ’lognormal’, and ’loglogistic’) and is the parametric model used.

aft.weib <- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat*surviv2$moun, dist="weibull")
summary(aft.weib)
anova(aft.weib)
aft.exp<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat*surviv2$moun, dist="exponential") 
 summary(aft.exp)
 anova(aft.exp)

aft.gau<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat*surviv2$moun, dist="gaussian") 
 summary(aft.gau)
 anova(aft.gau)


aft.log<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat*surviv2$moun, dist="logistic") 
 summary(aft.log)
 anova(aft.log)
 
 aft.logn<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat*surviv2$moun, dist="lognormal") 
 summary(aft.logn)
 anova(aft.logn)

  aft.logg<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat*surviv2$moun, dist="loglogistic") 
 summary(aft.logg)
 anova(aft.logg)

  res <- pairwise_survdiff(Surv(time2, binomial) ~ treat, data=surviv2,p.adjust.method = "holm")   #library survminer


#####  tratamiento en un modelo #######

######### survfit() #########

#Use a Surv object as the response on the left of the ~ operator and, if desired, terms separated by + operators on the right. One of the terms may be a strata object. For a single survival curve the right hand side should be ~ 1.

fit.surv <- survfit(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat)     #Survfit crea las curvas de supervivencia

plot(fit.surv , lty=c(1,2,3) ,main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival function")
summary(fit.surv )
 summary(fit.surv)$surv
#The median survival time is defined #to be the time t0.5 such that S(t0.5) = 0.5. Given an estimate of the survival function using Kaplan-Meier, this may be obtained graphically by drawing a horizontal line at 0.5.

print(survfit(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat), show.rmean=TRUE)

#########Tests for two or more samples###### 

  #Given two or more samples, is there a difference between the survival times? Setting up hypotheses for this problem, Then the test statistic is given by X2 = Z0ˆ-1Z,
  #which, under the null hypothesis, is distributed as a 2 distribution with n degrees of freedom
	#• H0 : h1(t) = h2(t) = · · · = hn(t) for all t.
	#• HA : hi(t0) 6= hj(t0) for at least one pair i, j and time t0.


  #The first argument is a survival object against a categorical covariate variable that is typically a variable designating which groups correspond to which survival times.
  #The second argument shown, rho, designates the weights. To give greater weight to the first part of the survival curves, use rho larger than 0. To give weight to the later part of the survival curves, use rho smaller than 0. 
  #The output of survdiff is relatively self-explanatory. A X2 statistic is computed along with a p-value.   #The test statistic is given by X2 = Z0ˆ-1Z, which, under the null hypothesis, is distributed as a X2 distribution with n degrees of freedom.

survdiff(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, rho=0)  
  
  #The test statistic is given by X2 = Z0ˆ-1Z, which, under the null hypothesis, is distributed as a X2 distribution with n degrees of freedom.

##############################################################
##########4.3 An accelerated failure-time (AFT) model##########

# The function survreg() is used for AFT modeling. 
# The argument dist has several options (’weibull’, ’exponential’, ’gaussian’, ’logistic’, ’lognormal’, and ’loglogistic’) and is the parametric model used.

aft.weib <- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, dist="weibull")
summary(aft.weib)
anova(aft.weib)
aft.exp<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, dist="exponential") 
 summary(aft.exp)
 anova(aft.exp)

aft.gau<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, dist="gaussian") 
 summary(aft.gau)
 anova(aft.gau)

aft.log<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, dist="logistic") 
 summary(aft.log)
 anova(aft.log)
 
 aft.logn<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, dist="lognormal") 
 summary(aft.logn)
 anova(aft.logn)

  aft.logg<- survreg(Surv(surviv2$time2, surviv2$binomial) ~ surviv2$treat, dist="loglogistic") 
 summary(aft.logg)
 anova(aft.logg)



