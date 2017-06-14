
#Es un script para hacer un modelo sencillo en el que el area de la planta depnde de la altura de la población (Variable "plot":bajo/alto) 
#y de la densidad de indivduos alrededor (variable continua, sería como tu peso de semilla)          

####en S_geneflow_fieldexperiment.R hay una parte con esta aplicación ya implementada!!
                
library(MuMIn)
options(na.action = "na.fail") #MuMIn
library(lme4)

#### Global model

#Primero definimos el modelo completo


area.lm= lm(log(area)~density*plot ,  data=area.sam)
anova(area.lm, test="Chi")
summary(area.lm)
par(mfrow=c(2,2))
plot(area.lm)

###Aic selection

#Luego hago AIC selection y model averaging con MuMIn. 


ms.area.lm<-dredge(area.lm,rank="AICc",extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))
# Te hace todos los modelos posibles y además te añade el estadistico R2
plot(ms.area.lm, xlab=c("dns","plt", "dns:plt")) # Gráfico útil para ver que variables salen significativas en cada modelo
delta2=subset(ms.area.lm, subset = delta < 2) #Te devuelve los modelos con delta < 2
summary(model.avg(ms.area.lm, subset = delta < 2)) 

#Predict with global model

#Expand.grid() crea data.frame con dos variables (plot y density). En plot se repetirá siempre "low" y density contiene la secuenciaque has denifido (0 a 6 que varía cada 0.1 )                 
plot.low=expand.grid(plot=("low"),density=seq(0,6, by=0.1))       #En tu caso density sería el peso de la semilla y haces le vector que te guste
plot.high=expand.grid(plot=("high"),density=seq(0,6, by=0.1))

area.lm.4= lm(log(area)~density+plot ,  data=area.sam)
area.lm.8= lm(log(area)~density*plot ,  data=area.sam)
avg=model.avg(area.lm.4,area.lm.8) #Vuelvo a ahcer el model averaging con los modelos que he seleccioando por AIC.

#Creo un objeto con los valores de Y para plot.high y plot.low utilizanod el objeto "avg" que he creado antes.

pred.high <- MuMIn:::predict.averaging(avg,newdata=plot.high,se.fit=TRUE,full = TRUE)  #usas full averaging y predices también con los SE.
pred.low <- MuMIn:::predict.averaging(avg, newdata=plot.low,se.fit=TRUE,full = TRUE)

#Draw de figure

#la variable X está en plot.high y la variable Y (predichos) está en el objeto pred.high. Sólo queda pintar la figura.

plot(plot.high$density,pred.high$fit, type="lines", xlim=c(0,10), ylim=c(-3.5,1.5),ylab="Crown area(log m^2)", xlab="Neighbourhood density (individuals)")
lines(plot.low$density,pred.low$fit, lty="dashed")
legend(7.5,1, inset=0.03,c("High", "Low" ), lty = c(1, 2), seg.len=2, bty="n", trace=TRUE, cex=0.8)
