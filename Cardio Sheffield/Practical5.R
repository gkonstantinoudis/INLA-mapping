### Practical5

remove(list=ls())
setwd("/Users/akeskonstantinoudes/Desktop/INLA course/Cardio Sheffield/")

#Remember to set the right directory


###################################################
### code chunk number 4: Practical5.Rnw:52-56
###################################################
#--- Import the data ---#
stroke <- read.csv("Stroke.csv")
head(stroke)
plot(stroke$stroke_exp)
table(stroke$Townsend.class)
table(stroke$NOx.class)
library(maptools)
sheffield.gen <- readShapePoly("Sheffield.shp")
plot(sheffield.gen)
nrow(stroke)

###################################################
### code chunk number 5: Practical5.Rnw:60-61
###################################################
head(stroke)


###################################################
### code chunk number 6: Practical5.Rnw:64-65
###################################################
names(sheffield.gen)


###################################################
### code chunk number 7: Practical5.Rnw:74-77 (eval = FALSE)
###################################################
library(spdep)
#make an adjacency matrix
nb2INLA("Sheffield.graph",poly2nb(sheffield.gen))
sheffield.adj <- paste(getwd(),"/sheffield.graph",sep="")

###################################################
### code chunk number 9: Practical5.Rnw:87-91
###################################################
stroke<- stroke[match(sheffield.gen$SP_ID,stroke$SP_ID),]
#Check for the first 10 areas
stroke$SP_ID[1:10]
sheffield.gen$SP_ID[1:10]


###################################################
### code chunk number 10: Practical5.Rnw:108-109
###################################################
stroke$ID<- seq(1,1030)


###################################################
### code chunk number 11: Practical5.Rnw:114-116
###################################################
formula.DM <- y ~  f(ID,model="bym", graph=sheffield.adj, 
                     hyper=list(prec.spatial=list(param=c(1,1))))


###################################################
### code chunk number 12: Practical5.Rnw:120-124
###################################################
library(INLA)
model.DM <- inla(formula.DM,family="binomial",
          data=stroke, offset=Offset, Ntrials=pop, 
          control.compute=list(dic=TRUE))
plot(model.DM)
names(model.DM)
class(model.DM$marginals.fitted.values)
plot(model.DM$marginals.fitted.values$fitted.Predictor.0008, type = "l")
###################################################
### code chunk number 15: Practical5.Rnw:142-146
###################################################
#Random effect
summary.random <- model.DM$summary.random$ID
dim(summary.random)
head(round(summary.random,3))


###################################################
### code chunk number 16: Practical5.Rnw:150-153
###################################################
#Hyperparameters
summary.hyper<- model.DM$summary.hyper
round(summary.hyper,3)


###################################################
### code chunk number 17: Practical5.Rnw:157-159
###################################################
marg.hyper <- inla.hyperpar.sample(100000,model.DM)
colnames(marg.hyper)
class(marg.hyper)
dim(marg.hyper)


###################################################
### code chunk number 18: Practical5.Rnw:164-166
###################################################
perc.var.u1 <- mean(marg.hyper[,1] / (marg.hyper[,1]+marg.hyper[,2]))
perc.var.u1


###################################################
### code chunk number 19: >
###################################################
#Natural scale
exp.uv <- lapply(model.DM$marginals.random$ID[1:1030], function(x) 
                inla.emarginal(exp,x))

exp.u <- lapply(model.DM$marginals.random$ID[1031:2060], function(x) 
                inla.emarginal(exp,x))
#Categories
cutoff=c(0.5,0.8,0.95,1.05,1.2,18.5)
exp.uv.cat <- cut(unlist(exp.uv),breaks=cutoff,include.lowest=TRUE)
exp.u.cat <- cut(unlist(exp.u),breaks=cutoff,include.lowest=TRUE)


###################################################
### code chunk number 20: Practical5.Rnw:197-208
###################################################
data.exp.BYM <- data.frame(SP_ID=stroke$SP_ID, 
                           exp.u=exp.u.cat, exp.uv=exp.uv.cat)
row.names(data.exp.BYM) <- seq(1,1030)

head(data.exp.BYM)

#Merge exp(v), exp(v+u) and the sheffield shapefile 
sheffield <- sheffield.gen
data.sheffield <- attr(sheffield, "data")
attr(sheffield, "data") <- merge(data.sheffield,
                                 data.exp.BYM, by="SP_ID")


###################################################
### code chunk number 21: figMapBYM
###################################################
library(RColorBrewer)
spplot(obj=sheffield, zcol=c("exp.u","exp.uv"),
       col.regions= brewer.pal(5, "BrBG"), main="")


###################################################
### code chunk number 22: Practical5.Rnw:232-236
###################################################
formula.hier <- y ~  f(ID,model="iid")
model.hier <- inla(formula.hier,family="binomial",
            data=stroke, offset=offset, Ntrials=pop, 
            control.compute=list(dic=TRUE))


###################################################
### code chunk number 23: figHyper
###################################################
plot(inla.tmarginal(function(x) 1/x, 
                   model.hier$marginals.hyper[[1]]), main="", type = "l", xlim = c(0.25,0.6),
     ylim = c(0,12))
lines(inla.tmarginal(function(x) 1/x,
                   model.DM$marginals.hyper[[1]]),col="red")

#the residual variance has reduced when ucluded the random effect.

###################################################
### code chunk number 24: Practical5.Rnw:257-260
###################################################
compareDIC=cbind(model.hier$dic[1:4],model.DM$dic[1:4])
colnames(compareDIC) = c("Hier", "DM")
compareDIC
#sort of validates the above results, since the DIC is smaller for the
#spatial model (some of the variability of the model was explained by v)
#the difference between the models is sufficient and we actually see that we
#need the random effect.

###################################################
### code chunk number 25: Practical5.Rnw:279-286
###################################################
formula.reg <- y ~  f(ID,model="bym", graph=sheffield.adj, 
                      hyper=list(prec.spatial=list(param=c(1,1)))) + 
                      Townsend + NOx
model.reg <- inla(formula.reg,family="binomial",
              data=stroke, offset=offset, Ntrials=pop, 
              control.compute=list(dic=TRUE), 
              control.fixed=list(prec=list(mean=0, prec=0.5)))


###################################################
### code chunk number 26: Practical5.Rnw:289-290
###################################################
round(model.reg$summary.fixed,3)


###################################################
### code chunk number 27: Practical5.Rnw:293-302
###################################################
names(model.reg$marginals.fixed)
log.Townsend <- model.reg$marginals.fixed[[2]]
log.NOx <- model.reg$marginals.fixed[[3]]

prob.Townsend <- 1-inla.pmarginal(0,log.Townsend)
prob.NOx <- 1-inla.pmarginal(0,log.NOx)

prob.Townsend
prob.NOx


###################################################
### code chunk number 28: Practical5.Rnw:307-327
###################################################
#Natural scale
exp.u.reg <- lapply(model.reg$marginals.random$ID[1031:2060], 
                    function(x) inla.emarginal(exp,x))

#Categories
cutoff=c(0.5,0.8,0.95,1.05,1.2,18.5)
exp.u.reg.cat <- cut(unlist(exp.u.reg),breaks=cutoff,
                     include.lowest=TRUE)

data.exp.BYM <- data.frame(SP_ID=stroke$SP_ID, 
                          exp.u=exp.u.cat, exp.u.reg=exp.u.reg.cat)
row.names(data.exp.BYM) <- seq(1,1030)

head(data.exp.BYM)

#Merge exp(u) for the two models with the sheffield shapefile 
sheffield <- sheffield.gen
data.sheffield <- attr(sheffield, "data")
attr(sheffield, "data") <- merge(data.sheffield,data.exp.BYM, 
                                 by="SP_ID")


###################################################
### code chunk number 29: figMapEcreg
###################################################
spplot(obj=sheffield, zcol=c("exp.u","exp.u.reg"),
       col.regions= brewer.pal(5, "BrBG"), main="")

###################################################
### code chunk number 30: Practical5.Rnw:346-349
###################################################
marg.hyper2 <- inla.hyperpar.sample(100000,model.reg)
perc.var.u2 <- mean(marg.hyper2[,1] /(marg.hyper2[,1]+marg.hyper2[,2]))
perc.var.u2


