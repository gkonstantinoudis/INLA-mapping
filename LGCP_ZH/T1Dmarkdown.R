## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----results='hide', message=FALSE, warning=FALSE------------------------

library(INLA)
library(raster)
library(maptools)
library(spatstat)
library(fields)
library(raster)
library(latticeExtra)

# load the datafiles
# cases
load("caseproces.RData")
# population
load("poppointproc.RData")
# shapefile
load("ZHshapefile.RData")

## ---- fig.align = "center"-----------------------------------------------
centroid.x <- c(681689, 696839, 700017)
centroid.y <- c(247620, 261850, 240222)

centroids <- data.frame(centroid.x = centroid.x, centroid.y = centroid.y)

radius <- 5000

n = 1000
pts = seq(0, 2 * pi, length.out = n)
ZH = cbind(centroid.x[1] + radius * sin(pts), centroid.y[1] + radius * cos(pts))
WH <- cbind(centroid.x[2] + radius * sin(pts), centroid.y[2] + radius * cos(pts))
GL <- cbind(centroid.x[3] + radius * sin(pts), centroid.y[3] + radius * cos(pts))
sl = SpatialPolygons(list(Polygons(list(Polygon(ZH)), 1), Polygons(list(Polygon(WH)), 2), 
                          Polygons(list(Polygon(GL)), 3))) 
crs(sl) <- CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")

plot(mun2017)
plot(sl, border = "red", add=T, lwd = 2)
points(caseproces, cex = .2, pch = 19)


## ----fig.height = 6, fig.width = 5, fig.align = "center", message=F, warning=F----

cantonZH <- aggregate(mun2017)
plot(cantonZH)
W <- as.owin(cantonZH)
boundary <-  inla.mesh.segment(cbind(W$bdry[[1]]$x, W$bdry[[1]]$y))

mesh <- inla.mesh.create.helper(boundary = boundary, 
                                cutoff=500, max.edge=c(1400, 10000), offset=c(.1,10000))

plot(mesh)


## ----fig.align = "center"------------------------------------------------
# load the function
source("VoronoiFun.R")
VoronoiTes <- inla.mesh.dual(mesh)
crs(VoronoiTes) <- cantonZH@proj4string
plot(VoronoiTes)

## ----fig.align = "center"------------------------------------------------
# First I will define the polygons inside and outside the canton. We will use this to assign NA population density.
temp.spp <- SpatialPoints(coordinates(VoronoiTes), cantonZH@proj4string)
ov_1 <- over(temp.spp, cantonZH)

newdata <- data.frame(ID = 1:length(VoronoiTes@polygons))
row.names(newdata) <- newdata$ID
VoronoiTes <- SpatialPolygonsDataFrame(VoronoiTes, newdata)

VoronoiTes$InOut <- "In"
VoronoiTes$InOut[is.na(ov_1)] <- "Out"

# Now count population per polygon
 
spat.points.temp <- SpatialPoints(poppointproc, proj4string = cantonZH@proj4string)
over.tes <- over(spat.points.temp, VoronoiTes)
over.tes$x <-  spat.points.temp@coords[,1]
over.tes$y <-  spat.points.temp@coords[,2]

# add 1 as counts to be aggregated
  
over.tes$count <- 1
temp.ag <- aggregate(over.tes$count, by = list(ID = over.tes$ID), sum)
  
temp.merge <- merge(VoronoiTes, temp.ag, by.x= "ID", by.y = "ID")
# everything flagged as out stays NA, everything in becomes 0
temp.merge$x[temp.merge$InOut %in% "In" & is.na(temp.merge$x)] <- 0

# and divide with the area to calculate the population density
temp.merge$dens <- temp.merge$x/unlist(lapply(temp.merge@polygons, function(X) slot(X, "area")))
spplot(temp.merge, "dens", col = "transparent", main = "Population density ZH")


## ------------------------------------------------------------------------

spde <- inla.spde2.pcmatern(mesh=mesh, alpha=2, 
                            prior.range=c(60000,0.5), # Pr(phi < 60000) = 0.5
                            prior.sigma=c(1,0.01) # Pr(sigma > 1) = 0.01
)

## ------------------------------------------------------------------------

# the case vector. It is NA on the mesh nodes outside the domain, 0 on the mesh nodes inside the domain, and 1 for the case locations
temp.y <- rep(NA, times = (mesh$n))
temp.y[(VoronoiTes$InOut %in% "In")] <- 0
y.pp <- c(temp.y, rep(1, times = nrow(caseproces)))

# the e.pp is a vector of the areas of the Voronoi polygons, and it is 0 on the cases (because they are points)
e.pp <- c(diag(spde$param.inla$M0), rep(0,nrow(caseproces)))

# projector matrix, required when the locations are not vertices
lmat <- inla.spde.make.A(mesh = mesh, loc = as.matrix(caseproces))
imat <- Diagonal(mesh$n, rep(1, mesh$n))

# the entire projector matrix is
A.pp <- rBind(imat, lmat)

# For the offset I will adjust for area and population
ForOffset <- c(temp.merge$dens, rep(0, times = nrow(caseproces)))

# create a stack
stk.pp <- inla.stack(data=list(y=y.pp, e=e.pp*ForOffset), A=list(A.pp, 1), 
                     tag='pp', effects=c(list(i=1:mesh$n), list(intercept=rep(1, length(y.pp)))))


## ---- message=F, warning=F-----------------------------------------------

lgcp_mod <- inla(y ~ -1 + intercept + f(i, model = spde), 
               family='poisson', data = inla.stack.data(stk.pp),
               control.predictor = list(A = inla.stack.A(stk.pp), link = 1),
               E = inla.stack.data(stk.pp)$e, verbose=F, control.compute=list(dic=TRUE,cpo=TRUE, config = TRUE), 
               control.inla=list(strategy="simplified.laplace", int.strategy="eb"))


## ----width = 2, fig.align = "center"-------------------------------------

local.plot.field = function(field, mesh, ...){
  stopifnot(length(field) == mesh$n)
  proj = inla.mesh.projector(mesh, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), ...)
}

local.plot.field(lgcp_mod$summary.fitted.values$mean[1:mesh$n], mesh)


## ----fig.align = "center"------------------------------------------------

# create grid

grdpts <- makegrid(cantonZH, cellsize = c(1000,1000))
spgrd <- SpatialPoints(grdpts, proj4string = CRS(proj4string(cantonZH)))
spgrdWithin <- SpatialPixels(spgrd[cantonZH,])

# and project on the centroids
pgrid0 <- inla.mesh.projector(mesh, loc = coordinates(spgrdWithin))
# the median random effect
prd0.m <- inla.mesh.project(pgrid0, lgcp_mod$summary.random$i$`0.5quant`)

# covert to polygon
pol <- as(spgrdWithin, "SpatialPolygons")
pol_df <- SpatialPolygonsDataFrame(pol, data.frame(median.ranef = prd0.m, row.names = names(pol)))

spplot(pol_df, "median.ranef", col = "transparent", main = "Median Random Effects")


## ----fig.align = "center"------------------------------------------------

marfitted <- lgcp_mod$marginals.fitted.values[1:mesh$n]
threshold <- 0.0075
exceed.prob <- lapply(X= marfitted, FUN = function(x) inla.pmarginal(marginal = x, threshold))
exceed.prob <- 1 - unlist(exceed.prob)

# and we can also project this on the 1km grid as before

exceed.prob <- inla.mesh.project(pgrid0, exceed.prob)
pol_df$exceed.prob <- exceed.prob
spplot(pol_df, "exceed.prob", col = "transparent", main = "Exceedance Probability")

## ----fig.align = "center"------------------------------------------------

temp.ex <-  pol_df[pol_df$exceed.prob >=0.95,]
temp.ex <- aggregate(temp.ex)

spplot(pol_df, "exceed.prob", col = "transparent", main = "Exceedance Probability", col.regions = colorRampPalette(c('white', "seagreen"))(16)) + layer(sp.polygons(temp.ex, col="red", lwd = 2)) +  layer(sp.polygons(sl, col="yellow", lwd = 2))

## ----fig.align = "center", height = 3------------------------------------
par(mfrow = c(1,2), mai = c(.8,.8,.2,.1))
prior.median.range <- 60000
tmp = inla.tmarginal(function(x) 1/x, lgcp_mod$marginals.hyperpar$`Range for i`) 
lambda = -log(.5)/(1/prior.median.range)
plot(tmp, type = "l", xlab = "inverse range nominal", ylab = "Density", xlim = c(0,0.0002), 
     ylim = c(0, lambda))
xvals = seq(0, 0.0004, length.out=1000)
lines(xvals, lambda*exp(-lambda*xvals), col = "blue")
grid()

prior.median.sd <- 1
tmp = inla.tmarginal(function(x) x, lgcp_mod$marginals.hyperpar$`Stdev for i`) 
plot(tmp, type = "l", xlab = expression(sigma[u]), ylab = "Density", xlim = c(0, 1))
xvals = seq(0, 1, length.out=1000)
lambda = -log(.01)/prior.median.sd; lines(xvals, lambda*exp(-lambda*xvals), col = "blue")
grid()


