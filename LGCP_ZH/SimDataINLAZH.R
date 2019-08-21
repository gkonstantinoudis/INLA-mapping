

# Simulate data for the LGCP ZH

# Created 17.09.2018


############################################################################################

#

load("O:/PhD/Projects/LGCP for Childhood cancers/Data/Simulation/CAR_offset.RData")
spplot(mun2017, "x")

library(sp)
ts <- spsample(mun2017, n = 1000000, type = "clustered")

s <- sapply(slot(mun2017, 'polygons'), function(i) spsample(i, n=100, type='hexagonal', offset=c(0,0)))

temp.2 <- mun2017[complete.cases(mun2017$x),]


t <- spsample(slot(mun2017, 'polygons')[[1]], n=temp.2$x[1], type='random', offset=c(0,0))
plot(t)

for(i in 2:170){
  print(i)
  t <- rbind(t,
  spsample(slot(mun2017, 'polygons')[[i]], n=temp.2$x[i], type='random', offset=c(0,0))
  )
}

class(t)
t.dat <- as.data.frame(t)
plot(t.dat, cex = .2)


plot(mun2017)
points(spsample(slot(mun2017, 'polygons')[[i]], n=temp.2$x[i], type='random'), cex = .1, pch = 19)

# try smth else

data_url <- "http://www.web.statistik.zh.ch/cms_vis/2017_HP_Hexbins/data/hex_bev_ausl_alter_2016_round_3_4.json"
data_file <- "O:/srt_ZH.json"
# for some reason, I can't read from the url directly, though the tutorial
# says I can
download.file(data_url, data_file)

library(rjson)

json_data <- fromJSON(file=data_file)

class(json_data)
length(json_data)
class(json_data$features)
length(json_data$features)
json_data$features[[2500]]

plot(
c(8.781361, 8.780026, 8.777431, 8.776173, 8.777508, 8.780103, 8.781361),
c(47.315618, 47.314106, 47.314136, 47.315679, 47.317192, 47.317161, 47.315618)
)

library(sp)

tt = Polygon(t(matrix(unlist(json_data$features[[i]]$geometry$coordinates), ncol = 7)))
ps = Polygons(list(tt),i)

sp <- list(Polygons(list(Polygon(ZH)), 1), Polygons(list(Polygon(WH)), 2), 
      Polygons(list(Polygon(GL)), 3))

plot(sps_1)

id <- einw <- numeric(8510)
pol_list <- list()

for(i in 1:8510){

  print(i)
  id[i] <-  json_data$features[[i]]$properties$GRID_ID
  einw[i] <- json_data$features[[i]]$properties$N_EINW
  pol_list[[i]] <- Polygons(list(Polygon(t(matrix(unlist(json_data$features[[i]]$geometry$coordinates), ncol = 7)))), 
                                 i)
    
}

sl <- SpatialPolygons(pol_list)
plot(sl)


centroids <- coordinates(sl)
x <- centroids[,1]
y <- centroids[,2]

ex_1.7 <- SpatialPolygonsDataFrame(sl, data=data.frame(x=x, y=y, n = einw, row.names=row.names(sl)))

ex_1.7@data
spplot(ex_1.7, "n", col = "transparent")


o <- readShapePoly("C:/Users/konstantinoudis.CAMPUS/Downloads/2017_HP_Hexbins/shp/Hexgrid_100000_m2_Gemeinde")
plot(o)
sum(centroids %in% coordinates(o) )
rm(centroids)

# https://statistik.zh.ch/internet/justiz_inneres/statistik/de/aktuell/mitteilungen/2017/bevoelkerung2017.html
# try to merge them
ex_1.7@data$id <- id
sum(o$GRID_ID %in% ex_1.7@data$id)
o$n <- NA
o$n[o$GRID_ID %in% ex_1.7@data$id] <- ex_1.7@data$n[ex_1.7@data$id %in% o$GRID_ID[o$GRID_ID %in% ex_1.7@data$id]]
spplot(o, "n", col = "transparent")

popdensZH <- o
save(popdensZH, file = "O:/popdensZH.RData")


###################################################################################################################

# simulate data for the inla course

load("G:/Back up 03.04.2018/PhD/lectures/181030_INLA_ZH/LGCP_ZH/popdensZH.RData")

# I will sample approximately 200000, to do so I will create weights and multiply with 200000

popdensZH$weights <- round((popdensZH$n/sum(popdensZH$n, na.rm = T))*200000)
popdensZH$weights[is.na(popdensZH$weights)] <- 0

sum(popdensZH$weights)

# now sample

len <- nrow(popdensZH)

start <- which(popdensZH$weights!=0)[1]
t <- spsample(slot(popdensZH, 'polygons')[[start]], n=popdensZH$weights[start], type='random', offset=c(0,0))
plot(t)

for(i in 2:len){
  
  print(i)
  
  if(popdensZH$weights[i]!=0){
    t <- rbind(t,
               spsample(slot(popdensZH, 'polygons')[[i]], n=popdensZH$weights[i], type='random', offset=c(0,0))
    )
  }

}

t.dat <- as.data.frame(t)

# t.dat will be the truth

poppointproc <- t.dat
save(poppointproc, file = "G:/Back up 03.04.2018/PhD/lectures/181030_INLA_ZH/LGCP_ZH/poppointproc.RData")

# now I will create a hypothetical scenario of risk consistent with the paper
library(crs)
name <- c("Zürich", "Winterthur","Gossau")

centroid.x <- c(2681689, 2696839, 2700017)
centroid.y <- c(1247620, 1261850, 1240222)

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

plot(sl, add = T, border = "green")

crs(t) <- CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")

t.ov <- over(t, sl)

t.dat$weights <- 1
t.dat$weights[!is.na(t.ov)] <- 3
t.dat$risk <- (t.dat$weights/sum(t.dat$weights))*1500

set.seed(852)
q <- runif(n = nrow(t.dat), min = 0, max = 1)

t.dat$case <- 0
t.dat$case[t.dat$risk > q] <- 1

sum(t.dat$case == 1)
plot(t.dat[t.dat$case == 1, c(1:2)])

caseproces <- t.dat
save(caseproces, file = "G:/Back up 03.04.2018/PhD/lectures/181030_INLA_ZH/LGCP_ZH/caseproces.RData")


##################################################################################################################
