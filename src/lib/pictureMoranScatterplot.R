pictureMoranScatterplot <- function(dfPooled, df, year, shp_comm=NA, shp_sll=NA, 
                                    scatterPlot=TRUE,
                                    map=FALSE,
                                    title="", 
                                    typeCities="metro"){
  
  #Moran scatteplot
  #qq <- quantile(df$log.Pop, c(0.25, 0.50, 0.75))
  cc <- c(0.50, 0.9, 0.99)
  qq <- quantile(df$log.Pop, cc)

  df <- df %>% mutate(log.Pop.b=ifelse(log.Pop<qq[1], 0.1,
                                       ifelse(log.Pop<qq[2], 0.3,
                                              ifelse(log.Pop<qq[3], 0.5, 1))))

  df <- df %>% mutate(cols=ifelse(log.Pop<qq[1], "#FFCC66",
                                  ifelse(log.Pop<qq[2], "#FF9900",
                                         ifelse(log.Pop<qq[3], "#CC6633", "#660000"))))
  
  
  #moran scatter plot
  moran.sm_Lag1 <- sm.regression(dfPooled$x,dfPooled$wx,se=TRUE,display="none")

  if (scatterPlot==TRUE){
    dev.new()
  plot(df$x, df$wx, xlab= "Population density (log)", ylab = "spatially lagged (LLA) population density (log)", main=title,  cex=df$log.Pop.b, xlim=c(0, 10), ylim=c(0,10), bg=24,  lwd=.4,  pch=19, col=df$cols)
  points(df$x[df$cities==1], df$wx[df$cities==1], col="#FFFF00", pch=19, cex=1)
  points(df$x[df$cities==1], df$wx[df$cities==1], col="black", pch=21, cex=1)
  abline(h = mean(dfPooled$wx), lty = 2)
  abline(v = mean(dfPooled$x), lty = 2)
  lines(moran.sm_Lag1$eval.points,moran.sm_Lag1$estimate, col="blue")
  lines(moran.sm_Lag1$eval.points,moran.sm_Lag1$estimate+2*moran.sm_Lag1$se, lty = 2, col="blue")
  lines(moran.sm_Lag1$eval.points,moran.sm_Lag1$estimate-2*moran.sm_Lag1$se, lty = 2, col="blue")
  lines(c(0,15), c(0,15))
  grid()
  text(df$x[df$cities==1], df$wx[df$cities==1], labels=df$name.com[df$cities==1], pos=3, cex=0.8, col="black")
  legend("topleft", legend=c(paste(0, "<=population<", cc[1]*100, "%", sep=""), paste(cc[1]*100, "<=population<", cc[2]*100, "%", sep=""), paste(cc[2]*100, "<=population<", cc[3]*100, "%", sep=""), paste(cc[3]*100, "<=population<", 100, "%", sep="")), col=c("#FFCC66", "#FF9900", "#CC6633", "#660000"), pt.cex=c(0.2, 0.3, 0.5, 1), pch=19, cex=1)

  dev.copy2eps(file=paste(title, ".eps", sep=""))
  }
  
  
  if (map==TRUE){
    
    #Map
    shp_comm@data <- shp_comm@data %>% left_join(df, by=c("PRO_COM"="geo"))
    shp_comm@data <- shp_comm@data %>% mutate(cities=ifelse(cities==1, 1, NA))
    shp_comm@data <- shp_comm@data %>% mutate(citiesHighPop=ifelse(citiesHighPop==1, 1, NA))
    shp_comm_sf <- st_as_sf(shp_comm)
    
    if (typeCities=="metro"){
      cities <- shp_comm_sf %>% filter(cities==1)
      cities <- st_intersection(cities, st_union(shp_comm_sf))
      cities <- st_cast(cities,"POLYGON")
    }
    if (typeCities=="highPop"){
      cities <- shp_comm_sf %>% filter(citiesHighPop==1)
      cities <- st_intersection(cities, st_union(shp_comm_sf))
      cities <- st_cast(cities,"POLYGON")
    }
    if (typeCities=="fua"){
      cities <- shp_comm_sf %>% filter(citiesFUA==1)
      cities <- st_intersection(cities, st_union(shp_comm_sf))
      cities <- st_cast(cities,"POLYGON")
    }
    

    qq <- quantile(dfPooled$x, seq(0, 1, by=0.2))

    # tm_obj <- tm_shape(shp_comm) + tm_fill("x", title="Population density (log)", breaks=qq, palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5)

    tm_obj <- tm_shape(shp_comm) + tm_fill("x", title="Population density (log)", breaks=qq, palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(cities) +  tm_dots(size = 0.05, col="red")   + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5)
    
    dev.new()
      tm_obj
    tmap_save(tm = tm_obj,  filename = paste("map_popDensSLL_", year, ".pdf", sep=""))
    
    # qq <- quantile(dfPooled$wx, seq(0, 1, by=0.1))
    # tm_obj <- tm_shape(shp_comm) + tm_fill("wx", title="Spatially lagged (LLA) population denstiy (log)", breaks=qq, palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(shp_comm_cities) + tm_fill("cities", palette=c("#FFFF00"), colorNA = "transparent",legend.show=FALSE) + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5) #+ tm_layout(legend.show=TRUE)
    # 
    # dev.new()
    # tm_obj
    # tmap_save(tm = tm_obj,  filename = paste("map_popDensWSLL_", year, ".pdf", sep=""))
  }
  
  
}