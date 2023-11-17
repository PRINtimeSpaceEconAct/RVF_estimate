dfForPicture <- function(df, years, fua){
    yearsObs <- df$time
    tt <- yearsObs%in%years
    df <- df[tt,]
    df <- df %>% dplyr::select(geo, Comune, time, logPopDensity, WlogPopDensity, cities, pop)  %>% dplyr::rename(name.com=Comune, x=logPopDensity, wx=WlogPopDensity) %>% mutate(log.Pop=log(pop))
    
    #Change names cities (in engligh)
    df <- df %>% mutate(name.com=ifelse(name.com=="Milano", "Milan", 
                                        ifelse(name.com=="Firenze", "Florence",
                                               ifelse(name.com=="Torino", "Turin",
                                                      ifelse(name.com=="Venezia", "Venice",
                                                             ifelse(name.com=="Napoli", "Naples",
                                                                    ifelse(name.com=="Genova", "Genoa",
                                                                           ifelse(name.com=="Roma", "Rome",name.com))))))))
    
    #Cities bigger that 100000 inhabitants
    df <- df %>% mutate(citiesHighPop=if_else(pop>=20000, 1, 0))
    
    #FUA as eurostat
    df <- df %>% left_join(fua, by=c("name.com"="NAME"))
    df <- df %>% mutate(citiesFUA=if_else(!is.na(CODE), 1, 0))
    df <- df %>% dplyr::select(-c(CODE))
    return(df)
}