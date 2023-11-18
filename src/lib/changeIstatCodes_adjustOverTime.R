changeIstatCodes_adjustOverTime <- function(df){
  
  #Change istat code
  df <- df %>% 
    mutate(geo=ifelse((nchar(geo)==1), paste("00000", geo, sep=""),  
                           ifelse((nchar(geo)==4), paste("00", geo, sep=""),  
                                  ifelse((nchar(geo)==5), paste("0", geo, sep=""), 
                                         ifelse((nchar(geo)==6), geo, NA )))))
  
  df <- df %>% arrange(geo)
  years <- sort(unique(df$time))
  
  #Some municipalities changed code (province) during the period
  check.codeMun <- df %>% dplyr::group_by(Comune) %>% dplyr::summarise(Ncodes=length(unique(geo))) %>% filter(Ncodes>1)
  #Some municipalities have same name but different codes because in different regions!
  check.yearCodeMun <- df %>% dplyr::group_by(Comune) %>% dplyr::summarise(Nyears=length(time)) %>% filter(Nyears>length(years))
  check.codeMun <- check.codeMun %>% filter(Comune%in%check.yearCodeMun$Comune!=TRUE)
  
  #Replace the code of all years with the code of the last year available (2019)
  for (i in 1:nrow(check.codeMun)){
    iii <- df$Comune%in%check.codeMun$Comune[i] 
    lastYear <- which.max(df$time[iii])
    last.code <- df$geo[iii][lastYear]
    df$geo[iii] <- rep(last.code, sum(iii))
  }
  
  return(df)
  
}




