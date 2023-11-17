decompositionOfVariance <- function(data, var.name, group.name){
    require(dplyr)

    df = data %>% dplyr::select(value=!!var.name,group=!!group.name)
    
    vars <- df %>% dplyr::summarize(variance=var(value))
    vars <- vars*(nrow(df)-1)
    
    meanGroups <- df %>% group_by(group) %>% dplyr::summarise(meanGroups=mean(value))  %>% as.data.frame() 
    df <- df %>% dplyr::left_join(meanGroups, by="group")
    
    df <- df %>% dplyr::mutate(dist2MeanGroups=(value-meanGroups)^2)
    
    varGroups <- df %>% dplyr::group_by(group) %>% dplyr::summarise(varGroups=sum(dist2MeanGroups)) %>% as.data.frame() 
    varGroupsWithin <- sum(varGroups$varGroups)
    
    nGroups <- df %>% dplyr::group_by(group) %>% dplyr::summarise(n=length((dist2MeanGroups)))
    
    meanGroups <- meanGroups %>% dplyr::left_join(nGroups, by="group")    
    meanGroups <- meanGroups %>% dplyr::mutate(meanValue=mean(df$value))
    meanGroups <- meanGroups %>% dplyr::mutate(dist2Mean=n*(meanGroups-meanValue)^2)
    
    varGroupsBetween <- sum(meanGroups$dist2Mean)
    
    varWithinPerc = as.numeric(varGroupsWithin/vars*100)
    varBetweenPerc = as.numeric(varGroupsBetween/vars*100)
    varTot = var(df$value)
    
    list(varWithinPerc=varWithinPerc, varBetweenPerc=varBetweenPerc, varTot=varTot)
    
    }













