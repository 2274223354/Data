# PSEM --------------------------------------------------------------------
rm(list=ls())
library(nlme)
library(car)
library(lmerTest)
library(piecewiseSEM)
library(data.table)
library(pROC)
library(dplyr)
library(plyr)
library(doParallel)
##===============================================================================
## HN_FD:NaturalForests ---------
##===============================================================================
rm(list=ls())
####log and scale####
data_sem <-read.csv( "E:/GEB data/HN_data_sem_NaturalForests.csv" )

data_sem %>%
  dplyr::count(AgeGroup)

data_sem=data_sem[,c( "PLT_CN",  "INVYR.x","FLDTYPCD","COUNTYCD.x" ,"ECO_ID", "LON", "LAT","STDAGE" ,
                      "Carbon_Mg_ha" , "cv_dbh",  "VPD" , "Elevation" ,  "MAT" ,"MAP",                   
                      "pH_mean","bdod_mean" , "clay_mean", "nitrogen_mean" , "FDis" , "PC1.CWM",  "PC2.CWM" ,  "AgeGroup"           )]

# PSEM ------------------------------------------------------------------
####PSEM_220_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ #Elevation +   MAT + MAP+VPD +  clay_mean + pH_mean+
      
      FD_All ,
    
    # STDAGE +PC1,
    
    random = ~1 | ECO_ID/COUNTYCD.x,  
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(
    PC1 ~ pH_mean #+
    #VPD + STDAGE 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,  
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + #FD_All +  
              PC1 +
              STDAGE# + 
            # 
            # VPD + 
            
            #clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1
)


summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_HN_Multi-FD_Age2_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_HN_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age3####

data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  nlme::lme(FD_All ~ #Elevation +
              clay_mean #+ MAP+ VPD
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(PC1 ~ 
              # clay_mean + 
              MAT,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme::lme(Carbon_Mg_ha ~#Elevation +   MAT + MAP+VPD +  clay_mean + pH_mean+
              
              FD_All +PC1,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2)
  ,FD_All %~~%  PC1
)
summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_HN_Multi-FD_Age3_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_HN_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



#################################################################################
#################################################################################
####PSEM_630_Age1####

data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  nlme::lme(
    DBH_CV ~ STDAGE +MAT,
    
    random = ~1 | ECO_ID/COUNTYCD.x,  
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV +# FD_All +  PC1+PC2+
              
              STDAGE #+ 
            
            #  MAP
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_HN_Multi-FD_Age1_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_HN_Multi-FD_Age1_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_630_Age2####

data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))

mod_630_FD <- psem(
  nlme::lme(
    DBH_CV ~ STDAGE + PC1 #+ FD_All
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,  
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme::lme(
    PC1 ~ VPD #+clay_mean
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,  
    
    data = data_sem2
  ),
  nlme::lme(Carbon_Mg_ha ~ DBH_CV  + 
              
              STDAGE + #PC1+PC2+
              
              VPD 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_HN_Multi-FD_Age2_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_HN_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


#################################################################################
#################################################################################
####PSEM_620_Age1####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))
#write.csv(data_sem2,"E:/GEB data/620_data_HN_Multi-FD_Age1_NaturalForests_线.R")

mod_620_FD <- psem(
  
  
  nlme::lme(
    PC1 ~clay_mean
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 200)
  ),
  nlme::lme(Carbon_Mg_ha ~DBH_CV + #FD_All + 
              MAT+PC1# +#PC2
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
)

summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_HN_Multi-FD_Age1_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_HN_Multi-FD_Age1_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~   MAT +
      
      # PC2+
      
      STDAGE +FD_All ,#+PC1,
    
    random = ~1 | ECO_ID/COUNTYCD.x,  
    
    data = data_sem2
  ),
  nlme::lme(
    FD_All ~  MAT
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,  
    
    data = data_sem2
  ),
  nlme::lme(
    PC2 ~ MAT+
      VPD 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,  
    
    data = data_sem2
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All +  
              
              STDAGE 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,  
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC2
)

summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_HN_Multi-FD_Age2_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD




num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_HN_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


##===============================================================================
## American_FD:NaturalForests ---------
##===============================================================================
rm(list=ls())
####log and scale####
data_sem <-read.csv( "E:/GEB data/FIA_data_sem_NaturalForests.csv" )

data_sem=data_sem[,c( "PLT_CN",  "INVYR.x","FLDTYPCD","COUNTYCD.x" ,"ECO_ID", "LON", "LAT","STDAGE" ,
                      "Carbon_Mg_ha" , "cv_dbh",  "VPD" , "Elevation" ,  "MAT" ,"MAP",                   
                      "pH_mean","bdod_mean" , "clay_mean", "nitrogen_mean" , "FDis" , "PC1.CWM",  "PC2.CWM" ,  "AgeGroup"           )]



# PSEM ------------------------------------------------------------------
####PSEM_220_Age1####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
library(dplyr)



data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))

mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2
  ),
  
  nlme::lme(FD_All ~ 
              MAP ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2)
  
)
summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age1_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age1_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age1_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)
####PSEM_220_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ clay_mean + FD_All + PC1 + STDAGE ,#+ Elevation +MAT,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ STDAGE + MAT + #MAP + VPD + pH_mean + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  
  
  nlme::lme(PC1 ~ STDAGE +
              MAT + #MAP + 
              # VPD + 
              clay_mean,#+pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  nlme:: lme(PC2 ~ STDAGE + MAT + 
               #MAP + VPD + 
               clay_mean,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + PC2 + 
              STDAGE + #MAT + #MAP + #VPD + 
              #Elevation + 
              pH_mean +
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age2_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ clay_mean +FD_All +#pH_mean +  
      #VPD +#MAP +
      PC2+
      PC1 #+ STDAGE
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ STDAGE + MAT +# MAP + 
              #VPD + #Elevation +#pH_mean + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme::lme(PC1 ~ STDAGE +
              MAT +MAP + 
              # VPD+#pH_mean+
              
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE +MAP ++#pH_mean +Elevation
               
               MAT + 
               VPD + clay_mean,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + PC2 + 
              STDAGE + 
              MAT + #MAP +
              #  VPD + 
              #Elevation + 
              pH_mean +
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age3_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)

####PSEM_220_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All + PC2+#clay_mean +
      PC1 + #Elevation# +
      #PC2#+ #STDAGE 
      +MAT
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ #Elevation + #STDAGE + 
              MAT + #MAP + 
              #VPD + #pH_mean + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, #correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2
            #,
            #control = nlme::lmeControl(opt = "optim", maxIter = 50)
  ),
  
  
  
  nlme::lme(PC1 ~ STDAGE +
              MAT + #MAP + 
              #VPD + 
              clay_mean#+pH_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 50)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE + 
               MAT + 
               # MAP + 
               VPD + clay_mean ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 50)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ STDAGE + DBH_CV + FD_All + PC1 + PC2 + #MAT+ #VPD+
              clay_mean
            
            
            #+  MAP + VPD +  Elevation  +  + pH_mean +
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age4_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age4_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~  FD_All + clay_mean +
      PC1 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ #Elevation + #STDAGE + 
              #MAT + #MAP + 
              #VPD + #pH_mean + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, #correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2
            #,
            #control = nlme::lmeControl(opt = "optim", maxIter = 50)
  ),
  
  
  
  nlme::lme(PC1 ~ #STDAGE +
              clay_mean +MAT 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 50)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE + 
               # MAT + 
               # MAP + 
               VPD + clay_mean,#pH_mean +
             #Elevation,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 50)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + PC2 + 
              #STDAGE + 
              #MAT #+# MAP 
              # + VPD 
              #Elevation  +
              #pH_mean +
              clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age5_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age5_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


#################################################################################
#################################################################################
####PSEM_630_Age1####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ clay_mean +FD_All ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV
            # STDAGE + MAT + 
            #  MAP + 
            #  VPD + clay_mean+
            # pH_mean +
            #   Elevation
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  #,FD_All%~~% PC1#,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_630_FD) 

AIC(mod_630_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/630_psem_American_Multi-FD_Age1_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_American_Multi-FD_Age1_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



####PSEM_630_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ STDAGE + FD_All + 
      PC1 +# MAT+#MAP+
      clay_mean
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ STDAGE + #MAT + 
              # MAP +
              #Elevation+ VPD +
              pH_mean +  clay_mean
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  
  
  nlme::lme(PC1 ~ STDAGE +
              MAT ,#+ #MAP + 
            #VPD + 
            #clay_mean+
            # pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  # nlme:: lme(PC2 ~ STDAGE + 
  #             
  #              #MAP + 
  #              VPD #+ clay_mean+#pH_mean 
  #            #+Elevation
  #            ,
  #           random = ~1 | ECO_ID/COUNTYCD.x, 
  #            data = data_sem2),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + #PC1 + PC2 + 
              STDAGE + MAP + 
              # Elevation + 
              pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1#,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_American_Multi-FD_Age2_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD




num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_American_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_630_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]

data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ PC2+FD_All +MAP
    #  PC1 + 
    # STDAGE+clay_mean  +
    # MAT#+Elevation 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 200)
  ),
  
  nlme::lme(FD_All ~ #STDAGE + 
              MAT + 
              MAP
            #clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2
  ),
  
  nlme::lme(PC1 ~ 
              #  MAT+
              MAP+
              #VPD + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE + 
               clay_mean + MAT+
               VPD + pH_mean ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All +  PC2 + 
              STDAGE + #MAT + #MAP +
              VPD + 
              # Elevation + 
              pH_mean,# +
            #clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_American_Multi-FD_Age3_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD


num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_American_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_630_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ PC2+FD_All + 
      PC1 + 
      #  STDAGE+#+#VPD+ #clay_mean  +
      MAT#+Elevation 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ #STDAGE + 
              MAT 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(PC1 ~ STDAGE +#Elevation + 
              
              VPD ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ STDAGE +
               VPD,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + VPD+ 
              STDAGE ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  ,FD_All%~~% PC2 ,FD_All%~~% PC1, PC1%~~% PC2
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_American_Multi-FD_Age4_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



num_correlations <- 3
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_American_Multi-FD_Age4_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



####PSEM_630_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ PC2+FD_All + 
      PC1 + 
      STDAGE+#+#VPD+ #clay_mean  +
      MAT#+Elevation 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  # nlme::lme(FD_All ~ STDAGE + 
  #              MAT + 
  #             
  #             VPD
  #           ,
  #           random = ~1 | ECO_ID/COUNTYCD.x, 
  #           data = data_sem2#,
  #           #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  # ),
  # 
  # nlme::lme(PC1 ~ STDAGE +#Elevation + 
  #             
  #             VPD + 
  #             #clay_mean+
  #             pH_mean,
  #           random = ~1 | ECO_ID/COUNTYCD.x, 
  #           data = data_sem2#,
  #           #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  # ),
  
  nlme:: lme(PC2 ~ pH_mean +MAP+#STDAGE + MAT+
               VPD,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + 
              PC2 + 
              STDAGE+
              #+ MAT + MAP +
              #  VPD+
              # Elevation + 
              # pH_mean +
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  ,FD_All%~~% PC1, FD_All%~~% PC2, PC1%~~% PC2
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_American_Multi-FD_Age5_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD

num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_American_Multi-FD_Age5_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



#################################################################################
#################################################################################
####PSEM_620_Age1####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))

mod_620_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All +PC1#+Elevation
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV,
            
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
  #,PC1%~~% PC2
)


summary(mod_620_FD) 

AIC(mod_620_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/620_psem_American_Multi-FD_Age1_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD




num_correlations <- 0
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_American_Multi-FD_Age1_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



####PSEM_620_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All +PC1+VPD+
      STDAGE 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ MAT+
              # MAP +
              #Elevation+ VPD + 
              pH_mean  #+ clay_mean
            
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  
  
  nlme::lme(PC1 ~ STDAGE +
              # Elevation + 
              MAT + #MAP + 
              VPD + 
              clay_mean+
              pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  nlme:: lme(PC2 ~ STDAGE + 
               #  MAT + 
               #MAP + 
               VPD #+ clay_mean
             +pH_mean 
             # +Elevation
             ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + #PC2 + 
              clay_mean +   pH_mean + STDAGE 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_American_Multi-FD_Age2_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD


num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_American_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ PC2+FD_All +VPD+ 
      PC1 +  #pH_mean + 
      STDAGE+MAP
    #  +Elevation 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ 
              MAT + #MAP+
              
              
              #Elevation +
              pH_mean + 
              clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(PC1 ~
              clay_mean +   MAT + MAP+VPD+
              # Elevation + 
              pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE +
               MAT + MAP+VPD+
               clay_mean+
               pH_mean ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All +# PC1 + 
              PC2+
              STDAGE +  #VPD+ 
              pH_mean +
              MAT,#+ Elevation, #+
            
            
            
            #clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, 
  PC1%~~% PC2
)


summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_American_Multi-FD_Age3_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD


num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_American_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ #PC2+MAP+
      FD_All + 
      PC1 + 
      STDAGE
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ #STDAGE + 
              #MAT+
              # MAP + 
              VPD  #Elevation +
            #pH_mean #+ 
            #clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(PC1 ~ #STDAGE +
              #MAP+
              MAT + #pH_mean + 
              VPD, #+ 
            #clay_mean+
            # ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE +
               MAT + #MAP+ 
               VPD+
               
               pH_mean ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + #FD_All + #PC1 +
              PC2 + 
              STDAGE +  #MAP +
              # VPD+ 
              # Elevation + 
              pH_mean,# +
            #clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, 
  PC1%~~% PC2
)


summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_American_Multi-FD_Age4_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD

num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_American_Multi-FD_Age4_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Broadleaf Forest" ) &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ #PC2+
      FD_All + MAP+
      # clay_mean + 
      PC1 +  STDAGE +MAT
    
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ STDAGE +#MAT + #VPD +
              #pH_mean +  
              clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x,
            data = data_sem2
  ),
  
  nlme::lme(PC1 ~ STDAGE +VPD+
              MAT+ MAP +clay_mean+ pH_mean,
            
            
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ STDAGE + #MAT +  MAP +
               VPD + clay_mean+
               pH_mean ,
             random = ~1 | ECO_ID/COUNTYCD.x,
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV +#FD_All+
              # PC1 +
              VPD+ 
              #  PC2 + 
              STDAGE# #MAT + MAP ,
            ,
            
            #pH_mean +
            #clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  ,FD_All%~~% PC1, PC1%~~% PC2#,FD_All%~~% PC2
)

summary(mod_620_FD) 


# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_American_Multi-FD_Age5_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD


num_correlations <- 2 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_American_Multi-FD_Age5_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



##===============================================================================
## Japan_FD:NaturalForests ---------
##===============================================================================
rm(list=ls())
####log and scale####
data_sem <-read.csv( "E:/GEB data/Japan_data_sem_NaturalForests.csv" )
colnames(data_sem )
data_sem=data_sem[,c( "PLT_CN",  "INVYR.x","FLDTYPCD","COUNTYCD.x" ,"ECO_ID", "LON", "LAT","STDAGE" ,
                      "Carbon_Mg_ha" , "cv_dbh",  "VPD" , "Elevation" ,  "MAT" ,"MAP",                   
                      "pH_mean","bdod_mean" , "clay_mean", "nitrogen_mean" , "FDis" , "PC1.CWM",  "PC2.CWM" ,  "AgeGroup"           )]


# PSEM ------------------------------------------------------------------
####PSEM_220_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    PC2 ~  MAP
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ PC2+DBH_CV
            
            
            
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
)


summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age3_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  nlme:: lme(DBH_CV ~FD_All+#Elevation +  STDAGE +
               MAT+
               VPD,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme:: lme(PC2 ~#Elevation +  STDAGE +
               MAT+
               VPD,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ PC2+DBH_CV+
              MAT+
              
              # STDAGE + MAT + MAP +
              #   VPD+ 
              #  Elevation +
              pH_mean, #+
            # clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  
)




summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age4_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age4_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~  FD_All + STDAGE #+PC2
    
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV +  PC2 +FD_All
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age5_NaturalForests.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD



mod_bor = mod_220_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age5_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


#################################################################################
#################################################################################
####PSEM_630_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]

data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )

data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ STDAGE ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(
    Carbon_Mg_ha ~ DBH_CV + PC2 #+
    # STDAGE + MAT + MAP +VPD+ Elevation+pH_mean+clay_mean
    
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
)


summary(mod_630_FD)


summary(mod_630_FD) 

AIC(mod_630_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/630_psem_Japan_Multi-FD_Age3_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD




num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_Japan_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_630_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )


data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ #FD_All +VPD+
      PC1 + 
      PC2+
      pH_mean
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2
  ),
  
  nlme::lme(FD_All ~ 
              VPD
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + PC2 ,#+
            # STDAGE + MAT + MAP +VPD+ Elevation+pH_mean+clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2)
  
)


summary(mod_630_FD) 

# AIC(mod_630_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/630_psem_Japan_Multi-FD_Age4_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD


num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_Japan_Multi-FD_Age4_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



####PSEM_630_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c("Conifer-Broadleaf Mixed Forest" ) &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )


data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_630_FD <- psem(
  
  
  nlme:: lme(PC2 ~ clay_mean   ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + PC2 + #FD_All+
              #STDAGE + MAT + MAP + Elevation+pH_mean+
              clay_mean+
              VPD,
            # Elevation , 
            #pH_mean,# +
            #clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  
)


summary(mod_630_FD) 

AIC(mod_630_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/630_psem_Japan_Multi-FD_Age5_NaturalForests.R")
save(mod_630_FD,file=fn)

mod_bor=mod_630_FD



mod_bor = mod_630_FD


num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/630_psem_Newmode_Japan_Multi-FD_Age5_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)



#################################################################################
#################################################################################
####PSEM_620_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest")
                   &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All  
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV 
            + PC2 +
              STDAGE+
              
              clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~%  PC2
)


summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_Japan_Multi-FD_Age2_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD




num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_Japan_Multi-FD_Age2_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest")
                   &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ FD_All 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(
    FD_All ~ pH_mean
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    
    data = data_sem2
  ),
  
  nlme::lme(
    PC2 ~ MAP
    ,random = ~1 | ECO_ID/COUNTYCD.x,   data = data_sem2  ),
  
  
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV +PC2
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~%  PC2 
)


summary(mod_620_FD) 

# AIC(mod_620_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/620_psem_Japan_Multi-FD_Age3_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD


num_correlations <- 1
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_Japan_Multi-FD_Age3_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest")
                   &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ PC2+FD_All + 
      PC1 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ 
              
              MAP 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(PC1 ~ 
              MAT +MAP , 
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ 
               MAT,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(Carbon_Mg_ha ~ VPD  +
              DBH_CV + FD_All,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,#FD_All%~~%  PC2, 
  PC1%~~% PC2
)


summary(mod_620_FD) 

AIC(mod_620_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/620_psem_Japan_Multi-FD_Age4_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD



mod_bor = mod_620_FD


num_correlations <- 2 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_Japan_Multi-FD_Age4_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_620_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest")
                   &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh+1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_620_FD <- psem(
  
  nlme::lme(
    DBH_CV ~ PC2 +VPD+MAT
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ 
              
              MAP 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme::lme(Carbon_Mg_ha ~ #VPD  +
              DBH_CV + FD_All,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  
)

summary(mod_620_FD) 


AIC(mod_620_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/620_psem_Japan_Multi-FD_Age5_NaturalForests.R")
save(mod_620_FD,file=fn)

mod_bor=mod_620_FD


num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
  
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/620_psem_Newmode_Japan_Multi-FD_Age5_NaturalForests.R"

save(newmod_bor, file = fn_new_bor)


###################################################################################
###################################################################################
##===============================================================================
## HN_FD:Plantations ---------
##===============================================================================
rm(list=ls())
####log and scale####
data_sem <-read.csv( "E:/GEB data/HN_data_sem_Plantations.csv" )

data_sem %>%
  dplyr::count(AgeGroup)


data_sem=data_sem[,c( "PLT_CN",  "INVYR.x","FLDTYPCD","COUNTYCD.x" ,"ECO_ID", "LON", "LAT","STDAGE" ,
                      "Carbon_Mg_ha" , "cv_dbh",  "VPD" , "Elevation" ,  "MAT" ,"MAP",                   
                      "pH_mean","bdod_mean" , "clay_mean", "nitrogen_mean" , "FDis" , "PC1.CWM",  "PC2.CWM" ,  "AgeGroup"           )]

# PSEM ------------------------------------------------------------------
####PSEM_220_Age1####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ #Elevation +   MAT + MAP+VPD +  clay_mean + pH_mean+
      
      
      
      STDAGE +FD_All +PC1,
    
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2
  ),
  
  nlme::lme(
    PC1 ~ #pH_mean +
      VPD #+ STDAGE 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + PC1,
            
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2
  )
  
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_HN_Multi-FD_Age1_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_HN_Multi-FD_Age1_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))



mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ #Elevation +   MAT + MAP+VPD +  clay_mean + pH_mean+
      
      
      
      STDAGE +FD_All +PC1,
    
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2#,
    #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(
    PC1 ~ pH_mean +
      VPD 
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    
    data = data_sem2
  ),
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV +# FD_All +  
              PC1 + PC2+ MAT+
              STDAGE  
            
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2#,
            # control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,PC2%~~% PC1
)

summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_HN_Multi-FD_Age2_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_HN_Multi-FD_Age2_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  nlme::lme(DBH_CV ~#Elevation +   MAT + MAP+VPD +  clay_mean + 
              # pH_mean+
              FD_All ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2),
  
  
  
  nlme::lme(PC1 ~   pH_mean ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme::lme(Carbon_Mg_ha ~#Elevation +   MAT + MAP+VPD +  clay_mean + pH_mean+
              DBH_CV+  
              PC1,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2)
  ,FD_All %~~%  PC1
)
summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_HN_Multi-FD_Age3_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_HN_Multi-FD_Age3_Plantations.R"

save(newmod_bor, file = fn_new_bor)



#################################################################################
#################################################################################
##===============================================================================
## Jpan_FD:Plantations ---------
##===============================================================================
rm(list=ls())
####log and scale####
data_sem <-read.csv( "E:/GEB data/Japan_data_sem_Plantations.csv" )

data_sem %>%
  dplyr::count(AgeGroup)

data_sem=data_sem[,c( "PLT_CN",  "INVYR.x","FLDTYPCD","COUNTYCD.x" ,"ECO_ID", "LON", "LAT","STDAGE" ,
                      "Carbon_Mg_ha" , "cv_dbh",  "VPD" , "Elevation" ,  "MAT" ,"MAP",                   
                      "pH_mean","bdod_mean" , "clay_mean", "nitrogen_mean" , "FDis" , "PC1.CWM",  "PC2.CWM" ,  "AgeGroup"           )]



# PSEM ------------------------------------------------------------------
####PSEM_220_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All
    ,
    random = ~1 | ECO_ID/COUNTYCD.x,
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + 
              # FD_All + 
              PC1+
              #PC2 + 
              STDAGE + 
              #MAT + MAP + VPD + 
              #Elevation  +
              # pH_mean +
              clay_mean
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  # ,FD_All%~~%  PC2#,PC1%~~%  PC2
)


summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age2_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


num_correlations <- 0 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age2_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All +#pH_mean +  
      #VPD+
      PC2+
      STDAGE
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ STDAGE + 
              VPD ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  
  nlme:: lme(PC2 ~ #STDAGE + MAT + 
               MAP,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ FD_All + 
              DBH_CV + 
              # PC1 + 
              PC2 + 
              STDAGE 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~%  PC2
)


summary(mod_220_FD)


fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age3_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


#pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
mod_bor = mod_220_FD




num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age3_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))

mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~STDAGE + 
      FD_All +#MAT +  
      
      PC2
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~
              VPD ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme:: lme(PC2 ~ #STDAGE + 
               MAT + 
               VPD,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ #FD_All+
              DBH_CV + PC1 + PC2 + VPD +
              STDAGE 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~%  PC2,FD_All%~~%  PC1, PC1%~~% PC2
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age4_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age4_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age5####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 5) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))

mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~STDAGE + 
      
      FD_All #+#MAT +  
    
    #PC2
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme:: lme(PC2 ~ STDAGE + 
               MAT + 
               VPD,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ #FD_All+
              DBH_CV +# PC1 +
              PC2 #+ VPD +
            #STDAGE 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~%  PC2
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_Japan_Multi-FD_Age5_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_Japan_Multi-FD_Age5_Plantations.R"

save(newmod_bor, file = fn_new_bor)
#################################################################################
#################################################################################
##===============================================================================
## American_FD:Plantations ---------
##===============================================================================
rm(list=ls())
####log and scale####
data_sem <-read.csv( "E:/GEB data/FIA_data_sem_Plantations.csv" )

data_sem %>%
  dplyr::count(AgeGroup)
data_sem=data_sem[,c( "PLT_CN",  "INVYR.x","FLDTYPCD","COUNTYCD.x" ,"ECO_ID", "LON", "LAT","STDAGE" ,
                      "Carbon_Mg_ha" , "cv_dbh",  "VPD" , "Elevation" ,  "MAT" ,"MAP",                   
                      "pH_mean","bdod_mean" , "clay_mean", "nitrogen_mean" , "FDis" , "PC1.CWM",  "PC2.CWM" ,  "AgeGroup"           )]



# PSEM ------------------------------------------------------------------
####PSEM_220_Age1####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 1) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <-psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All + PC1,# + PC2,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2#,
    #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ STDAGE + VPD,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  nlme::lme(PC1 ~ STDAGE +
              #MAT + #MAP + 
              VPD + 
              clay_mean,#+pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2),
  
  # nlme:: lme(PC2 ~ STDAGE + #MAT + 
  #               VPD + clay_mean,
  #            random = ~1 | ECO_ID/COUNTYCD.x, 
  #            data = data_sem2),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All +#PC1 +# PC2 + 
              STDAGE,# +
            # MAT + MAP + VPD + 
            # Elevation + 
            # pH_mean +
            # clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  ,FD_All%~~% PC1#,FD_All%~~%  PC2, PC1%~~% PC2
)

summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age1_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age1_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age2####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 2) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ #pH_mean +
      FD_All +# PC1+
      # PC2 +
      STDAGE +clay_mean #+MAT
    +VPD,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme::lme(FD_All ~ STDAGE + #MAT + 
              #MAP + 
              VPD + 
              # pH_mean + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)),
  
  nlme::lme(PC1 ~  STDAGE +
              MAT + #MAP + 
              #VPD + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)),
  
  nlme:: lme(PC2 ~ STDAGE + MAT +# MAP+
               VPD + 
               clay_mean ,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2,
             control = nlme::lmeControl(opt = "optim", maxIter = 100)),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + #FD_All +
              PC1 + PC2 + 
              STDAGE + MAT + #MAP + #VPD + 
              
              #pH_mean +
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100))
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_220_FD) 

AIC(mod_220_FD, AIC.type = "dsep",
    aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age2_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age2_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age3####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 3) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~clay_mean + FD_All +#pH_mean +  
      #MAT +
      VPD+#MAP +
      PC2#+
    #PC1 #+ STDAGE#+MAT
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2,
    control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(FD_All ~ STDAGE +#MAT + 
              VPD+
              # pH_mean + 
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  
  nlme::lme(PC1 ~ STDAGE +
              MAT +#MAP + 
              #VPD+
              
              clay_mean+pH_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  nlme:: lme(PC2 ~ #STDAGE +
               MAT + 
               #MAP +
               VPD + clay_mean,
             random = ~1 | ECO_ID/COUNTYCD.x, 
             data = data_sem2#,
             #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV + FD_All + PC1 + PC2 + 
              STDAGE + 
              MAT +# MAP +
              # VPD + 
              # Elevation + 
              # pH_mean +
              clay_mean,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 100)
  )
  ,FD_All%~~% PC1,FD_All%~~%  PC2, PC1%~~% PC2
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age3_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD


#pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
mod_bor = mod_220_FD




num_correlations <- 3 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age3_Plantations.R"

save(newmod_bor, file = fn_new_bor)


####PSEM_220_Age4####
data_sem2=data_sem[data_sem$FLDTYPCD %in%  c( "Conifer Forest" ) &data_sem$AgeGroup %in%  c( 4) ,]
data_sem2 <- data_sem2 %>%
  mutate(
    DBH_CV               = scale(log(cv_dbh + 1)),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),
    PC1= scale( PC1.CWM ),
    PC2= scale( PC2.CWM ) , 
    FD_All=scale(FDis),
    STDAGE=scale(STDAGE), 
    Carbon_Mg_ha=log(Carbon_Mg_ha+1  )
    
  )
data_sem2 <- as.data.frame(lapply(data_sem2, function(x) {
  if(is.matrix(x)) return(as.numeric(x)) else return(x)
}))


mod_220_FD <- psem(
  
  
  nlme::lme(
    DBH_CV ~ FD_All +PC1+
      VPD
    ,
    random = ~1 | ECO_ID/COUNTYCD.x, 
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_sem2, control = nlme::lmeControl(opt = "optim", maxIter = 200)
  ),
  
  
  nlme::lme(FD_All ~ MAT 
            ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(PC1 ~ 
              MAT ,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            data = data_sem2#,
            #control = nlme::lmeControl(opt = "optim", maxIter = 100)
  ),
  
  
  nlme::lme(Carbon_Mg_ha ~ DBH_CV+VPD,
            random = ~1 | ECO_ID/COUNTYCD.x, 
            correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE, fixed = FALSE),
            data = data_sem2,
            control = nlme::lmeControl(opt = "optim", maxIter = 200)
  )
  ,FD_All%~~%  PC1
)


summary(mod_220_FD) 

# AIC(mod_220_FD, AIC.type = "dsep",
#     aicc = T )

fn<-paste("E:/GEB data/220_psem_American_Multi-FD_Age4_Plantations.R")
save(mod_220_FD,file=fn)

mod_bor=mod_220_FD

num_correlations <- 1 
model_range <- 1:(length(mod_bor) - 1 - num_correlations)


pz_bor <- lapply(mod_bor[model_range], function(x) { 
  
  p_table <- summary(x)$tTable
  p_col <- grep('P-value', colnames(p_table))
  
 
  which(p_table[, p_col] > 0.05) 
})


print(pz_bor)

newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "E:/GEB data/220_psem_Newmode_American_Multi-FD_Age4_Plantations.R"

save(newmod_bor, file = fn_new_bor)


#################################################################################
#################################################################################
