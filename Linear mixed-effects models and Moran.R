##===============================================================================
## American_NaturalForests: LMM and Moran's I ---------
##===============================================================================
rm(list=ls())
####220NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/FIA_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_NaturalForests_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_220_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_220_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_NaturalForests_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_220_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_220_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
####630NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/FIA_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer-Broadleaf Mixed Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 12)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_NaturalForests_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_630_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_630_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        #correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        # 将 correlation 部分修改为：
        correlation = nlme::corExp(value = 0.1, form = ~ LON + LAT, nugget = TRUE, fixed = FALSE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_NaturalForests_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_630_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_630_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
####620NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/FIA_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}


#########FD: Moran's I####


age_groups <- c( "1", "2", "3", "4", "5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 45) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 6)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_NaturalForests_620_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_620_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_620_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", 
            "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 45) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_NaturalForests_620_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_620_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_620_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
##===============================================================================
## Japan_NaturalForests: LMM and Moran's I ---------
##===============================================================================
rm(list=ls())
####220NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/Japan_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}
#########FD: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_Japan_NaturalForests_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_220_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_220_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_Japan_NaturalForests_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_220_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_220_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
####630NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/Japan_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer-Broadleaf Mixed Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####

age_groups <- c("1","2","3", "4", "5")

TRAITS <- c(
  "Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
  "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass"
)


library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()

output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data)
  
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c(
    "Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation",
    "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x"
  )
  
  for (tr in TRAITS) {
    
    target_trait <- paste0(tr, ".FDis")  
    
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    cat(sprintf(
      "\n>>> 正在全量计算: AgeGroup %s | 目标性状: %s (样本量: %d) <<<\n",
      ag, target_trait, nrow(temp_data)
    ))
    
    # 样本量阈值
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0(
      "Carbon_Mg_ha ~ ",
      target_trait,
      " + cv_dbh + MAT + MAP + VPD + pH_mean + clay_mean"
    )
    
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      
      # =========================================================
      # A. 先尝试带空间相关结构
      # =========================================================
      m <- tryCatch({
        
        nlme::lme(
          fixed = fixed_formula,
          random = ~1 | ECO_ID/COUNTYCD.x,
          correlation = nlme::corExp(
            form = ~ LON + LAT,
            nugget = TRUE
          ),
          data = temp_data,
          control = nlme::lmeControl(
            opt = "optim",
            maxIter = 100,
            msMaxIter = 100
          )
        )
        
      }, error = function(e1) {
        
        cat("  [警告] 空间相关结构失败，自动移除 correlation 重新运行\n")
        cat("  [原始错误]: ", e1$message, "\n")
        
        # =====================================================
        # B. 如果失败，则删除 correlation 后重新运行
        # =====================================================
        nlme::lme(
          fixed = fixed_formula,
          random = ~1 | ECO_ID/COUNTYCD.x,
          data = temp_data,
          control = nlme::lmeControl(
            opt = "optim",
            maxIter = 100,
            msMaxIter = 100
          )
        )
      })
      
      # --- B. 提取系数和 R2 ---
      at <- as.data.frame(summary(m)$tTable)
      
      r2_val <- MuMIn::r.squaredGLMM(m)
      
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor <- rownames(at)
      at$N_obs <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(
        as.matrix(coords_clean),
        k = 8
      )
      
      nb <- spdep::knn2nb(knn)
      
      lw <- spdep::nb2listw(
        nb,
        style = "W",
        zero.policy = TRUE
      )
      
      moran_res <- spdep::moran.test(
        res_clean,
        lw
      )
      
      # --- D. ---
      img_file <- paste0(
        figure_path,
        "Variogram_Japan_NaturalForests_630_",
        ag,
        "_",
        tr,
        ".png"
      )
      
      png(
        img_file,
        width = 800,
        height = 600,
        res = 120
      )
      
      vg <- nlme::Variogram(
        m,
        form = ~ LON + LAT,
        resType = "normalized"
      )
      
      print(
        plot(
          vg,
          main = paste0(
            "Variogram: AgeGroup ",
            ag,
            " | Trait: ",
            tr
          )
        )
      )
      
      dev.off()
      
      # --- E. 存入 Moran's I 简表 ---
      trait_p <- at[
        grep(target_trait, rownames(at)),
        "p-value"
      ][1]
      
      cv_p <- at[
        grep("cv_dbh", rownames(at)),
        "p-value"
      ][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag,
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      
      
      final_results <<- rbind(
        final_results,
        data.frame(
          AgeGroup = ag,
          Trait = tr,
          Status = paste("Failed:", e$message)
        )
      )
    })
    
    rm(temp_data)
    gc()
  }
}
write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_630_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_630_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####

age_groups <- c("1","2","3", "4", "5")

TRAITS <- c(
  "Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
  "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass"
)


library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()

output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data)
  
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c(
    "Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation",
    "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x"
  )
  
  for (tr in TRAITS) {
    
    target_trait <- paste0(tr, ".CWM")  
    
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    cat(sprintf(
      "\n>>> 正在全量计算: AgeGroup %s | 目标性状: %s (样本量: %d) <<<\n",
      ag, target_trait, nrow(temp_data)
    ))
    
    # 样本量阈值
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0(
      "Carbon_Mg_ha ~ ",
      target_trait,
      " + cv_dbh + MAT + MAP + VPD + pH_mean + clay_mean"
    )
    
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      
      # =========================================================
      # A. 先尝试带空间相关结构
      # =========================================================
      m <- tryCatch({
        
        nlme::lme(
          fixed = fixed_formula,
          random = ~1 | ECO_ID/COUNTYCD.x,
          correlation = nlme::corExp(
            form = ~ LON + LAT,
            nugget = TRUE
          ),
          data = temp_data,
          control = nlme::lmeControl(
            opt = "optim",
            maxIter = 100,
            msMaxIter = 100
          )
        )
        
      }, error = function(e1) {
        
        cat("  [警告] 空间相关结构失败，自动移除 correlation 重新运行\n")
        cat("  [原始错误]: ", e1$message, "\n")
        
        # =====================================================
        # B. 如果失败，则删除 correlation 后重新运行
        # =====================================================
        nlme::lme(
          fixed = fixed_formula,
          random = ~1 | ECO_ID/COUNTYCD.x,
          data = temp_data,
          control = nlme::lmeControl(
            opt = "optim",
            maxIter = 100,
            msMaxIter = 100
          )
        )
      })
      
      # --- B. 提取系数和 R2 ---
      at <- as.data.frame(summary(m)$tTable)
      
      r2_val <- MuMIn::r.squaredGLMM(m)
      
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor <- rownames(at)
      at$N_obs <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(
        as.matrix(coords_clean),
        k = 8
      )
      
      nb <- spdep::knn2nb(knn)
      
      lw <- spdep::nb2listw(
        nb,
        style = "W",
        zero.policy = TRUE
      )
      
      moran_res <- spdep::moran.test(
        res_clean,
        lw
      )
      
      # --- D. ---
      img_file <- paste0(
        figure_path,
        "Variogram_Japan_NaturalForests_630_",
        ag,
        "_",
        tr,
        ".png"
      )
      
      png(
        img_file,
        width = 800,
        height = 600,
        res = 120
      )
      
      vg <- nlme::Variogram(
        m,
        form = ~ LON + LAT,
        resType = "normalized"
      )
      
      print(
        plot(
          vg,
          main = paste0(
            "Variogram: AgeGroup ",
            ag,
            " | Trait: ",
            tr
          )
        )
      )
      
      dev.off()
      
      # --- E. 存入 Moran's I 简表 ---
      trait_p <- at[
        grep(target_trait, rownames(at)),
        "p-value"
      ][1]
      
      cv_p <- at[
        grep("cv_dbh", rownames(at)),
        "p-value"
      ][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag,
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      
      
      final_results <<- rbind(
        final_results,
        data.frame(
          AgeGroup = ag,
          Trait = tr,
          Status = paste("Failed:", e$message)
        )
      )
    })
    
    rm(temp_data)
    gc()
  }
}

write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_630_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_630_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
####620NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))

data_sem=read.csv("E:/GEB data/Japan_data_sem_NaturalForests.csv" )

data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest") , ]
fdis_cols <- c(
  "Shade.tolerance.FDis", "Drought.tolerance.FDis", "WD_mean_add.FDis",
  "LMA.FDis", "Nmass.FDis", "Pmass.FDis", 
  "Leaf.longevity.FDis"
)


data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()





log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}


#########FD: Moran's I####


age_groups <- c("2" ,"3", "4", "5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_Japan_NaturalForests_620_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_620_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_620_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####

age_groups <- c("2","3", "4", "5")

TRAITS <- c(
  "Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
  "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass"
)


library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()

output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data)
  
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c(
    "Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation",
    "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x"
  )
  
  for (tr in TRAITS) {
    
    target_trait <- paste0(tr, ".CWM")  
    
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    cat(sprintf(
      "\n>>> 正在全量计算: AgeGroup %s | 目标性状: %s (样本量: %d) <<<\n",
      ag, target_trait, nrow(temp_data)
    ))
    
    # 样本量阈值
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0(
      "Carbon_Mg_ha ~ ",
      target_trait,
      " + cv_dbh + MAT + MAP + VPD + pH_mean + clay_mean"
    )
    
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      
      # =========================================================
      # A. 先尝试带空间相关结构
      # =========================================================
      m <- tryCatch({
        
        nlme::lme(
          fixed = fixed_formula,
          random = ~1 | COUNTYCD.x,
          correlation = nlme::corExp(
            form = ~ LON + LAT,
            nugget = TRUE
          ),
          data = temp_data,
          control = nlme::lmeControl(
            opt = "optim",
            maxIter = 100,
            msMaxIter = 100
          )
        )
        
      }, error = function(e1) {
        
        cat("  [警告] 空间相关结构失败，自动移除 correlation 重新运行\n")
        cat("  [原始错误]: ", e1$message, "\n")
        
        # =====================================================
        # B. 如果失败，则删除 correlation 后重新运行
        # =====================================================
        nlme::lme(
          fixed = fixed_formula,
          random = ~1 | COUNTYCD.x,
          data = temp_data,
          control = nlme::lmeControl(
            opt = "optim",
            maxIter = 100,
            msMaxIter = 100
          )
        )
      })
      
      # --- B. 提取系数和 R2 ---
      at <- as.data.frame(summary(m)$tTable)
      
      r2_val <- MuMIn::r.squaredGLMM(m)
      
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor <- rownames(at)
      at$N_obs <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(
        as.matrix(coords_clean),
        k = 8
      )
      
      nb <- spdep::knn2nb(knn)
      
      lw <- spdep::nb2listw(
        nb,
        style = "W",
        zero.policy = TRUE
      )
      
      moran_res <- spdep::moran.test(
        res_clean,
        lw
      )
      
      # --- D. ---
      img_file <- paste0(
        figure_path,
        "Variogram_Japan_NaturalForests_630_",
        ag,
        "_",
        tr,
        ".png"
      )
      
      png(
        img_file,
        width = 800,
        height = 600,
        res = 120
      )
      
      vg <- nlme::Variogram(
        m,
        form = ~ LON + LAT,
        resType = "normalized"
      )
      
      print(
        plot(
          vg,
          main = paste0(
            "Variogram: AgeGroup ",
            ag,
            " | Trait: ",
            tr
          )
        )
      )
      
      dev.off()
      
      # --- E. 存入 Moran's I 简表 ---
      trait_p <- at[
        grep(target_trait, rownames(at)),
        "p-value"
      ][1]
      
      cv_p <- at[
        grep("cv_dbh", rownames(at)),
        "p-value"
      ][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag,
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      
      
      final_results <<- rbind(
        final_results,
        data.frame(
          AgeGroup = ag,
          Trait = tr,
          Status = paste("Failed:", e$message)
        )
      )
    })
    
    rm(temp_data)
    gc()
  }
}

write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_620_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_620_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
##===============================================================================
## HN_NaturalForests: LMM and Moran's I ---------
##===============================================================================
rm(list=ls())
####220NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/HN_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_220_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_220_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_220_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_220_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
####630NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/HN_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer-Broadleaf Mixed Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_630_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_630_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_630_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_630_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_630_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_630_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
####620NaturalForests####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/HN_data_sem_NaturalForests.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Broadleaf Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean),FDis=scale(FDis),STDAGE=scale(STDAGE)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}


#########FD: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 20) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_620_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_620_NaturalForests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_620_NaturalForests.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", 
            "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 30) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_620_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_620_NaturalForests_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_620_NaturalForests_CWM.csv"), row.names = FALSE)

#################################################################################
#################################################################################
#################################################################################
#################################################################################

data_orig <- read.csv("E:/GEB data/FIA_Plot_AnalysesPre.csv")

data_sem_remaining <- data_orig[!(data_orig$FLDTYPCD %in% c(103, 161, 141, 142, 162) & 
                                    data_orig$STDORGCD == 1 & 
                                    data_orig$STDAGE >= 11 & 
                                    data_orig$STDAGE <= 30), ]

data_2_ids <- read.csv("E:/GEB data/FIA_Plot_Plantations220_Age2_稀疏.csv")

data_2_full_rows <- data_orig[data_orig$PLT_CN %in% data_2_ids$PLT_CN, ]

data_final <- rbind(data_sem_remaining, data_2_full_rows)
write.csv( data_final, "E:/GEB data/FIA_Plot_AnalysesPre_稀疏.csv" )

##===============================================================================
## American_Plantations: LMM and Moran's I ---------
##===============================================================================
rm(list=ls())
####220Plantations####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/FIA_data_sem_Plantations.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}
#########FD: Moran's I####


age_groups <- c( "1", "2", "3", "4")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_Plantations_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_220_Plantations.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_220_Plantations.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3", "4")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_American_Plantations_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_American_220_Plantations_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_American_220_Plantations_CWM.csv"), row.names = FALSE)

##===============================================================================
## Japan_Plantations:LMM and  Moran's I ---------
##===============================================================================
rm(list=ls())
####220Plantations####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/Japan_data_sem_Plantations.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_Japan_Plantations_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_220_Plantations.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_220_Plantations.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3", "4","5")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data#,
        #control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_Japan_Plantations_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_Japan_220_Plantations_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_Japan_220_Plantations_CWM.csv"), row.names = FALSE)

#################################################################################
##===============================================================================
## HN_Plantations: LMM and Moran's I ---------
##===============================================================================
rm(list=ls())
####220Plantations####
TRAITS<-c(
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,
  "WD_mean_add",
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"
  
)

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))


data_sem=read.csv("E:/GEB data/HN_data_sem_Plantations.csv" )
data_sem <- data_sem[data_sem$FLDTYPCD %in% c("Conifer Forest") , ]

data_sem <- data_sem %>%
  mutate(
    cv_dbh               = scale(log(cv_dbh+1 )),
    MAT                  = scale(MAT),
    MAP                  = scale(MAP),
    VPD                  = scale(VPD),
    Elevation            = scale(Elevation),
    pH_mean              = scale(pH_mean),
    bdod_mean            = scale(bdod_mean),
    clay_mean            = scale(clay_mean),
    nitrogen_mean        = scale(nitrogen_mean)
  )
sum(duplicated(data_sem[, c("LON","LAT")]))

library(dplyr)
data_sem$PLT_CN=as.character(data_sem$PLT_CN )

data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()


data_sem$Carbon_Mg_ha=log(data_sem$Carbon_Mg_ha+1  )
data_sem <- data_sem %>% na.omit()

log.norm.scale <- function(df) {
  
  log_df <-df#log(df+1)
  # return((log_df - min(log_df)) / (max(log_df) - min(log_df)))
  return(scale(log_df))
}
###

for (c in c("Shade.tolerance.FDis",
            "Drought.tolerance.FDis",
            "WD_mean_add.FDis",
            "LMA.FDis",
            "Nmass.FDis",
            "Pmass.FDis",
            "Leaf.longevity.FDis",
            "Max.height.FDis",
            "Leaf.FDis",
            "Shade.tolerance.CWM",
            "Drought.tolerance.CWM",
            "WD_mean_add.CWM",
            "LMA.CWM",
            "Nmass.CWM",
            "Pmass.CWM",
            "Leaf.longevity.CWM",
            
            "Max.height.CWM",
            "Leaf.CWM"
            
)) {
  data_sem[, c] <- log.norm.scale(data_sem[, c])
}

#########FD: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".FDis")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | ECO_ID/COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_NaturalForests_220_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}



write.csv(meta.table, paste0(output_path, "LME_Detail_HN_220_Plantations.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_220_Plantations.csv"), row.names = FALSE)

#########CWM: Moran's I####


age_groups <- c( "1", "2", "3")
TRAITS <- c("Leaf", "Max.height", "Shade.tolerance", "Drought.tolerance", 
            "WD_mean_add", "LMA", "Leaf.longevity", "Nmass", "Pmass") 



library(nlme)
library(MuMIn)
library(spdep)

meta.table <- data.frame()
final_results <- data.frame()
output_path <- "E:/GEB data/"
figure_path <- "E:/GEB data/Plot_Moran/"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(figure_path)) dir.create(figure_path, recursive = TRUE)


for (ag in age_groups) {
  
  
  if(exists("sub_data")) rm(sub_data) 
  sub_data_ag <- data_sem[data_sem$AgeGroup == ag, ]
  
  
  needed_cols <- c("Carbon_Mg_ha", "cv_dbh", "MAT", "MAP", "VPD", "Elevation", 
                   "pH_mean", "clay_mean", "LON", "LAT", "ECO_ID", "COUNTYCD.x")
  
  for (tr in TRAITS) {
    target_trait <- paste0(tr, ".CWM")
    if (!(target_trait %in% colnames(sub_data_ag))) next
    
    
    temp_data <- sub_data_ag[, c(needed_cols, target_trait)]
    temp_data <- na.omit(temp_data)
    
    
    
    
    if(nrow(temp_data) < 35) {
      cat("[Skipped] Insufficient sample size")
      next
    }
    
    
    formula_str <- paste0("Carbon_Mg_ha ~ ", target_trait, " + cv_dbh + MAT + MAP + VPD  + pH_mean + clay_mean")
    fixed_formula <- as.formula(formula_str)
    
    tryCatch({
      # --- A.  ---
      m <- nlme::lme(
        fixed = fixed_formula,
        random = ~1 | COUNTYCD.x, 
        correlation = nlme::corExp(form = ~ LON + LAT, nugget = TRUE),
        data = temp_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      
      # --- B.  ---
      at <- as.data.frame(summary(m)$tTable)
      r2_val <- MuMIn::r.squaredGLMM(m)
      at$R2m <- r2_val[1, "R2m"]
      at$R2c <- r2_val[1, "R2c"]
      at$factor   <- rownames(at)
      at$N_obs    <- nrow(temp_data)
      at$TRAIT_ID <- tr
      at$AgeGroup <- ag
      
      
      meta.table <<- rbind(meta.table, at)
      
      # --- C. Moran's I ---
      res_norm <- resid(m, type = "normalized")
      coords <- temp_data[, c("LON", "LAT")]
      
      keep <- complete.cases(coords, res_norm)
      res_clean <- res_norm[keep]
      coords_clean <- coords[keep, ]
      
      knn <- spdep::knearneigh(as.matrix(coords_clean), k = 8)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_res <- spdep::moran.test(res_clean, lw)
      
      # --- D. ---
      img_file <- paste0(figure_path, "Variogram_HN_Plantations_220_CWM_", ag, "_", tr, ".png")
      png(img_file, width = 800, height = 600, res = 120)
      vg <- nlme::Variogram(m, form = ~ LON + LAT, resType = "normalized")
      print(plot(vg, main = paste0("Variogram: AgeGroup ", ag, " | Trait: ", tr)))
      dev.off()
      
      # --- E.  ---
      
      trait_p <- at[grep(target_trait, rownames(at)), "p-value"][1]
      cv_p    <- at[grep("cv_dbh", rownames(at)), "p-value"][1]
      
      tmp_df <- data.frame(
        AgeGroup = ag, 
        Trait = tr,
        Trait_P = trait_p,
        CV_DBH_P = cv_p,
        Moran_I = moran_res$estimate[1],
        Moran_P = moran_res$p.value,
        Status = "Success"
      )
      final_results <<- rbind(final_results, tmp_df)
      
      
      
    }, error = function(e) {
      
      final_results <<- rbind(final_results, data.frame(AgeGroup=ag, Trait=tr, Status=paste("Failed:", e$message)))
    })
    
    rm(temp_data)
    gc()
  }
}


write.csv(meta.table, paste0(output_path, "LME_Detail_HN_220_Plantations_CWM.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_path, "Moran_Summary_HN_220_Plantations_CWM.csv"), row.names = FALSE)


##===============================================================================
## American:NaturalForests ---------
##===============================================================================
rm(list=ls())
library(MuMIn)
#### meta data ####
meta.data_220_FD <-read.csv("E:/GEB data/LME_Detail_American_220_NaturalForests.csv")
meta.data_220_CWM <-read.csv("E:/GEB data/LME_Detail_American_220_NaturalForests_CWM.csv")
meta_data_220_FD <- meta.data_220_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_220_CWM <- meta.data_220_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_220=rbind(meta_data_220_FD ,meta_data_220_CWM )

meta.data_630_FD <-read.csv("E:/GEB data/LME_Detail_American_630_NaturalForests.csv")
meta.data_630_CWM <-read.csv("E:/GEB data/LME_Detail_American_630_NaturalForests_CWM.csv")
meta_data_630_FD <- meta.data_630_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_630_CWM <- meta.data_630_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_630=rbind(meta_data_630_FD ,meta_data_630_CWM )

meta.data_620_FD <-read.csv("E:/GEB data/LME_Detail_American_620_NaturalForests.csv")
meta.data_620_CWM <-read.csv("E:/GEB data/LME_Detail_American_620_NaturalForests_CWM.csv")
meta_data_620_FD <- meta.data_620_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_620_CWM <- meta.data_620_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_620=rbind(meta_data_620_FD ,meta_data_620_CWM )


meta.data_220$Types <-"220_General"
meta.data_630$Types <-"630_General"
meta.data_620$Types <-"620_General"


meta.data=rbind(meta.data_220,meta.data_630,meta.data_620
                
)

str(meta.data)
meta.data$RESPONSE="Carbon_Mg_ha"
meta.data$trait=meta.data$factor
meta.data$CAT <- ifelse(grepl("\\.CWM$", meta.data$factor), "CWM_A", 
                        ifelse(grepl("\\.FDis$", meta.data$factor), "FD_A", NA))
meta.data =meta.data [,c( "Value",     "Std.Error" ,"DF" ,       "t.value" ,  "p.value" ,  "R2m" ,
                          "R2c"  ,     "factor",    "N_obs" ,  "RESPONSE", "CAT" , "TRAIT_ID" , "AgeGroup" , "Types"   )]

colnames(meta.data )=c( "Estimate" ,  "Std..Error", "df","t.value" , "Pr...t.." ,  
                        "R2m" ,"R2c",  "factor" ,"Nptag",  "RESPONSE", "CAT" , "TRAIT" , "AgeGroup" , "Types" )




t.cor<-c(  
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,"LDMC"  ,  
  "WD_mean_add",          
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"  
  
)

#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate,
                                meta.data$Estimate)

### subset data for analyses
dat<-droplevels(subset(meta.data, RESPONSE=="Carbon_Mg_ha")) ## accumulated stand volume
dat$var=rownames(dat)


unique(dat)
dat <- dat %>%
  dplyr::mutate(
    AGE = case_when(
      Types == "220_General" & AgeGroup == 1 ~ "8",  
      Types == "220_General" & AgeGroup == 2 ~ "21", 
      Types == "220_General" & AgeGroup == 3 ~ "40",
      Types == "220_General" & AgeGroup == 4 ~ "60",
      Types == "220_General" & AgeGroup == 5 ~ "82", 
      
      Types == "630_General" & AgeGroup == 1 ~ "8",  
      Types == "630_General" & AgeGroup == 2 ~ "21", 
      Types == "630_General" & AgeGroup == 3 ~ "41", 
      Types == "630_General" & AgeGroup == 4 ~ "60", 
      Types == "630_General" & AgeGroup == 5 ~ "82", 
      
      
      Types == "620_General" & AgeGroup == 1 ~ "8",  
      Types == "620_General" & AgeGroup == 2 ~ "22", 
      Types == "620_General" & AgeGroup == 3 ~ "41", Types == "620_General" & AgeGroup == 4 ~ "61", 
      Types == "620_General" & AgeGroup == 5 ~ "84",
      TRUE ~ as.character(NA)  # For other rows, assign NA or other default value
    )
  )


## center AGE 
dat$AGE2 <- scale(as.numeric(dat$AGE)  )
#dat$AGE2 <- dat$AGE- mean(dat$AGE, na.rm = TRUE)  
dat$AGE2.sq<-dat$AGE2^2
##########
##General
dat_220_General=dat[dat$Types == "220_General", ]
# m2_220_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "220_General", ])

dat_220_General_FD_A <- subset(dat_220_General, CAT == "FD_A")
dat_220_General_CWM_A <- subset(dat_220_General, CAT == "CWM_A")

library(mgcv)
m2_220_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_FD_A)
summary(m2_220_General_FD_A)
r.squaredGLMM(m2_220_General_FD_A)
m2_220_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_CWM_A)
summary(m2_220_General_CWM_A )
r.squaredGLMM(m2_220_General_CWM_A)

dat_220_General_FD_A$predict.Est.corr.m2 <- predict(m2_220_General_FD_A)
dat_220_General_CWM_A$predict.Est.corr.m2 <- predict(m2_220_General_CWM_A)



dat_220_General_pred <- rbind(
  dat_220_General_FD_A, dat_220_General_CWM_A )

dat_220_General_FD_A$predict.Est.corr.m2<-predict(m2_220_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=mean)
sd_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=sd)
n_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=length)
sd_220_General_FD_A$Estimate.corr.se<-sd_220_General_FD_A$Estimate.corr/sqrt(n_220_General_FD_A$Estimate.corr)
msd_220_General_FD_A<-merge(m_220_General_FD_A,sd_220_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

dat_220_General_CWM_A$predict.Est.corr.m2<-predict(m2_220_General_CWM_A)



# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=mean)
sd_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=sd)
n_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=length)
sd_220_General_CWM_A$Estimate.corr.se<-sd_220_General_CWM_A$Estimate.corr/sqrt(n_220_General_CWM_A$Estimate.corr)
msd_220_General_CWM_A<-merge(m_220_General_CWM_A,sd_220_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_220_General=rbind(
  msd_220_General_FD_A, msd_220_General_CWM_A
)

dat_630_General=dat[dat$Types == "630_General", ]


dat_630_General_FD_A <- subset(dat_630_General, CAT == "FD_A")
dat_630_General_CWM_A <- subset(dat_630_General, CAT == "CWM_A")




m2_630_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_FD_A)
summary(m2_630_General_FD_A)
r.squaredGLMM(m2_630_General_FD_A)
m2_630_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_CWM_A)
summary(m2_630_General_CWM_A )
r.squaredGLMM(m2_630_General_CWM_A)


dat_630_General_FD_A$predict.Est.corr.m2 <- predict(m2_630_General_FD_A)
dat_630_General_CWM_A$predict.Est.corr.m2 <- predict(m2_630_General_CWM_A)

dat_630_General_pred <- rbind(
  dat_630_General_FD_A, dat_630_General_CWM_A )

dat_630_General_FD_A$predict.Est.corr.m2<-predict(m2_630_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=mean)
sd_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=sd)
n_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=length)
sd_630_General_FD_A$Estimate.corr.se<-sd_630_General_FD_A$Estimate.corr/sqrt(n_630_General_FD_A$Estimate.corr)
msd_630_General_FD_A<-merge(m_630_General_FD_A,sd_630_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_630_General_CWM_A$predict.Est.corr.m2<-predict(m2_630_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=mean)
sd_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=sd)
n_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=length)
sd_630_General_CWM_A$Estimate.corr.se<-sd_630_General_CWM_A$Estimate.corr/sqrt(n_630_General_CWM_A$Estimate.corr)
msd_630_General_CWM_A<-merge(m_630_General_CWM_A,sd_630_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_630_General=rbind(
  msd_630_General_FD_A, msd_630_General_CWM_A
)

##General
dat_620_General=dat[dat$Types == "620_General", ]
# m2_620_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "620_General", ])

dat_620_General_FD_A <- subset(dat_620_General, CAT == "FD_A")
dat_620_General_CWM_A <- subset(dat_620_General, CAT == "CWM_A")

m2_620_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_FD_A)
summary(m2_620_General_FD_A)
r.squaredGLMM(m2_620_General_FD_A)
m2_620_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_CWM_A)
summary(m2_620_General_CWM_A )
r.squaredGLMM(m2_620_General_CWM_A)


dat_620_General_FD_A$predict.Est.corr.m2 <- predict(m2_620_General_FD_A)
dat_620_General_CWM_A$predict.Est.corr.m2 <- predict(m2_620_General_CWM_A)

dat_620_General_pred <- rbind(
  dat_620_General_FD_A, dat_620_General_CWM_A )

dat_620_General_FD_A$predict.Est.corr.m2<-predict(m2_620_General_FD_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=mean)
sd_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=sd)
n_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=length)
sd_620_General_FD_A$Estimate.corr.se<-sd_620_General_FD_A$Estimate.corr/sqrt(n_620_General_FD_A$Estimate.corr)
msd_620_General_FD_A<-merge(m_620_General_FD_A,sd_620_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_620_General_CWM_A$predict.Est.corr.m2<-predict(m2_620_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=mean)
sd_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=sd)
n_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=length)
sd_620_General_CWM_A$Estimate.corr.se<-sd_620_General_CWM_A$Estimate.corr/sqrt(n_620_General_CWM_A$Estimate.corr)
msd_620_General_CWM_A<-merge(m_620_General_CWM_A,sd_620_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_620_General=rbind(
  msd_620_General_FD_A, msd_620_General_CWM_A
)



msd_220_General$Types="220_General"
msd_630_General$Types="630_General"
msd_620_General$Types="620_General"

msd=rbind(msd_220_General,msd_630_General,  msd_620_General     )


####save data####
msd$Nation="American"
write.csv(msd,"E:/GEB data/msd_American_NaturalForests.csv"  )
msd_5=read.csv("E:/GEB data/msd_American_NaturalForests.csv")
###### paired-t-test for CAT difference per year
msd$Types[which(msd$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd$Types[which(msd$Types %in% c( "630_General" ))] = c("CBMF_General")
msd$Types[which(msd$Types %in% c("620_General"))] = c("BMF_General")



msd$Types <- factor(msd$Types, levels = c("CP_General","CBMF_General","BMF_General"
))
msd_5$Types[which(msd_5$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd_5$Types[which(msd_5$Types %in% c( "630_General" ))] = c("CBMF_General")
msd_5$Types[which(msd_5$Types %in% c("620_General"))] = c("BMF_General")


msd_5$Types <- factor(msd_5$Types, levels = c("CP_General","CBMF_General","BMF_General" ))

write.csv(msd,"E:/GEB data/Meta_dat1_American_NaturalForests.csv" )

write.csv(dat,"E:/GEB data/Meta_dat2_American_NaturalForests.csv" )

##===============================================================================
## Japan:NaturalForests ---------
##===============================================================================
rm(list=ls())
#### meta data ####
meta.data_220_FD <-read.csv("E:/GEB data/LME_Detail_Japan_220_NaturalForests.csv")
meta.data_220_CWM <-read.csv("E:/GEB data/LME_Detail_Japan_220_NaturalForests_CWM.csv")
meta_data_220_FD <- meta.data_220_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_220_CWM <- meta.data_220_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_220=rbind(meta_data_220_FD ,meta_data_220_CWM )

meta.data_630_FD <-read.csv("E:/GEB data/LME_Detail_Japan_630_NaturalForests.csv")
meta.data_630_CWM <-read.csv("E:/GEB data/LME_Detail_Japan_630_NaturalForests_CWM.csv")
meta_data_630_FD <- meta.data_630_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_630_CWM <- meta.data_630_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_630=rbind(meta_data_630_FD ,meta_data_630_CWM )

meta.data_620_FD <-read.csv("E:/GEB data/LME_Detail_Japan_620_NaturalForests.csv")
meta.data_620_CWM <-read.csv("E:/GEB data/LME_Detail_Japan_620_NaturalForests_CWM.csv")
meta_data_620_FD <- meta.data_620_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_620_CWM <- meta.data_620_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_620=rbind(meta_data_620_FD ,meta_data_620_CWM )


meta.data_220$Types <-"220_General"
meta.data_630$Types <-"630_General"
meta.data_620$Types <-"620_General"


meta.data=rbind(meta.data_220,meta.data_630,meta.data_620
                
)

str(meta.data)
meta.data$RESPONSE="Carbon_Mg_ha"
meta.data$trait=meta.data$factor
meta.data$CAT <- ifelse(grepl("\\.CWM$", meta.data$factor), "CWM_A", 
                        ifelse(grepl("\\.FDis$", meta.data$factor), "FD_A", NA))
meta.data =meta.data [,c( "Value",     "Std.Error" ,"DF" ,       "t.value" ,  "p.value" ,  "R2m" ,
                          "R2c"  ,     "factor",    "N_obs" ,  "RESPONSE", "CAT" , "TRAIT_ID" , "AgeGroup" , "Types"   )]

colnames(meta.data )=c( "Estimate" ,  "Std..Error", "df","t.value" , "Pr...t.." ,  
                        "R2m" ,"R2c",  "factor" ,"Nptag",  "RESPONSE", "CAT" , "TRAIT" , "AgeGroup" , "Types" )




t.cor<-c(  
  "Leaf", 
  "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,"LDMC"  ,  
  "WD_mean_add",          
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"  
  
)

#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate,
                                meta.data$Estimate)

### subset data for analyses
dat<-droplevels(subset(meta.data, RESPONSE=="Carbon_Mg_ha")) ## accumulated stand volume
dat$var=rownames(dat)


unique(dat)
dat <- dat %>%
  dplyr::mutate(
    AGE = case_when(
      Types == "220_General" & AgeGroup == 1 ~ "7",  
      Types == "220_General" & AgeGroup == 2 ~ "23", 
      Types == "220_General" & AgeGroup == 3 ~ "43",
      Types == "220_General" & AgeGroup == 4 ~ "61",
      Types == "220_General" & AgeGroup == 5 ~ "83", 
      
      Types == "630_General" & AgeGroup == 1 ~ "8",  
      Types == "630_General" & AgeGroup == 2 ~ "24", 
      Types == "630_General" & AgeGroup == 3 ~ "41", 
      Types == "630_General" & AgeGroup == 4 ~ "60", 
      Types == "630_General" & AgeGroup == 5 ~ "86", 
      
      
      Types == "620_General" & AgeGroup == 1 ~ "8",  
      Types == "620_General" & AgeGroup == 2 ~ "25", 
      Types == "620_General" & AgeGroup == 3 ~ "43", Types == "620_General" & AgeGroup == 4 ~ "60", 
      Types == "620_General" & AgeGroup == 5 ~ "83",
      TRUE ~ as.character(NA)  # For other rows, assign NA or other default value
    )
  )


## center AGE 
dat$AGE2 <- scale(as.numeric(dat$AGE)  )
#dat$AGE2 <- dat$AGE- mean(dat$AGE, na.rm = TRUE)  
dat$AGE2.sq<-as.numeric(dat$AGE)^2
##########
##General
dat_220_General=dat[dat$Types == "220_General", ]
# m2_220_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "220_General", ])

dat_220_General_FD_A <- subset(dat_220_General, CAT == "FD_A")
dat_220_General_CWM_A <- subset(dat_220_General, CAT == "CWM_A")

library(mgcv)

m2_220_General_FD_A <- lm(
  Estimate.corr ~ AGE2 ,
  data = dat_220_General_FD_A
)
summary(m2_220_General_FD_A)
r.squaredGLMM(m2_220_General_FD_A)
m2_220_General_CWM_A <-lm(
  Estimate.corr ~ AGE2 , data = dat_220_General_CWM_A)
summary(m2_220_General_CWM_A )
r.squaredGLMM(m2_220_General_CWM_A)

dat_220_General_FD_A$predict.Est.corr.m2 <- predict(m2_220_General_FD_A)
dat_220_General_CWM_A$predict.Est.corr.m2 <- predict(m2_220_General_CWM_A)



dat_220_General_pred <- rbind(
  dat_220_General_FD_A, dat_220_General_CWM_A )

dat_220_General_FD_A$predict.Est.corr.m2<-predict(m2_220_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=mean)
sd_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=sd)
n_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=length)
sd_220_General_FD_A$Estimate.corr.se<-sd_220_General_FD_A$Estimate.corr/sqrt(n_220_General_FD_A$Estimate.corr)
msd_220_General_FD_A<-merge(m_220_General_FD_A,sd_220_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

dat_220_General_CWM_A$predict.Est.corr.m2<-predict(m2_220_General_CWM_A)



# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=mean)
sd_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=sd)
n_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=length)
sd_220_General_CWM_A$Estimate.corr.se<-sd_220_General_CWM_A$Estimate.corr/sqrt(n_220_General_CWM_A$Estimate.corr)
msd_220_General_CWM_A<-merge(m_220_General_CWM_A,sd_220_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_220_General=rbind(
  msd_220_General_FD_A, msd_220_General_CWM_A
)

dat_630_General=dat[dat$Types == "630_General", ]


dat_630_General_FD_A <- subset(dat_630_General, CAT == "FD_A")
dat_630_General_CWM_A <- subset(dat_630_General, CAT == "CWM_A")




m2_630_General_FD_A <- lm(
  Estimate.corr ~ AGE2 ,
  data = dat_630_General_FD_A
)
summary(m2_630_General_FD_A)
r.squaredGLMM(m2_630_General_FD_A)
m2_630_General_CWM_A <-lm(
  Estimate.corr ~ AGE2 , data = dat_630_General_CWM_A)
summary(m2_630_General_CWM_A )
r.squaredGLMM(m2_630_General_CWM_A)


dat_630_General_FD_A$predict.Est.corr.m2 <- predict(m2_630_General_FD_A)
dat_630_General_CWM_A$predict.Est.corr.m2 <- predict(m2_630_General_CWM_A)

dat_630_General_pred <- rbind(
  dat_630_General_FD_A, dat_630_General_CWM_A )

dat_630_General_FD_A$predict.Est.corr.m2<-predict(m2_630_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=mean)
sd_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=sd)
n_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=length)
sd_630_General_FD_A$Estimate.corr.se<-sd_630_General_FD_A$Estimate.corr/sqrt(n_630_General_FD_A$Estimate.corr)
msd_630_General_FD_A<-merge(m_630_General_FD_A,sd_630_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_630_General_CWM_A$predict.Est.corr.m2<-predict(m2_630_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=mean)
sd_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=sd)
n_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=length)
sd_630_General_CWM_A$Estimate.corr.se<-sd_630_General_CWM_A$Estimate.corr/sqrt(n_630_General_CWM_A$Estimate.corr)
msd_630_General_CWM_A<-merge(m_630_General_CWM_A,sd_630_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_630_General=rbind(
  msd_630_General_FD_A, msd_630_General_CWM_A
)

##General
dat_620_General=dat[dat$Types == "620_General", ]
# m2_620_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "620_General", ])

dat_620_General_FD_A <- subset(dat_620_General, CAT == "FD_A")
dat_620_General_CWM_A <- subset(dat_620_General, CAT == "CWM_A")

m2_620_General_FD_A <- lm(
  Estimate.corr ~ AGE2 ,
  data = dat_620_General_FD_A
)
summary(m2_620_General_FD_A)
r.squaredGLMM(m2_620_General_FD_A)
m2_620_General_CWM_A <-lm(
  Estimate.corr ~ AGE2 , data = dat_620_General_CWM_A)
summary(m2_620_General_CWM_A )
r.squaredGLMM(m2_620_General_CWM_A)


dat_620_General_FD_A$predict.Est.corr.m2 <- predict(m2_620_General_FD_A)
dat_620_General_CWM_A$predict.Est.corr.m2 <- predict(m2_620_General_CWM_A)

dat_620_General_pred <- rbind(
  dat_620_General_FD_A, dat_620_General_CWM_A )

dat_620_General_FD_A$predict.Est.corr.m2<-predict(m2_620_General_FD_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=mean)
sd_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=sd)
n_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=length)
sd_620_General_FD_A$Estimate.corr.se<-sd_620_General_FD_A$Estimate.corr/sqrt(n_620_General_FD_A$Estimate.corr)
msd_620_General_FD_A<-merge(m_620_General_FD_A,sd_620_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_620_General_CWM_A$predict.Est.corr.m2<-predict(m2_620_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=mean)
sd_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=sd)
n_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=length)
sd_620_General_CWM_A$Estimate.corr.se<-sd_620_General_CWM_A$Estimate.corr/sqrt(n_620_General_CWM_A$Estimate.corr)
msd_620_General_CWM_A<-merge(m_620_General_CWM_A,sd_620_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_620_General=rbind(
  msd_620_General_FD_A, msd_620_General_CWM_A
)



msd_220_General$Types="220_General"
msd_630_General$Types="630_General"
msd_620_General$Types="620_General"

msd=rbind(msd_220_General,msd_630_General,  msd_620_General     )


####save data####
msd$Nation="Japan"
write.csv(msd,"E:/GEB data/msd_Japan_NaturalForests.csv"  )
msd_5=read.csv("E:/GEB data/msd_Japan_NaturalForests.csv")
###### paired-t-test for CAT difference per year
msd$Types[which(msd$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd$Types[which(msd$Types %in% c( "630_General" ))] = c("CBMF_General")
msd$Types[which(msd$Types %in% c("620_General"))] = c("BMF_General")



msd$Types <- factor(msd$Types, levels = c("CP_General","CBMF_General","BMF_General"
))
msd_5$Types[which(msd_5$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd_5$Types[which(msd_5$Types %in% c( "630_General" ))] = c("CBMF_General")
msd_5$Types[which(msd_5$Types %in% c("620_General"))] = c("BMF_General")


msd_5$Types <- factor(msd_5$Types, levels = c("CP_General","CBMF_General","BMF_General" ))

write.csv(msd,"E:/GEB data/Meta_dat1_Japan_NaturalForests.csv" )

write.csv(dat,"E:/GEB data/Meta_dat2_Japan_NaturalForests.csv" )

##===============================================================================
## HN:NaturalForests ---------
##===============================================================================
rm(list=ls())
#### meta data ####
meta.data_220_FD <-read.csv("E:/GEB data/LME_Detail_HN_220_NaturalForests.csv")
meta.data_220_CWM <-read.csv("E:/GEB data/LME_Detail_HN_220_NaturalForests_CWM.csv")
meta_data_220_FD <- meta.data_220_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_220_CWM <- meta.data_220_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_220=rbind(meta_data_220_FD ,meta_data_220_CWM )

meta.data_630_FD <-read.csv("E:/GEB data/LME_Detail_HN_630_NaturalForests.csv")
meta.data_630_CWM <-read.csv("E:/GEB data/LME_Detail_HN_630_NaturalForests_CWM.csv")
meta_data_630_FD <- meta.data_630_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_630_CWM <- meta.data_630_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_630=rbind(meta_data_630_FD ,meta_data_630_CWM )

meta.data_620_FD <-read.csv("E:/GEB data/LME_Detail_HN_620_NaturalForests.csv")
meta.data_620_CWM <-read.csv("E:/GEB data/LME_Detail_HN_620_NaturalForests_CWM.csv")
meta_data_620_FD <- meta.data_620_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_620_CWM <- meta.data_620_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_620=rbind(meta_data_620_FD ,meta_data_620_CWM )


meta.data_220$Types <-"220_General"
meta.data_630$Types <-"630_General"
meta.data_620$Types <-"620_General"


meta.data=rbind(meta.data_220,meta.data_630,meta.data_620
                
)

str(meta.data)
meta.data$RESPONSE="Carbon_Mg_ha"
meta.data$trait=meta.data$factor
meta.data$CAT <- ifelse(grepl("\\.CWM$", meta.data$factor), "CWM_A", 
                        ifelse(grepl("\\.FDis$", meta.data$factor), "FD_A", NA))
meta.data =meta.data [,c( "Value",     "Std.Error" ,"DF" ,       "t.value" ,  "p.value" ,  "R2m" ,
                          "R2c"  ,     "factor",    "N_obs" ,  "RESPONSE", "CAT" , "TRAIT_ID" , "AgeGroup" , "Types"   )]

colnames(meta.data )=c( "Estimate" ,  "Std..Error", "df","t.value" , "Pr...t.." ,  
                        "R2m" ,"R2c",  "factor" ,"Nptag",  "RESPONSE", "CAT" , "TRAIT" , "AgeGroup" , "Types" )




t.cor<-c(  
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,"LDMC"  ,  
  "WD_mean_add",          
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"  
  
)

#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate,
                                meta.data$Estimate)

### subset data for analyses
dat<-droplevels(subset(meta.data, RESPONSE=="Carbon_Mg_ha")) ## accumulated stand volume
dat$var=rownames(dat)


unique(dat)
dat <- dat %>%
  dplyr::mutate(
    AGE = case_when(
      Types == "220_General" & AgeGroup == 1 ~ "9",  
      Types == "220_General" & AgeGroup == 2 ~ "23", 
      Types == "220_General" & AgeGroup == 3 ~ "37",
      Types == "220_General" & AgeGroup == 4 ~ "55",
      Types == "220_General" & AgeGroup == 5 ~ "85", 
      
      Types == "630_General" & AgeGroup == 1 ~ "8",  
      Types == "630_General" & AgeGroup == 2 ~ "20", 
      Types == "630_General" & AgeGroup == 3 ~ "36", 
      Types == "630_General" & AgeGroup == 4 ~ "60", 
      Types == "630_General" & AgeGroup == 5 ~ "83", 
      
      
      Types == "620_General" & AgeGroup == 1 ~ "8",  
      Types == "620_General" & AgeGroup == 2 ~ "19", 
      Types == "620_General" & AgeGroup == 3 ~ "39", Types == "620_General" & AgeGroup == 4 ~ "61", 
      Types == "620_General" & AgeGroup == 5 ~ "80",
      TRUE ~ as.character(NA)  # For other rows, assign NA or other default value
    )
  )


## center AGE 
dat$AGE2 <- scale(as.numeric(dat$AGE)  )
#dat$AGE2 <- dat$AGE- mean(dat$AGE, na.rm = TRUE)  
dat$AGE2.sq<-dat$AGE2^2
##########
##General
dat_220_General=dat[dat$Types == "220_General", ]
# m2_220_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "220_General", ])

dat_220_General_FD_A <- subset(dat_220_General, CAT == "FD_A")
dat_220_General_CWM_A <- subset(dat_220_General, CAT == "CWM_A")

library(mgcv)
m2_220_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_FD_A)
summary(m2_220_General_FD_A)
r.squaredGLMM(m2_220_General_FD_A)
m2_220_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_CWM_A)
summary(m2_220_General_CWM_A )
r.squaredGLMM(m2_220_General_CWM_A)

dat_220_General_FD_A$predict.Est.corr.m2 <- predict(m2_220_General_FD_A)
dat_220_General_CWM_A$predict.Est.corr.m2 <- predict(m2_220_General_CWM_A)



dat_220_General_pred <- rbind(
  dat_220_General_FD_A, dat_220_General_CWM_A )

dat_220_General_FD_A$predict.Est.corr.m2<-predict(m2_220_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=mean)
sd_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=sd)
n_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=length)
sd_220_General_FD_A$Estimate.corr.se<-sd_220_General_FD_A$Estimate.corr/sqrt(n_220_General_FD_A$Estimate.corr)
msd_220_General_FD_A<-merge(m_220_General_FD_A,sd_220_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

dat_220_General_CWM_A$predict.Est.corr.m2<-predict(m2_220_General_CWM_A)



# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=mean)
sd_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=sd)
n_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=length)
sd_220_General_CWM_A$Estimate.corr.se<-sd_220_General_CWM_A$Estimate.corr/sqrt(n_220_General_CWM_A$Estimate.corr)
msd_220_General_CWM_A<-merge(m_220_General_CWM_A,sd_220_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_220_General=rbind(
  msd_220_General_FD_A, msd_220_General_CWM_A
)

dat_630_General=dat[dat$Types == "630_General", ]


dat_630_General_FD_A <- subset(dat_630_General, CAT == "FD_A")
dat_630_General_CWM_A <- subset(dat_630_General, CAT == "CWM_A")




m2_630_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_FD_A)
summary(m2_630_General_FD_A)
r.squaredGLMM(m2_630_General_FD_A)
m2_630_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_CWM_A)
summary(m2_630_General_CWM_A )
r.squaredGLMM(m2_630_General_CWM_A)


dat_630_General_FD_A$predict.Est.corr.m2 <- predict(m2_630_General_FD_A)
dat_630_General_CWM_A$predict.Est.corr.m2 <- predict(m2_630_General_CWM_A)

dat_630_General_pred <- rbind(
  dat_630_General_FD_A, dat_630_General_CWM_A )

dat_630_General_FD_A$predict.Est.corr.m2<-predict(m2_630_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=mean)
sd_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=sd)
n_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=length)
sd_630_General_FD_A$Estimate.corr.se<-sd_630_General_FD_A$Estimate.corr/sqrt(n_630_General_FD_A$Estimate.corr)
msd_630_General_FD_A<-merge(m_630_General_FD_A,sd_630_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_630_General_CWM_A$predict.Est.corr.m2<-predict(m2_630_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=mean)
sd_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=sd)
n_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=length)
sd_630_General_CWM_A$Estimate.corr.se<-sd_630_General_CWM_A$Estimate.corr/sqrt(n_630_General_CWM_A$Estimate.corr)
msd_630_General_CWM_A<-merge(m_630_General_CWM_A,sd_630_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_630_General=rbind(
  msd_630_General_FD_A, msd_630_General_CWM_A
)

##General
dat_620_General=dat[dat$Types == "620_General", ]
# m2_620_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "620_General", ])

dat_620_General_FD_A <- subset(dat_620_General, CAT == "FD_A")
dat_620_General_CWM_A <- subset(dat_620_General, CAT == "CWM_A")

m2_620_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_FD_A)
summary(m2_620_General_FD_A)
r.squaredGLMM(m2_620_General_FD_A)
m2_620_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_CWM_A)
summary(m2_620_General_CWM_A )
r.squaredGLMM(m2_620_General_CWM_A)


dat_620_General_FD_A$predict.Est.corr.m2 <- predict(m2_620_General_FD_A)
dat_620_General_CWM_A$predict.Est.corr.m2 <- predict(m2_620_General_CWM_A)

dat_620_General_pred <- rbind(
  dat_620_General_FD_A, dat_620_General_CWM_A )

dat_620_General_FD_A$predict.Est.corr.m2<-predict(m2_620_General_FD_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=mean)
sd_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=sd)
n_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=length)
sd_620_General_FD_A$Estimate.corr.se<-sd_620_General_FD_A$Estimate.corr/sqrt(n_620_General_FD_A$Estimate.corr)
msd_620_General_FD_A<-merge(m_620_General_FD_A,sd_620_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_620_General_CWM_A$predict.Est.corr.m2<-predict(m2_620_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=mean)
sd_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=sd)
n_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=length)
sd_620_General_CWM_A$Estimate.corr.se<-sd_620_General_CWM_A$Estimate.corr/sqrt(n_620_General_CWM_A$Estimate.corr)
msd_620_General_CWM_A<-merge(m_620_General_CWM_A,sd_620_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_620_General=rbind(
  msd_620_General_FD_A, msd_620_General_CWM_A
)



msd_220_General$Types="220_General"
msd_630_General$Types="630_General"
msd_620_General$Types="620_General"

msd=rbind(msd_220_General,msd_630_General,  msd_620_General     )


####save data####
msd$Nation="HN"
write.csv(msd,"E:/GEB data/msd_HN_NaturalForests.csv"  )
msd_5=read.csv("E:/GEB data/msd_HN_NaturalForests.csv")
###### paired-t-test for CAT difference per year
msd$Types[which(msd$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd$Types[which(msd$Types %in% c( "630_General" ))] = c("CBMF_General")
msd$Types[which(msd$Types %in% c("620_General"))] = c("BMF_General")



msd$Types <- factor(msd$Types, levels = c("CP_General","CBMF_General","BMF_General"
))
msd_5$Types[which(msd_5$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd_5$Types[which(msd_5$Types %in% c( "630_General" ))] = c("CBMF_General")
msd_5$Types[which(msd_5$Types %in% c("620_General"))] = c("BMF_General")


msd_5$Types <- factor(msd_5$Types, levels = c("CP_General","CBMF_General","BMF_General" ))

write.csv(msd,"E:/GEB data/Meta_dat1_HN_NaturalForests.csv" )

write.csv(dat,"E:/GEB data/Meta_dat2_HN_NaturalForests.csv" )


##===============================================================================
##  Plot1:Succession ---------
##===============================================================================
rm(list=ls())
library(ggplot2)
####American####
msd=read.csv("E:/GEB data/Meta_dat1_American_NaturalForests.csv" )
unique(msd$Types)
N1=read.csv("E:/GEB data/FIA_Plot_N220_NaturalForests.csv" )
N2=read.csv("E:/GEB data/FIA_Plot_N630_NaturalForests.csv" )
N3=read.csv("E:/GEB data/FIA_Plot_N620_NaturalForests.csv" )
N1$Types="CP_General"
N2$Types="CBMF_General"
N3$Types="BMF_General"
N=rbind(N1,N2,N3)

msd=merge(msd,N,by.x= c("AGE", "Types") ,by.y= c("STDAGE_mean", "Types"))

# Define a mapping between original and new type labels
type_mapping <- list(
  "CP_General" = "Conifer",
  "CBMF_General" = "Conifer–broadleaf mixed",
  "BMF_General" = "Broadleaf"
  
)

# Apply the mapping to the Types column
msd$Types <- sapply(msd$Types, function(x) type_mapping[[x]])


msd$Types = factor(msd$Types, levels = c("Conifer", "Conifer–broadleaf mixed",
                                         "Broadleaf")
                   
)

head(msd)

msd$Significant <- ifelse(msd$Estimate.corr - msd$Estimate.corr.se <= 0 & msd$Estimate.corr + 
                            msd$Estimate.corr.se >= 0, "NS", "Significant")

msd<- subset(msd, Types  %in%  c( "Conifer","Conifer–broadleaf mixed","Broadleaf" ))
msd$CAT = factor(msd$CAT, levels = c("FD_C", "FD_B", "FD_A", 
                                     "CWM_C", "CWM_B", "CWM_A"),
                 labels = c("FD Conifer", "FD Broadleaf",  "FD",
                            "CWM Conifer", "CWM Broadleaf",  "CWM")) 

msd <- msd %>%
  subset(CAT %in% c("FD", "CWM")) %>%
  dplyr::mutate(CAT = dplyr::recode(CAT,
                                    FD = "FD",
                                    CWM = "FI"))







msd <- msd %>%
  mutate(n_discrete = case_when(
    n <= 100              ~ "1 \u2264 n \u2264 100",   
    n > 100 & n <= 500    ~ "100 < n \u2264 500",
    n > 500 & n <= 2000   ~ "500 < n \u2264 2,000",
    n > 2000 & n <= 4000  ~ "2,000 < n \u2264 4,000"
  ))


n_levels_math <- c(
  "1 \u2264 n \u2264 100", 
  "100 < n \u2264 500", 
  "500 < n \u2264 2,000", 
  "2,000 < n \u2264 4,000"
)

msd$n_discrete <- factor(msd$n_discrete, levels = n_levels_math)

msd <- msd %>%
  add_count(Types, CAT, name = "n_points")

head(msd)
p1= ggplot(msd, aes(x = AGE, y = Estimate.corr, color = CAT, fill=CAT)) +
  
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = TRUE, 
              fill = "#dcdcd6", 
              alpha = 0.2) +
  geom_errorbar(data = msd, aes(x = AGE, ymin = Estimate.corr - Estimate.corr.se,
                                ymax = Estimate.corr + Estimate.corr.se, col = CAT),
                size = 1) + 
  # geom_point(data = msd, aes(x = AGE, y = Estimate.corr, shape = Significant, col = CAT),
  #            size = 3.5, alpha = 1) +
  geom_point(aes(size = n_discrete, shape = Significant), 
             alpha = 0.7, stroke = 0.8) +
  
  scale_size_manual(
    values = c(
      "1 \u2264 n \u2264 100"    = 2,
      "100 < n \u2264 500"       = 3,
      "500 < n \u2264 2,000"     = 4,
      "2,000 < n \u2264 4,000"   = 5
    ),
    name = "Sample size (n)",
    drop = FALSE
  ) +
  
  scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) + # Solid for 'sign', hollow for 'nosign'
  
  #facet_wrap(~ Types, ncol = 1, scales = "fixed", axes = "all_x")+
  facet_wrap(~ Types, ncol = 1, scales = "free", strip.position = "right", axes = "all_x") +
  
  scale_color_manual(values = c("FD" = "#D37C5E", 
                                "FI" = "#9DC7C6") )+ 
  
  scale_y_continuous(limits = c(-0.1, 0.2), breaks = seq(-0.1, 0.2, by = 0.1)) +
  
  
  theme_bw( )+
  
  theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank())+
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  ggtitle("(a) Southeastern US")+
  #scale_y_continuous(limits = c(-0.2, 0.4)#, breaks = seq(0, 100, by = 40)
  # ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # 设置 x 轴为1位小数
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  
  # #220 General
  
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 40, y = -0.1, label = "*"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 60, y = -0.1, label = "**"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 82, y = -0.1, label = "**"), size = 5, color = "black")+
  
  # #630 General
  geom_text(data = msd[msd$Types == "Conifer–broadleaf mixed", ],
            aes(x = 21, y = -0.1, label = "**"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer–broadleaf mixed", ],
            aes(x = 41, y = -0.1, label = "***"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer–broadleaf mixed", ],
            aes(x = 60, y = -0.1, label = "***"), size = 5, color = "black")+
  
  
  # #620 General
  
  geom_text(data = msd[msd$Types == "Broadleaf", ],
            aes(x = 22, y = -0.1, label = "**"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Broadleaf", ],
            aes(x = 41, y = -0.1, label = "**"), size = 5, color = "black")+
  
  theme(
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    plot.title = element_text(size=16, color="black",family="serif"),
    # legend.title = element_text(size = 12,family="serif"), 
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16, family = "serif", color="black")  ,
    axis.line = element_line(size = 0.5), 
    axis.ticks = element_line(size = 0.5),
    
    axis.title.y = element_text(size = 16,family="serif", color="black"), 
    axis.text.y = element_text(size = 16,family="serif", color="black"),
    axis.title.x = element_text(size = 16,family="serif", color="black"), 
    axis.text.x = element_text(size = 16,family="serif", color="black"),
    
    strip.background = element_rect(fill = "white",
                                    color = "NA"),  
    strip.text =element_blank(),  #element_text(size = 16, face = "bold",family="serif", color="black"),
    
    legend.position = "right" ) 


p1



#################################################################################
library(scales)
####Japan####

msd=read.csv("E:/GEB data/Meta_dat1_Japan_NaturalForests.csv" )
unique(msd$Types)
N1=read.csv("E:/GEB data/Japan_Plot_N220_NaturalForests.csv" )
N2=read.csv("E:/GEB data/Japan_Plot_N630_NaturalForests.csv" )
N3=read.csv("E:/GEB data/Japan_Plot_N620_NaturalForests.csv" )
N1$Types="CP_General"
N2$Types="CBMF_General"
N3$Types="BMF_General"
N=rbind(N1,N2,N3)

msd=merge(msd,N,by.x= c("AGE", "Types") ,by.y= c("STDAGE_mean", "Types"),all.x=T)
# Define a mapping between original and new type labels
type_mapping <- list(
  "CP_General" = "Conifer",
  "CBMF_General" = "Conifer–broadleaf mixed",
  "BMF_General" = "Broadleaf"
  
)

# Apply the mapping to the Types column
msd$Types <- sapply(msd$Types, function(x) type_mapping[[x]])


msd$Types = factor(msd$Types, levels = c("Conifer", "Conifer–broadleaf mixed",
                                         "Broadleaf")
                   
)


head(msd)

msd$Significant <- ifelse(msd$Estimate.corr - msd$Estimate.corr.se <= 0 & msd$Estimate.corr + 
                            msd$Estimate.corr.se >= 0, "NS", "Significant")

msd<- subset(msd, Types  %in%  c( "Conifer","Conifer–broadleaf mixed","Broadleaf" ))
msd$CAT = factor(msd$CAT, levels = c("FD_C", "FD_B", "FD_A", 
                                     "CWM_C", "CWM_B", "CWM_A"),
                 labels = c("FD Conifer", "FD Broadleaf",  "FD",
                            "CWM Conifer", "CWM Broadleaf",  "CWM")) 

msd <- msd %>%
  subset(CAT %in% c("FD", "CWM")) %>%
  dplyr::mutate(CAT = dplyr::recode(CAT,
                                    FD = "FD",
                                    CWM = "FI"))


msd <- msd %>%
  mutate(n_discrete = case_when(
    n <= 100              ~ "1 \u2264 n \u2264 100",   
    n > 100 & n <= 500    ~ "100 < n \u2264 500",
    n > 500 & n <= 2000   ~ "500 < n \u2264 2,000",
    n > 2000 & n <= 4000  ~ "2,000 < n \u2264 4,000"
  ))


n_levels_math <- c(
  "1 \u2264 n \u2264 100", 
  "100 < n \u2264 500", 
  "500 < n \u2264 2,000", 
  "2,000 < n \u2264 4,000"
)

msd$n_discrete <- factor(msd$n_discrete, levels = n_levels_math)


msd <- msd %>%
  add_count(Types, CAT, name = "n_points")

head(msd)

p2= ggplot(msd, aes(x = AGE, y = Estimate.corr, color = CAT, fill=CAT)) +
  
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = TRUE, 
              fill = "#dcdcd6", 
              alpha = 0.2) +
  geom_errorbar(data = msd, aes(x = AGE, ymin = Estimate.corr - Estimate.corr.se,
                                ymax = Estimate.corr + Estimate.corr.se, col = CAT),
                size = 1) + 
  # geom_point(data = msd, aes(x = AGE, y = Estimate.corr, shape = Significant, col = CAT),
  #            size = 3.5, alpha = 1) +
  
  geom_point(aes(size = n_discrete, shape = Significant), 
             alpha = 0.7, stroke = 0.8) +
  
  scale_size_manual(
    values = c(
      "1 \u2264 n \u2264 100"    = 2,
      "100 < n \u2264 500"       = 3,
      "500 < n \u2264 2,000"     = 4,
      "2,000 < n \u2264 4,000"   = 5
    ),
    name = "Sample size (n)",
    drop = FALSE
  ) +
  
  
  
  scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) + # Solid for 'sign', hollow for 'nosign'
  
  #facet_wrap(~ Types, ncol = 1, scales = "free"#, axes = "all_x" )+
  
  facet_wrap(~ Types, ncol = 1, scales = "free", strip.position = "right", axes = "all_x") +
  
  
  scale_color_manual(values = c("FD" = "#D37C5E", "FI" = "#9DC7C6") )+ 
  
  
  
  scale_x_continuous(limits = c(20, 90)) +
  scale_y_continuous(
    n.breaks = 4,                                  
    labels = label_number(accuracy = 0.1)          
  )+
  
  theme_bw( )+
  
  theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank())+
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  ggtitle("(d) Japan")+
  #scale_y_continuous(limits = c(-0.2, 0.4)#, breaks = seq(0, 100, by = 40)
  # ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # 设置 x 轴为1位小数
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  
  # #220 General
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 43, y = -1, label = "*"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 43, y = 1, label = ""), size = 5, color = "black")+
  theme(
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    plot.title = element_text(size=16, color="black",family="serif"),
    # legend.title = element_text(size = 12,family="serif"), 
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16, family = "serif", color="black")  ,
    axis.line = element_line(size = 0.5), 
    axis.ticks = element_line(size = 0.5),
    
    axis.title.y = element_text(size = 16,family="serif", color="black"), 
    axis.text.y = element_text(size = 16,family="serif", color="black"),
    axis.title.x = element_text(size = 16,family="serif", color="black"), 
    axis.text.x = element_text(size = 16,family="serif", color="black"),
    
    strip.background = element_rect(fill = "white",
                                    color = "NA"),  
    strip.text = element_blank(), #element_text(size = 16, face = "bold",family="serif", color="black"),
    
    legend.position = "bottom" ) 


p2


#################################################################################
####HN####

msd=read.csv("E:/GEB data/Meta_dat1_HN_NaturalForests.csv" )
# Define a mapping between original and new type labels
unique(msd$Types)
N1=read.csv("E:/GEB data/HN_Plot_N220_NaturalForests.csv" )
N2=read.csv("E:/GEB data/HN_Plot_N630_NaturalForests.csv" )
N3=read.csv("E:/GEB data/HN_Plot_N620_NaturalForests.csv" )
N1$Types="CP_General"
N2$Types="CBMF_General"
N3$Types="BMF_General"
N=rbind(N1,N2,N3)

msd=merge(msd,N,by.x= c("AGE", "Types") ,by.y= c("STDAGE_mean", "Types"),all.x=T)




type_mapping <- list(
  "CP_General" = "Conifer",
  "CBMF_General" = "Conifer–broadleaf mixed",
  "BMF_General" = "Broadleaf"
  
)

# Apply the mapping to the Types column
msd$Types <- sapply(msd$Types, function(x) type_mapping[[x]])


msd$Types = factor(msd$Types, levels = c("Conifer", "Conifer–broadleaf mixed",
                                         "Broadleaf")
                   
)


head(msd)

msd$Significant <- ifelse(msd$Estimate.corr - msd$Estimate.corr.se <= 0 & msd$Estimate.corr + 
                            msd$Estimate.corr.se >= 0, "NS", "Significant")

msd<- subset(msd, Types  %in%  c( "Conifer","Conifer–broadleaf mixed","Broadleaf" ))
msd$CAT = factor(msd$CAT, levels = c("FD_C", "FD_B", "FD_A", 
                                     "CWM_C", "CWM_B", "CWM_A"),
                 labels = c("FD Conifer", "FD Broadleaf",  "FD",
                            "CWM Conifer", "CWM Broadleaf",  "CWM")) 

msd <- msd %>%
  subset(CAT %in% c("FD", "CWM")) %>%
  dplyr::mutate(CAT = dplyr::recode(CAT,
                                    FD = "FD",
                                    CWM = "FI"))



msd <- msd %>%
  mutate(n_discrete = case_when(
    n <= 100              ~ "1 \u2264 n \u2264 100",   
    n > 100 & n <= 500    ~ "100 < n \u2264 500",
    n > 500 & n <= 2000   ~ "500 < n \u2264 2,000",
    n > 2000 & n <= 4000  ~ "2,000 < n \u2264 4,000"
  ))


n_levels_math <- c(
  "1 \u2264 n \u2264 100", 
  "100 < n \u2264 500", 
  "500 < n \u2264 2,000", 
  "2,000 < n \u2264 4,000"
)

msd$n_discrete <- factor(msd$n_discrete, levels = n_levels_math)

msd <- msd %>%
  add_count(Types, CAT, name = "n_points")

head(msd)
p3= ggplot(msd, aes(x = AGE, y = Estimate.corr, color = CAT, fill=CAT)) +
  
  geom_smooth(data = subset(msd, n_points >= 3), 
              method = "lm", 
              formula = y ~ x, 
              se = TRUE, 
              fill = "#dcdcd6", 
              alpha = 0.2) +
  geom_errorbar(data = msd, aes(x = AGE, ymin = Estimate.corr - Estimate.corr.se,
                                ymax = Estimate.corr + Estimate.corr.se, col = CAT),
                size = 1) + 
  # geom_point(data = msd, aes(x = AGE, y = Estimate.corr, shape = Significant, col = CAT),
  #            size = 3.5, alpha = 1) +
  facet_wrap(~ Types, ncol = 1, scales = "free", strip.position = "right", axes = "all_x") +
  
  geom_point(aes(size = n_discrete, shape = Significant), 
             alpha = 0.7, stroke = 0.8) +
  
  scale_size_manual(
    values = c(
      "1 \u2264 n \u2264 100"    = 2,
      "100 < n \u2264 500"       = 3,
      "500 < n \u2264 2,000"     = 4,
      "2,000 < n \u2264 4,000"   = 5
    ),
    name = "Sample size (n)",
    drop = FALSE
  ) +
  
  scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) + # Solid for 'sign', hollow for 'nosign'
  
  # facet_wrap(~ Types, ncol = 1, scales = "free", axes = "all_x")+
  
  
  scale_color_manual(values = c("FD" = "#D37C5E",   "FI" = "#9DC7C6") )+ 
  
  
  scale_y_continuous(
    n.breaks = 4,                                  
    labels = label_number(accuracy = 0.1)          
  )+
  
  
  theme_bw( )+
  
  theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank())+
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  ggtitle("(g) Hunan, China")+
  
  scale_x_continuous(limits = c(6, 40)   ) +
  
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # 设置 x 轴为1位小数
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  
  # #220 General
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 37, y = -0.1, label = "*"), size = 5, color = "black")+
  
  # #630 General
  geom_text(data = msd[msd$Types == "Conifer–broadleaf mixed", ], 
            aes(x = 20, y = -0.5, label = "**"), size = 5, color = "black")+
  
  
  # # #620 General
  geom_text(data = msd[msd$Types == "Broadleaf", ], 
            aes(x = 39, y = -0.1, label = "*"), size = 5, color = "black")+
  # geom_text(data = msd[msd$Types == "Broadleaf", ], 
  #           aes(x = 40, y = -0.8, label = "*"), size = 5, color = "black")+
  # # geom_text(data = msd[msd$Types == "Broadleaf", ], 
  # #           aes(x = 21, y = -0.25, label = "***"), size = 5, color = "black")+
  # # 
  
  theme(
    panel.spacing.y = unit(2, "lines"),
    
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    plot.title = element_text(size=16, color="black",family="serif"),
    # legend.title = element_text(size = 12,family="serif"), 
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16, family = "serif", color="black")  ,
    axis.line = element_line(size = 0.5), 
    axis.ticks = element_line(size = 0.5),
    
    axis.title.y = element_text(size = 16,family="serif", color="black"), 
    axis.text.y = element_text(size = 16,family="serif", color="black"),
    axis.title.x = element_text(size = 16,family="serif", color="black"), 
    axis.text.x = element_text(size = 16,family="serif", color="black"),
    
    strip.background = element_rect(fill = "white",
                                    color = "NA"),  
    strip.text =element_blank(),  #element_text(size = 16, face = "bold",family="serif", color="black"),
    
    legend.position = "bottom" ) 


p3


#####Plot merge#############
library(stats)
library(ggpubr)


p11=p1+theme(axis.title.y = element_blank())
p21=p2+theme(axis.title.y = element_blank())
p31=p3+theme(axis.title.y = element_blank())



fig1 = ggarrange(p11, p21, p31,
                 ncol = 3, nrow = 1,
                 common.legend = TRUE, 
                 legend = "bottom",
                 align = "hv")
fig1

fig1_final <- annotate_figure(fig1,
                              left = text_grob(expression(paste( "Trait effect on AGC (slope)")), 
                                               rot = 90, vjust = 1, size = 16, family = "serif"))


fig1_final <- fig1_final + theme(plot.background = element_rect(fill = "white", color = NA))

print(fig1_final)



ggsave("Succession_NaturalForests_R.png", path = "E:/GEB data/figure",width =8, height =9,
       units = "in",dpi=600, plot=fig1_final)

##===============================================================================
## American:Plantations ---------
##===============================================================================
rm(list=ls())
#### meta data ####
meta.data_220_FD <-read.csv("E:/GEB data/LME_Detail_American_220_Plantations.csv")
meta.data_220_CWM <-read.csv("E:/GEB data/LME_Detail_American_220_Plantations_CWM.csv")
meta_data_220_FD <- meta.data_220_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_220_CWM <- meta.data_220_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_220=rbind(meta_data_220_FD ,meta_data_220_CWM )

meta.data_630_FD <-read.csv("E:/GEB data/LME_Detail_American_630_Plantations.csv")
meta.data_630_CWM <-read.csv("E:/GEB data/LME_Detail_American_630_Plantations_CWM.csv")
meta_data_630_FD <- meta.data_630_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_630_CWM <- meta.data_630_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_630=rbind(meta_data_630_FD ,meta_data_630_CWM )

meta.data_620_FD <-read.csv("E:/GEB data/LME_Detail_American_620_Plantations.csv")
meta.data_620_CWM <-read.csv("E:/GEB data/LME_Detail_American_620_Plantations_CWM.csv")
meta_data_620_FD <- meta.data_620_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_620_CWM <- meta.data_620_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_620=rbind(meta_data_620_FD ,meta_data_620_CWM )


meta.data_220$Types <-"220_General"
meta.data_630$Types <-"630_General"
meta.data_620$Types <-"620_General"


meta.data=rbind(meta.data_220,meta.data_630,meta.data_620
                
)

str(meta.data)
meta.data$RESPONSE="Carbon_Mg_ha"
meta.data$trait=meta.data$factor
meta.data$CAT <- ifelse(grepl("\\.CWM$", meta.data$factor), "CWM_A", 
                        ifelse(grepl("\\.FDis$", meta.data$factor), "FD_A", NA))
meta.data =meta.data [,c( "Value",     "Std.Error" ,"DF" ,       "t.value" ,  "p.value" ,  "R2m" ,
                          "R2c"  ,     "factor",    "N_obs" ,  "RESPONSE", "CAT" , "TRAIT_ID" , "AgeGroup" , "Types"   )]

colnames(meta.data )=c( "Estimate" ,  "Std..Error", "df","t.value" , "Pr...t.." ,  
                        "R2m" ,"R2c",  "factor" ,"Nptag",  "RESPONSE", "CAT" , "TRAIT" , "AgeGroup" , "Types" )




t.cor<-c(  
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,"LDMC"  ,  
  "WD_mean_add",          
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"  
  
)

#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate,
                                meta.data$Estimate)

### subset data for analyses
dat<-droplevels(subset(meta.data, RESPONSE=="Carbon_Mg_ha")) ## accumulated stand volume
dat$var=rownames(dat)


unique(dat)
dat <- dat %>%
  dplyr::mutate(
    AGE = case_when(
      Types == "220_General" & AgeGroup == 1 ~ "8",  
      Types == "220_General" & AgeGroup == 2 ~ "21", 
      Types == "220_General" & AgeGroup == 3 ~ "37",
      Types == "220_General" & AgeGroup == 4 ~ "59",
      Types == "220_General" & AgeGroup == 5 ~ "79", 
      
      Types == "630_General" & AgeGroup == 1 ~ "6",  
      Types == "630_General" & AgeGroup == 2 ~ "19", 
      Types == "630_General" & AgeGroup == 3 ~ "39", 
      Types == "630_General" & AgeGroup == 4 ~ "58", 
      Types == "630_General" & AgeGroup == 5 ~ "75", 
      
      
      Types == "620_General" & AgeGroup == 1 ~ "7",  
      Types == "620_General" & AgeGroup == 2 ~ "18", 
      Types == "620_General" & AgeGroup == 3 ~ "39", Types == "620_General" & AgeGroup == 4 ~ "59", 
      Types == "620_General" & AgeGroup == 5 ~ "79",
      TRUE ~ as.character(NA)  # For other rows, assign NA or other default value
    )
  )


## center AGE 
dat$AGE2 <- scale(as.numeric(dat$AGE)  )
#dat$AGE2 <- dat$AGE- mean(dat$AGE, na.rm = TRUE)  
dat$AGE2.sq<-dat$AGE2^2
##########
##General
dat_220_General=dat[dat$Types == "220_General", ]
# m2_220_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "220_General", ])

dat_220_General_FD_A <- subset(dat_220_General, CAT == "FD_A")
dat_220_General_CWM_A <- subset(dat_220_General, CAT == "CWM_A")

library(mgcv)
m2_220_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_FD_A)
summary(m2_220_General_FD_A)
r.squaredGLMM(m2_220_General_FD_A)
m2_220_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_CWM_A)
summary(m2_220_General_CWM_A )
r.squaredGLMM(m2_220_General_CWM_A)

dat_220_General_FD_A$predict.Est.corr.m2 <- predict(m2_220_General_FD_A)
dat_220_General_CWM_A$predict.Est.corr.m2 <- predict(m2_220_General_CWM_A)



dat_220_General_pred <- rbind(
  dat_220_General_FD_A, dat_220_General_CWM_A )

dat_220_General_FD_A$predict.Est.corr.m2<-predict(m2_220_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=mean)
sd_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=sd)
n_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=length)
sd_220_General_FD_A$Estimate.corr.se<-sd_220_General_FD_A$Estimate.corr/sqrt(n_220_General_FD_A$Estimate.corr)
msd_220_General_FD_A<-merge(m_220_General_FD_A,sd_220_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

dat_220_General_CWM_A$predict.Est.corr.m2<-predict(m2_220_General_CWM_A)



# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=mean)
sd_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=sd)
n_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=length)
sd_220_General_CWM_A$Estimate.corr.se<-sd_220_General_CWM_A$Estimate.corr/sqrt(n_220_General_CWM_A$Estimate.corr)
msd_220_General_CWM_A<-merge(m_220_General_CWM_A,sd_220_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_220_General=rbind(
  msd_220_General_FD_A, msd_220_General_CWM_A
)

dat_630_General=dat[dat$Types == "630_General", ]


dat_630_General_FD_A <- subset(dat_630_General, CAT == "FD_A")
dat_630_General_CWM_A <- subset(dat_630_General, CAT == "CWM_A")




m2_630_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_FD_A)
summary(m2_630_General_FD_A)
r.squaredGLMM(m2_630_General_FD_A)
m2_630_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_CWM_A)
summary(m2_630_General_CWM_A )
r.squaredGLMM(m2_630_General_CWM_A)


dat_630_General_FD_A$predict.Est.corr.m2 <- predict(m2_630_General_FD_A)
dat_630_General_CWM_A$predict.Est.corr.m2 <- predict(m2_630_General_CWM_A)

dat_630_General_pred <- rbind(
  dat_630_General_FD_A, dat_630_General_CWM_A )

dat_630_General_FD_A$predict.Est.corr.m2<-predict(m2_630_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=mean)
sd_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=sd)
n_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=length)
sd_630_General_FD_A$Estimate.corr.se<-sd_630_General_FD_A$Estimate.corr/sqrt(n_630_General_FD_A$Estimate.corr)
msd_630_General_FD_A<-merge(m_630_General_FD_A,sd_630_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_630_General_CWM_A$predict.Est.corr.m2<-predict(m2_630_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=mean)
sd_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=sd)
n_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=length)
sd_630_General_CWM_A$Estimate.corr.se<-sd_630_General_CWM_A$Estimate.corr/sqrt(n_630_General_CWM_A$Estimate.corr)
msd_630_General_CWM_A<-merge(m_630_General_CWM_A,sd_630_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_630_General=rbind(
  msd_630_General_FD_A, msd_630_General_CWM_A
)

##General
dat_620_General=dat[dat$Types == "620_General", ]
# m2_620_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "620_General", ])

dat_620_General_FD_A <- subset(dat_620_General, CAT == "FD_A")
dat_620_General_CWM_A <- subset(dat_620_General, CAT == "CWM_A")

m2_620_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_FD_A)
summary(m2_620_General_FD_A)
r.squaredGLMM(m2_620_General_FD_A)
m2_620_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_CWM_A)
summary(m2_620_General_CWM_A )
r.squaredGLMM(m2_620_General_CWM_A)


dat_620_General_FD_A$predict.Est.corr.m2 <- predict(m2_620_General_FD_A)
dat_620_General_CWM_A$predict.Est.corr.m2 <- predict(m2_620_General_CWM_A)

dat_620_General_pred <- rbind(
  dat_620_General_FD_A, dat_620_General_CWM_A )

dat_620_General_FD_A$predict.Est.corr.m2<-predict(m2_620_General_FD_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=mean)
sd_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=sd)
n_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=length)
sd_620_General_FD_A$Estimate.corr.se<-sd_620_General_FD_A$Estimate.corr/sqrt(n_620_General_FD_A$Estimate.corr)
msd_620_General_FD_A<-merge(m_620_General_FD_A,sd_620_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_620_General_CWM_A$predict.Est.corr.m2<-predict(m2_620_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=mean)
sd_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=sd)
n_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=length)
sd_620_General_CWM_A$Estimate.corr.se<-sd_620_General_CWM_A$Estimate.corr/sqrt(n_620_General_CWM_A$Estimate.corr)
msd_620_General_CWM_A<-merge(m_620_General_CWM_A,sd_620_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_620_General=rbind(
  msd_620_General_FD_A, msd_620_General_CWM_A
)



msd_220_General$Types="220_General"
msd_630_General$Types="630_General"
msd_620_General$Types="620_General"

msd=rbind(msd_220_General,msd_630_General,  msd_620_General     )


####save data####
msd$Nation="American"
write.csv(msd,"E:/GEB data/msd_American_Plantations.csv"  )
msd_5=read.csv("E:/GEB data/msd_American_Plantations.csv")
###### paired-t-test for CAT difference per year
msd$Types[which(msd$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd$Types[which(msd$Types %in% c( "630_General" ))] = c("CBMF_General")
msd$Types[which(msd$Types %in% c("620_General"))] = c("BMF_General")



msd$Types <- factor(msd$Types, levels = c("CP_General","CBMF_General","BMF_General"
))
msd_5$Types[which(msd_5$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd_5$Types[which(msd_5$Types %in% c( "630_General" ))] = c("CBMF_General")
msd_5$Types[which(msd_5$Types %in% c("620_General"))] = c("BMF_General")


msd_5$Types <- factor(msd_5$Types, levels = c("CP_General","CBMF_General","BMF_General" ))

write.csv(msd,"E:/GEB data/Meta_dat1_American_Plantations.csv" )

write.csv(dat,"E:/GEB data/Meta_dat2_American_Plantations.csv" )

##===============================================================================
## Japan:Plantations ---------
##===============================================================================
rm(list=ls())
#### meta data ####
meta.data_220_FD <-read.csv("E:/GEB data/LME_Detail_Japan_220_Plantations.csv")
meta.data_220_CWM <-read.csv("E:/GEB data/LME_Detail_Japan_220_Plantations_CWM.csv")
meta_data_220_FD <- meta.data_220_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_220_CWM <- meta.data_220_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_220=rbind(meta_data_220_FD ,meta_data_220_CWM )

meta.data_630_FD <-read.csv("E:/GEB data/LME_Detail_Japan_630_Plantations.csv")
meta.data_630_CWM <-read.csv("E:/GEB data/LME_Detail_Japan_630_Plantations_CWM.csv")
meta_data_630_FD <- meta.data_630_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_630_CWM <- meta.data_630_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_630=rbind(meta_data_630_FD ,meta_data_630_CWM )

meta.data_620_FD <-read.csv("E:/GEB data/LME_Detail_Japan_620_Plantations.csv")
meta.data_620_CWM <-read.csv("E:/GEB data/LME_Detail_Japan_620_Plantations_CWM.csv")
meta_data_620_FD <- meta.data_620_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_620_CWM <- meta.data_620_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_620=rbind(meta_data_620_FD ,meta_data_620_CWM )


meta.data_220$Types <-"220_General"
meta.data_630$Types <-"630_General"
meta.data_620$Types <-"620_General"


meta.data=rbind(meta.data_220,meta.data_630,meta.data_620
                
)

str(meta.data)
meta.data$RESPONSE="Carbon_Mg_ha"
meta.data$trait=meta.data$factor
meta.data$CAT <- ifelse(grepl("\\.CWM$", meta.data$factor), "CWM_A", 
                        ifelse(grepl("\\.FDis$", meta.data$factor), "FD_A", NA))
meta.data =meta.data [,c( "Value",     "Std.Error" ,"DF" ,       "t.value" ,  "p.value" ,  "R2m" ,
                          "R2c"  ,     "factor",    "N_obs" ,  "RESPONSE", "CAT" , "TRAIT_ID" , "AgeGroup" , "Types"   )]

colnames(meta.data )=c( "Estimate" ,  "Std..Error", "df","t.value" , "Pr...t.." ,  
                        "R2m" ,"R2c",  "factor" ,"Nptag",  "RESPONSE", "CAT" , "TRAIT" , "AgeGroup" , "Types" )




t.cor<-c(  
  "Leaf", 
  "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,"LDMC"  ,  
  "WD_mean_add",          
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"  
  
)

#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate,
                                meta.data$Estimate)

### subset data for analyses
dat<-droplevels(subset(meta.data, RESPONSE=="Carbon_Mg_ha")) ## accumulated stand volume
dat$var=rownames(dat)


unique(dat)
dat <- dat %>%
  dplyr::mutate(
    AGE = case_when(
      Types == "220_General" & AgeGroup == 1 ~ "7",  
      Types == "220_General" & AgeGroup == 2 ~ "24", 
      Types == "220_General" & AgeGroup == 3 ~ "43",
      Types == "220_General" & AgeGroup == 4 ~ "58",
      Types == "220_General" & AgeGroup == 5 ~ "84", 
      
      Types == "630_General" & AgeGroup == 1 ~ "9",  
      Types == "630_General" & AgeGroup == 2 ~ "24", 
      Types == "630_General" & AgeGroup == 3 ~ "42", 
      Types == "630_General" & AgeGroup == 4 ~ "59", 
      Types == "630_General" & AgeGroup == 5 ~ "85", 
      
      
      Types == "620_General" & AgeGroup == 1 ~ "7",  
      Types == "620_General" & AgeGroup == 2 ~ "24", 
      Types == "620_General" & AgeGroup == 3 ~ "42", Types == "620_General" & AgeGroup == 4 ~ "57", 
      Types == "620_General" & AgeGroup == 5 ~ "82",
      TRUE ~ as.character(NA)  # For other rows, assign NA or other default value
    )
  )


## center AGE 
dat$AGE2 <- scale(as.numeric(dat$AGE)  )
#dat$AGE2 <- dat$AGE- mean(dat$AGE, na.rm = TRUE)  
dat$AGE2.sq<-as.numeric(dat$AGE)^2
##########
##General
dat_220_General=dat[dat$Types == "220_General", ]
# m2_220_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "220_General", ])

dat_220_General_FD_A <- subset(dat_220_General, CAT == "FD_A")
dat_220_General_CWM_A <- subset(dat_220_General, CAT == "CWM_A")

library(mgcv)

m2_220_General_FD_A <- lm(
  Estimate.corr ~ AGE2 ,
  data = dat_220_General_FD_A
)
summary(m2_220_General_FD_A)
r.squaredGLMM(m2_220_General_FD_A)
m2_220_General_CWM_A <-lm(
  Estimate.corr ~ AGE2 , data = dat_220_General_CWM_A)
summary(m2_220_General_CWM_A )
r.squaredGLMM(m2_220_General_CWM_A)

dat_220_General_FD_A$predict.Est.corr.m2 <- predict(m2_220_General_FD_A)
dat_220_General_CWM_A$predict.Est.corr.m2 <- predict(m2_220_General_CWM_A)



dat_220_General_pred <- rbind(
  dat_220_General_FD_A, dat_220_General_CWM_A )

dat_220_General_FD_A$predict.Est.corr.m2<-predict(m2_220_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=mean)
sd_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=sd)
n_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=length)
sd_220_General_FD_A$Estimate.corr.se<-sd_220_General_FD_A$Estimate.corr/sqrt(n_220_General_FD_A$Estimate.corr)
msd_220_General_FD_A<-merge(m_220_General_FD_A,sd_220_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

dat_220_General_CWM_A$predict.Est.corr.m2<-predict(m2_220_General_CWM_A)



# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=mean)
sd_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=sd)
n_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=length)
sd_220_General_CWM_A$Estimate.corr.se<-sd_220_General_CWM_A$Estimate.corr/sqrt(n_220_General_CWM_A$Estimate.corr)
msd_220_General_CWM_A<-merge(m_220_General_CWM_A,sd_220_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_220_General=rbind(
  msd_220_General_FD_A, msd_220_General_CWM_A
)

dat_630_General=dat[dat$Types == "630_General", ]


dat_630_General_FD_A <- subset(dat_630_General, CAT == "FD_A")
dat_630_General_CWM_A <- subset(dat_630_General, CAT == "CWM_A")




m2_630_General_FD_A <- lm(
  Estimate.corr ~ AGE2 ,
  data = dat_630_General_FD_A
)
summary(m2_630_General_FD_A)
r.squaredGLMM(m2_630_General_FD_A)
m2_630_General_CWM_A <-lm(
  Estimate.corr ~ AGE2 , data = dat_630_General_CWM_A)
summary(m2_630_General_CWM_A )
r.squaredGLMM(m2_630_General_CWM_A)


dat_630_General_FD_A$predict.Est.corr.m2 <- predict(m2_630_General_FD_A)
dat_630_General_CWM_A$predict.Est.corr.m2 <- predict(m2_630_General_CWM_A)

dat_630_General_pred <- rbind(
  dat_630_General_FD_A, dat_630_General_CWM_A )

dat_630_General_FD_A$predict.Est.corr.m2<-predict(m2_630_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=mean)
sd_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=sd)
n_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=length)
sd_630_General_FD_A$Estimate.corr.se<-sd_630_General_FD_A$Estimate.corr/sqrt(n_630_General_FD_A$Estimate.corr)
msd_630_General_FD_A<-merge(m_630_General_FD_A,sd_630_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_630_General_CWM_A$predict.Est.corr.m2<-predict(m2_630_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=mean)
sd_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=sd)
n_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=length)
sd_630_General_CWM_A$Estimate.corr.se<-sd_630_General_CWM_A$Estimate.corr/sqrt(n_630_General_CWM_A$Estimate.corr)
msd_630_General_CWM_A<-merge(m_630_General_CWM_A,sd_630_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_630_General=rbind(
  msd_630_General_FD_A, msd_630_General_CWM_A
)

##General
dat_620_General=dat[dat$Types == "620_General", ]
# m2_620_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "620_General", ])

dat_620_General_FD_A <- subset(dat_620_General, CAT == "FD_A")
dat_620_General_CWM_A <- subset(dat_620_General, CAT == "CWM_A")

m2_620_General_FD_A <- lm(
  Estimate.corr ~ AGE2 ,
  data = dat_620_General_FD_A
)
summary(m2_620_General_FD_A)
r.squaredGLMM(m2_620_General_FD_A)
m2_620_General_CWM_A <-lm(
  Estimate.corr ~ AGE2 , data = dat_620_General_CWM_A)
summary(m2_620_General_CWM_A )
r.squaredGLMM(m2_620_General_CWM_A)


dat_620_General_FD_A$predict.Est.corr.m2 <- predict(m2_620_General_FD_A)
dat_620_General_CWM_A$predict.Est.corr.m2 <- predict(m2_620_General_CWM_A)

dat_620_General_pred <- rbind(
  dat_620_General_FD_A, dat_620_General_CWM_A )

dat_620_General_FD_A$predict.Est.corr.m2<-predict(m2_620_General_FD_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=mean)
sd_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=sd)
n_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=length)
sd_620_General_FD_A$Estimate.corr.se<-sd_620_General_FD_A$Estimate.corr/sqrt(n_620_General_FD_A$Estimate.corr)
msd_620_General_FD_A<-merge(m_620_General_FD_A,sd_620_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_620_General_CWM_A$predict.Est.corr.m2<-predict(m2_620_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=mean)
sd_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=sd)
n_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=length)
sd_620_General_CWM_A$Estimate.corr.se<-sd_620_General_CWM_A$Estimate.corr/sqrt(n_620_General_CWM_A$Estimate.corr)
msd_620_General_CWM_A<-merge(m_620_General_CWM_A,sd_620_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_620_General=rbind(
  msd_620_General_FD_A, msd_620_General_CWM_A
)



msd_220_General$Types="220_General"
msd_630_General$Types="630_General"
msd_620_General$Types="620_General"

msd=rbind(msd_220_General,msd_630_General,  msd_620_General     )


####save data####
msd$Nation="Japan"
write.csv(msd,"E:/GEB data/msd_Japan_Plantations.csv"  )
msd_5=read.csv("E:/GEB data/msd_Japan_Plantations.csv")
###### paired-t-test for CAT difference per year
msd$Types[which(msd$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd$Types[which(msd$Types %in% c( "630_General" ))] = c("CBMF_General")
msd$Types[which(msd$Types %in% c("620_General"))] = c("BMF_General")



msd$Types <- factor(msd$Types, levels = c("CP_General","CBMF_General","BMF_General"
))
msd_5$Types[which(msd_5$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd_5$Types[which(msd_5$Types %in% c( "630_General" ))] = c("CBMF_General")
msd_5$Types[which(msd_5$Types %in% c("620_General"))] = c("BMF_General")


msd_5$Types <- factor(msd_5$Types, levels = c("CP_General","CBMF_General","BMF_General" ))

write.csv(msd,"E:/GEB data/Meta_dat1_Japan_Plantations.csv" )

write.csv(dat,"E:/GEB data/Meta_dat2_Japan_Plantations.csv" )

##===============================================================================
## HN:Plantations ---------
##===============================================================================
rm(list=ls())
#### meta data ####
meta.data_220_FD <-read.csv("E:/GEB data/LME_Detail_HN_220_Plantations.csv")
meta.data_220_CWM <-read.csv("E:/GEB data/LME_Detail_HN_220_Plantations_CWM.csv")
meta_data_220_FD <- meta.data_220_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_220_CWM <- meta.data_220_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_220=rbind(meta_data_220_FD ,meta_data_220_CWM )

meta.data_630_FD <-read.csv("E:/GEB data/LME_Detail_HN_630_Plantations.csv")
meta.data_630_CWM <-read.csv("E:/GEB data/LME_Detail_HN_630_Plantations_CWM.csv")
meta_data_630_FD <- meta.data_630_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_630_CWM <- meta.data_630_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_630=rbind(meta_data_630_FD ,meta_data_630_CWM )

meta.data_620_FD <-read.csv("E:/GEB data/LME_Detail_HN_620_Plantations.csv")
meta.data_620_CWM <-read.csv("E:/GEB data/LME_Detail_HN_620_Plantations_CWM.csv")
meta_data_620_FD <- meta.data_620_FD %>%
  filter(stringr::str_detect(factor, "\\.FDis")) %>%
  mutate(factor = droplevels(as.factor(factor)))
meta_data_620_CWM <- meta.data_620_CWM %>%
  filter(stringr::str_detect(factor, "\\.CWM")) %>%
  mutate(factor = droplevels(as.factor(factor)))

meta.data_620=rbind(meta_data_620_FD ,meta_data_620_CWM )


meta.data_220$Types <-"220_General"
meta.data_630$Types <-"630_General"
meta.data_620$Types <-"620_General"


meta.data=rbind(meta.data_220,meta.data_630,meta.data_620
                
)

str(meta.data)
meta.data$RESPONSE="Carbon_Mg_ha"
meta.data$trait=meta.data$factor
meta.data$CAT <- ifelse(grepl("\\.CWM$", meta.data$factor), "CWM_A", 
                        ifelse(grepl("\\.FDis$", meta.data$factor), "FD_A", NA))
meta.data =meta.data [,c( "Value",     "Std.Error" ,"DF" ,       "t.value" ,  "p.value" ,  "R2m" ,
                          "R2c"  ,     "factor",    "N_obs" ,  "RESPONSE", "CAT" , "TRAIT_ID" , "AgeGroup" , "Types"   )]

colnames(meta.data )=c( "Estimate" ,  "Std..Error", "df","t.value" , "Pr...t.." ,  
                        "R2m" ,"R2c",  "factor" ,"Nptag",  "RESPONSE", "CAT" , "TRAIT" , "AgeGroup" , "Types" )




t.cor<-c(  
  "Leaf", "Max.height","Shade.tolerance" ,"Drought.tolerance"  ,"LDMC"  ,  
  "WD_mean_add",          
  "LMA", "Leaf.longevity",
  
  "Nmass" , "Pmass"  
  
)

#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate,
                                meta.data$Estimate)

### subset data for analyses
dat<-droplevels(subset(meta.data, RESPONSE=="Carbon_Mg_ha")) ## accumulated stand volume
dat$var=rownames(dat)


unique(dat)
dat <- dat %>%
  dplyr::mutate(
    AGE = case_when(
      Types == "220_General" & AgeGroup == 1 ~ "7",  
      Types == "220_General" & AgeGroup == 2 ~ "20", 
      Types == "220_General" & AgeGroup == 3 ~ "37",
      Types == "220_General" & AgeGroup == 4 ~ "52",
      Types == "220_General" & AgeGroup == 5 ~ "85", 
      
      Types == "630_General" & AgeGroup == 1 ~ "7",  
      Types == "630_General" & AgeGroup == 2 ~ "19", 
      Types == "630_General" & AgeGroup == 3 ~ "35", 
      Types == "630_General" & AgeGroup == 4 ~ "60", 
      Types == "630_General" & AgeGroup == 5 ~ "83", 
      
      
      Types == "620_General" & AgeGroup == 1 ~ "7",  
      Types == "620_General" & AgeGroup == 2 ~ "19", 
      Types == "620_General" & AgeGroup == 3 ~ "34", Types == "620_General" & AgeGroup == 4 ~ "61", 
      Types == "620_General" & AgeGroup == 5 ~ "85",
      TRUE ~ as.character(NA)  # For other rows, assign NA or other default value
    )
  )




## center AGE 
dat$AGE2 <- scale(as.numeric(dat$AGE)  )
#dat$AGE2 <- dat$AGE- mean(dat$AGE, na.rm = TRUE)  
dat$AGE2.sq<-dat$AGE2^2
##########
##General
dat_220_General=dat[dat$Types == "220_General", ]
# m2_220_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "220_General", ])

dat_220_General_FD_A <- subset(dat_220_General, CAT == "FD_A")
dat_220_General_CWM_A <- subset(dat_220_General, CAT == "CWM_A")

library(mgcv)
m2_220_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_FD_A)
summary(m2_220_General_FD_A)
r.squaredGLMM(m2_220_General_FD_A)
m2_220_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_220_General_CWM_A)
summary(m2_220_General_CWM_A )
r.squaredGLMM(m2_220_General_CWM_A)

dat_220_General_FD_A$predict.Est.corr.m2 <- predict(m2_220_General_FD_A)
dat_220_General_CWM_A$predict.Est.corr.m2 <- predict(m2_220_General_CWM_A)



dat_220_General_pred <- rbind(
  dat_220_General_FD_A, dat_220_General_CWM_A )

dat_220_General_FD_A$predict.Est.corr.m2<-predict(m2_220_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=mean)
sd_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=sd)
n_220_General_FD_A<-aggregate(dat_220_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_FD_A$AGE,CAT=dat_220_General_FD_A$CAT),FUN=length)
sd_220_General_FD_A$Estimate.corr.se<-sd_220_General_FD_A$Estimate.corr/sqrt(n_220_General_FD_A$Estimate.corr)
msd_220_General_FD_A<-merge(m_220_General_FD_A,sd_220_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

dat_220_General_CWM_A$predict.Est.corr.m2<-predict(m2_220_General_CWM_A)



# calculate the mean value to include in the figure as point and get values for in the text
m_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=mean)
sd_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=sd)
n_220_General_CWM_A<-aggregate(dat_220_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_220_General_CWM_A$AGE,CAT=dat_220_General_CWM_A$CAT),FUN=length)
sd_220_General_CWM_A$Estimate.corr.se<-sd_220_General_CWM_A$Estimate.corr/sqrt(n_220_General_CWM_A$Estimate.corr)
msd_220_General_CWM_A<-merge(m_220_General_CWM_A,sd_220_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_220_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_220_General=rbind(
  msd_220_General_FD_A, msd_220_General_CWM_A
)

dat_630_General=dat[dat$Types == "630_General", ]


dat_630_General_FD_A <- subset(dat_630_General, CAT == "FD_A")
dat_630_General_CWM_A <- subset(dat_630_General, CAT == "CWM_A")




m2_630_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_FD_A)
summary(m2_630_General_FD_A)
r.squaredGLMM(m2_630_General_FD_A)
m2_630_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_630_General_CWM_A)
summary(m2_630_General_CWM_A )
r.squaredGLMM(m2_630_General_CWM_A)


dat_630_General_FD_A$predict.Est.corr.m2 <- predict(m2_630_General_FD_A)
dat_630_General_CWM_A$predict.Est.corr.m2 <- predict(m2_630_General_CWM_A)

dat_630_General_pred <- rbind(
  dat_630_General_FD_A, dat_630_General_CWM_A )

dat_630_General_FD_A$predict.Est.corr.m2<-predict(m2_630_General_FD_A)

# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=mean)
sd_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=sd)
n_630_General_FD_A<-aggregate(dat_630_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_FD_A$AGE,CAT=dat_630_General_FD_A$CAT),FUN=length)
sd_630_General_FD_A$Estimate.corr.se<-sd_630_General_FD_A$Estimate.corr/sqrt(n_630_General_FD_A$Estimate.corr)
msd_630_General_FD_A<-merge(m_630_General_FD_A,sd_630_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_630_General_CWM_A$predict.Est.corr.m2<-predict(m2_630_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=mean)
sd_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=sd)
n_630_General_CWM_A<-aggregate(dat_630_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_630_General_CWM_A$AGE,CAT=dat_630_General_CWM_A$CAT),FUN=length)
sd_630_General_CWM_A$Estimate.corr.se<-sd_630_General_CWM_A$Estimate.corr/sqrt(n_630_General_CWM_A$Estimate.corr)
msd_630_General_CWM_A<-merge(m_630_General_CWM_A,sd_630_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_630_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_630_General=rbind(
  msd_630_General_FD_A, msd_630_General_CWM_A
)

##General
dat_620_General=dat[dat$Types == "620_General", ]
# m2_620_General<-lmer(Estimate.corr~CAT*(AGE2#+AGE2.sq)
#                      +(1|CAT:TRAIT), data = dat[dat$Types == "620_General", ])

dat_620_General_FD_A <- subset(dat_620_General, CAT == "FD_A")
dat_620_General_CWM_A <- subset(dat_620_General, CAT == "CWM_A")

m2_620_General_FD_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_FD_A)
summary(m2_620_General_FD_A)
r.squaredGLMM(m2_620_General_FD_A)
m2_620_General_CWM_A <- lm(Estimate.corr ~ AGE2, data = dat_620_General_CWM_A)
summary(m2_620_General_CWM_A )
r.squaredGLMM(m2_620_General_CWM_A)


dat_620_General_FD_A$predict.Est.corr.m2 <- predict(m2_620_General_FD_A)
dat_620_General_CWM_A$predict.Est.corr.m2 <- predict(m2_620_General_CWM_A)

dat_620_General_pred <- rbind(
  dat_620_General_FD_A, dat_620_General_CWM_A )

dat_620_General_FD_A$predict.Est.corr.m2<-predict(m2_620_General_FD_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=mean)
sd_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=sd)
n_620_General_FD_A<-aggregate(dat_620_General_FD_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_FD_A$AGE,CAT=dat_620_General_FD_A$CAT),FUN=length)
sd_620_General_FD_A$Estimate.corr.se<-sd_620_General_FD_A$Estimate.corr/sqrt(n_620_General_FD_A$Estimate.corr)
msd_620_General_FD_A<-merge(m_620_General_FD_A,sd_620_General_FD_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_FD_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


dat_620_General_CWM_A$predict.Est.corr.m2<-predict(m2_620_General_CWM_A)


# calculate the mean value to include in the figure as point and get values for in the text
m_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=mean)
sd_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=sd)
n_620_General_CWM_A<-aggregate(dat_620_General_CWM_A[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat_620_General_CWM_A$AGE,CAT=dat_620_General_CWM_A$CAT),FUN=length)
sd_620_General_CWM_A$Estimate.corr.se<-sd_620_General_CWM_A$Estimate.corr/sqrt(n_620_General_CWM_A$Estimate.corr)
msd_620_General_CWM_A<-merge(m_620_General_CWM_A,sd_620_General_CWM_A[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd_620_General_CWM_A$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes


msd_620_General=rbind(
  msd_620_General_FD_A, msd_620_General_CWM_A
)



msd_220_General$Types="220_General"
msd_630_General$Types="630_General"
msd_620_General$Types="620_General"

msd=rbind(msd_220_General,msd_630_General,  msd_620_General     )


####save data####
msd$Nation="HN"
write.csv(msd,"E:/GEB data/msd_HN_Plantations.csv"  )
msd_5=read.csv("E:/GEB data/msd_HN_Plantations.csv")
###### paired-t-test for CAT difference per year
msd$Types[which(msd$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd$Types[which(msd$Types %in% c( "630_General" ))] = c("CBMF_General")
msd$Types[which(msd$Types %in% c("620_General"))] = c("BMF_General")



msd$Types <- factor(msd$Types, levels = c("CP_General","CBMF_General","BMF_General"
))
msd_5$Types[which(msd_5$Types %in%   c( "220_General" ) )] =  c("CP_General")
msd_5$Types[which(msd_5$Types %in% c( "630_General" ))] = c("CBMF_General")
msd_5$Types[which(msd_5$Types %in% c("620_General"))] = c("BMF_General")


msd_5$Types <- factor(msd_5$Types, levels = c("CP_General","CBMF_General","BMF_General" ))

write.csv(msd,"E:/GEB data/Meta_dat1_HN_Plantations.csv" )

write.csv(dat,"E:/GEB data/Meta_dat2_HN_Plantations.csv" )


##===============================================================================
##  Plot1:Succession ---------
##===============================================================================
rm(list=ls())
library(ggplot2)
####American####

msd=read.csv("E:/GEB data/Meta_dat1_American_Plantations.csv" )
unique(msd$Types)
N=read.csv("E:/GEB data/FIA_Plot_N220_Plantations.csv" )

N$Types="CP_General"


msd=merge(msd,N,by.x= c("AGE", "Types") ,by.y= c("STDAGE_mean", "Types"))

# Define a mapping between original and new type labels
type_mapping <- list(
  "CP_General" = "Conifer",
  "CBMF_General" = "Conifer–broadleaf mixed",
  "BMF_General" = "Broadleaf"
  
)

# Apply the mapping to the Types column
msd$Types <- sapply(msd$Types, function(x) type_mapping[[x]])


msd$Types = factor(msd$Types, levels = c("Conifer", "Conifer–broadleaf mixed",
                                         "Broadleaf")
                   
)


head(msd)

msd$Significant <- ifelse(msd$Estimate.corr - msd$Estimate.corr.se <= 0 & msd$Estimate.corr + 
                            msd$Estimate.corr.se >= 0, "NS", "Significant")

msd<- subset(msd, Types  %in%  c( "Conifer","Conifer–broadleaf mixed","Broadleaf" ))
msd$CAT = factor(msd$CAT, levels = c("FD_C", "FD_B", "FD_A", 
                                     "CWM_C", "CWM_B", "CWM_A"),
                 labels = c("FD Conifer", "FD Broadleaf",  "FD",
                            "CWM Conifer", "CWM Broadleaf",  "CWM")) 

msd <- msd %>%
  subset(CAT %in% c("FD", "CWM")) %>%
  dplyr::mutate(CAT = dplyr::recode(CAT,
                                    FD = "FD",
                                    CWM = "FI"))







msd <- msd %>%
  mutate(n_discrete = case_when(
    n <= 100              ~ "1 \u2264 n \u2264 100",   
    n > 100 & n <= 500    ~ "100 < n \u2264 500",
    n > 500 & n <= 2000   ~ "500 < n \u2264 2,000",
    n > 2000 & n <= 4000  ~ "2,000 < n \u2264 4,000"
  ))


n_levels_math <- c(
  "1 \u2264 n \u2264 100", 
  "100 < n \u2264 500", 
  "500 < n \u2264 2,000", 
  "2,000 < n \u2264 4,000"
)

msd$n_discrete <- factor(msd$n_discrete, levels = n_levels_math)

msd <- msd %>%
  add_count(Types, CAT, name = "n_points")

head(msd)
p1= ggplot(msd, aes(x = AGE, y = Estimate.corr, color = CAT, fill=CAT)) +
  
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = TRUE, 
              fill = "#dcdcd6", 
              alpha = 0.2) +
  geom_errorbar(data = msd, aes(x = AGE, ymin = Estimate.corr - Estimate.corr.se,
                                ymax = Estimate.corr + Estimate.corr.se, col = CAT),
                size = 1,width=2) + 
  # geom_point(data = msd, aes(x = AGE, y = Estimate.corr, shape = Significant, col = CAT),
  #            size = 3.5, alpha = 1) +
  geom_point(aes(size = n_discrete, shape = Significant), 
             alpha = 0.7, stroke = 0.8) +
  
  scale_size_manual(
    values = c(
      "1 \u2264 n \u2264 100"    = 2,
      "100 < n \u2264 500"       = 3,
      "500 < n \u2264 2,000"     = 4,
      "2,000 < n \u2264 4,000"   = 5
    ),
    name = "Sample size (n)",
    drop = FALSE
  ) +
  
  scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) + # Solid for 'sign', hollow for 'nosign'
  
  #facet_wrap(~ Types, ncol = 1, scales = "fixed", axes = "all_x")+
  facet_wrap(~ Types, ncol = 1, scales = "free", strip.position = "right", axes = "all_x") +
  
  scale_color_manual(values = c(
    "FD" = "#D37C5E", 
    "FI" = "#9DC7C6") )+ 
  
  
  
  
  scale_y_continuous(
    n.breaks = 4#,                                  
    # labels = label_number(accuracy = 0.1)          
  )+
  scale_x_continuous(limits = c(0, 70)) +
  
  theme_bw( )+
  
  theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank())+
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  ggtitle("(a) Southeastern US")+
  #scale_y_continuous(limits = c(-0.2, 0.4)#, breaks = seq(0, 100, by = 40)
  # ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # 设置 x 轴为1位小数
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  
  # #220 General
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 8, y = -0.25, label = "**"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 21, y = -0.25, label = "*"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 59, y = -0.25, label = "**"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ],
            aes(x = 8, y = 0.25, label = ""), size = 5, color = "black")+
  
  
  theme(
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    plot.title = element_text(size=16, color="black",family="serif"),
    # legend.title = element_text(size = 12,family="serif"), 
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16, family = "serif", color="black")  ,
    axis.line = element_line(size = 0.5), 
    axis.ticks = element_line(size = 0.5),
    
    axis.title.y = element_text(size = 16,family="serif", color="black"), 
    axis.text.y = element_text(size = 16,family="serif", color="black"),
    axis.title.x = element_text(size = 16,family="serif", color="black"), 
    axis.text.x = element_text(size = 16,family="serif", color="black"),
    
    strip.background = element_rect(fill = "white",
                                    color = "NA"),  
    strip.text = element_blank(),
    
    legend.position = "bottom" ) 

p1


#################################################################################
####Japan####

msd=read.csv("E:/GEB data/Meta_dat1_Japan_Plantations.csv" )
msd=msd %>%filter(Types == "CP_General")
N=read.csv("E:/GEB data/Japan_Plot_N220_Plantations.csv" )
N$Types="CP_General"

msd=merge(msd,N,by.x= c("AGE", "Types") ,by.y= c("STDAGE_mean", "Types"),all.x=T)
# Define a mapping between original and new type labels
type_mapping <- list(
  "CP_General" = "Conifer",
  "CBMF_General" = "Conifer–broadleaf mixed",
  "BMF_General" = "Broadleaf"
  
)

# Apply the mapping to the Types column
msd$Types <- sapply(msd$Types, function(x) type_mapping[[x]])


msd$Types = factor(msd$Types, levels = c("Conifer", "Conifer–broadleaf mixed",
                                         "Broadleaf")
                   
)


head(msd)

msd$Significant <- ifelse(msd$Estimate.corr - msd$Estimate.corr.se <= 0 & msd$Estimate.corr + 
                            msd$Estimate.corr.se >= 0, "NS", "Significant")

msd<- subset(msd, Types  %in%  c( "Conifer","Conifer–broadleaf mixed","Broadleaf" ))
msd$CAT = factor(msd$CAT, levels = c("FD_C", "FD_B", "FD_A", 
                                     "CWM_C", "CWM_B", "CWM_A"),
                 labels = c("FD Conifer", "FD Broadleaf",  "FD",
                            "CWM Conifer", "CWM Broadleaf",  "CWM")) 

msd <- msd %>%
  subset(CAT %in% c("FD", "CWM")) %>%
  dplyr::mutate(CAT = dplyr::recode(CAT,
                                    FD = "FD",
                                    CWM = "FI"))


msd <- msd %>%
  mutate(n_discrete = case_when(
    n <= 100              ~ "1 \u2264 n \u2264 100",   
    n > 100 & n <= 500    ~ "100 < n \u2264 500",
    n > 500 & n <= 2000   ~ "500 < n \u2264 2,000",
    n > 2000 & n <= 4000  ~ "2,000 < n \u2264 4,000"
  ))


n_levels_math <- c(
  "1 \u2264 n \u2264 100", 
  "100 < n \u2264 500", 
  "500 < n \u2264 2,000", 
  "2,000 < n \u2264 4,000"
)

msd$n_discrete <- factor(msd$n_discrete, levels = n_levels_math)


msd <- msd %>%
  add_count(Types, CAT, name = "n_points")

head(msd)

p2= ggplot(msd, aes(x = AGE, y = Estimate.corr, color = CAT, fill=CAT)) +
  
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = TRUE, 
              fill = "#dcdcd6", 
              alpha = 0.2) +
  geom_errorbar(data = msd, aes(x = AGE, ymin = Estimate.corr - Estimate.corr.se,
                                ymax = Estimate.corr + Estimate.corr.se, col = CAT),
                size = 1,width=2) + 
  # geom_point(data = msd, aes(x = AGE, y = Estimate.corr, shape = Significant, col = CAT),
  #            size = 3.5, alpha = 1) +
  
  geom_point(aes(size = n_discrete, shape = Significant), 
             alpha = 0.7, stroke = 0.8) +
  
  scale_size_manual(
    values = c(
      "1 \u2264 n \u2264 100"    = 2,
      "100 < n \u2264 500"       = 3,
      "500 < n \u2264 2,000"     = 4,
      "2,000 < n \u2264 4,000"   = 5
    ),
    name = "Sample size (n)",
    drop = FALSE
  ) +
  
  
  
  scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) + # Solid for 'sign', hollow for 'nosign'
  
  #facet_wrap(~ Types, ncol = 1, scales = "free"#, axes = "all_x" )+
  
  facet_wrap(~ Types, ncol = 1, scales = "free", strip.position = "right", axes = "all_x") +
  
  
  scale_color_manual(values = c(
    
    "FD" = "#D37C5E", 
    "FI" = "#9DC7C6") )+ 
  
  
  
  
  scale_x_continuous(limits = c(16, 92)) +
  
  theme_bw( )+
  
  theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank())+
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  ggtitle("(b) Japan")+
  #scale_y_continuous(limits = c(-0.2, 0.4)#, breaks = seq(0, 100, by = 40)
  # ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # 设置 x 轴为1位小数
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  
  # #220 General
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 24, y = -0.25, label = "*"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 43, y = -0.25, label = "*"), size = 5, color = "black")+
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 58, y = -0.25, label = "*"), size = 5, color = "black")+
  
  geom_text(data = msd[msd$Types == "Conifer", ], 
            aes(x = 58, y = 0.1, label = ""), size = 5, color = "black")+
  
  theme(
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    plot.title = element_text(size=16, color="black",family="serif"),
    # legend.title = element_text(size = 12,family="serif"), 
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16, family = "serif", color="black")  ,
    axis.line = element_line(size = 0.5), 
    axis.ticks = element_line(size = 0.5),
    
    axis.title.y = element_text(size = 16,family="serif", color="black"), 
    axis.text.y = element_text(size = 16,family="serif", color="black"),
    axis.title.x = element_text(size = 16,family="serif", color="black"), 
    axis.text.x = element_text(size = 16,family="serif", color="black"),
    
    strip.background = element_rect(fill = "white",
                                    color = "NA"),  
    strip.text = element_blank(),
    
    legend.position = "bottom" ) 

p2


#################################################################################
####HN####

msd=read.csv("E:/GEB data/Meta_dat1_HN_Plantations.csv" )
msd=msd %>%filter(Types == "CP_General")
N=read.csv("E:/GEB data/HN_Plot_N220_Plantations.csv" )

N$Types="CP_General"

msd=merge(msd,N,by.x= c("AGE", "Types") ,by.y= c("STDAGE_mean", "Types"),all.x=T)




type_mapping <- list(
  "CP_General" = "Conifer",
  "CBMF_General" = "Conifer–broadleaf mixed",
  "BMF_General" = "Broadleaf"
  
)

# Apply the mapping to the Types column
msd$Types <- sapply(msd$Types, function(x) type_mapping[[x]])


msd$Types = factor(msd$Types, levels = c("Conifer", "Conifer–broadleaf mixed",
                                         "Broadleaf")
                   
)


head(msd)

msd$Significant <- ifelse(msd$Estimate.corr - msd$Estimate.corr.se <= 0 & msd$Estimate.corr + 
                            msd$Estimate.corr.se >= 0, "NS", "Significant")

msd<- subset(msd, Types  %in%  c( "Conifer","Conifer–broadleaf mixed","Broadleaf" ))
msd$CAT = factor(msd$CAT, levels = c("FD_C", "FD_B", "FD_A", 
                                     "CWM_C", "CWM_B", "CWM_A"),
                 labels = c("FD Conifer", "FD Broadleaf",  "FD",
                            "CWM Conifer", "CWM Broadleaf",  "CWM")) 

msd <- msd %>%
  subset(CAT %in% c("FD", "CWM")) %>%
  dplyr::mutate(CAT = dplyr::recode(CAT,
                                    FD = "FD",
                                    CWM = "FI"))



msd <- msd %>%
  mutate(n_discrete = case_when(
    n <= 100              ~ "1 \u2264 n \u2264 100",   
    n > 100 & n <= 500    ~ "100 < n \u2264 500",
    n > 500 & n <= 2000   ~ "500 < n \u2264 2,000",
    n > 2000 & n <= 4000  ~ "2,000 < n \u2264 4,000"
  ))


n_levels_math <- c(
  "1 \u2264 n \u2264 100", 
  "100 < n \u2264 500", 
  "500 < n \u2264 2,000", 
  "2,000 < n \u2264 4,000"
)

msd$n_discrete <- factor(msd$n_discrete, levels = n_levels_math)

msd <- msd %>%
  add_count(Types, CAT, name = "n_points")

head(msd)
p3= ggplot(msd, aes(x = AGE, y = Estimate.corr, color = CAT, fill=CAT)) +
  
  geom_smooth(data = subset(msd, n_points >= 3), 
              method = "lm", 
              formula = y ~ x, 
              se = TRUE, 
              fill = "#dcdcd6", 
              alpha = 0.2) +
  geom_errorbar(data = msd, aes(x = AGE, ymin = Estimate.corr - Estimate.corr.se,
                                ymax = Estimate.corr + Estimate.corr.se, col = CAT),
                size = 1,width=2) + 
  # geom_point(data = msd, aes(x = AGE, y = Estimate.corr, shape = Significant, col = CAT),
  #            size = 3.5, alpha = 1) +
  facet_wrap(~ Types, ncol = 1, scales = "free", strip.position = "right", axes = "all_x") +
  
  geom_point(aes(size = n_discrete, shape = Significant), 
             alpha = 0.7, stroke = 0.8) +
  
  scale_size_manual(
    values = c(
      "1 \u2264 n \u2264 100"    = 2,
      "100 < n \u2264 500"       = 3,
      "500 < n \u2264 2,000"     = 4,
      "2,000 < n \u2264 4,000"   = 5
    ),
    name = "Sample size (n)",
    drop = FALSE
  ) +
  
  scale_shape_manual(values = c("Significant" = 16, "NS" = 1)) + # Solid for 'sign', hollow for 'nosign'
  
  # facet_wrap(~ Types, ncol = 1, scales = "free", axes = "all_x")+
  
  
  scale_color_manual(values = c( "FD" = "#D37C5E",  "FI" = "#9DC7C6") )+ 
  
  
  scale_y_continuous(
    n.breaks = 4#,                                  
    # labels = label_number(accuracy = 0.1)          
  )+
  
  
  theme_bw( )+
  
  theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank())+
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  
  
  labs(x="Stand age (year)", y = "Trait effect on AGC (slope)" )+
  
  ggtitle("(c) Hunan")+
  
  
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  
  theme(
    panel.spacing.y = unit(2, "lines"),
    
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    plot.title = element_text(size=16, color="black",family="serif"),
    # legend.title = element_text(size = 12,family="serif"), 
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16, family = "serif", color="black")  ,
    axis.line = element_line(size = 0.5), 
    axis.ticks = element_line(size = 0.5),
    
    axis.title.y = element_text(size = 16,family="serif", color="black"), 
    axis.text.y = element_text(size = 16,family="serif", color="black"),
    axis.title.x = element_text(size = 16,family="serif", color="black"), 
    axis.text.x = element_text(size = 16,family="serif", color="black"),
    
    strip.background = element_rect(fill = "white",
                                    color = "NA"),  
    strip.text = element_blank(),
    
    legend.position = "bottom" ) 

p3


#####Plot merge#############
library(stats)
library(ggpubr)


p11=p1+theme(axis.title.y = element_blank())
p21=p2+theme(axis.title.y = element_blank())
p31=p3+theme(axis.title.y = element_blank())



fig1 = ggarrange(p11, p21, p31,
                 
                 ncol = 3, nrow =3,
                 common.legend = TRUE, 
                 legend = "bottom",
                 align = "hv")
fig1

fig1_final <- annotate_figure(fig1,
                              left = text_grob(expression(paste( "Trait effect on AGC (slope)")), 
                                               rot = 90, vjust = 1, size = 16, family = "serif"))


fig1_final <- fig1_final + theme(plot.background = element_rect(fill = "white", color = NA))

print(fig1_final)



ggsave("Succession_Plantations_R.png", path = "E:/GEB data/figure",width =8, height =9,
       units = "in",dpi=600, plot=fig1_final)
