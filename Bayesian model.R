#####bayes_model#####
bayes_sem_model <- "
model {
  
 
  for(i in 1:N){
   
    mu_CV[i] <- alpha_CV[bin[i]] +
      delta1[bin[i]] * trait[i] +
      delta2[bin[i]] * trait_CWM[i] +
      gamma_CV[1] * MAT[i] +
      gamma_CV[2] * MAP[i] +
      gamma_CV[3] * VPD[i] + 
      gamma_CV[5] * pH_mean[i] +
      gamma_CV[6] * clay_mean[i] +
      gamma_CV[7] * STDAGE[i]    
    
    DBH_CV[i] ~ dnorm(mu_CV[i], tau_CV)
  }
  
  # bin-level random effects for CV
  for(j in 1:N_bin){
 
    alpha_CV_mu[j] <- gamma0_CV + 
      gamma_CV_bin[1] * MAT_bin[j] + 
      gamma_CV_bin[2] * MAP_bin[j] +
      gamma_CV_bin[3] * VPD_bin[j] + 
      gamma_CV_bin[5] * pH_mean_bin[j] +
      gamma_CV_bin[6] * clay_mean_bin[j] +
      gamma_CV_bin[7] * STDAGE_bin[j]    
    
    alpha_CV[j] ~ dnorm(alpha_CV_mu[j], tau_alpha_CV)
    
    delta1[j] ~ dnorm(0, tau_delta1)
    delta2[j] ~ dnorm(0, tau_delta2)
  }
  
 
  for(i in 1:N){
   
    mu_AGB[i] <- alpha0 + alpha_AGB[bin[i]] +   
      beta3[bin[i]] * trait[i] +
      beta4[bin[i]] * trait_CWM[i] +
      beta5[bin[i]] * DBH_CV[i] +
      gamma_AGB[1] * MAT[i] + 
      gamma_AGB[2] * MAP[i] +
      gamma_AGB[3] * VPD[i] + 
      gamma_AGB[5] * pH_mean[i] +
      gamma_AGB[6] * clay_mean[i] +
      gamma_AGB[7] * STDAGE[i]    
    
    Carbon_Mg_ha[i] ~ dnorm(mu_AGB[i], tau_AGB)
  }
  
  # bin-level random effects for AGB
  for(j in 1:N_bin){
    alpha_AGB_mu[j] <- 0   
    alpha_AGB[j] ~ dnorm(alpha_AGB_mu[j], tau_alpha_AGB)
    
    beta3[j] ~ dnorm(0, tau_beta3)
    beta4[j] ~ dnorm(0, tau_beta4)
    
   
    beta5_mu[j] <- beta5_0 + 
      beta5_1 * FD_bin[j] + 
      beta5_2 * FI_bin[j] + 
      beta5_3 * MAT_bin[j] + 
      beta5_4 * MAP_bin[j] +
      beta5_5 * VPD_bin[j] + 
      beta5_7 * pH_mean_bin[j] +
      beta5_8 * clay_mean_bin[j] +
      beta5_9 * STDAGE_bin[j]    
    beta5[j] ~ dnorm(beta5_mu[j], tau_beta5)
  }
  
  alpha0 ~ dnorm(0, 1.0E-4)   
  gamma0_CV ~ dnorm(0, 1.0E-4) 

  tau_CV <- pow(sigma_CV, -2)
  sigma_CV ~ dunif(0, 5)        
  
  tau_AGB <- pow(sigma_AGB, -2)
  sigma_AGB ~ dunif(0, 5)
  
  tau_alpha_CV <- pow(sigma_alpha_CV, -2)
  sigma_alpha_CV ~ dunif(0, 5)
  
  tau_alpha_AGB <- pow(sigma_alpha_AGB, -2)
  sigma_alpha_AGB ~ dunif(0, 5)
  
  tau_delta1 <- pow(sigma_delta1, -2)
  sigma_delta1 ~ dunif(0, 5)
  
  tau_delta2 <- pow(sigma_delta2, -2)
  sigma_delta2 ~ dunif(0, 5)
  
  tau_beta3 <- pow(sigma_beta3, -2)
  sigma_beta3 ~ dunif(0, 5)
  
  tau_beta4 <- pow(sigma_beta4, -2)
  sigma_beta4 ~ dunif(0, 5)
  
  tau_beta5 <- pow(sigma_beta5, -2)
  sigma_beta5 ~ dunif(0, 5)
  
 
  for(k in c(1,2,3,5,6,7)) {
    gamma_CV[k] ~ dnorm(0, 1.0) 
    gamma_CV_bin[k] ~ dnorm(0, 0.25) 
    gamma_AGB[k] ~ dnorm(0, 1.0)
  }
  
  
  beta5_0 ~ dnorm(0, 1.0)
  beta5_1 ~ dnorm(0, 0.25)
  beta5_2 ~ dnorm(0, 0.25)
  beta5_3 ~ dnorm(0, 0.25)
  beta5_4 ~ dnorm(0, 0.25)
  beta5_5 ~ dnorm(0, 0.25)
  beta5_7 ~ dnorm(0, 0.25)
  beta5_8 ~ dnorm(0, 0.25)
  beta5_9 ~ dnorm(0, 0.25)
}"

write(bayes_sem_model, "E:/GEB data/bayes_sem_model.txt")


##===============================================================================
## Figure 4  ---------
##===============================================================================
library(R2jags)
##===============================================================================
## FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "FDis", "PC1.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    FDis          = as.numeric(scale(FDis)),
    PC1.CWM       = as.numeric(scale(PC1.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$FDis       
trait_CWM <- data_all$PC1.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_FDis_620_Age3_NaturalForests.csv" )  

prop_FDis=read.csv("E:/GEB data/American_FDis_620_Age3_NaturalForests.csv"   )
unique(prop_FDis$path)
prop_FDis_filtered <- prop_FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_FDis_filtered$path <- factor(
  prop_FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_FDis_filtered$effect_type <- factor(
  prop_FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_FDis <- ggplot(prop_FDis_filtered,
                   aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_FDis 


##===============================================================================
## Max.height.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Max.height.FDis", "Max.height.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Max.height.FDis          = as.numeric(scale(Max.height.FDis)),
    Max.height.CWM       = as.numeric(scale(Max.height.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Max.height.FDis       
trait_CWM <- data_all$Max.height.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Max.height.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Max.height.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Max.height.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Max.height.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Max.height.FDis_620_Age3_NaturalForests.csv" )  

prop_Max.height.FDis=read.csv("E:/GEB data/American_Max.height.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Max.height.FDis$path)
prop_Max.height.FDis_filtered <- prop_Max.height.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Max.height.FDis_filtered$path <- factor(
  prop_Max.height.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Max.height.FDis_filtered$effect_type <- factor(
  prop_Max.height.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Max.height.FDis <- ggplot(prop_Max.height.FDis_filtered,
                              aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Max.height.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Max.height.FDis 


##===============================================================================
## WD_mean_add.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "WD_mean_add.FDis", "WD_mean_add.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    WD_mean_add.FDis          = as.numeric(scale(WD_mean_add.FDis)),
    WD_mean_add.CWM       = as.numeric(scale(WD_mean_add.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$WD_mean_add.FDis       
trait_CWM <- data_all$WD_mean_add.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_WD_mean_add.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_WD_mean_add.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_WD_mean_add.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max WD_mean_add.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_WD_mean_add.FDis_620_Age3_NaturalForests.csv" )  

prop_WD_mean_add.FDis=read.csv("E:/GEB data/American_WD_mean_add.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_WD_mean_add.FDis$path)
prop_WD_mean_add.FDis_filtered <- prop_WD_mean_add.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_WD_mean_add.FDis_filtered$path <- factor(
  prop_WD_mean_add.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_WD_mean_add.FDis_filtered$effect_type <- factor(
  prop_WD_mean_add.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_WD_mean_add.FDis <- ggplot(prop_WD_mean_add.FDis_filtered,
                               aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max WD_mean_add.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_WD_mean_add.FDis 


##===============================================================================
## LMA.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "LMA.FDis", "LMA.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    LMA.FDis          = as.numeric(scale(LMA.FDis)),
    LMA.CWM       = as.numeric(scale(LMA.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$LMA.FDis       
trait_CWM <- data_all$LMA.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_LMA.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_LMA.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_LMA.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max LMA.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_LMA.FDis_620_Age3_NaturalForests.csv" )  

prop_LMA.FDis=read.csv("E:/GEB data/American_LMA.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_LMA.FDis$path)
prop_LMA.FDis_filtered <- prop_LMA.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_LMA.FDis_filtered$path <- factor(
  prop_LMA.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_LMA.FDis_filtered$effect_type <- factor(
  prop_LMA.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_LMA.FDis <- ggplot(prop_LMA.FDis_filtered,
                       aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max LMA.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_LMA.FDis 


##===============================================================================
## Leaf.longevity.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Leaf.longevity.FDis", "Leaf.longevity.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Leaf.longevity.FDis          = as.numeric(scale(Leaf.longevity.FDis)),
    Leaf.longevity.CWM       = as.numeric(scale(Leaf.longevity.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Leaf.longevity.FDis       
trait_CWM <- data_all$Leaf.longevity.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Leaf.longevity.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Leaf.longevity.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Leaf.longevity.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Leaf.longevity.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Leaf.longevity.FDis_620_Age3_NaturalForests.csv" )  

prop_Leaf.longevity.FDis=read.csv("E:/GEB data/American_Leaf.longevity.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Leaf.longevity.FDis$path)
prop_Leaf.longevity.FDis_filtered <- prop_Leaf.longevity.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Leaf.longevity.FDis_filtered$path <- factor(
  prop_Leaf.longevity.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Leaf.longevity.FDis_filtered$effect_type <- factor(
  prop_Leaf.longevity.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Leaf.longevity.FDis <- ggplot(prop_Leaf.longevity.FDis_filtered,
                                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Leaf.longevity.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Leaf.longevity.FDis 


##===============================================================================
## Nmass.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Nmass.FDis", "Nmass.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Nmass.FDis          = as.numeric(scale(Nmass.FDis)),
    Nmass.CWM       = as.numeric(scale(Nmass.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Nmass.FDis       
trait_CWM <- data_all$Nmass.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Nmass.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Nmass.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Nmass.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Nmass.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Nmass.FDis_620_Age3_NaturalForests.csv" )  

prop_Nmass.FDis=read.csv("E:/GEB data/American_Nmass.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Nmass.FDis$path)
prop_Nmass.FDis_filtered <- prop_Nmass.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Nmass.FDis_filtered$path <- factor(
  prop_Nmass.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Nmass.FDis_filtered$effect_type <- factor(
  prop_Nmass.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Nmass.FDis <- ggplot(prop_Nmass.FDis_filtered,
                         aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Nmass.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Nmass.FDis 


##===============================================================================
## Pmass.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Pmass.FDis", "Pmass.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Pmass.FDis          = as.numeric(scale(Pmass.FDis)),
    Pmass.CWM       = as.numeric(scale(Pmass.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Pmass.FDis       
trait_CWM <- data_all$Pmass.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Pmass.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Pmass.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Pmass.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Pmass.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Pmass.FDis_620_Age3_NaturalForests.csv" )  

prop_Pmass.FDis=read.csv("E:/GEB data/American_Pmass.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Pmass.FDis$path)
prop_Pmass.FDis_filtered <- prop_Pmass.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Pmass.FDis_filtered$path <- factor(
  prop_Pmass.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Pmass.FDis_filtered$effect_type <- factor(
  prop_Pmass.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Pmass.FDis <- ggplot(prop_Pmass.FDis_filtered,
                         aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Pmass.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Pmass.FDis 


##===============================================================================
## Shade.tolerance.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Shade.tolerance.FDis", "Shade.tolerance.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Shade.tolerance.FDis          = as.numeric(scale(Shade.tolerance.FDis)),
    Shade.tolerance.CWM       = as.numeric(scale(Shade.tolerance.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Shade.tolerance.FDis       
trait_CWM <- data_all$Shade.tolerance.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Shade.tolerance.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Shade.tolerance.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Shade.tolerance.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Shade.tolerance.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Shade.tolerance.FDis_620_Age3_NaturalForests.csv" )  

prop_Shade.tolerance.FDis=read.csv("E:/GEB data/American_Shade.tolerance.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Shade.tolerance.FDis$path)
prop_Shade.tolerance.FDis_filtered <- prop_Shade.tolerance.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Shade.tolerance.FDis_filtered$path <- factor(
  prop_Shade.tolerance.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Shade.tolerance.FDis_filtered$effect_type <- factor(
  prop_Shade.tolerance.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Shade.tolerance.FDis <- ggplot(prop_Shade.tolerance.FDis_filtered,
                                   aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Shade.tolerance.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Shade.tolerance.FDis 


##===============================================================================
## Drought.tolerance.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Drought.tolerance.FDis", "Drought.tolerance.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Drought.tolerance.FDis          = as.numeric(scale(Drought.tolerance.FDis)),
    Drought.tolerance.CWM       = as.numeric(scale(Drought.tolerance.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Drought.tolerance.FDis       
trait_CWM <- data_all$Drought.tolerance.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Drought.tolerance.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Drought.tolerance.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Drought.tolerance.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Drought.tolerance.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Drought.tolerance.FDis_620_Age3_NaturalForests.csv" )  

prop_Drought.tolerance.FDis=read.csv("E:/GEB data/American_Drought.tolerance.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Drought.tolerance.FDis$path)
prop_Drought.tolerance.FDis_filtered <- prop_Drought.tolerance.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Drought.tolerance.FDis_filtered$path <- factor(
  prop_Drought.tolerance.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Drought.tolerance.FDis_filtered$effect_type <- factor(
  prop_Drought.tolerance.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Drought.tolerance.FDis <- ggplot(prop_Drought.tolerance.FDis_filtered,
                                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Drought.tolerance.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Drought.tolerance.FDis 


##===============================================================================
## Leaf.FDis：Age3 ---------
##===============================================================================
rm(list=ls())

C_datatotal_all <- fread("E:/GEB data/FIA_data_sem_NaturalForests.csv")

C_bin <- read.csv("E:/GEB data/Class_Climate_NaturalForests.csv")
C_bin <- C_bin[, c("MAT", "MAP", "MAT_MAP_bin")]
data_sem <- merge(C_datatotal_all, C_bin, by = c("MAT", "MAP"))


data_sem <- data_sem[, c("PLT_CN", "INVYR.x", "FLDTYPCD", "LON", "LAT", "STDAGE",
                         "Carbon_Mg_ha", "cv_dbh", "VPD", "Elevation", "MAT", "MAP", "MAT_MAP_bin",               
                         "pH_mean", "bdod_mean", "clay_mean", "nitrogen_mean", "Leaf.FDis", "Leaf.CWM", "AgeGroup")]

data_sem$PLT_CN <- as.character(data_sem$PLT_CN)


data_sem <- data_sem %>%
  group_by(LON, LAT) %>%
  slice(1) %>%
  ungroup()

data_sem <- na.omit(data_sem)



data_all <- data_sem %>%
  filter(FLDTYPCD == "Broadleaf Forest" & AgeGroup == 3) %>%
  mutate(
    DBH_CV        = as.numeric(scale(log(cv_dbh + 1))),
    MAT           = as.numeric(scale(MAT)),
    MAP           = as.numeric(scale(MAP)),
    VPD           = as.numeric(scale(VPD)),
    pH_mean       = as.numeric(scale(pH_mean)),
    clay_mean     = as.numeric(scale(clay_mean)),
    STDAGE        = as.numeric(scale(STDAGE)), 
    Leaf.FDis          = as.numeric(scale(Leaf.FDis)),
    Leaf.CWM       = as.numeric(scale(Leaf.CWM)),
    Carbon_Mg_ha  = log(Carbon_Mg_ha + 1) 
  )


y_CV      <- data_all$DBH_CV
y_AGB     <- data_all$Carbon_Mg_ha
trait     <- data_all$Leaf.FDis       
trait_CWM <- data_all$Leaf.CWM    
MAT       <- data_all$MAT
MAP       <- data_all$MAP
VPD       <- data_all$VPD 
pH_mean   <- data_all$pH_mean
clay_mean <- data_all$clay_mean
STDAGE    <- data_all$STDAGE 


bin   <- as.numeric(factor(data_all$MAT_MAP_bin))
N     <- nrow(data_all)
N_bin <- max(bin)


FD_bin        <- as.vector(tapply(trait, bin, mean))
FI_bin        <- as.vector(tapply(trait_CWM, bin, mean))
MAT_bin       <- as.vector(tapply(MAT, bin, mean))
MAP_bin       <- as.vector(tapply(MAP, bin, mean))
VPD_bin       <- as.vector(tapply(VPD, bin, mean))
pH_mean_bin   <- as.vector(tapply(pH_mean, bin, mean))
clay_mean_bin <- as.vector(tapply(clay_mean, bin, mean))
STDAGE_bin    <- as.vector(tapply(STDAGE, bin, mean))


jags_data <- list(
  N             = N,
  N_bin         = N_bin,
  DBH_CV        = y_CV,
  Carbon_Mg_ha  = y_AGB,
  trait         = trait,
  trait_CWM     = trait_CWM,
  MAT           = MAT,
  MAP           = MAP,
  bin           = bin,
  FD_bin        = FD_bin,
  FI_bin        = FI_bin,
  MAT_bin       = MAT_bin,
  MAP_bin       = MAP_bin,
  VPD           = VPD,
  pH_mean       = pH_mean,
  clay_mean     = clay_mean,
  STDAGE        = STDAGE,
  VPD_bin       = VPD_bin,
  pH_mean_bin   = pH_mean_bin,
  clay_mean_bin = clay_mean_bin,
  STDAGE_bin    = STDAGE_bin
)


jags_data <- lapply(jags_data, as.vector)




params <- c(
  "delta1", "delta2",
  "beta3", "beta4", "beta5",
  "beta5_0", "beta5_1", "beta5_2", "beta5_3", "beta5_4", "beta5_5", "beta5_7", "beta5_8", "beta5_9", 
  "gamma0_CV", "gamma_CV", "gamma_AGB",   
  "alpha0", "alpha_CV", "alpha_AGB",                     
  "sigma_CV", "sigma_AGB", "sigma_alpha_CV", "sigma_alpha_AGB"
)
fit_sem_jags <- jags(
  model.file         = "E:/GEB data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 30000,    
  n.burnin           = 15000,    
  n.thin             = 15        
)


post <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
rhats <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]


bad_rhats <- rhats[rhats > 1.1]
if(length(bad_rhats) > 0) {
  print("Warning: There are still non-converged parameters:")
  print(bad_rhats)
} else {
  print("Convergence achieved: All parameters successfully converged (Rhat < 1.1).")
}
write.csv(post, "E:/GEB data/posterior_American_Leaf.FDis_620_Age3_NaturalForests.csv", row.names = FALSE)

saveRDS(fit_sem_jags, "E:/GEB data/beiyesi_American_Leaf.FDis_620_Age3_NaturalForests.rds")

write.csv(fit_sem_jags$BUGSoutput$summary, "E:/GEB data/Beiyrsi_Summary_American_Leaf.FDis_620_Age3_NaturalForests.csv")

get_effect_prop <- function(post, param_MAPfix, path_label){
  
  summary_mat <- apply(
    post[, grep(paste0("^", param_MAPfix, "\\."), names(post))],
    2,
    function(x){
      c(
        mean  = mean(x),
        lower = quantile(x, 0.015),
        upper = quantile(x, 0.975)
      )
    }
  )
  
  df <- as.data.frame(t(summary_mat))
  
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0,
    "positive significant",
    ifelse(
      df$lower < 0 & df$upper < 0,
      "negative significant",
      ifelse(
        df$mean > 0,
        "positive non-significant",
        "negative non-significant"
      )
    )
  )
  
  dplyr::count(df, effect_type) %>%
    mutate(
      percent = n / sum(n) * 100,
      path = path_label
    )
}


all_prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV")
)

all_prop$effect_type <- factor(
  all_prop$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)

all_prop$path <- factor(
  all_prop$path,
  levels = c(
    "FD → AGB",
    "CWM → AGB",
    "DBHCV → AGB",
    "FD → DBHCV",
    "CWM → DBHCV"
  )
)


# delta1: FD -> DBHCV
# beta5: DBHCV -> AGB



delta1_samples <- post[, grep("^delta1\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta1_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)
library(scales)

library(ggplot2)
indirect_plot <- ggplot(indirect_prop, aes(x = "Indirect effect: FD → DBHCV → AGB", y = percent, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#92C5DE",
      "positive non-significant" = "#F4A582",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of bins (%)",
    fill = "Effect type"
  ) +
  theme_classic(base_size = 13)

indirect_plot
indirect_prop$path="FD → DBHCV → AGB"

head(indirect_prop)
head(all_prop)

prop=rbind(all_prop,indirect_prop   )


fig_all <- ggplot(prop,
                  aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#2166AC",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Leaf.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_all







# delta2: CWM -> DBHCV
# beta5: DBHCV -> AGB



delta2_samples <- post[, grep("^delta2\\.", names(post))]
beta5_samples  <- post[, grep("^beta5\\.", names(post))]


indirect_samples <- delta2_samples * beta5_samples  # matrix: iterations x bins


indirect_summary <- apply(indirect_samples, 2, function(x){
  c(
    mean = mean(x),
    lower = quantile(x, 0.015),
    upper = quantile(x, 0.975)
  )
})

indirect_df <- as.data.frame(t(indirect_summary))
indirect_df$bin <- rownames(indirect_df)
rownames(indirect_df) <- NULL


indirect_df$effect_type <- ifelse(
  indirect_df$lower > 0 & indirect_df$upper > 0,
  "positive significant",
  ifelse(
    indirect_df$lower < 0 & indirect_df$upper < 0,
    "negative significant",
    ifelse(
      indirect_df$mean > 0,
      "positive non-significant",
      "negative non-significant"
    )
  )
)

indirect_df$effect_type <- factor(
  indirect_df$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)


indirect_prop_CWM <- indirect_df %>%
  dplyr::count(`effect_type`) %>%   
  mutate(percent = n / sum(n) * 100)

indirect_prop_CWM$path="CWM → DBHCV → AGB"




prop=rbind(all_prop,indirect_prop, indirect_prop_CWM  )
write.csv(prop, "E:/GEB data/American_Leaf.FDis_620_Age3_NaturalForests.csv" )  

prop_Leaf.FDis=read.csv("E:/GEB data/American_Leaf.FDis_620_Age3_NaturalForests.csv"   )
unique(prop_Leaf.FDis$path)
prop_Leaf.FDis_filtered <- prop_Leaf.FDis %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))

prop_Leaf.FDis_filtered$path <- factor(
  prop_Leaf.FDis_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  )
)
prop_Leaf.FDis_filtered$effect_type <- factor(
  prop_Leaf.FDis_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  )
)    



fig_Leaf.FDis <- ggplot(prop_Leaf.FDis_filtered,
                        aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     ="#F4A582"
      
    ),
    drop = FALSE
  ) +
  labs(
    x = "Proportion of bins (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  ggtitle("Max Leaf.FDis") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

fig_Leaf.FDis 


#####Plot#################################################################

prop_620_FD3=read.csv("E:/GEB data/American_FDis_620_Age3_NaturalForests.csv"  )

prop_620_Leaf3=read.csv("E:/GEB data/American_Leaf.FDis_620_Age3_NaturalForests.csv"  )

prop_620_LMA3=read.csv("E:/GEB data/American_LMA.FDis_620_Age3_NaturalForests.csv"  )

prop_620_WD3=read.csv("E:/GEB data/American_WD_mean_add.FDis_620_Age3_NaturalForests.csv"  )


prop_620_Nmass3=read.csv("E:/GEB data/American_Nmass.FDis_620_Age3_NaturalForests.csv"  )

prop_620_Pmass3=read.csv("E:/GEB data/American_Pmass.FDis_620_Age3_NaturalForests.csv"  )

prop_620_ST3=read.csv("E:/GEB data/American_Shade.tolerance.FDis_620_Age3_NaturalForests.csv"  )

prop_620_DT3=read.csv("E:/GEB data/American_Drought.tolerance.FDis_620_Age3_NaturalForests.csv"  )

prop_620_LL3=read.csv("E:/GEB data/American_Leaf.longevity.FDis_620_Age3_NaturalForests.csv"  )

prop_620_height3=read.csv("E:/GEB data/American_Max.height.FDis_620_Age3_NaturalForests.csv"  )

add_age_trait <- function(df, age, trait){
  df$Age   <- age
  df$Trait <- trait
  return(df)
}

prop_620_all <- rbind(
  
  add_age_trait(prop_620_FD3, "Age3", "FD"),
  
  add_age_trait(prop_620_LMA3, "Age3", "LMA"),
 
  add_age_trait(prop_620_WD3, "Age3", "WD"),
 
  add_age_trait(prop_620_Nmass3, "Age3", "Nmass"),
 
  add_age_trait(prop_620_Pmass3, "Age3", "Pmass"),
  
  add_age_trait(prop_620_ST3, "Age3", "ST"),
  
  add_age_trait(prop_620_DT3, "Age3", "DT"),
 
  add_age_trait(prop_620_LL3, "Age3", "LL"),
 
  add_age_trait(prop_620_height3, "Age3", "Height"),
  
  add_age_trait(prop_620_Leaf3, "Age3", "Leaf")
)

prop_620_all




prop_height_filtered <- prop_620_all %>%
  filter(!path %in% c("DBHCV → AGB", "FD → DBHCV", "CWM → DBHCV"))%>%
  filter(!Age %in% c("Age1","Age2", "Age4", "Age5"))

unique(prop_620_all$Trait)

prop_height_filtered$path <- factor(
  prop_height_filtered$path,
  levels = c(
    "FD → AGB","FD → DBHCV → AGB",
    "CWM → AGB","CWM → DBHCV → AGB"  ),
  labels = c(
    "FD → AGC","FD → DBH variation → AGC",
    "FI → AGC","FI → DBH variation → AGC"
  )
)
prop_height_filtered$effect_type <- factor(
  prop_height_filtered$effect_type,
  levels = c(
    "negative significant",
    "negative non-significant",
    "positive non-significant",
    "positive significant"
  ),
  labels = c(
    "Significant negative",
    "Non-significant negative",
    "Non-significant positive",
    "Significant positive"
  )
)    

head(prop_height_filtered)

prop_height_filtered$Age <- factor(
  prop_height_filtered$Age,
  levels = c("Age1", "Age2", "Age3", "Age4")
)
prop_height_filtered$Trait <- factor(
  prop_height_filtered$Trait,
  levels = c("FD","WD","Height", "LMA", "LL", "Leaf", "Nmass", "Pmass", "ST", "DT"),
  ,
  labels = c(
    "Multi-trait",
    "Wood density","Max height","LMA", "Leaf longevity","Leaf structure",
    "Nmass", "Pmass",
    "Shade tolerance", "Drought tolerance" 
  )
  
  
)
prop_height_filtered <- prop_height_filtered %>%
  arrange(Trait, Age) %>%
  mutate(
    combo_effect = factor(
      paste(Age, Trait, sep = " | "),
      levels = unique(paste(Age, Trait, sep = " | "))
    )
  )

unique(prop_height_filtered$combo_effect)

#library(lemon)

fig_620 <- ggplot(
  prop_height_filtered,
  aes(
    x = percent,
    y = Trait,
    fill = effect_type
  )
) +
  geom_bar(stat = "identity", width = 0.65) +
  
  facet_wrap(~  path, ncol = 2) +
  
  scale_x_continuous(
    #labels = scales::percent_format(scale = 1),
    expand = c(0, 2.5)
  ) +
  
  scale_fill_manual(
    values = c(
      
      
      "Significant negative"  = "#99cc33",
      "Non-significant negative" = "#ebffac",
      "Non-significant positive"= "#fadfd8",
      "Significant positive" = "#e990ab"
      
      
      
    ),
    drop = FALSE
  ) +
  
  labs(
    x = "Proportion of positive and negative coefficient\n estimates in the climatic units (%)",
    y = NULL,
    fill = "Effect type"
  ) +
  
  theme_classic(base_size = 13) +
 
  theme(legend.position="no",
        
        legend.text = element_text(size=10, family = "serif"),
        legend.title = element_blank(),
        legend.text.align = 1,
        legend.key.size = unit(0.6, "cm"),
        panel.grid = element_blank() ,
        
        axis.text.x=element_text(colour='black',family = "serif",
                                 size=14),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line.y = element_blank( ), 
        axis.ticks.y = element_blank( ),
        
        axis.line.x = element_line(size = 0.5), 
        axis.ticks.x = element_line(size = 0.5),
        
        plot.title = element_text(size=4, family = "serif"),
        plot.subtitle = element_text(hjust = 0, size = 14, family = "serif", color = "black"),
        axis.text.y = element_text(size = 14, family = "serif", color = "black"),
        axis.title = element_text(size = 14, family = "serif", color = "black"),
        strip.text = element_text(size = 14, family="serif"),
        
        panel.background = element_blank(),
        panel.border = element_blank()
  )+
  

  theme(panel.spacing.y = unit(1, "lines")) +
  theme(strip.background = element_blank())+ 
  scale_y_discrete(limits = rev)+
  guides(color = guide_legend(
    nrow = 2,
    direction = "horizontal"  
  ))

fig_620
