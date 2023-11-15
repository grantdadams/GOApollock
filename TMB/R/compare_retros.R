## Load functions and packages ----
library(tidyverse)
theme_set(theme_bw())
source("TMB/R/prepare_tmb_objects_2023.R")
source("TMB/R/proj_calcs/run_projections.R")


## Load data/runs ----
# load("TMB/Selectivity_runs.RData")
peels <- readRDS("TMB/Output/2023_peels.RDS")
model_names <- c(
  'Mod 0: Constant',
  'Mod 1: ParDevs',
  'Mod 2: Log-AR1-Age',
  'Mod 3: Log-AR1-Yr',
  'Mod 4: Log-2D-AR1',
  'Mod 5: Age-specific',
  'Mod 6: AR1-Yr',
  'Mod 7: 2D-AR1',
  "Mod 8: 3D-AR1cond"
  # , 'Mod 9: 3D-AR1mar'
  # "Mod 10: 2D-AR1condCoh"
)


## Calculate MSE and relative error (Mohn's rho) ----
mse_mat_proj <- mse_mat_avg5 <- re_mat_proj <- re_mat_avg5 <- matrix(0, nrow = length(peels), ncol = dat$nages + 1) # Age- specific sel
mse_brp_proj <- mse_brp_avg5 <- re_brp_proj <- re_brp_avg5 <- matrix(0, nrow = length(peels), ncol = 4) # B0, B40, OFL, ABC
re_mat_ssb <- rep(0, length(peels)) # SSB rho
re_mat_rec <- rep(0, length(peels)) # Rec rho
aic_mat <- matrix(0, nrow = length(peels), ncol = length(peels[[1]])) # AIC



# - Loop through models
ind <- 1
for(i in 1:length(peels)){ # - Loop models
  for(j in 1:(length(peels[[i]])-1)){ # - Loop peels (first peel is base model)

    ## AIC ----
    aic_mat[i,1] <- as.numeric(peels[[i]][[1]]$opt$AIC)
    aic_mat[i,j+1] <- as.numeric(peels[[i]][[j+1]]$opt$AIC)

    ## Selectivity ----
    # Terminal year of peel
    nyrs_peel <- peels[[i]][[j+1]]$obj$env$data$nyrs
    rec_yrs <- peels[[i]][[j+1]]$obj$env$data$styr:(peels[[i]][[j+1]]$obj$env$data$endyr -1)
    rec_ind <- which(rec_yrs >= 1978)

    # - Get selectivity in projection year of peel
    slctfsh_base <- peels[[i]][[1]]$rep$slctfsh[nyrs_peel+1,]
    slctfsh_proj <- peels[[i]][[j+1]]$rep$slctfsh[nyrs_peel+1,]
    slctfsh_avg5 <- colMeans(peels[[i]][[j+1]]$rep$slctfsh[(nyrs_peel-4):(nyrs_peel),]) # Average of 5 previous years

    # - Calculate metrics
    # -- Mean squared error
    # - Age specific
    mse_mat_proj[i,1:dat$nages] <- mse_mat_proj[i,1:dat$nages] + (slctfsh_base - slctfsh_proj)^2/(length(peels[[i]])-1)
    mse_mat_avg5[i,1:dat$nages] <- mse_mat_avg5[i,1:dat$nages] + (slctfsh_base - slctfsh_avg5)^2/(length(peels[[i]])-1)

    # - Across all ages
    mse_mat_proj[i,(dat$nages + 1)] <- mse_mat_proj[i,(dat$nages + 1)] + mean(mse_mat_proj[i,2:dat$nages])/(length(peels[[i]])-1)
    mse_mat_avg5[i,(dat$nages + 1)] <- mse_mat_avg5[i,(dat$nages + 1)] + mean(mse_mat_avg5[i,2:dat$nages])/(length(peels[[i]])-1)

    # -- Mean relative error (Rho)
    # - Age specific
    re_mat_proj[i,1:dat$nages] <- re_mat_proj[i,1:dat$nages] + (slctfsh_proj - slctfsh_base)/slctfsh_base/(length(peels[[i]])-1)
    re_mat_avg5[i,1:dat$nages] <- re_mat_avg5[i,1:dat$nages] + (slctfsh_avg5 - slctfsh_base)/slctfsh_base/(length(peels[[i]])-1)

    # - Across all ages
    re_mat_proj[i,(dat$nages + 1)] <- re_mat_proj[i,(dat$nages + 1)] + mean(re_mat_proj[i,2:dat$nages])/(length(peels[[i]])-1)
    re_mat_avg5[i,(dat$nages + 1)] <- re_mat_avg5[i,(dat$nages + 1)] + mean(re_mat_avg5[i,2:dat$nages])/(length(peels[[i]])-1)


    ##  Terminal SSB Rho ----
    # - Get SSB in terminal year of peel
    ssb_base <- peels[[i]][[1]]$rep$Espawnbio[nyrs_peel]
    ssb_peel <- peels[[i]][[j+1]]$rep$Espawnbio[nyrs_peel]
    re_mat_ssb[i] <- re_mat_ssb[i] + (ssb_peel - ssb_base)/ssb_base/(length(peels[[i]])-1)

    ##  Terminal REC Rho ----
    # - Get rec in terminal year of peel
    rec_base <- mean(peels[[i]][[1]]$rep$recruit[rec_ind])
    rec_peel <- mean(peels[[i]][[j+1]]$rep$recruit[rec_ind])
    re_mat_rec[i] <- re_mat_rec[i] + (rec_peel - rec_base)/rec_base/(length(peels[[i]])-1)

    ## BRPs ----
    abc_base <- abc_calc_tmb(replist = peels[[i]][[1]]$rep, datlist = dat, peel = j, proj_sel = TRUE)$proj_scens %>% # Dynamic RP for peel year
      filter(Alternative == 1)
    abc_proj_sel <- abc_calc_tmb(replist = peels[[i]][[j+1]]$rep, datlist = peels[[i]][[j+1]]$obj$env$data, peel = 0, proj_sel = TRUE)$proj_scens %>% # Dynamic RP for peel year
      filter(Alternative == 1)
    abc_mean_sel <- abc_calc_tmb(replist = peels[[i]][[j+1]]$rep, datlist = peels[[i]][[j+1]]$obj$env$data, peel = 0, proj_sel = FALSE)$proj_scens %>% # Dynamic RP for peel year
      filter(Alternative == 1)

    # - Mean squared error
    # - B0, B40, OFL, ABC
    # -- Projected sel
    mse_brp_proj[i,1] <- mse_brp_proj[i,1] + (abc_base$B0[1]/1e6 - abc_proj_sel$B0[1]/1e6)^2/(length(peels[[i]])-1) # B0
    mse_brp_proj[i,2] <- mse_brp_proj[i,2] + (abc_base$B40[1]/1e6 - abc_proj_sel$B40[1]/1e6)^2/(length(peels[[i]])-1) #B40
    mse_brp_proj[i,3] <- mse_brp_proj[i,3] + (abc_base$OFL[2]/1e6 - abc_proj_sel$OFL[2]/1e6)^2/(length(peels[[i]])-1) # OFL
    mse_brp_proj[i,4] <- mse_brp_proj[i,4] + (abc_base$ABC[2]/1e6 - abc_proj_sel$ABC[2]/1e6)^2/(length(peels[[i]])-1) # ABC

    # -- Mean sel
    mse_brp_avg5[i,1] <- mse_brp_avg5[i,1] + (abc_base$B0[1]/1e6 - abc_mean_sel$B0[1]/1e6)^2/(length(peels[[i]])-1) # B0
    mse_brp_avg5[i,2] <- mse_brp_avg5[i,2] + (abc_base$B40[1]/1e6 - abc_mean_sel$B40[1]/1e6)^2/(length(peels[[i]])-1) #B40
    mse_brp_avg5[i,3] <- mse_brp_avg5[i,3] + (abc_base$OFL[2]/1e6 - abc_mean_sel$OFL[2]/1e6)^2/(length(peels[[i]])-1) # OFL
    mse_brp_avg5[i,4] <- mse_brp_avg5[i,4] + (abc_base$ABC[2]/1e6 - abc_mean_sel$ABC[2]/1e6)^2/(length(peels[[i]])-1) # ABC

    # - Mean relative error (Rho)
    # -- Projected sel
    re_brp_proj[i,1] <- re_brp_proj[i,1] + (abc_base$B0[1] - abc_proj_sel$B0[1])/abc_base$B0[1]/(length(peels[[i]])-1) # B0
    re_brp_proj[i,2] <- re_brp_proj[i,2] + (abc_base$B40[1] - abc_proj_sel$B40[1])/abc_base$B40[1]/(length(peels[[i]])-1) #B40
    re_brp_proj[i,3] <- re_brp_proj[i,3] + (abc_base$OFL[2] - abc_proj_sel$OFL[2])/abc_base$OFL[2]/(length(peels[[i]])-1) # OFL
    re_brp_proj[i,4] <- re_brp_proj[i,4] + (abc_base$ABC[2] - abc_proj_sel$ABC[2])/abc_base$ABC[2]/(length(peels[[i]])-1) # ABC

    # -- Mean sel
    re_brp_avg5[i,1] <- re_brp_avg5[i,1] + (abc_base$B0[1] - abc_mean_sel$B0[1])/abc_base$B0[1]/(length(peels[[i]])-1) # B0
    re_brp_avg5[i,2] <- re_brp_avg5[i,2] + (abc_base$B40[1] - abc_mean_sel$B40[1])/abc_base$B40[1]/(length(peels[[i]])-1) #B40
    re_brp_avg5[i,3] <- re_brp_avg5[i,3] + (abc_base$OFL[2] - abc_mean_sel$OFL[2])/abc_base$OFL[2]/(length(peels[[i]])-1) # OFL
    re_brp_avg5[i,4] <- re_brp_avg5[i,4] + (abc_base$ABC[2] - abc_mean_sel$ABC[2])/abc_base$ABC[2]/(length(peels[[i]])-1) # ABC
  }
}


## Reformat ----
# -- SEL
mse_mat_proj <- as.data.frame(mse_mat_proj)
mse_mat_avg5 <- as.data.frame(mse_mat_avg5)

re_mat_proj <- as.data.frame(re_mat_proj)
re_mat_avg5 <- as.data.frame(re_mat_avg5)
colnames(mse_mat_proj) <- colnames(mse_mat_avg5) <- colnames(re_mat_proj) <- colnames(re_mat_avg5) <- c(paste0("Age", 1:dat$nages), "Combined")

re_mat_ssb <- data.frame(Mohns = re_mat_ssb, MohnsRec = re_mat_rec)

aic_mat <- as.data.frame(aic_mat)
colnames(aic_mat) <- paste0("Peel", 0:7)

# -- BRPS
mse_brp_proj <- as.data.frame(mse_brp_proj)
mse_brp_avg5 <- as.data.frame(mse_brp_avg5)
re_brp_proj <- as.data.frame(re_brp_proj)
re_brp_avg5 <- as.data.frame(re_brp_avg5)
colnames(mse_brp_proj) <- colnames(mse_brp_avg5) <- colnames(re_brp_proj) <- colnames(re_brp_avg5) <- c("B0", "B40", "OFL", "ABC")

# - Metric used
mse_mat_proj$Metric <- "MSE"
mse_mat_avg5$Metric <- "MSE"

re_mat_proj$Metric <- "RE"
re_mat_avg5$Metric <- "RE"

re_mat_ssb$Metric <- "RE"

mse_brp_proj$Metric <- "MSE"
mse_brp_avg5$Metric <- "MSE"

re_brp_proj$Metric <- "RE"
re_brp_avg5$Metric <- "RE"

# - Selectivity projection
mse_mat_proj$Selectivity <- "Projected"
mse_mat_avg5$Selectivity <- "Average"

re_mat_proj$Selectivity <- "Projected"
re_mat_avg5$Selectivity <- "Average"

mse_brp_proj$Selectivity <- "Projected"
mse_brp_avg5$Selectivity <- "Average"

re_brp_proj$Selectivity <- "Projected"
re_brp_avg5$Selectivity <- "Average"

# - Model
aic_mat$Model <- sapply(peels, function(x) x[[1]]$version)
mse_mat_proj$Model <- sapply(peels, function(x) x[[1]]$version)
mse_mat_avg5$Model <- sapply(peels, function(x) x[[1]]$version)

re_mat_proj$Model <- sapply(peels, function(x) x[[1]]$version)
re_mat_avg5$Model <- sapply(peels, function(x) x[[1]]$version)

re_mat_ssb$Model <- sapply(peels, function(x) x[[1]]$version)

mse_brp_proj$Model <- sapply(peels, function(x) x[[1]]$version)
mse_brp_avg5$Model <- sapply(peels, function(x) x[[1]]$version)

re_brp_proj$Model <- sapply(peels, function(x) x[[1]]$version)
re_brp_avg5$Model <- sapply(peels, function(x) x[[1]]$version)


##  Combine and save ----
# - Sel
results <- rbind(mse_mat_proj, mse_mat_avg5, re_mat_proj, re_mat_avg5) %>%
  arrange(Model, Metric, Selectivity) %>%
  select(Model, Metric, Selectivity, paste0("Age", 1:dat$nages), Combined)

write.csv(results, file = "TMB/Results/Table3_Retrospective_selectivity.csv")

# - Brp
brp_results <- rbind(mse_brp_proj, mse_brp_avg5, re_brp_proj, re_brp_avg5) %>%
  arrange(Model, Metric, Selectivity) %>%
  select(Model, Metric, Selectivity, "B0", "B40", "OFL", "ABC")
write.csv(brp_results, file = "TMB/Results/Table3_Retrospective_BRPs.csv")

# - AIC
write.csv(aic_mat, file = "TMB/Results/Table5_Retrospective_AIC.csv")
aic_mat_keep <- aic_mat[c(2,8,9),]
delta_aic <- as.data.frame(apply(aic_mat_keep[,1:8], 2, function(x) x - min(x)))
delta_aic$Model <- aic_mat_keep$Model
write.csv(delta_aic, file = "TMB/Results/Table5_Retrospective_dAIC_main.csv")


## Plots ----

# * Get ssb ----
ssb <- data.frame(model=model_names[1], Mohns = re_mat_ssb$Mohns[1], peel = 0,
                  SSB = peels[[1]][[1]]$rep$Espawnbio,
                  SE = peels[[1]][[1]]$sd$se[which(peels[[1]][[1]]$sd$name == "Espawnbio")],
                  Year = peels[[1]][[1]]$obj$env$data$styr:peels[[1]][[1]]$obj$env$data$endyr,
                  SSB_Pct_Diff = 100*(peels[[1]][[1]]$rep$Espawnbio-peels[[1]][[1]]$rep$Espawnbio)/peels[[1]][[1]]$rep$Espawnbio)

for(i in 1:length(peels)){
  for(j in 1:length(peels[[i]])){
    ssb <- rbind(ssb,
                 data.frame(model=model_names[i], Mohns = re_mat_ssb$Mohns[i], peel = j-1,
                            SSB = peels[[i]][[j]]$rep$Espawnbio,
                            SE = peels[[i]][[j]]$sd$se[which(peels[[i]][[j]]$sd$name == "Espawnbio")],
                            Year = peels[[i]][[j]]$obj$env$data$styr:peels[[i]][[j]]$obj$env$data$endyr,

                            SSB_Pct_Diff = 100*(peels[[i]][[j]]$rep$Espawnbio-peels[[i]][[1]]$rep$Espawnbio[1:length(peels[[i]][[j]]$obj$env$data$styr:peels[[i]][[j]]$obj$env$data$endyr)])/peels[[i]][[1]]$rep$Espawnbio[1:length(peels[[i]][[j]]$obj$env$data$styr:peels[[i]][[j]]$obj$env$data$endyr)])
    )
  }
}
ssb <- ssb[-1,]


# * Plot SSB ----
ssb %>%
  filter(grepl(x=model, paste0("Mod ", c(0,1,7,8),":", collapse = "|"))) %>%
  filter(Year >= 2010) %>%
  mutate(label = paste0(model, " (Mohn's Rho =", round(Mohns,3), ")"),
         peel = as.factor(peel)) %>%
  ggplot(aes(Year, SSB, group = peel, color=peel)) +
  #geom_ribbon(alpha=.5) +
  geom_line(lwd=1) +
  facet_wrap(~label)


# * Cole's plot ----

## Plot it
ggw <- 5.5; ggh <- 4
ssb <- ssb %>%
  mutate(
    UPR = SSB + 1.92 * SE,
    LWR = SSB - 1.92 * SE
  )

thisyear <- 2023
for(i in c(1,2,8,9)){
  tmp1 <- ssb %>% filter(model == model_names[i])
  g1 <- ggplot(tmp1 %>% filter(peel == 0),
               aes(Year, SSB, group=NULL, fill=NULL, color=NULL,
                   ymin=LWR, ymax=UPR)) + geom_ribbon(alpha=.25)

  g1 <- g1 + geom_line(data = tmp1,
                       aes(Year, SSB, group=peel, lwr=NULL,
                           upr=NULL, fill=NULL, color=factor(peel)))
  tmp2 <- filter(tmp1, thisyear-peel==Year)

  g1 <- g1 +
    geom_point(data=tmp2, aes(y=SSB, color=factor(peel)), size=2) +
    theme(legend.position='none') +
    annotate('label', x=2008,y=.05, label=paste0(tmp2$model[1], " (Mohn's Rho =", round(tmp2$Mohns[1],3), ")")) +
    labs(x=NULL, y='Spawning Biomass (million t)')


  # Rel difference
  g2 <- ggplot(tmp1, aes(Year, SSB_Pct_Diff, group=peel, color=factor(peel))) + geom_line() +
    labs(x=NULL, y='Percent difference from peel 0', color=NULL)

  g2 <- g2 + geom_point(data=tmp1 %>% filter(thisyear-peel==Year), size=2)+
    theme(legend.position=c(.35,.15)) +
    guides(color=guide_legend(nrow=2))

  g <- cowplot::plot_grid(g1,g2, nrow=2)

  ggsave(paste0("TMB/Results/Figure10 Model ", i-1, " retros.png"), g, width=ggw, height=ggh*1.5, dpi=200)
}



## Plot sel Mohn's by age ----
age_sel <- results %>%
  select(-Combined) %>%
  pivot_longer(!c(Model, Metric, Selectivity), names_to = "Age", values_to = "Bias") %>%
  mutate(Age = as.numeric(gsub("Age","", Age)),
         Model = as.factor(Model)) %>%
  mutate(Model_name = case_when(
    Model == "mod8" ~ model_names[9],
    Model == "mod1" ~ model_names[2],
    TRUE ~ Model
  ),
  )

g1 <- ggplot(age_sel %>% filter(Metric == "RE" & Model %in% c("mod1", "mod8")), aes(x = Age, y = Bias, fill=Model_name)) +
  geom_hline(yintercept=0) +
  geom_point(aes(shape=Selectivity), size=4) +
  scale_shape_manual(values=c(21,24)) +
  labs(x="Age", y="Mohn's Rho") +
  facet_wrap(~ Model_name) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = NA)),
         colour = guide_legend(override.aes = list(shape = 21)))

g2 <- ggplot(age_sel %>% filter(Metric == "MSE" & Model %in% c("mod1", "mod8")), aes(x = Age, y = Bias, fill=Model_name )) +
  geom_point(aes(shape=Selectivity), size=4) +
  scale_shape_manual(values=c(21,24)) +
  labs(x="Age", y="MSE") +
  facet_wrap(~ Model_name) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = NA)),
         colour = guide_legend(override.aes = list(shape = 21)))

g <- cowplot::plot_grid(g1,g2, nrow=2)
g

ggsave(paste0("TMB/Results/Figure11 Model age-selex retros.png"), g, width=ggw*2, height=ggh*1.5, dpi=200)

