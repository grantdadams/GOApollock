## Load functions and packages ----
library(tidyverse)
theme_set(theme_bw())
source("TMB/R/prepare_tmb_objects_2023.R")


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
  "Mod 8: 3D-AR1cond",
  # , 'Mod 9: 3D-AR1mar'
  "Mod 10: 2D-AR1condCoh"
)


## Calculate MSE and relative error (Mohn's rho) ----
mse_mat_proj <- mse_mat_avg5 <- re_mat_proj <- re_mat_avg5 <- matrix(0, nrow = length(peels), ncol = dat$nages + 1)
re_mat_ssb <- rep(0, length(peels))


# - Loop through models
ind <- 1
for(i in 1:length(peels)){ # - Loop models
  for(j in 1:(length(peels[[i]])-1)){ # - Loop peels (first peel is base model)

    # Terminal year of peel
    nyrs_peel <- peels[[i]][[j+1]]$obj$env$data$nyrs

    # - Get selectivity in projection year of peel
    slctfsh_base <- peels[[i]][[1]]$rep$slctfsh[nyrs_peel+1,]
    slctfsh_proj <- peels[[i]][[j+1]]$rep$slctfsh[nyrs_peel+1,]
    slctfsh_avg5 <- colMeans(peels[[i]][[j+1]]$rep$slctfsh[(nyrs_peel-4):(nyrs_peel),]) # Average of 5 previous years

    # - Get SSB in terminal year of peel
    ssb_base <- peels[[i]][[1]]$rep$Espawnbio[nyrs_peel]
    ssb_peel <- peels[[i]][[j+1]]$rep$Espawnbio[nyrs_peel]

    # - Calculate metrics
    # -- Mean squared error
    # - Age specific
    mse_mat_proj[i,1:dat$nages] <- mse_mat_proj[i,1:dat$nages] + (slctfsh_base - slctfsh_proj)^2/(length(peels[[i]])-1)
    mse_mat_avg5[i,1:dat$nages] <- mse_mat_avg5[i,1:dat$nages] + (slctfsh_base - slctfsh_avg5)^2/(length(peels[[i]])-1)

    # - Across all ages
    mse_mat_proj[i,(dat$nages + 1)] <- mse_mat_proj[i,(dat$nages + 1)] + mean(mse_mat_proj[i,1:dat$nages])/(length(peels[[i]])-1)
    mse_mat_avg5[i,(dat$nages + 1)] <- mse_mat_avg5[i,(dat$nages + 1)] + mean(mse_mat_avg5[i,1:dat$nages])/(length(peels[[i]])-1)

    # -- Mean relative error
    # - Age specific
    re_mat_proj[i,1:dat$nages] <- re_mat_proj[i,1:dat$nages] + (slctfsh_proj - slctfsh_base)/slctfsh_base/(length(peels[[i]])-1)
    re_mat_avg5[i,1:dat$nages] <- re_mat_avg5[i,1:dat$nages] + (slctfsh_avg5 - slctfsh_base)/slctfsh_base/(length(peels[[i]])-1)

    # - Across all ages
    re_mat_proj[i,(dat$nages + 1)] <- re_mat_proj[i,(dat$nages + 1)] + mean(re_mat_proj[i,1:dat$nages])/(length(peels[[i]])-1)
    re_mat_avg5[i,(dat$nages + 1)] <- re_mat_avg5[i,(dat$nages + 1)] + mean(re_mat_avg5[i,1:dat$nages])/(length(peels[[i]])-1)

    # - Terminal SSB
    re_mat_ssb[i] <- re_mat_ssb[i] + (ssb_peel - ssb_base)/ssb_base/(length(peels[[i]])-1)
  }
}


## Reformat ----
mse_mat_proj <- as.data.frame(mse_mat_proj)
mse_mat_avg5 <- as.data.frame(mse_mat_avg5)

re_mat_proj <- as.data.frame(re_mat_proj)
re_mat_avg5 <- as.data.frame(re_mat_avg5)

re_mat_ssb <- data.frame(Mohns = re_mat_ssb)

colnames(mse_mat_proj) <- colnames(mse_mat_avg5) <- colnames(re_mat_proj) <- colnames(re_mat_avg5) <- c(paste0("Age", 1:dat$nages), "Combined")

# - Metric used
mse_mat_proj$Metric <- "MSE"
mse_mat_avg5$Metric <- "MSE"

re_mat_proj$Metric <- "RE"
re_mat_avg5$Metric <- "RE"

re_mat_ssb$Metric <- "RE"

# - Selectivity projection
mse_mat_proj$Selectivity <- "Projected"
mse_mat_avg5$Selectivity <- "Average"

re_mat_proj$Selectivity <- "Projected"
re_mat_avg5$Selectivity <- "Average"

# - Model
mse_mat_proj$Model <- sapply(peels, function(x) x[[1]]$version)
mse_mat_avg5$Model <- sapply(peels, function(x) x[[1]]$version)

re_mat_proj$Model <- sapply(peels, function(x) x[[1]]$version)
re_mat_avg5$Model <- sapply(peels, function(x) x[[1]]$version)

re_mat_ssb$Model <- sapply(peels, function(x) x[[1]]$version)


##  Combine ----
results <- rbind(mse_mat_proj, mse_mat_avg5, re_mat_proj, re_mat_avg5) %>%
  arrange(Model, Metric, Selectivity) %>%
  select(Model, Metric, Selectivity, paste0("Age", 1:dat$nages), Combined)


write.csv(results, file = "TMB/Results/Table3_Retrospective_selectivity.csv")


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
    annotate('label', x=2010,y=.05, label=paste0(tmp2$model[1], " (Mohn's Rho =", 100 * round(tmp2$Mohns[1],3), ")")) +
    labs(x=NULL, y='Spawning Biomass (million t)')


  # Rel difference
  g2 <- ggplot(tmp1, aes(Year, SSB_Pct_Diff, group=peel, color=factor(peel))) + geom_line() +
    labs(x=NULL, y='Percent difference from peel 0', color=NULL)

  g2 <- g2 + geom_point(data=tmp1 %>% filter(thisyear-peel==Year), size=2)+
    theme(legend.position=c(.35,.15)) +
    guides(color=guide_legend(nrow=2))

  g <- cowplot::plot_grid(g1,g2, nrow=2)

  ggsave(paste0("TMB/Results/Figure10_Model ", i-1, "retros.png"), g, width=ggw, height=ggh*1.5, dpi=200)
}

