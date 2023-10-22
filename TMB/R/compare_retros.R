## Load functions and packages ----
library(tidyverse)
theme_set(theme_bw())
source("TMB/R/prepare_tmb_objects_2023.R")


## Load data/runs ----
# load("TMB/Selectivity_runs.RData")
peels <- readRDS("TMB/Output/2023_peels.RDS")


## Calculate MSE and relative error ----
mse_mat_proj <- mse_mat_avg5 <- re_mat_proj <- re_mat_avg5 <- matrix(0, nrow = length(peels), ncol = dat$nages + 1)


# - Loop through models
ind <- 1
for(i in 1:length(peels)){ # - Loop models
  for(j in 1:(length(peels[[i]])-1)){ # - Loop peels

    # - Get selectivity
    nyrs <- peels[[i]][[j]]$obj$env$data$nyrs
    slctfsh <- peels[[i]][[j]]$rep$slctfsh[nyrs,]
    slctfsh_proj <- peels[[i]][[j+1]]$rep$slctfsh[nyrs,]
    slctfsh_avg5 <- colMeans(peels[[i]][[j+1]]$rep$slctfsh[(nyrs-5):(nyrs-1),]) # Average of 5 previous years

    # - Calculate metrics
    # -- Mean squared error
    # - Age specific
    mse_mat_proj[i,1:dat$nages] <- mse_mat_proj[i,1:dat$nages] + (slctfsh - slctfsh_proj)^2/(length(peels[[i]])-1)
    mse_mat_avg5[i,1:dat$nages] <- mse_mat_avg5[i,1:dat$nages] + (slctfsh - slctfsh_avg5)^2/(length(peels[[i]])-1)

    # - Across all ages
    mse_mat_proj[i,(dat$nages + 1)] <- mse_mat_proj[i,(dat$nages + 1)] + mean(mse_mat_proj[i,1:dat$nages])/(length(peels[[i]])-1)
    mse_mat_avg5[i,(dat$nages + 1)] <- mse_mat_avg5[i,(dat$nages + 1)] + mean(mse_mat_avg5[i,1:dat$nages])/(length(peels[[i]])-1)

    # -- Mean relative error
    # - Age specific
    re_mat_proj[i,1:dat$nages] <- re_mat_proj[i,1:dat$nages] + (slctfsh_proj - slctfsh)/slctfsh/(length(peels[[i]])-1)
    re_mat_avg5[i,1:dat$nages] <- re_mat_avg5[i,1:dat$nages] + (slctfsh_avg5 - slctfsh)/slctfsh/(length(peels[[i]])-1)

    # - Across all ages
    re_mat_proj[i,(dat$nages + 1)] <- re_mat_proj[i,(dat$nages + 1)] + mean(re_mat_proj[i,1:dat$nages])/(length(peels[[i]])-1)
    re_mat_avg5[i,(dat$nages + 1)] <- re_mat_avg5[i,(dat$nages + 1)] + mean(re_mat_avg5[i,1:dat$nages])/(length(peels[[i]])-1)
  }
}

## Reformat ----
mse_mat_proj <- as.data.frame(mse_mat_proj)
mse_mat_avg5 <- as.data.frame(mse_mat_avg5)

re_mat_proj <- as.data.frame(re_mat_proj)
re_mat_avg5 <- as.data.frame(re_mat_avg5)

colnames(mse_mat_proj) <- colnames(mse_mat_avg5) <- colnames(re_mat_proj) <- colnames(re_mat_avg5) <- c(paste0("Age", 1:dat$nages), "Combined")

# - Metric used
mse_mat_proj$Metric <- "MSE"
mse_mat_avg5$Metric <- "MSE"

re_mat_proj$Metric <- "RE"
re_mat_avg5$Metric <- "RE"

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

##  Combine
results <- rbind(mse_mat_proj, mse_mat_avg5, re_mat_proj, re_mat_avg5) %>%
  arrange(Model, Metric, Selectivity) %>%
  select(Model, Metric, Selectivity, paste0("Age", 1:dat$nages), Combined)
