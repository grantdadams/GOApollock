## LIBRARIES AND DATA ----
library(GOApollock)
source("TMB/R/prepare_tmb_objects_2023.R")
source("TMB/R/Functions/phaser.R")
source("TMB/R/Functions/create_bounds.R")
source("~/GitHub/GOApollock/TMB/R/Functions/retro functions.R")

# map <- lapply(map, function(x) as.factor(as.numeric(x)*NA))

# trace = 1, print fixed effects vector?
# map off rho?


## COMPILE AND BUILD TMB ----
compile("TMB/src/goa_pk_tmb.cpp")
dyn.load('TMB/src/goa_pk_tmb.dll')


## MODEL 0 ----
# - Double logistic
# - Turn on slp and inf parameters
map_mod0 <- map
dat$seltype <- 0

# -- Fixed effect pars
map_mod0$inf1_fsh_mean <- as.factor(1)
map_mod0$inf2_fsh_mean <- as.factor(1)
map_mod0$log_slp1_fsh_mean <- as.factor(1)
map_mod0$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 1

# - Run model 0
input_mod0 <- list(version="mod0", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars, map=map_mod0, random=NULL)
fit_mod0 <- fit_pk(input=input_mod0, getsd=TRUE, control=control, use_bounds = TRUE)
peels_mod0 <- fit_pk_retros(fit_mod0, peels = 1:5)
peels_mod0 <- c(list(fit_mod0), peels_mod0)


## MODEL 1 ----
# - Double logistic with random effects on ascending portion
# - Turn on sel sd and slp and inf parameters
map_mod1 <- map

# -- Random effect pars
map_mod1$slp1_fsh_dev <- as.factor(1:length(pars$slp1_fsh_dev))
map_mod1$inf1_fsh_dev <- as.factor(1:length(pars$inf1_fsh_dev))

# -- Fixed effect pars
map_mod1$ln_sel_sd <- as.factor(1)
map_mod1$inf2_fsh_mean <- as.factor(1)
map_mod1$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 1

# - Build model 1
random <- c("slp1_fsh_dev", "inf1_fsh_dev", "slp2_fsh_dev", "inf2_fsh_dev")

# - Run model 1
input_mod1 <- list(version="mod1", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars, map=map_mod1, random=random)
fit_mod1 <- fit_pk(input=input_mod1, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=3, silent = TRUE)
peels_mod1 <- fit_pk_retros(fit_mod1, peels = 1:5)
peels_mod1 <- c(list(fit_mod1), peels_mod1)


## MODEL 2 ----
# - Double logistic with AR1 age random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho
map_mod2 <- map

# -- Random effect pars
map_mod2$selpars_re <- as.factor(1:length(pars$selpars_re))

# -- Fixed effect pars
map_mod2$ln_sel_sd <- as.factor(1)
map_mod2$sel_rho_a <- as.factor(1)

map_mod2$inf1_fsh_mean <- as.factor(1)
map_mod2$inf2_fsh_mean <- as.factor(1)
map_mod2$log_slp1_fsh_mean <- as.factor(1)
map_mod2$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 2

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod2 <- list(version="mod2", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars, map=map_mod2, random=random)
fit_mod2 <- fit_pk(input=input_mod2, getsd=TRUE, control=control, use_bounds = TRUE)
peels_mod2 <- fit_pk_retros(fit_mod2, peels = 1:5)
peels_mod2 <- c(list(fit_mod2), peels_mod2)


## MODEL 3 ----
# - Double logistic with AR1 year random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho y
map_mod3 <- map
pars_mod3 <- pars

# -- Random effect pars
pars_mod3$selpars_re <- matrix(0, nrow = 1, ncol = dat$nyrs + dat$projfsh_nyrs)
map_mod3$selpars_re <- as.factor(1:length(pars_mod3$selpars_re))

# -- Fixed effect pars
map_mod3$ln_sel_sd <- as.factor(1)
map_mod3$sel_rho_y <- as.factor(1)

map_mod3$inf1_fsh_mean <- as.factor(1)
map_mod3$inf2_fsh_mean <- as.factor(1)
map_mod3$log_slp1_fsh_mean <- as.factor(1)
map_mod3$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 3

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod3 <- list(version="mod3", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod3, map=map_mod3, random=random)
fit_mod3 <- fit_pk(input=input_mod3, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=2)
peels_mod3 <- fit_pk_retros(fit_mod3, peels = 1:5)
peels_mod3 <- c(list(fit_mod3), peels_mod3)


## MODEL 4 ----
# - Double logistic with 2D-AR1 age, year random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho
map_mod4 <- map
pars_mod4 <- pars

# -- Random effect pars
pars_mod4$selpars_re <- matrix(0, nrow = dat$trmage, ncol = dat$nyrs + dat$projfsh_nyrs)
map_mod4$selpars_re <- as.factor(1:length(pars_mod4$selpars_re))

# -- Fixed effect pars
map_mod4$ln_sel_sd <- as.factor(1)
map_mod4$sel_rho_a <- as.factor(1)
map_mod4$sel_rho_y <- as.factor(1)

map_mod4$inf1_fsh_mean <- as.factor(1)
map_mod4$inf2_fsh_mean <- as.factor(1)
map_mod4$log_slp1_fsh_mean <- as.factor(1)
map_mod4$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 4

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod4 <- list(version="mod4", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod4, map=map_mod4, random=random)
fit_mod4 <- fit_pk(input=input_mod4, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=2)
peels_mod4 <- fit_pk_retros(fit_mod4, peels = 1:5)
peels_mod4 <- c(list(fit_mod4), peels_mod4)


## MODEL 5 ----
# - Non-parametric by age (No random effects)
# - Turn on mean_sel
map_mod5 <- map

# -- Fixed effect pars
map_mod5$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 5

# - Build model

# - Run model 1
input_mod5 <- list(version="mod5", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars, map=map_mod5, random=random)
fit_mod5 <- fit_pk(input=input_mod5, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=2)
peels_mod5 <- fit_pk_retros(fit_mod5, peels = 1:5)
peels_mod5 <- c(list(fit_mod5), peels_mod5)


## MODEL 6 ----
# - Non-parametric 1D-AR1 year random effects
# - Turn on mean_sel, rho_y, sd, and ranef vector
map_mod6 <- map
pars_mod6 <- pars

# -- Random effect pars
pars_mod6$selpars_re <- matrix(0, nrow = 1, ncol = dat$nyrs + dat$projfsh_nyrs)
map_mod6$selpars_re <- as.factor(1:length(pars_mod6$selpars_re))

# -- Fixed effect pars
map_mod6$ln_sel_sd <- as.factor(1)
map_mod6$sel_rho_y <- as.factor(1)
map_mod6$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 6

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod6 <- list(version="mod6", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod6, map=map_mod6, random=random)
fit_mod6 <- fit_pk(input=input_mod6, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=3)
peels_mod6 <- fit_pk_retros(fit_mod6, peels = 1:5)
peels_mod6 <- c(list(fit_mod6), peels_mod6)


## MODEL 7 ----
# - Non-parametric 2D-AR1 age, year random effects
# - Turn on mean_sel, rho, sd, and ranef vector
map_mod7 <- map
pars_mod7 <- pars

# -- Random effect pars
pars_mod7$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map_mod7$selpars_re <- as.factor(1:length(pars_mod7$selpars_re))

# -- Fixed effect pars
map_mod7$ln_sel_sd <- as.factor(1)
map_mod7$sel_rho_a <- as.factor(1)
map_mod7$sel_rho_y <- as.factor(1)
map_mod7$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 7

# - Build model 7
random <- c("selpars_re")

# - Run model 1
input_mod7 <- list(version="mod7", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod7, map=map_mod7, random=random)
fit_mod7 <- fit_pk(input=input_mod7, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=2)
peels_mod7 <- fit_pk_retros(fit_mod7, peels = 1:5)
peels_mod7 <- c(list(fit_mod7), peels_mod7)


## MODEL 8 ----
# - Non-parametric 3D-AR1 age, year, cohort random effects using conditional var
# - Turn on mean_sel, rho (a,y,c), sd, and ranef vector
map_mod8 <- map
pars_mod8 <- pars

# -- Random effect pars
pars_mod8$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map_mod8$selpars_re <- as.factor(1:length(pars_mod8$selpars_re))

# -- Fixed effect pars
map_mod8$ln_sel_sd <- as.factor(1)
map_mod8$sel_rho_a <- as.factor(1)
map_mod8$sel_rho_y <- as.factor(1)
map_mod8$sel_rho_c <- as.factor(1)
map_mod8$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 8
dat$sel_vartype <- 0

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod8 <- list(version="mod8", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod8, map=map_mod8, random=random)
fit_mod8 <- fit_pk(input=input_mod8, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=2)
peels_mod8 <- fit_pk_retros(fit_mod8, peels = 1:5)
peels_mod8 <- c(list(fit_mod8), peels_mod8)


## MODEL 9 ----
# - Non-parametric 3D-AR1 age, year, cohort random effects using marginal var
# - Turn on mean_sel, rho (a,y,c), sd, and ranef vector
map_mod9 <- map
pars_mod9 <- pars

# -- Random effect pars
pars_mod9$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map_mod9$selpars_re <- as.factor(1:length(pars_mod9$selpars_re))

# -- Fixed effect pars
map_mod9$ln_sel_sd <- as.factor(1)
map_mod9$sel_rho_a <- as.factor(1)
map_mod9$sel_rho_y <- as.factor(1)
map_mod9$sel_rho_c <- as.factor(1)
map_mod9$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 8
dat$sel_vartype <- 1

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod9 <- list(version="mod9", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod9, map=map_mod9, random=random)
fit_mod9 <- fit_pk(input=input_mod9, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=3)
peels_mod9 <- fit_pk_retros(fit_mod9, peels = 1:5)
peels_mod9 <- c(list(fit_mod9), peels_mod9)


## MODEL 10 ----
# - Non-parametric 3D-AR1 age and cohort random effects using conditional var (no year effect)
# - Turn on mean_sel, rho (a,c), sd, and ranef vector
map_mod10 <- map
pars_mod10 <- pars

# -- Random effect pars
pars_mod10$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map_mod10$selpars_re <- as.factor(1:length(pars_mod10$selpars_re))

# -- Fixed effect pars
map_mod10$ln_sel_sd <- as.factor(1)
map_mod10$sel_rho_a <- as.factor(1)
map_mod10$sel_rho_c <- as.factor(1)
map_mod10$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 8
dat$sel_vartype <- 0

# - Build model
random <- c("selpars_re")

# - Run model 1
input_mod10 <- list(version="mod10", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars_mod10, map=map_mod10, random=random)
fit_mod10 <- fit_pk(input=input_mod10, getsd=TRUE, control=control, use_bounds = TRUE, newtonsteps=2)
peels_mod10 <- fit_pk_retros(fit_mod10, peels = 1:5)
peels_mod10 <- c(list(fit_mod10), peels_mod10)


## Save ----
# - Models
fits <- list(fit_mod0=fit_mod0,
             fit_mod1=fit_mod1,
             fit_mod2=fit_mod2,
             fit_mod3=fit_mod3,
             fit_mod4=fit_mod4,
             fit_mod5=fit_mod5,
             fit_mod6=fit_mod6,
             fit_mod7=fit_mod7,
             fit_mod8=fit_mod8,
             fit_mod9=fit_mod9,
             fit_mod10=fit_mod10)
saveRDS(fits, 'TMB/Output/2023_fits.RDS')


# - Peels
peels <- list(peels_mod0=peels_mod0,
             peels_mod1=peels_mod1,
             peels_mod2=peels_mod2,
             peels_mod3=peels_mod3,
             peels_mod4=peels_mod4,
             peels_mod5=peels_mod5,
             peels_mod6=peels_mod6,
             peels_mod7=peels_mod7,
             peels_mod8=peels_mod8,
             peels_mod9=peels_mod9,
             peels_mod10=peels_mod10)
saveRDS(peels, 'TMB/Output/2023_peels.RDS')

