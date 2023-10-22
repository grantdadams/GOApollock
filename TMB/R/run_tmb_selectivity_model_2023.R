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

# - Build model 1
input_mod0 <- list(version="mod0", path="TMB/src/",
                   modfile="goa_pk_tmb",
                   dat=dat, pars=pars, map=map_mod0, random=NULL)
fit_mod0 <- fit_pk(input=input_mod0, getsd=TRUE, control=control, use_bounds = TRUE)
check <- fit_peel(fit_mod0, 1)


# -- Build model obj
lwr <- get_bounds(obj_mod0)$lwr
upr <- get_bounds(obj_mod0)$upr
obj_mod0 <- MakeADFun(data=dat2, parameters=pars2, map=map2, random = NULL, silent=FALSE)

# - Optimize
# opt_mod0 <- with(obj_mod0, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod0 <- TMBhelper::fit_tmb(obj_mod0, control=control, newtonsteps=1, upp)
fit_mod0$obj <- obj_mod0
fit_mod0$path <- "TMB/src/"
fit_mod0$modfile <- "goa_pk_tmb"
fit_peel(fit_mod0, 5)


## MODEL 1 ----
# - Double logistic with random effects on ascending portion
# - Turn on sel sd and slp and inf parameters
map_mod1 <- map

# -- Random effect pars
map_mod1$slp1_fsh_dev <- as.factor(1:length(pars$slp1_fsh_dev))

map_mod1$inf1_fsh_dev <- as.factor(1:length(pars$inf1_fsh_dev))
#map_mod1$slp2_fsh_dev <- as.factor(1:length(pars$slp2_fsh_dev))
#map_mod1$inf2_fsh_dev <- as.factor(1:length(pars$inf2_fsh_dev))

# -- Fixed effect pars
map_mod1$ln_sel_sd <- as.factor(1)
# map_mod1$inf1_fsh_mean <- as.factor(1) # removing mean because no zero constraint
map_mod1$inf2_fsh_mean <- as.factor(1)
# map_mod1$log_slp1_fsh_mean <- as.factor(1) # removing mean because no zero constraint
map_mod1$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 1

# - Build model 1
random <- c("slp1_fsh_dev", "inf1_fsh_dev", "slp2_fsh_dev", "inf2_fsh_dev")

# -- Phase
# phased_pars_mod1 <- TMBphase(
#   data = dat,
#   parameters = pars,
#   map = map_mod1,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod1 <- MakeADFun(data=dat, parameters=pars, map=map_mod1, DLL = "goa_pk_tmb", random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod1)$lwr
upr <- get_bounds(obj_mod1)$upr
# opt_mod1 <- with(obj_mod1, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod1 <- TMBhelper::fit_tmb(obj_mod1, control=control, newtonsteps=3)
fit_mod1$obj <- obj_mod1
fit_mod1$path <- "TMB/src/"
fit_mod1$modfile <- "goa_pk_tmb"
fit_peel(fit_mod1, 5)


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

# -- Phase
# phased_pars_mod2 <- TMBphase(
#   data = dat,
#   parameters = pars,
#   map = map_mod2,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod2 <- MakeADFun(data=dat, parameters=pars, map=map_mod2, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod2)$lwr
upr <- get_bounds(obj_mod2)$upr
# opt_mod2 <- with(obj_mod2, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod2 <- TMBhelper::fit_tmb(obj_mod2, control=control, newtonsteps=1)
fit_mod2$obj <- obj_mod2


## MODEL 3 ----
# - Double logistic with AR1 year random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho y
map_mod3 <- map

# -- Random effect pars
pars_mod3 <- pars
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

# -- Phase
# phased_pars_mod3 <- TMBphase(
#   data = dat,
#   parameters = pars_mod3,
#   map = map_mod3,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod3 <- MakeADFun(data=dat, parameters=pars_mod3, map=map_mod3, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod3)$lwr
upr <- get_bounds(obj_mod3)$upr
# opt_mod3 <- with(obj_mod3, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod3 <- TMBhelper::fit_tmb(obj_mod3, loopnum=3, control=control, newtonsteps=2)
fit_mod3$obj <- obj_mod3

## MODEL 4 ----
# - Double logistic with 2D-AR1 age, year random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho
map_mod4 <- map
pars_mod4 <- pars

# -- Random effect pars
pars_mod4$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
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

# -- Phase
# phased_pars_mod4 <- TMBphase(
#   data = dat,
#   parameters = pars_mod4,
#   map = map_mod4,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod4 <- MakeADFun(data=dat, parameters=pars_mod4, map=map_mod4, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod4)$lwr
upr <- get_bounds(obj_mod4)$upr
# opt_mod4 <- with(obj_mod4, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod4 <- TMBhelper::fit_tmb(obj_mod4, loopnum=3, control=control, newtonsteps=2)
fit_mod4$obj <- obj_mod4


## MODEL 5 ----
# - Non-parametric by age (No random effects)
# - Turn on mean_sel
map_mod5 <- map

# -- Fixed effect pars
map_mod5$mean_sel <- as.factor(1:dat$nages)

# -- Data switch
dat$seltype <- 5

# - Build model

# -- Phase
# phased_pars_mod5 <- TMBphase(
#   data = dat,
#   parameters = pars,
#   map = map_mod5,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod5 <- MakeADFun(data=dat, parameters=pars, map=map_mod5, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod5)$lwr
upr <- get_bounds(obj_mod5)$upr
# opt_mod5 <- with(obj_mod5, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod5 <- TMBhelper::fit_tmb(obj_mod5, loopnum=3, control=control, newtonsteps=2)
fit_mod5$obj <- obj_mod5


## MODEL 6 ----
# - Non-parametric 1D-AR1 year random effects
# - Turn on mean_sel, rho_y, sd, and ranef vector
map_mod6 <- map

# -- Random effect pars
pars_mod6 <- pars
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

# -- Phase
# phased_pars_mod6 <- TMBphase(
#   data = dat,
#   parameters = pars_mod6,
#   map = map_mod6,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod6 <- MakeADFun(data=dat, parameters=pars_mod6, map=map_mod6, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod6)$lwr
upr <- get_bounds(obj_mod6)$upr
# opt_mod6 <- with(obj_mod6, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod6 <- TMBhelper::fit_tmb(obj_mod6, loopnum=3, control=control, newtonsteps=3)
fit_mod6$obj <- obj_mod6


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

# -- Phase
# phased_pars_mod7 <- TMBphase(
#   data = dat,
#   parameters = pars_mod7,
#   map = map_mod7,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod7 <- MakeADFun(data=dat, parameters=pars_mod7, map=map_mod7, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod7)$lwr
upr <- get_bounds(obj_mod7)$upr
# opt_mod7 <- with(obj_mod7, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod7 <- TMBhelper::fit_tmb(obj_mod7, loopnum=3, control=control, newtonsteps=2)
fit_mod7$obj <- obj_mod7


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

# -- Phase
# phased_pars_mod8 <- TMBphase(
#   data = dat,
#   parameters = pars_mod8,
#   map = map_mod8,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod8 <- MakeADFun(data=dat, parameters=pars_mod8, map=map_mod8, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod8)$lwr
upr <- get_bounds(obj_mod8)$upr
#opt_mod8 <- with(obj_mod8, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod8 <- fit_tmb(obj_mod8, loopnum=3, control=control, newtonsteps=2)
fit_mod8$obj <- obj_mod8


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

# -- Phase
# phased_pars_mod9 <- TMBphase(
#   data = dat,
#   parameters = pars_mod9,
#   map = map_mod9,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod9 <- MakeADFun(data=dat, parameters=pars_mod9, map=map_mod9, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod9)$lwr
upr <- get_bounds(obj_mod9)$upr
# opt_mod9 <- with(obj_mod9, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod9 <- fit_tmb(obj_mod9, loopnum=3, control=control, newtonsteps=0, getSD=FALSE)
fit_mod9$obj <- obj_mod9


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

# -- Phase
# phased_pars_mod10 <- TMBphase(
#   data = dat,
#   parameters = pars_mod10,
#   map = map_mod10,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod10 <- MakeADFun(data=dat, parameters=pars_mod10, map=map_mod10, random=random, silent=TRUE)

# - Optimize
lwr <- get_bounds(obj_mod10)$lwr
upr <- get_bounds(obj_mod10)$upr
#opt_mod10 <- with(obj_mod10, nlminb(par,fn, gr, control = control, lower=lwr, upper=upr))
fit_mod10 <- fit_tmb(obj_mod10, loopnum=3, control=control, newtonsteps=2)
fit_mod10$obj <- obj_mod10


## Save Image ----
# save.image(file='TMB/Selectivity_runs.RData')
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

