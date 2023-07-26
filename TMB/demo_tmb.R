library(GOApollock)
library(TMBhelper)
library(ggplot2)
library(TMB)


## The dat and rep files from the 2022 final model
dat <- readRDS('TMB/datfile.RDS')
replist <- readRDS("TMB/repfile.RDS")
stdadmb <- readRDS("TMB/stdfile.RDS")
pars <- readRDS('TMB/pars.RDS')
map <- readRDS('TMB/map.RDS')

# Map with parameters turned on
map0 <- lapply(map, function(x) as.numeric(x))
map0 <- lapply(map0, function(x) replace(x, values = 1:length(x)))
map0$natMscalar <- NA
map0$sigmaR <- NA
map0 <- lapply(map0, function(x) as.factor(x))
map <- map0

## Fit the model and compare to ADMB. Should be identical
compile("source/goa_pk_tmb.cpp", flag='-w')
dyn.load('source/goa_pk_tmb.dll')
obj <- MakeADFun(data=dat, parameters=pars, map=map, random=NULL, silent=FALSE)

## optimize and compare
opt <- with(obj, nlminb(par,fn,gr))
opt <- TMBhelper::fit_tmb(obj, control=list(trace=50))
sdrep <- sdreport(obj)
stdtmb <- with(sdrep, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup

ssb <- rbind(cbind(model='TMB',filter(stdtmb, par=='Espawnbio')),
             cbind(model='ADBM',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))

filter(ssb, year>2015) %>% arrange(year)

ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

