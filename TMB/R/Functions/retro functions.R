
#' Fit a sequence of models to peeled data set for retrospective
#' calculations
#'
#' @param obj A fitted model from \code{fit_pk}.
#' @param peels Vector of peels to fit, with 0 being the original
#'   data set
#' @param getsd Whether to run sdreport
#' @param ... Additional arguments to pass to \link{\code{fit_pk}}
#' @return A list of model fits of class 'pkfit'
#' @details This function fits a series of models to peels based
#'   off an original fit. The data, parameter, and map lists are
#'   modified based on the peel in a way that the original MLEs
#'   are preserved so that it runs faster. Pass an unfitted obj
#'   if this is undesired behavior.
#' @export
fit_pk_retros <- function(fit, peels=0:7, getsd=TRUE, ...){
  if(class(fit)[1]!='pkfit') stop("fit argument is not a fitted model")
  retros <- lapply(peels, function(i) fit_peel(fit, peel=i, getsd=getsd, ...))
  return(retros)
}

#' Modify parameter list to match peel
peel_pars <- function(pars, dat,  peel){
  p <- pars
  stopifnot(peel>=0)
  nyrs <- length(p$dev_log_recruit)
  ind <- 1:(nyrs-peel)

  # Parameter vectors
  x <- c("dev_log_recruit", "dev_log_F", "log_q1_dev",
         "log_q2_dev", "log_q3_dev")
  for(i in x) p[[i]] <- p[[i]][ind]

  # Fish par (has projection years)
  # - vectors
  nyrs_fsh <- length(p$inf1_fsh_dev)
  ind_fsh <- 1:(nyrs_fsh-peel)
  x_fsh <- c("slp1_fsh_dev", "inf1_fsh_dev",
             "slp2_fsh_dev", "inf2_fsh_dev")
  for(i in x_fsh) p[[i]] <- p[[i]][ind_fsh]

  # - matrices
  if(dat$seltype %in% c(4, 7, 8)){
    x = "selpars_re"
    for(i in x) p[[i]] <- p[[i]][,ind_fsh]
  }

  if(peel==0) stopifnot(all.equal(p,pars))

  return(p)
}

#' Modify map to match peel
peel_map <- function(map, pars){
  ## tricky part is map elements may or may not be specified so
  ## do this generically
  m <- map
  p <- pars
  stopifnot(is.list(map))
  ## ind <- 1:length(yrs)
  ## nyrs <- length(ind)
  ## for(i in 1:length(m)){
  ##   if(length(m[[i]])>nyrs) m[[i]] <- m[[i]][ind]
  ## }
  ## if(any(sapply(m, NROW)>nyrs))
  ##   stop("Some pars too long in peel ",peel)
  ##  if(peel==0) stopifnot(all.equal(m,map))
  for(i in names(m)){
    #print(i)
    #if(i=='log_q2_dev') browser()
    m[[i]] <- as.factor(
      as.numeric(m[[i]][1:length(p[[i]])]
      )
    )
  }
  return(m)
}

#' Modify data list to match peel
#'
peel_data <- function(dat, peel){
  # Checks
  stopifnot(peel>=0)
  stopifnot(is.list(dat))
  d <- dat

  # Years
  d$endyr <- as.integer(dat$endyr-peel)
  d$yrs <- dat$styr:d$endyr
  d$nyrs <- length(d$yrs)

  # Fishery data
  # - Catch
  i0 <- d$yrs-dat$styr+1
  d$cattot <- d$cattot[i0]
  d$cattot_log_sd <- d$cattot_log_sd[i0]
  d$wt_fsh <- d$wt_fsh[i0,]

  # - Catch-at-age
  i1 <- which(d$fshyrs<=d$endyr)
  d$fshyrs <- d$fshyrs[i1]                # age comp yrs
  d$nyrs_fsh <- length(i1)
  d$multN_fsh <- d$multN_fsh[i1]
  d$ac_yng_fsh <- d$ac_yng_fsh[i1]
  d$ac_old_fsh <- d$ac_old_fsh[i1]
  d$catp <- d$catp[i1,]

  # - Catch-at-length
  i2 <- which(d$fshlenyrs<=d$endyr)
  d$fshlenyrs <- d$fshlenyrs[i2]
  d$multNlen_fsh <- d$multNlen_fsh[i2]
  d$rwlk_sd <- d$rwlk_sd[i0] # Unused
  d$lenp <- d$lenp[i2,]

  # Survey 1  (Acoustic) EK500 (biosonic deleted in 2022)
  # - Index
  i3 <- which(d$srvyrs1 <=d$endyr)
  d$srvyrs1 <- d$srvyrs1[i3]
  d$nyrs_srv1 <- length(d$srvyrs1)
  d$indxsurv1 <- d$indxsurv1[i3]
  d$indxsurv_log_sd1 <- d$indxsurv_log_sd1[i3]
  d$q1_rwlk_sd <- d$q1_rwlk_sd[i0]
  d$yrfrct_srv1 <- d$yrfrct_srv1[i0]
  d$wt_srv1 <- d$wt_srv1[i0,]

  # - Age comp
  i4 <- which(d$srv_acyrs1<=d$endyr)
  d$srv_acyrs1 <- d$srv_acyrs1[i4]
  d$nyrsac_srv1 <- length(d$srv_acyrs1)
  d$multN_srv1 <- d$multN_srv1[i4]
  d$ac_yng_srv1 <- d$ac_yng_srv1[i4]
  d$ac_old_srv1 <- d$ac_old_srv1[i4]
  d$est_q1_dev <- d$yrs %in% do.call(seq, as.list(range(d$srv_acyrs1)))

  # - Length comp
  i5 <- which(d$srv_lenyrs1<=d$endyr)
  d$srv_lenyrs1 <- d$srv_lenyrs1[i5]
  d$nyrslen_srv1 <- length(d$srv_lenyrs1)
  d$multNlen_srv1 <- d$multNlen_srv1[i5]
  d$srvp1 <- d$srvp1[i4,]
  d$srvlenp1 <- d$srvlenp1[i5,]

  # Survey 2 (Bottom trawl)
  # - Index
  i3 <- which(d$srvyrs2 <=d$endyr)
  d$srvyrs2 <- d$srvyrs2[i3]
  d$nyrs_srv2 <- length(d$srvyrs2)
  d$indxsurv2 <- d$indxsurv2[i3]
  d$indxsurv_log_sd2 <- d$indxsurv_log_sd2[i3]
  d$q2_rwlk_sd <- d$q2_rwlk_sd[i0]
  d$yrfrct_srv2 <- d$yrfrct_srv2[i0]
  d$wt_srv2 <- d$wt_srv2[i0,]

  # - Age comp
  i4 <- which(d$srv_acyrs2<=d$endyr)
  d$srv_acyrs2 <- d$srv_acyrs2[i4]
  d$nyrsac_srv2 <- length(d$srv_acyrs2)
  d$multN_srv2 <- d$multN_srv2[i4]
  d$ac_yng_srv2 <- d$ac_yng_srv2[i4]
  d$ac_old_srv2 <- d$ac_old_srv2[i4]
  d$srvp2 <- d$srvp2[i4,]
  d$est_q2_dev <- d$yrs %in% do.call(seq, as.list(range(d$srv_acyrs2)))

  # - Length comp
  i5 <- which(d$srv_lenyrs2<=d$endyr)
  d$srv_lenyrs2 <- d$srv_lenyrs2[i5]
  d$nyrslen_srv2 <- length(d$srv_lenyrs2)
  d$multNlen_srv2 <- d$multNlen_srv2[i5]
  d$srvlenp2 <- d$srvlenp2[i5,]

  # Survey 3 (ADFG coastal survey)
  # - Index
  i3 <- which(d$srvyrs3 <= d$endyr)
  d$srvyrs3 <- d$srvyrs3[i3]
  d$nyrs_srv3 <- length(d$srvyrs3)
  d$indxsurv3 <- d$indxsurv3[i3]
  d$indxsurv_log_sd3 <- d$indxsurv_log_sd3[i3]
  d$q3_rwlk_sd <- d$q3_rwlk_sd[i0]
  d$yrfrct_srv3 <- d$yrfrct_srv3[i0]
  d$wt_srv3 <- d$wt_srv3[i0,]

  # - Age comp
  i4 <- which(d$srv_acyrs3<=d$endyr)
  d$srv_acyrs3 <- d$srv_acyrs3[i4]
  d$nyrsac_srv3 <- length(d$srv_acyrs3)
  d$multN_srv3 <- d$multN_srv3[i4]
  d$ac_yng_srv3 <- d$ac_yng_srv3[i4]
  d$ac_old_srv3 <- d$ac_old_srv3[i4]
  d$srvp3 <- d$srvp3[i4,]
  d$est_q3_dev <- d$yrs %in% do.call(seq, as.list(range(d$srv_acyrs3)))

  # - Length comp
  i5 <- which(d$srv_lenyrs3<=d$endyr)
  d$srv_lenyrs3 <- d$srv_lenyrs3[i5]
  d$nyrslen_srv3 <- length(d$srv_lenyrs3)
  d$multNlen_srv3 <- d$multNlen_srv3[i5]
  d$srvlenp3 <- d$srvlenp3[i5,]

  # Survey 6 (Summer acoustic)
  # - Index
  i3 <- which(d$srvyrs6 <=d$endyr)
  d$srvyrs6 <- d$srvyrs6[i3]
  d$nyrs_srv6 <- length(d$srvyrs6)
  d$indxsurv6 <- d$indxsurv6[i3]
  d$indxsurv_log_sd6 <- d$indxsurv_log_sd6[i3]
  d$q6_rwlk_sd <- d$q6_rwlk_sd[i0]
  d$yrfrct_srv6 <- d$yrfrct_srv6[i0]
  d$wt_srv6 <- d$wt_srv6[i0,]

  # - Age comp
  i4 <- which(d$srv_acyrs6<=d$endyr)
  d$srv_acyrs6 <- d$srv_acyrs6[i4]
  d$nyrsac_srv6 <- length(d$srv_acyrs6)
  d$multN_srv6 <- d$multN_srv6[i4]
  d$ac_yng_srv6 <- d$ac_yng_srv6[i4]
  d$ac_old_srv6 <- d$ac_old_srv6[i4]
  d$srvp6 <- d$srvp6[i4,]

  # - Length comp
  i5 <- which(d$srv_lenyrs6<=d$endyr)
  d$srv_lenyrs6 <- d$srv_lenyrs6[i5]
  d$nyrslen_srv6 <- length(d$srv_lenyrs6)
  d$multNlen_srv6 <- d$multNlen_srv6[i5]
  d$srvlenp6 <- d$srvlenp6[i5,, drop=FALSE]

  # Survey 4 (Age 1 acoustic)
  i3 <- which(d$srvyrs4 <=d$endyr)
  d$srvyrs4 <- d$srvyrs4[i3]
  d$nyrs_srv4 <- length(d$srvyrs4)
  d$indxsurv4 <- d$indxsurv4[i3]
  d$indxsurv_log_sd4 <- d$indxsurv_log_sd4[i3]

  # Survey 5 (Age 2 acoustic)
  i3 <- which(d$srvyrs5 <=d$endyr)
  d$srvyrs5 <- d$srvyrs5[i3]
  d$nyrs_srv5 <- length(d$srvyrs5)
  d$indxsurv5 <- d$indxsurv5[i3]
  d$indxsurv_log_sd5 <- d$indxsurv_log_sd5[i3]



  d$ay_Index <- as.matrix(expand.grid("age" = seq_len(d$nages),
                                      "year" = seq_len(d$nyrs + d$projfsh_nyrs) ))

  if(peel==0) stopifnot(all.equal(d,dat))
  ## this breaks makeadfun for some reason if not done??
  ##   for(i in 1:length(dat)){
  ## ##     browser()
  ##     print(cbind(i=i, names(dat)[i], class(dat[[i]])[1], class(d[[i]])[1]))
  ##     ## if(class(dat[[i]])[1]=='integer'){
  ##     ##   message("converting ", names(dat)[i])
  ##     ##   class(d[[i]]) <- 'integer'
  ##     ## }
  ##     class(d[[i]])[1] <- class(dat[[i]])[1]
  ##     print(cbind(i=i, names(dat)[i], class(dat[[i]])[1], class(d[[i]])[1]))
  ##   }
  ## attributes(d) <- attributes(dat)
  ## attributes(d) <- attributes(dat)
  ## attributes(d)$checks.passed <- NULL
  return(d)
}

#' Internal wrapper function to fit a single peel
#'
fit_peel <- function(fit, peel, getsd=TRUE, ...){
  control <- list(eval.max=10000, iter.max=10000, trace=0)
  stopifnot(peel>=0)
  dat2 <- peel_data(dat=fit$obj$env$data, peel)
  attributes(dat2) <- attributes(fit$obj$env$data)
  attributes(dat2)$check.passed <- NULL
  pars2 <- peel_pars(pars=fit$obj$env$parList(), dat=dat2, peel)
  map2 <- peel_map(map=fit$obj$env$map, pars2)
  input2 <- list(version=paste0('peel',peel), path=fit$path,
                 modfile=fit$modfile,
                 dat=dat2, pars=pars2, map=map2, random=fit$obj$env$random)
  message("Starting optimization for peel=",peel)
  fit <- fit_pk(input=input2, getsd=getsd, control=control, ...)
  return(fit)
}

