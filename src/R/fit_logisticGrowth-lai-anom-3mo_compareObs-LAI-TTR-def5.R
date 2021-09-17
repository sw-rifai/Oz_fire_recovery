library(tidyverse);
library(stars);
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table 

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
sdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")

# some smoothed lai values got slightly negative, so resetting to ~0.01
dat[,`:=`(slai=ifelse(slai<0.01,0.01,slai))]
dat[,`:=`(slai_12mo=slai_anom+malai)]

#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5_lai)==FALSE][date_fire1<ymd('2015-03-01')][ttr5_lai>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_fire1,ttr5_lai)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>%
  mutate(recovery_date = date_fire1+days(ttr5_lai)) %>% 
  group_by(x,y,id) %>% 
  filter(date >= date_fire1) %>% 
  filter(date <= recovery_date + years(3)) %>% 
  ungroup() %>%
  mutate(post_days = as.double(date - date_fire1)) %>% 
  as.data.table() 
mdat[,`:=`(slai_12mo = slai_anom_12mo+malai)]
mdat[,`:=`(slai_1mo = slai_anom+malai)]
mdat[,`:=`(slai_1mo = ifelse(slai_1mo<0.01,0.01,slai_1mo))]

mdat[,`:=`(slai_3mo = slai_anom_3mo+malai)]
mdat[,`:=`(slai_3mo = ifelse(slai_3mo<0.01,0.01,slai_3mo))]

rm(dat); gc(full=TRUE)



fn_logistic_growth_v0 <- function(din){
  start_day <- din[post_days <= 366][slai_1mo == min(slai_1mo)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+365/2 ))]
  upper_K <- din$malai[1]+1*din$lai_yr_sd[1]
  lower_K <- din$malai[1]-1*din$lai_yr_sd[1]
  lower_K <- max(c(0.4,lower_K),na.rm=TRUE)
  min_slai_anom <- din[post_days<=366][slai_anom==min(slai_anom,na.rm=T)]$slai_anom
  date_min_slai_anom <- din[post_days<=366][slai_anom==min(slai_anom,na.rm=T)]$date
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  upper_L0 <- din$malai[1]
  lower_L0 <- 0.01
  
  try(fit <- nls_multstart(slai_1mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
                           data=din,
                           # iter=1,
                           iter=20,
                           supp_errors = 'Y',
                           start_lower = c(K=0.1*lower_K, L0=0.01, r=0),
                           start_upper = c(K=0.9*upper_K, L0=0.9*upper_K, r=0.001), 
                           lower= c(K=lower_K, L0=lower_L0, r=0.0001), 
                           upper = c(K=upper_K, 
                                     L0=upper_L0, 
                                     r=0.2))
      ,silent = TRUE)
  if(exists('fit')==FALSE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
    out$start_day <- start_day
    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai_1mo, 
                                      estimate = predict(fit))
    out$rmse <- yardstick::rmse_vec(truth = din$slai_1mo, 
                                    estimate = predict(fit))
    out$min_slai_anom <- min_slai_anom
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
    
  },silent=TRUE)
  out$nobs_til_recovery <- nrow(din)
  return(out)
}

# mdat <- arrow::read_parquet("/home/sami/scratch/mdat.parquet")

get_rmse <- function(fit){
  response <- fit$m$lhs()
  pred <- fit$m$predict()
  rmse <- yardstick::rmse_vec(truth=response, 
    estimate=pred, na_rm=T)
  return(rmse)
}

alt_nls <- function (formula, data = parent.frame(), iter, start_lower, 
    start_upper, supp_errors = c("Y", "N"), convergence_count = 100, 
    control, modelweights, ...) 
{
    if (missing(supp_errors)) {
        supp_errors <- "N"
    }
    formula <- stats::as.formula(formula)
    if (any(class(data) %in% c("data.frame", "list", "environment"))) {
        params_ind <- all.vars(formula[[3]])[all.vars(formula[[3]]) %in% 
            names(data)]
        params_est <- all.vars(formula[[3]])[!all.vars(formula[[3]]) %in% 
            names(data)]
        params_dep <- all.vars(formula[[2]])
    }
    else {
        stop("data should be a data.frame, list or an environment")
    }
    if (missing(start_lower) || missing(start_upper)) {
        cat("No boundaries specified for the starting values of sought parameters. \n", 
            "Default values of +/- 1e+10 will be used. This is likely \n", 
            "to slow the process of finding the best model. \n")
        r <- readline("Continue with default values [y/n]? ")
        if (tolower(r) == "n") {
            stop("Please enter upper and lower parameter boundaries as start_lower and start_upper in function argument.")
        }
    }
    if (missing(start_lower)) 
        start_lower <- rep(-10^10, length(params_est))
    if (missing(start_upper)) 
        start_upper <- rep(10^10, length(params_est))
    if (length(start_lower) != length(params_est) || length(start_upper) != 
        length(params_est)) {
        stop("There must be as many parameter starting bounds as there are parameters")
    }
    params_bds <- data.frame(param = params_est, low.bds = unlist(start_lower), 
        high.bds = unlist(start_upper), stringsAsFactors = FALSE)
    if (missing(control)) {
        control <- minpack.lm::nls.lm.control(maxiter = 1000, 
            ftol = .Machine$double.eps, ptol = .Machine$double.eps)
    }
    silent <- ifelse(supp_errors == "Y", TRUE, FALSE)
    if (silent == TRUE) {
        oo <- options(warn = -1)
        on.exit(options(oo))
    }
    if ("modelweights" %in% all.vars(formula)) {
        stop(paste0("The variable name 'modelweights' is reserved for model weights. Please change the name\n", 
            "of this variable"))
    }
    if (missing(modelweights)) {
        data$modelweights <- rep(1, length(data[[params_dep]]))
    }
    else {
        data$modelweights <- eval(substitute(modelweights), data)
    }
    if (length(iter) == 1) {
        multistart_type <- "shotgun"
    }
    else if (length(iter) == nrow(params_bds)) {
        multistart_type <- "gridstart"
    }
    else {
        stop(paste0("iter should be of length 1 for shotgun approach and of the same length as the\n", 
            "number of parameters for the gridstart approach."))
    }
    if (multistart_type == "gridstart" && convergence_count != 
        FALSE) {
        warning("A gridstart approach cannot be applied with convergence_count. Convergence count will be set to FALSE")
        convergence_count <- FALSE
    }
    if (multistart_type == "shotgun") {
        strt <- purrr::map2(params_bds$low.bds, params_bds$high.bds, 
            ~runif(iter, .x, .y))
        names(strt) <- params_bds$param
        strt <- dplyr::bind_rows(strt)
    }
    if (multistart_type == "gridstart") {
        params_bds$iter <- iter
        strt <- purrr::pmap(as.list(params_bds[, -1]), function(low.bds, 
            high.bds, iter) seq(from = low.bds, to = high.bds, 
            length.out = iter))
        names(strt) <- params_bds$param
        strt <- tibble::as_tibble(expand.grid(strt))
    }
    if (multistart_type == "shotgun" && convergence_count != 
        FALSE) {
        count <- 0
        stored_AIC <- Inf
        fit <- NULL
        fit_best <- NULL
        for (j in 1:nrow(strt)) {
            start.vals <- as.list(strt[j, ])
            try(fit <- minpack.lm::nlsLM(formula, start = start.vals, 
                control = control, data = data, weights = modelweights, 
                ...), silent = silent)
            if (is.null(fit) && is.null(fit_best)) {
                count <- 0
            }
            else {
                count <- ifelse(stored_AIC <= get_rmse(fit), 
                  count + 1, 0)
            }
            if (count == convergence_count) 
                break
            if (!is.null(fit) && stored_AIC > get_rmse(fit)) {
                stored_AIC <- get_rmse(fit)
                fit_best <- fit
            }
        }
    }
    
    if (convergence_count == FALSE) {
        strt$iteration <- 1:nrow(strt)
        allfits <- tidyr::nest(strt, startpars = -iteration)
        fit <- NULL
        fit_best <- NULL
        fit_aic <- function(startpars) {
            start.vals <- as.list(startpars[[1]])
            try(fit <- minpack.lm::nlsLM(formula, start = start.vals, 
                control = control, data = data, weights = modelweights, 
                ...), silent = silent)
            AICval <- ifelse(!is.null(fit), stats::get_rmse(fit), 
                Inf)
            return(AICval)
        }
        allfits <- dplyr::group_by(allfits, iteration)
        allfits <- dplyr::mutate(allfits, AICval = purrr::map_dbl(startpars, 
            ~fit_aic(startpars)))
        allfits <- dplyr::arrange(allfits, AICval)
        fit_best <- minpack.lm::nlsLM(formula, start = allfits$startpars[[1]], 
            control = control, data = data, weights = modelweights, 
            ...)
    }
    return(fit_best)
}

fn_logistic_growth <- function(din){
  # notes: 
  # K must be >= malai
  start_day <- din[post_days <= 366][slai_3mo == min(slai_3mo,na.rm=T)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+183))]
  min_slai_anom <- din[post_days<=366][near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_anom_3mo
  min_slai <- din[post_days<=366][near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_3mo
  malai <- din$malai[1]
  lai_yr_sd <- din$lai_yr_sd[1]
  offset <- abs(min_slai_anom + malai)
  date_min_slai_anom <- din[post_days<=366][slai_anom_3mo==min(slai_anom_3mo,na.rm=T)]$date
  day_ttr_offset <- as.numeric(date_min_slai_anom - first(din$date_fire1))
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  # upper_K <- din$malai[1] + 2*din$lai_yr_sd[1]
  # lower_K <- din$malai[1] - 2*din$lai_yr_sd[1]
  # lower_K <- ifelse(lower_K < 0.4, 0.4, lower_K)
  # lower_K <- max(c(0.4,lower_K),na.rm=TRUE)
  upper_K <- din$malai[1]+0.25*din$lai_yr_sd[1]
  lower_K <- din$malai[1]
  lower_K <- ifelse(upper_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  lower_K <- ifelse(lower_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  # upper_L0 <- din$malai[1]
  # lower_L0 <- 0.001
  # lower_L0 <- malai+min_slai_anom-0.25*din$lai_yr_sd[1]
  # lower_L0 <- max(c(0.001, lower_L0))
  # upper_L0 <- min_slai+1*din$lai_yr_sd[1]
  upper_L0 <- min_slai
  lower_L0 <- min_slai
  # din[,slai_1mo := slai_1mo+offset+lower_L0]
  b_coefs <- coef(fastLm(X=cbind(1,din$post_days),y=din$slai_3mo))

  try(fit <- alt_nls(slai_3mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
                           data=din,
                           # iter=1,
                           iter=1000,# was 20
                           # convergence_count = 100,
                           # control=nls.lm.control(),
                           supp_errors = 'Y',
                           start_lower = c(K=0.1*lower_K, L0=0.01, r=0),
                           start_upper = c(K=0.9*upper_K, L0=0.9*upper_K, r=0.001), 
                           lower= c(K=lower_K, 
                             L0=lower_L0, 
                             r=0.000025), 
                           upper = c(K=upper_K, 
                                     L0=upper_L0, 
                                     r=0.5))
      ,silent = TRUE)
  
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)

    if(exists('fit')==FALSE){
    out <- data.table(K=NA_real_,
      L0=NA_real_,
      r=NA_real_,
      isConv=FALSE,
      start_day=NA_real_,
      r2=NA_real_,
      rmse=NA_real_)
    }
  
  
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
    out$start_day <- start_day
    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai_1mo, 
                                      estimate = predict(fit))
    out$rmse <- yardstick::rmse_vec(truth = din$slai_1mo, 
                                    estimate = predict(fit))
  },silent=TRUE)
    out$min_slai_anom <- min_slai_anom
    out$min_slai <- min_slai
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
    out$lai_ma <- malai
  out$nobs_til_recovery <- nrow(din)
  out <- out[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
  out$b0 <- b_coefs[1]
  out$b1 <- b_coefs[2]
  return(out)
}


# # data.table approach -----------------------------------
# grpn <- uniqueN(mdat$id)
# pb <- txtProgressBar(min = 0, max = grpn, style = 3)
# out <- mdat[,{setTxtProgressBar(pb, .GRP); fn_w(.SD)}, by=.(x,y,id)]
# close(pb)
# arrow::write_parquet(merge(out, sdat, by=c("x","y","id")), 
#                      sink=paste0("outputs/weibull_fits_1burn_2001-2014fires_",Sys.time(),".parquet"))

# furrr approach -----------------------------------------------
vec_ids <- unique(mdat$id)
vec1 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[1]]
vec2 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[2]]
vec3 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[3]]
vec4 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[4]]


gc(full=TRUE)
plan(multisession, workers=20)
system.time(out1 <- mdat[id%in%vec1] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),
                .progress = TRUE, 
                .options = furrr_options(seed=333L)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out1)
out1[,`:=`(id=as.integer(id))]
plan(sequential)
plan(multisession, workers=20)
system.time(out2 <- mdat[id%in%vec2] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),
                .progress = TRUE, 
                .options = furrr_options(seed=333L)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out2)
out2[,`:=`(id=as.integer(id))]
plan(sequential)
plan(multisession, workers=20)
system.time(out3 <- mdat[id%in%vec3] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),
                .progress = TRUE, 
                .options = furrr_options(seed=333L)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out3)
out3[,`:=`(id=as.integer(id))]
plan(sequential)
plan(multisession, workers=20)
system.time(out4 <- mdat[id%in%vec4] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),
                .progress = TRUE, 
                .options = furrr_options(seed=333L)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out4)
out4[,`:=`(id=as.integer(id))]
plan(sequential)
fits <- rbindlist(list(out1,out2,out3,out4),use.names = TRUE)
gc(full=T)
arrow::write_parquet(merge(fits, sdat, by=c("id")), 
                     sink=paste0("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_",Sys.time(),".parquet"))
# END ****************************************************************

# # Load if not refitting
# # out <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-04-05 18:32:01.parquet")
# # out2 <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-03-31 12:48:06.parquet")
# out$r2 %>% hist
# out$Drop %>% hist
# out$lrc %>% hist
# 
# sdat <- sdat[,`:=`(id=as.integer(id))]
# out <- out[,`:=`(id=as.integer(id))]
# test <- merge(out, sdat, by=c("id")) %>% as.data.table
# test2 <-   test[Drop>0][r2>0.6][between(date_fire1,ymd("2002-10-01"),ymd("2003-03-01"))]#[sample(.N,5000)]
# test3 <- expand_grid(
#   test2,
#   # test %>% filter(between(r2,0.9,0.91)) %>% tb,
#   # test %>% filter(lrc <= -200) %>% tb,
#   post_days=floor(seq(1,3000,length.out=300))) %>% 
#   mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
#   mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
#   arrange(Drop) %>% 
#   mutate(recovered = ifelse(pred >= -0.1,1,0))
# 
# 
# test[isConv==TRUE][Drop>0][r2>0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
#   .[sample(.N,10000)] %>% 
#   expand_grid(
#     .,
#     post_days=floor(seq(1,3000,length.out=300))) %>% 
#   mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>%  
#   # mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
#   arrange(Drop) %>% 
#   mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days, na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% #pull(ttr_w) %>% summary
#   filter(ttr_w >= 365) %>% 
#   ggplot(data=.,aes(ttr5_kn, ttr_w))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_abline()+
#   geom_smooth(col='#CF0000',method='lm')+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='TTR Def 5', 
#        y='Weibull: Time to Recover (days)',
#        title=' Bushfires')+
#   theme_linedraw()
# 
# test[isConv==TRUE][Drop>0][r2>0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
#   .[sample(.N,1000)] %>% 
#   expand_grid(
#     .,
#     post_days=floor(seq(1,3000,length.out=300))) %>% 
#   mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>%  
#   # mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
#   # arrange(Drop) %>% 
#   # mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
#   # filter(recovered==1) %>% 
#   # group_by(id) %>% 
#   # filter(post_days == min(post_days, na.rm=TRUE)) %>% 
#   # ungroup() %>% 
#   # mutate(ttr_w = post_days) %>% #pull(ttr_w) %>% summary
#   # filter(ttr_w >= 365) %>% 
#   ggplot(data=.,aes(post_days,pred,group=id,color=ttr5_kn))+
#   geom_line(lwd=0.1)+
#   # ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='#CF0000')+
#   # geom_abline()+
#   # geom_smooth(col='#CF0000',method='lm')+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='post_days', 
#        y='NDVI anom',
#        title=' Bushfires')+
#   theme_linedraw()+
#   facet_wrap(~cut_number(lrc,4))
# 
# 
# vec_post_days <- sort(unique(mdat$post_days))
# test[isConv==TRUE][Drop>0][r2<0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
#   .[sample(.N,10)] %>% 
#   expand_grid(.,
#     post_days=vec_post_days) %>% 
#   mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
#   left_join(., mdat, by=c('id','post_days')) %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(post_days,pred,group=id,color=id))+
#   geom_point(aes(post_days,kn_anom,color=id,group=id),inherit.aes = F)+
#   geom_line(lwd=1)+
#   geom_hline(aes(yintercept=0),col='#CF0000')+
#   # geom_abline()+
#   # geom_smooth(col='#CF0000',method='lm')+
#   scico::scale_color_scico(end=0.9,palette = 'batlow')+
#   scale_x_continuous(limits=c(0,2500))+
#   labs(x='post_days', 
#        y='NDVI anom',
#        title=' Bushfires')+
#   theme_linedraw()+
#   facet_wrap(~cut_interval(r2,4))
# 
# 
# sum(test$lrc < -10,na.rm=TRUE)/dim(test)[1] # fraction with long recovery
# sum(test$r2 >0.2,na.rm=TRUE)/dim(test)[1] # fraction with reasonable fit
# 
# 
# out[,.(Asym,Drop,lrc,pwr,r2,nobs_til_recovery)][sample(.N,1000)] %>% 
#   GGally::ggpairs()
# 
# 
# test3 %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(ttr_w, ttr))+
#   geom_point()+
#   geom_smooth(method='lm')+
#   geom_abline(col='red')
#   
# scico::scico_palette_show()
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(pre_fire_vi_36mo, ttr_w))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                                select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 36 month NDVI anomaly', 
#        y='Weibull: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_Weibull-TTR-2002-fires_preFire36moVI-anom',Sys.time(),".png"))
# 
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(pre_fire_vi_12mo, ttr_w))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                                select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 12 month NDVI anomaly', 
#        y='Weibull: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_Weibull-TTR-2002-fires_preFire12moVI-anom',Sys.time(),".png"))
# 
# 
# lt <- arrow::read_parquet("outputs/linear_ttr_multiBurns_2001-2020_2021-01-20.parquet")
# library(mgcv)
# lt[is.na(ttr)==F][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))] %>% 
#   ggplot(data=.,aes(pre_fire_vi_36mo, ttr))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                           select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 36 month NDVI anomaly', 
#        y='Linear: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_linear-TTR-2002-fires_preFire36moVI-anom',Sys.time(),".png"))
# 
# 
# lt[is.na(ttr)==F][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))] %>% 
#   ggplot(data=.,aes(pre_fire_vi_12mo, ttr))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                                select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 12 month NDVI anomaly', 
#        y='Linear: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_linear-TTR-2002-fires_preFire12moVI-anom',Sys.time(),".png"))
# 
# 
# 
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.01,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(ttr, ttr_w))+
#   geom_point()+  
#   # geom_abline(col='red')+
#   geom_smooth()
# 
# 
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   lm(ttr_w~scale(Drop)+scale(pre_fire_vi_36mo), data=.) %>% 
#   summary
# 
# 
# test3 %>% 
#   filter(id %in% sample(unique(test3$id),100)) %>% 
#   # filter(lrc < -30) %>% 
#   ggplot(data=.,aes(post_days, pred, color=lrc, group=id))+
#   geom_line(lwd=0.05)+
#   geom_vline(aes(xintercept=365))+
#   scale_color_viridis_c(option='B', end=0.9, limits=c(-100,-20),oob=scales::squish, direction = -1)+
#   scale_x_continuous(limits=c(0,500))
#   # scale_y_continuous(limits=c(0,0.01))
# 
# test3 %>% 
#   ggplot(data=.,aes(post_days, p_diff, color=Drop, group=id))+
#   # geom_point(data=mdat[id%in%test2$id],
#   #            aes(post_days,kn_anom),inherit.aes = F)+
#   geom_line()+
#   scale_color_viridis_c(end=0.9)+
#   scale_y_continuous(limits=c(0,0.01))
# 
# 
# expand_grid(
#   test2,
#   # test %>% filter(between(r2,0.9,0.91)) %>% tb,
#   # test %>% filter(lrc <= -200) %>% tb,
#   post_days=floor(seq(1,3000,length.out=100))) %>% 
#   mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
#   ggplot(data=.,aes(post_days, pred, color=lrc, group=id))+
#   geom_point(data=mdat[id%in%test2$id],
#              aes(post_days,kn_anom),inherit.aes = F)+
#   geom_line()+
#   facet_wrap(~id)+
#   scale_color_viridis_c(end=0.9)
# 
# 
# 
# test %>% sample_n(1000) %>% 
#   ggplot(data=.,aes(ttr,lrc))+
#   geom_point()
# 
# 
# 
# 
