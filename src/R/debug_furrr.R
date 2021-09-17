# arrow::write_parquet(mdat, sink="/home/sami/scratch/mdat.parquet",
#   compression='snappy')
pacman::p_load(tidyverse, data.table, furrr, lubridate)
mdat <- arrow::read_parquet("/home/sami/scratch/mdat.parquet")

library(carrier)
crt1 <- crate(function(din){
  # notes: 
  # K must be >= malai
  library(data.table)
  # library(rethinking)
  # library(RcppArmadillo)
  first <- data.table::first
  `%>%` <- dplyr::'%>%'
  
  start_day <- din[post_days <= 366][slai_3mo == min(slai_3mo,na.rm=T)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+183))]
  min_slai_anom <- din[post_days<=366][dplyr::near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_anom_3mo
  min_slai <- din[post_days<=366][dplyr::near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_3mo
  malai <- din$malai[1]
  lai_yr_sd <- din$lai_yr_sd[1]
  offset <- abs(min_slai_anom + malai)
  date_min_slai_anom <- din[post_days<=366][slai_anom_3mo==min(slai_anom_3mo,na.rm=T)]$date
  day_ttr_offset <- as.numeric(date_min_slai_anom - first(din$date_fire1))
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  upper_K <- din$malai[1]+0.25*din$lai_yr_sd[1]
  lower_K <- din$malai[1]
  lower_K <- ifelse(upper_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  lower_K <- ifelse(lower_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  upper_L0 <- min_slai
  lower_L0 <- min_slai
  b_coefs <- stats::coef(RcppArmadillo::fastLm(X=cbind(1,din$post_days),y=din$slai_3mo))

  flist <- alist(
  # slai_3mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days))
  slai_3mo ~ dnorm(mu, sigma),
  mu <- K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
  K~ dnorm(upper_K,1),
  L0 ~ dnorm(min_slai,1),
  r <- dnorm(b_coefs[2],0.001),
  # r <- dgamma(b_shape=b_coefs[2],scale=1.05),
  # r ~ dgamma(b_coefs[2],b_coefs[2])
  sigma ~ dexp(1)
)
  # set.seed(333)
  suppressMessages(
  try(fit <- rethinking::quap(flist, 
    start=list(K=lower_K,
      L0=lower_L0,
      r=0.001,
      sigma=0.25),
    lower=list(K=lower_K,
      L0=lower_L0,
      r=0, 
      sigma=0.001),
    data=din, 
    hessian = F),
    silent = TRUE)
  )
    
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table::data.table(K=NA_real_,L0=NA_real_,r=NA_real_,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)

    if(exists('fit')==FALSE){
    out <- data.table::data.table(K=NA_real_,
      L0=NA_real_,
      r=NA_real_,
      start_day=NA_real_,
      r2=NA_real_,
      rmse=NA_real_)
    }


  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <-  as.data.frame(t(fit@coef)) %>% data.table::as.data.table()
    # out$isConv <- fit$convInfo$isConv
    out$start_day <- start_day
    test_df <- tidyr::expand_grid(out, post_days=din$post_days) %>%
    dplyr::mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days)))

    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai_3mo,
                                        estimate = test_df$pred)
    out$rmse <- yardstick::rmse_vec(truth = din$slai_3mo,
                                        estimate = test_df$pred)
  },silent=TRUE)
  # 
  
  
    out$lai_ma <- malai
    out$min_slai_anom <- min_slai_anom
    out$min_slai <- min_slai
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
  out$nobs_til_recovery <- nrow(din)
  out <- out[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
  out$b0 <- b_coefs[1]
  out$b1 <- b_coefs[2]
  # return(dgamma2(1, 0.001, scale=1))
  return(out)
})

junk <- sample(unique(mdat$id), 100)
din <- mdat[id==junk[1]]
crt1(mdat[id==junk[10]])



junk <- sample(unique(mdat$id),100)

gc(full=TRUE)
# plan(multisession, workers=20)
plan(sequential)
system.time(out1 <- mdat[id%in%junk] %>% 
              split(.$id) %>%
              future_map(~crt1(.x),
                .progress = TRUE, 
                .options = furrr_options(packages=c(
                  # "carrier",
                  "dplyr",
                  "data.table",
                  "lubridate"
                  # "rethinking"
                  ), 
                    seed=T,
                  globals=list(dgamma2 = rethinking::dgamma2, 
                               crt1 = crt1),
                  lazy=F)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
out1
gc(full=T)
setDT(out1)
out1[,`:=`(id=as.integer(id))]
plan(sequential)
head(out1)


fn_quap_logistic_growth(mdat[id==junk[3]])
crt1(mdat[id==junk[3]])

din <- mdat[id==junk[2]]

out1$r2 %>% summary

out1 %>% ggplot(data=.,aes(L0/K,r))+geom_point()+geom_smooth()
out1 %>% ggplot(data=.,aes(min_slai/lai_ma,r))+geom_point()+geom_smooth()



library(foreach)
library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl) # Register parallel backend for foreach
numCores <- getDoParWorkers()
print(numCores)

#! TESTING
df_out <- foreach(i = 1:length(junk),
  # .combine=rbind, 
  .packages=c('dplyr','data.table','rethinking')) %dopar%{
    crt1(mdat[id==junk[i]])
  }

df_out




junk <- sample(unique(mdat$id),100)
# data.table approach -----------------------------------
grpn <- uniqueN(mdat$id)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out <- mdat[id%in%junk][,{setTxtProgressBar(pb, .GRP); crt1(.SD)}, by=.(x,y,id)]
close(pb)
out

vec_ids <- unique(mdat$id)
vec1 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[1]]
vec2 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[2]]
vec3 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[3]]
vec4 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[4]]

grpn <- uniqueN(mdat$id)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out1 <- mdat[id%in%vec1][,{setTxtProgressBar(pb, .GRP); crt1(.SD)}, by=.(x,y,id)]
close(pb)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out2 <- mdat[id%in%vec2][,{setTxtProgressBar(pb, .GRP); crt1(.SD)}, by=.(x,y,id)]
close(pb)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out3 <- mdat[id%in%vec3][,{setTxtProgressBar(pb, .GRP); crt1(.SD)}, by=.(x,y,id)]
close(pb)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out4 <- mdat[id%in%vec4][,{setTxtProgressBar(pb, .GRP); crt1(.SD)}, by=.(x,y,id)]
close(pb)


out1[b1>0]$r2 %>% summary
out1[L0<K][b1>0][r2>0.3] %>% ggplot(data=.,aes(L0/K,r))+geom_point()+geom_smooth()
out1 %>% ggplot(data=.,aes(min_slai/lai_ma,r))+geom_point()+geom_smooth()





fn_test <- function(x){
  z <- x$z
  lower_L0 <- rnorm(1,sd=0.25)
  lower_K <- z
  upper_K <- z+1
  post_days <- 1:1e5
  K <- z
  L0 <- 0.1*z
  r <- (K-L0)/length(post_days)
  slai_3mo <- K/(1 + ((K-L0)/L0)*exp(-r*post_days)) + rnorm(100, 0.1, 0.1)
  din <- data.frame(slai_3mo)
  
  flist <- alist(
  slai_3mo ~ dnorm(mu, sigma),
  mu <- K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
  K~ dnorm(upper_K,1),
  L0 ~ dnorm(lower_K,1),
  r <- dnorm(r,0.001),
  sigma ~ dexp(1)
 )

  suppressMessages(
  try(fit <- rethinking::quap(flist, 
    start=list(K=lower_K,
      L0=lower_L0,
      r=0.001,
      sigma=0.25),
    lower=list(K=lower_K,
      L0=lower_L0,
      r=0, 
      sigma=0.001),
    data=din,
    hessian = F),
    silent = TRUE)
  )

  if(exists('fit')==F){
    # fout <- data.frame(K=NA_real_,L0=NA_real_,r=NA_real_,sigma=NA_real_)
    fout <- data.frame(K=0.0,L0=0.0,r=0.0,sigma=0.0)
  }else{
  fout <- data.frame(t(fit@coef))
  }
 return(fout)
}
system.time(fn_test(data.frame(z=5)))

bo <- data.table(id = 1:10, 
  z=rnorm(10, mean=6))

plan(multisession)
system.time(j1 <- bo %>% 
        split(.$id) %>%
        future_map(~fn_test(.x),
          .progress = TRUE, 
          .options = furrr_options(packages=c(
            "rethinking"), 
              seed=1)) %>% 
        future_map_dfr(~ as_tibble(.), .id='id')
)
j1

for(i in 1:nrow(bo)){
  print(fn_test(bo[i,]$z))
}


plan(sequential)
j1 <- bo %>% 
        split(.$id) %>%
        future_map(~mean(.x$z),
          .progress = TRUE, 
          .options = furrr_options(packages=c(
            "rethinking"), 
              seed=T,
            lazy=F)) %>% 
        future_map_dfr(~ as_tibble(.), .id='id')
j1






plan(multisession)
j1 <- bo %>% 
      split(.$id) %>%
      future_map(~fn_test(.x),
        .progress = TRUE, 
        .options = furrr_options(packages=c(
          "rethinking"), 
            seed=1)) %>% 
        future_map_dfr(~ as_tibble(.), .id='id')
j1



## Data
x <- rnorm(100)
y <- 2 * x + 0.2 + rnorm(100)
w <- 1 + x ^ 2
fitA %<-% lm(y ~ x, weights = w)      ## with offset
fitB %<-% lm(y ~ x - 1, weights = w)  ## without offset
fitC %<-% {
  w <- 1 + abs(x)
  lm(y ~ x, weights = w)
}
print(fitA)
print(fitB)
print(fitC)



library(bbmle)

fn_mle2 <- function(din){
  start_day <- din[post_days <= 366][slai_3mo == min(slai_3mo,na.rm=T)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+183))]
  min_slai_anom <- din[post_days<=366][dplyr::near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_anom_3mo
  min_slai <- din[post_days<=366][dplyr::near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_3mo
  malai <- din$malai[1]
  lai_yr_sd <- din$lai_yr_sd[1]
  offset <- abs(min_slai_anom + malai)
  date_min_slai_anom <- din[post_days<=366][slai_anom_3mo==min(slai_anom_3mo,na.rm=T)]$date
  day_ttr_offset <- as.numeric(date_min_slai_anom - first(din$date_fire1))
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  upper_K <- din$malai[1]+0.25*din$lai_yr_sd[1]
  lower_K <- din$malai[1]
  lower_K <- ifelse(upper_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  lower_K <- ifelse(lower_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  upper_L0 <- min_slai
  lower_L0 <- min_slai
  b_coefs <- stats::coef(RcppArmadillo::fastLm(X=cbind(1,din$post_days),y=din$slai_3mo))
  start_r <- b_coefs[2]/(nrow(din)*30.4)

  LL <- function(K,L0,r,sigma){
    -sum(stats::dnorm(slai_3mo, mean = (K/(1 + ((K-L0)/L0)*exp(-r*post_days))), sd=sigma,log=T))}
  fit <- mle2(LL,
    parameters=list(K~1,L0~1,r~0.1,sigma~1),
    data=din,
    method="L-BFGS-B",
    start=list(K=malai,
               L0=lower_L0,
               r=0.01,
               sigma=0.5),
    lower=c(K=lower_L0,
               L0=0.01,
               r=0.00001,
               sigma=0.1), 
        upper=c(K=10,
               L0=10,
               r=0.3,
               sigma=5))
  out <- data.frame(t(coef(fit)))
  return(out)
}

fn_mle2(mdat[id==1530])

plan(multisession)
system.time(mdat[id%in%junk[1:100]] %>% 
  split(.$id) %>% 
  future_map(~fn_mle2(.)))


plan(sequential)
system.time(
  out <- mdat[id%in%junk[1:100]] %>% 
  split(.$id) %>% 
  future_map(~fn_mle2(.), 
    .progress = T, 
    .options = furrr_options(seed=333)) 
  )

out <- out %>% 
        future_map_dfr(~ as_tibble(.), .id='id') %>% 
  mutate(id=as.numeric(id))

out %>% 
  filter(L0<K) %>% 
  filter(r<0.055) %>% 
  ggplot(data=.,aes(L0/K,r))+
  geom_point()



tibble(x = 1:10, 
       id = 1:10) %>% 
  split(.$id) %>% 
  future_map(~fn_mle2(.))



x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)

## in general it is best practice to use the `data' argument,
##  but variables can also be drawn from the global environment
LL <- function(ymax=15, xhalf=6)
    -sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE))
## uses default parameters of LL
(fit <- mle2(LL))
fit1F <- mle2(LL, fixed=list(xhalf=6))
coef(fit1F)
coef(fit1F,exclude.fixed=TRUE)
