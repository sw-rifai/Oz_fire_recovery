library(tidyverse); 
library(nls.multstart)

# Data generating function
fn_lg <- function(time, K,L0,r){
  K/(1+((K-L0)/L0)*exp(-r*time))
}

# solving for time at which L=0.5*K (midpoint): -log(1.0*L0/(K - L0))/r

# simulate data
sim <- expand_grid(time = seq(0,500, length.out=100), 
            K= 5, 
            L0=0.1, 
            r=0.0175) %>% 
  mutate(L = fn_lg(time,K,L0,r))

# plot simulated data
sim %>% 
  ggplot(data=., aes(time, L))+
  geom_line()+
  geom_vline( # draw midpoint
    aes(xintercept = -log(1.0*L0/(K - L0))/r), 
    lty=3
  )+
  geom_hline(aes(yintercept=0.5*K),lty=3)

fit_lf <- nls_multstart(L ~ fn_lg(time,K,L0,r),
                 data=sim %>% select(L,time), 
                 iter = 1000, # number of times nls will be fit
                 start_lower = c(K=0, # lower range of starting param guesses
                                 L0=0,
                                 r=0), 
                 start_upper = c(K=6, # upper range of starting param guesses
                                 L0=6, 
                                 r=0.2),
                 lower = c(K=0, # lower bound of params
                           L0=0, 
                           r=0),
                 upper = c(K=10, # upper bound params
                           L0=5,
                           r=0.2),
                 supp_errors = 'Y')
summary(fit_lf) # exact fit, too easy


# add noise to simulate observations
fn_noise <- function(x) rgamma(1, shape=10*x, scale=(1/10)) # ?rgamma: mean = shape*scale
sim <- sim %>% rowwise() %>% mutate(L_obs = fn_noise(L))
sim$L_obs %>% plot


fit_lf <- nls_multstart(L_obs ~ fn_lg(time,K,L0,r),
                        data=sim %>% select(L_obs,time), 
                        iter = 1000, # number of times nls will be fit
                        start_lower = c(K=0, # lower range of starting param guesses
                                        L0=0,
                                        r=0), 
                        start_upper = c(K=6, # upper range of starting param guesses
                                        L0=6, 
                                        r=0.2),
                        lower = c(K=0, # lower bound of params
                                  L0=0, 
                                  r=0),
                        upper = c(K=10, # upper bound params
                                  L0=5,
                                  r=0.2),
                        supp_errors = 'Y')
summary(fit_lf) # Still pretty close to true values

sim %>% 
  mutate(L_pred = fn_lg(time=time, 
                        K=coef(fit_lf)['K'], 
                        L0=coef(fit_lf)['L0'],
                        r=coef(fit_lf)['r'])) %>% 
  ggplot(data=.,aes(time, L_pred))+
  geom_point(aes(time,L_obs))+
  geom_line(col='blue')+
  geom_line(aes(time,L),col='red')+
  geom_line(col='blue')+
  geom_vline(
    aes(xintercept = -log(1.0*coef(fit_lf)['L0']/(coef(fit_lf)['K'] - coef(fit_lf)['L0']))/coef(fit_lf)['r']), 
    lty=3,col='blue'
  )+
  geom_hline(aes(yintercept=0.5*coef(fit_lf)['K']),lty=3,col='blue')





