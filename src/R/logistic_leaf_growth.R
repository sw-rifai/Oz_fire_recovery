library(tidyverse); 
library(brms)
library(deSolve)
options(mc.cores = parallel::detectCores())

# Logistic growth model -------------------------------------------------------
# dN/dt = r*N(1-N/K)
# dL/dt = lue*L(1-L/lai_cc)
# lai_cc ~ f(precip, vpd, pet)
# lue ~ f(precip,vpd,pet)
## How? R/GPP ~ Ta

time <- seq(from=0, to=10, by = 0.01)
pars1 <- c(lue = 1, lai_cc = 5)
pars2 <- c(lue = 1, lai_cc = 5)
pars3 <- c(lue = 1, lai_cc = 5)
state <- c(L = 0.1)
fn <- function(t, state, parameters){
  with(
    as.list(c(state, parameters)),{
      dL <- lue*L*(1-L/lai_cc)
      return(list(dL))
    }
  )
}

out1 <- ode(y = c(L=0.1), times = time, func = fn, parms = pars1)
out2 <- ode(y = c(L=0.5), times = time, func = fn, parms = pars2)
out3 <- ode(y = c(L=1), times = time, func = fn, parms = pars3)
plot(L~time,data=out1,col='red')
points(L~time,out2,col='black')
points(L~time,out3,col='blue')

curve(5 / 
        (1+((5-0.5)/0.5)*exp(-1*x)), 0,10,add=T,col='purple',lwd=3)

# Verhulst model ----------------------------------------------------------
# dN/dt = r*N - alpha*N**2 
# N* = r/alpha # N* is carrying capacity


time <- seq(from=0, to=10, by = 0.01)
pars1 <- c(lue = 0.75, lai_cc = 5, alpha = 0.75/5)

verhulst <- function(t, state, parameters){
  with(
    as.list(c(state, parameters)),{
      dL <- lue*L - alpha*L**2
      return(list(dL))
    }
  )
}

out1 <- ode(y = state, times = time, func = verhulst, parms = pars1)
out2 <- ode(y = state, times = time, func = fn, parms = pars2)
out3 <- ode(y = state, times = time, func = fn, parms = pars3)
plot(L~time,data=out1,col='red')
points(L~time,out2,col='black')
points(L~time,out3,col='blue')



curve(0.5 + 5*tanh(x/2), 0,20,col='blue'); abline(h=5+0.5,col='blue')
curve(0.25+4.5*tanh(x/3),add=T,col='black'); abline(h=4.5+0.25)
curve(0+4*tanh(x/4), add=T,col='red'); abline(h=4,col='red')






# Logistic growth model -------------------------------------------------------
# dN/dt = r*N(1-N/K)
# dL/dt = lue*L(1-L/lai_cc)
# lai_cc ~ f(precip, vpd, pet)
# lue ~ f(precip,vpd,pet)

time <- seq(from=0, to=10, by = 0.01)
state <- c(L = 0.1)
pars1 <- c(lue = 1, lai_cc = 5)
fn <- function(t, state, parameters){
  with(
    as.list(c(state, parameters)),{
      dL <- lue*L*(1-L/lai_cc)
      return(list(dL))
    }
  )
}

out1 <- ode(y = c(L=0.1), times = time, func = fn, parms = pars1)
out1 <- as.matrix(out1) 
out1 <- as.data.frame(out1)
out1 <- out1 %>% 
  as.data.frame() %>% 
  rowwise %>%
  mutate(Lobs = L+abs(rnorm(1, mean=0, sd=L**0.1)))
out1 %>% ggplot(data=.,aes(time,Lobs))+
  geom_point()+
  geom_smooth()
plot(L~time,data=out1,col='red')


library(brms)

f <- bf(Lobs ~ K/(1 + ((K-L0)/L0)*exp(-r*time)), 
   K ~ 1, 
   L0 ~ 1, 
   r ~ 1, 
   nl=TRUE)
bprior <- prior(nlpar=K, normal(5,2),lb=1)+
  prior(nlpar=L0, normal(0,1),lb=0.001)+
  prior(nlpar=r, normal(1,1),lb=0.5)
make_stancode(f, 
              data=out1, 
              family=gaussian(), 
              prior = bprior)

fit <- brm(f,
            data = out1, 
            prior = bprior, 
            # sample_prior = 'only')
            algorithm = 'sampling',
           # backend="cmdstanr", 
           chains=3,
           iter = 1000
            # control = list(adapt_delta=0.999)
)

plot(fit)
prior_summary(fit)


# Simulate model with multiple groups ------------------------------------------
# time <- seq(from=0,to=10,by=0.1)
time <- seq(from=0, to=10, by = 0.01)
fn <- function(t, state, parameters){
  with(
    as.list(c(state, parameters)),{
      dL <- r*L*(1-L/K)
      return(list(dL))
    }
  )
}
out2 <- tibble(K=rnorm(10,mean=5),
       r = rnorm(10,mean=1,sd=0.1),
       id=1:10)
fn_ode <- function(K,r){
  out <- ode(y=c(L=0.1), time=time, func=fn, parms=c("K"=K,"r"=r)) %>% 
    as.data.frame()
  return(out)
}
fn_ode(5,1) %>% str
sim2 <- out2 %>% #split(.$id) %>%
  pmap_dfr( ~fn_ode(.,.),.id = 'id') #%>%
  # map_dfr(~as_tibble(.), .id='id') 

sim2 %>% ggplot(data=.,aes(time,L,group=id,color=id))+
  geom_line()


tmp <- map2(out2$K, out2$r, ~fn_ode(.x,.y))
bind_rows(tmp,.id='id') %>% 
  ggplot(data=.,aes(time,L,group=id,color=id))+
  geom_line()

sim2 <- bind_rows(tmp,.id='id') %>% 
  as.data.frame() %>% 
  rowwise %>%
  mutate(Lobs = L+abs(rnorm(1, mean=0, sd=L**0.05)))
sim2 %>% ggplot(data=.,aes(time,Lobs,color=id))+
  geom_point()



f2 <- bf(Lobs ~ K/(1 + ((K-L0)/L0)*exp(-r*time)), 
        K ~ id, 
        L0 ~ id, 
        r ~ id, 
        nl=TRUE)
bprior2 <- prior(nlpar=K, normal(5,2),lb=1)+
  prior(nlpar=L0, normal(0,1),lb=0.001)+
  prior(nlpar=r, normal(1,1),lb=0.5)
make_stancode(f2, 
              data=sim2, 
              family=gaussian(), 
              prior = bprior2)

fit <- brm(f2,
           data = sim2, 
           prior = bprior2, 
           # sample_prior = 'only')
           algorithm = 'sampling',
           # backend="cmdstanr", 
           chains=3,
           iter = 1000
           # control = list(adapt_delta=0.999)
)
fit
plot(fit)
prior_summary(fit)
pp_check(fit,nsamples = 100)



# tmp1 <- tmp[[1]]
# tmp2 <- tmp[[2]]
# bind_rows(tmp[[1]],tmp[[2]],.id = 'id')$id %>% table
# tmp[[1]] %>% class
# tmp[[2]] %>% class
# 
# out2[1,] %>% as.list()
#   pmap( ~fn_ode(K,r))
# 
# 
# 
# 
# map2_df(out2$K, out2$r, fn_ode)
# 
# tibble(K=rnorm(10,mean=5),
#        r = rnorm(10,mean=0,sd=0.1),
#        id=1:10) %>% split(.$id) %>% 
#   map(~fn_ode(.x))
# 
# 
# 
# tibble(K=rnorm(10,mean=5),
#        r = rnorm(10,mean=0,sd=0.1),
#        id=1:10) %>% 
#   split(.$id) %>% 
#   pmap(., fn_ode(K,r))
# 
# # system.time(out <- mdat %>% 
# #               split(.$id) %>%
# #               future_map(~fn_w4(.x),.progress = TRUE) %>% 
# #               future_map_dfr(~ as_tibble(.), .id='id')
# # )
# 
# out2 <- expand_grid(tibble(time=time), 
#             tibble(K=rnorm(10,mean=5),
#                    r = rnorm(10,mean=0,sd=0.1),
#                    id=1:10))
# dim(out2)
# 
# 
# out2 %>% 
#   rowwise() %>% 
#   mutate(y = ode(y=c(L=0.1), 
#                  times=time, 
#                  func=fn, 
#                  parms=c(.K,.r))[,"L"])
# 
# out2 %>% 
#   split(.$id) %>% 
#   map(~ ode(y=c(L=0.1), 
#                 times=time, 
#                 func=fn, 
#                 parms=c(K,r))[,"L"])
# ode(y = c(L=0.1), times = time, func = fn, parms = pars1)
# 
# 
# 
# 
# 
# x <- list(1, 1, 1)
# y <- list(10, 20, 30)
# z <- list(100, 200, 300)
# l <- list(x,y,z)
# 
# pmap(l, function(c, b, a) (a + c) * b)
