library(data.table)
data.table::setDTthreads(1)
`%>%` <- magrittr::`%>%`
library(arrow)

# Data load and prep --------------------------
dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-05-24 16:20:29.parquet")

# landscape covars ------------
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")

# species ----------------------
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]


# pre & post_fire climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y","date",
                                                "map","mapet","mappet","mavpd15",
                                                "matmin","matmax",
                                                "precip_anom",
                                                "precip_anom_12mo",
                                                "post_precip_anom_12mo",
                                                "vpd15",
                                                "vpd15_u",
                                                "vpd15_anom",
                                                "vpd15_anom_12mo",
                                                "tmax_anom_12mo",
                                                "tmax",
                                                "tmax_u",
                                                "tmin",
                                                "tmin_u",
                                                "post_vpd15_anom_12mo", 
                                                "post_tmax_anom_12mo",
                                                "post_pet_anom_12mo")) %>% 
  as.data.table()

# next step breaks tune_grid???
post_clim <- post_clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                                       n = 3,fill = NA,align='right'), by=.(x,y)]



# data
library(AmesHousing)

# data cleaning
library(janitor)

# data prep
library(dplyr)

# tidymodels
library(rsample)
library(recipes)
library(parsnip)
library(tune)
library(dials)
library(workflows)
library(yardstick)

# speed up computation with parrallel processing (optional)
# library(doParallel)
# all_cores <- parallel::detectCores(logical = FALSE)
# registerDoParallel(cores = all_cores)


# set the random seed so we can reproduce any simulated results.
set.seed(1234)

# load the housing data and clean names
ames_data <- make_ames() %>%
  janitor::clean_names()


# split into training and testing datasets. Stratify by Sale price 
ames_split <- rsample::initial_split(
  ames_data, 
  prop = 0.8, 
  strata = sale_price
)

# preprocessing "recipe"
preprocessing_recipe <- 
  recipes::recipe(sale_price ~ ., data = training(ames_split)) %>%
  # convert categorical variables to factors
  recipes::step_string2factor(all_nominal()) %>%
  # combine low frequency factor levels
  recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information 
  recipes::step_nzv(all_nominal()) %>%
  prep()


ames_cv_folds <- 
  recipes::bake(
    preprocessing_recipe, 
    new_data = training(ames_split)
  ) %>%  
  rsample::vfold_cv(v = 5)


# XGBoost model specification
xgboost_model <- 
  parsnip::boost_tree(
    mode = "regression",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost", objective = "reg:squarederror")


# grid specification
xgboost_params <- 
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )

xgboost_grid <- 
  dials::grid_max_entropy(
    xgboost_params, 
    size = 6
  )

knitr::kable(head(xgboost_grid))


xgboost_wf <- 
  workflows::workflow() %>%
  add_model(xgboost_model) %>% 
  add_formula(sale_price ~ .)


library(doParallel)
all_cores <- parallel::detectCores(logical = FALSE)
registerDoParallel(cores = all_cores/2)

# hyperparameter tuning
xgboost_tuned <- tune::tune_grid(
  object = xgboost_wf,
  resamples = ames_cv_folds,
  grid = xgboost_grid,
  metrics = yardstick::metric_set(rmse, rsq, mae),
  control = tune::control_grid(verbose = TRUE)
)
