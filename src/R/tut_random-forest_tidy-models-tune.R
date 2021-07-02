library(tidymodels)
data(ames)
ames <- mutate(ames, Sale_Price = log10(Sale_Price))

set.seed(123)
ames_split <- initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train <- training(ames_split)
ames_test  <-  testing(ames_split)

# ames_rec <- 
#   recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
#            Latitude + Longitude, data = ames_train) %>%
#   step_log(Gr_Liv_Area, base = 10) %>% 
#   step_other(Neighborhood, threshold = 0.01) %>% 
#   step_dummy(all_nominal_predictors()) %>% 
#   step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") ) %>% 
#   step_ns(Latitude, Longitude, deg_free = 20)
# 
# lm_model <- linear_reg() %>% set_engine("lm")
# 
# lm_wflow <- 
#   workflow() %>% 
#   add_model(lm_model) %>% 
#   add_recipe(ames_rec)
# 
# lm_fit <- fit(lm_wflow, ames_train)
# 
# ames_test_res <- predict(lm_fit, new_data = ames_test %>% select(-Sale_Price))
# ames_test_res <- bind_cols(ames_test_res, ames_test %>% select(Sale_Price))
# ames_test_res
# 
# ggplot(ames_test_res, aes(x = Sale_Price, y = .pred)) + 
#   # Create a diagonal line:
#   geom_abline(lty = 2) + 
#   geom_point(alpha = 0.5) + 
#   labs(y = "Predicted Sale Price (log10)", x = "Sale Price (log10)") +
#   # Scale and size the x- and y-axis uniformly:
#   coord_obs_pred()
# 
# ames_metrics <- metric_set(rmse, rsq, mae)
# ames_metrics(ames_test_res, truth = Sale_Price, estimate = .pred)

# Random forest
rf_model <- 
  rand_forest(trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("regression")

rf_wflow <- 
  workflow() %>% 
  add_formula(
    Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
      Latitude + Longitude) %>% 
  add_model(rf_model) 

set.seed(55)
ames_folds <- vfold_cv(ames_train, v = 10)

keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)


# ??? _------------------------
rf_res <- rf_wflow %>% fit_resamples(resamples = ames_folds, control = keep_pred)


# hyperparameter tuning -----------------------------
rf_spec <- # This is the model object
  rand_forest(mtry = tune(),
              trees = tune(), 
              min_n = tune(),
              mode = 'regression') %>% 
  set_engine("ranger", 
             regularization.factor = tune("regularization"))

rf_param <- parameters(rf_spec)
rf_param

pca_rec <- 
  recipe(Sale_Price ~ ., data = ames_train) %>% 
  # Select the square-footage predictors and extract their PCA components:
  step_normalize(contains("SF")) %>% 
  # Select the number of components needed to capture 95% of
  # the variance in the predictors. 
  step_pca(contains("SF"), threshold = .95)

updated_param <- 
  workflow() %>% 
  add_model(rf_spec) %>% 
  add_recipe(pca_rec) %>% 
  parameters() %>% 
  finalize(ames_train) # key part to specify upper boundary for mtry
updated_param

updated_param %>% pull_dials_object("mtry") # verify it set the upper bound

# set up the grid of tuning parameters
rf_grid <- grid_latin_hypercube(finalize(mtry(), ames_train),
                                trees(),
                                min_n(),
                                regularization = regularization_factor(),
                                size=10)
rf_grid

# specify the 'workflow'
tune_wf <- workflow() %>% 
  add_recipe(pca_rec) %>% 
  add_model(rf_spec) 

# set up parallel processing
doParallel::registerDoParallel()
set.seed(345)
# actual grid tuning step
tune_res <- tune_grid(
  tune_wf,
  resamples = ames_folds,
  grid = rf_grid
)


# Plot the model fits to the tuning parameters
tune_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "R2")


# set up the 2nd grid of tuning parameters
rf_grid2 <- grid_regular(finalize(mtry(), ames_train),
                                trees(range = c(800,1000)),
                                min_n(),
                                regularization = regularization_factor(),
                                levels=c(trees=1, 
                                         mtry=3,
                                         min_n=3, 
                                         regularization=3))
rf_grid2

# Fit the revised tuning grid
tune_res2 <- tune_grid(
  tune_wf,
  resamples = ames_folds,
  grid = rf_grid2
)


# Plot the GOF of the 2nd gen of tuning model fits
tune_res2 %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

tune_res2 %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "R2")



# set up the 2nd grid of tuning parameters
rf_grid3 <- grid_regular(finalize(mtry(), ames_train),
                         trees(range = c(800,800)),
                         min_n(range= c(21,21)),
                         regularization = regularization_factor(range = c(0.25,1)),
                         levels=c(trees=1, 
                                  mtry=5,
                                  min_n=1, 
                                  regularization=5))
rf_grid3

# Fit the revised tuning grid
tune_res3 <- tune_grid(
  tune_wf,
  resamples = ames_folds,
  grid = rf_grid3
)


# Plot the GOF of the 2nd gen of tuning model fits
tune_res3 %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

tune_res3 %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "R2")


best_rmse <- select_best(tune_res3, metric = 'rmse')
best_rmse

final_rf <- finalize_model(rf_spec, parameters = best_rmse)
final_rf


# Examine variable importance ------------------
library(vip)
final_rf %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Sale_Price ~ .,
      data = ames_train
  ) %>%
  vip(geom = "point")


# Fit Final Model -------------------------
final_wf <- workflow() %>%
  add_model(final_rf) %>% 
  add_recipe(pca_rec)

final_res <- final_wf %>%
  last_fit(ames_split)

final_res %>%
  collect_metrics()


final_res %>% 
  collect_predictions() %>% 
  ggplot(data=.,aes(.pred, Sale_Price))+
  geom_point()+
  geom_abline(col='red')
