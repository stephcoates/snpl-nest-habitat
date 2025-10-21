# GLM covariate term testing ----
# this function tries each term in a list of variables as linear and quadratic, returns a list sorted by AIC value
run_glmer_models <- function(variables, data) {
  results <- data.frame(Variable = character(), Model = character(), AIC = numeric(), stringsAsFactors = FALSE)
  
  for (var in variables) {
    if (!var %in% names(data)) {
      warning(paste("Variable", var, "not found in the dataset. Skipping..."))
      next
    }
    
    # Linear model
    mod1 <- glmmTMB(as.formula(paste("Outcome ~ ", var)), data = data, family = binomial)
    aic1 <- AIC(mod1)
    
    # Quadratic model
    mod2 <- glmmTMB(as.formula(paste("Outcome ~ ", var, "+ I(", var, "^2)")), data = data, family = binomial)
    aic2 <- AIC(mod2)
    
    # Store results
    results <- rbind(results, 
                     data.frame(Variable = var, Model = "Linear", AIC = aic1),
                     data.frame(Variable = var, Model = "Quadratic", AIC = aic2))
  }
  
  # Sort results by AIC (ascending)
  results <- results[order(results$AIC), ]
  
  return(results)
}

# Function to check collinearity models in a list ----
check_all_collinearity <- function(model_list) {
  results <- lapply(model_list, function(m) {
    check_collinearity(m)
  })
  names(results) <- names(model_list)
  return(results)
}

# Function to check goodness-of-fit diagnostics a model ----
run_dharma_diagnostics <- function(model) {
  # Simulate residuals
  res <- simulateResiduals(fittedModel = model)
  
  # Plot residual diagnostics
  plot(res)
  
  # Run a KS test on the simulated residuals
  uniformity_test <- testUniformity(res)
  
  # Return the test results
  return(uniformity_test)
}

# Function to check goodness-of-fit for a list of models ----
# returns summary table w/ p values
run_all_dharma_diagnostics <- function(model_list, show_plots = FALSE) {
  
  # Helper function for a single model
  run_dharma_diagnostics <- function(model, show_plot = FALSE) {
    res <- simulateResiduals(fittedModel = model)
    
    if (show_plot) {
      plot(res)
    }
    
    uniformity_test <- testUniformity(res)
    return(list(res = res, test = uniformity_test))
  }
  
  # Run diagnostics for each model
  results <- map(model_list, ~ run_dharma_diagnostics(.x, show_plots))
  
  # Create a summary data frame of test results
  summary_df <- map_dfr(
    names(results),
    ~ data.frame(
      model = .x,
      statistic = results[[.x]]$test$statistic,
      p_value = results[[.x]]$test$p.value
    )
  )
  
  # Return both the raw results and the summary table
  return(list(
    diagnostics = results,
    summary = summary_df
  ))
}

# Create candidate model table from list of models ----
create_candidate_model_table <- function(models_list) {
  model_info <- lapply(names(models_list), function(name) {
    model <- models_list[[name]]
    
    # Extract key metrics
    model_aic <- AIC(model)
    model_dev <- deviance(model)
    model_df  <- df.residual(model)
    
    # Correctly count parameters (k)
    if ("glmmTMB" %in% class(model)) {
      k <- length(fixef(model)$cond)  # fixed effects only
    } else {
      k <- length(coef(model))        # works for glm and similar
    }
    
    # Convert formula to readable string
    formula_text <- paste(deparse(formula(model)), collapse = " ")
    
    data.frame(
      Model = name,
      Formula = formula_text,
      K = k,
      AIC = model_aic,
      DeltaAIC = NA,  # placeholder for later calculation
      Deviance = model_dev,
      df = model_df,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine and calculate ΔAIC
  model_table <- bind_rows(model_info) %>%
    mutate(
      DeltaAIC = AIC - min(AIC, na.rm = TRUE)
    ) %>%
    select(Model, Formula, K, AIC, DeltaAIC, Deviance, df) %>%
    arrange(AIC)
  
  return(model_table)
}
