


# Function to make 4PL formula with beta distribution for Fv/Fm data
# Simply adds response_var as a user option rather than hardcoding n_surv | n_trials
make_4pl_formula <- function(random_effects = NULL, 
                                  lower = 0, 
                                  upper = 1, 
                                  family = brms::Beta(link = "logit"),
                                  response_var = "fvfm_beta") {
  
  b <- compute_4pl_bounds(lower, upper)
  
  low_expr <- sprintf("(%.6f + inv_logit(lowraw) * %.6f)", b$low_min, b$low_w)
  up_expr <- sprintf("(%.6f + inv_logit(upraw)  * %.6f)", b$up_min, b$up_w)
  
  main_rhs <- sprintf("%s + (%s - %s) / (1 + exp(exp(logk) * (logd - mid)))", 
                      low_expr, up_expr, low_expr)
  
  mid_terms <- c("temp_c", tdt_format_random_effects(random_effects))
  mid_rhs <- paste(mid_terms, collapse = " + ")
  
  # Check if family is beta
  family_name <- if (inherits(family, "brmsfamily")) family$family else family$family
  
  if (identical(family_name, "beta")) {
    # Beta regression: no trials() needed
    response_formula <- sprintf("%s ~ logit(%s)", response_var, main_rhs)
  } else {
    # Binomial/beta-binomial: use trials()
    response_formula <- sprintf("n_surv | trials(n_total) ~ %s", main_rhs)
  }
  
  brms::bf(
    stats::as.formula(response_formula),
    stats::as.formula("lowraw ~ temp_c"),
    stats::as.formula("upraw  ~ temp_c"),
    stats::as.formula("logk   ~ temp_c"),
    stats::as.formula(paste("mid    ~", mid_rhs)),
    stats::as.formula("phi ~ 1"), 
    nl = TRUE,
    family = family
  )
}

# Then 


## Reloading functions that didn't appear when loading bayesTLS

# Compute 4PL pounds
compute_4pl_bounds <- function(lower = 0, upper = 1,
                               pad = 0.001, gap = 0.002) {
  if (upper <= lower)
    stop("upper must be strictly greater than lower.", call. = FALSE)
  if (2 * pad + gap >= (upper - lower))
    stop("pad and gap leave no room for asymptote intervals; ",
         "reduce pad/gap or widen lower/upper.", call. = FALSE)
  
  midpoint <- (lower + upper) / 2
  low_min  <- lower + pad
  low_max  <- midpoint - gap / 2
  up_min   <- midpoint + gap / 2
  up_max   <- upper - pad
  
  list(low_min  = low_min,
       low_max  = low_max,
       low_w    = low_max - low_min,
       up_min   = up_min,
       up_max   = up_max,
       up_w     = up_max - up_min,
       midpoint = midpoint)
}

# Format random effect specification
tdt_format_random_effects <- function(random_effects = NULL) {
  if (is.null(random_effects) || length(random_effects) == 0) return(character())
  out <- vapply(random_effects, function(term) {
    term <- trimws(term)
    if (grepl("^\\(", term)) term else paste0("(1 | ", term, ")")
  }, character(1))
  unname(out)
}
