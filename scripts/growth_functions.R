#' Corrects for background 
#' 
#' This function corrects a vector of cell counts or absorbance readings for the minimum value. 
#' Typically applied to grouped or nested datasets prior to fitting the logistic curve.
#' @inheritParams fit_logistic_curve
correct_for_min_background <- function(data_n) {
  if (all(is.na(data_n))) return(NA_real_)
  return(data_n - min(data_n, na.rm = TRUE))
}

#' Corrects for background 
#' 
#' Correct for background by subtracting a constant.
#' @inheritParams fit_logistic_curve
correct_for_fixed_background <- function(data_n, bgrd = 0) {
  return(data_n - bgrd)
}

#' Corrects for background 
#' 
#' Correct for background using a condition
#' @inheritParams fit_logistic_curve
correct_for_fixed_background <- function(data_n, cond) {
  stopifnot(!missing(data_n))
  stopifnot(!missing(cond))
  stopifnot(is.logical(cond))
  stopifnot(length(cond) == 1 || length(cond) == length(data_n))
  bgrd <- data_n[cond]
  if (all(is.na(bgrd))) return(data_n)
  return(data_n - mean(bgrd, na.rm = TRUE))
}

#' Mark data points after the population maximum
#' 
#' Function operates on a data frame holding growth curve(s) and marks data points that occur AFTER the peak in population as belonging to the death phase. Either pre-group or use the \code{group_by} parameter to specify what constitutes a single growth curve.
#' @inheritParams fit_logistic_curve
#' @return data frame with new column \code{death_phase}
mark_death_phase <- function(df, time, N, group_by = NULL, quiet = FALSE) {
  
  # safety checks
  if (missing(time)) stop("specify a 'time' column", call. = FALSE)
  if (missing(N)) stop("specify a population size 'N' column", call. = FALSE)
  
  # get expressions
  time_col <- rlang::enexpr(time)
  N_col <- rlang::enexpr(N)
  
  # grouping
  group_by_expr <- rlang::enexpr(group_by)
  if (!rlang::is_null(group_by_expr)) {
    if (rlang::is_call(group_by_expr) && rlang::call_name(group_by_expr) == "c") {
      group_by_cols <- rlang::call_args(group_by_expr)
    } else {
      group_by_cols <- list(group_by_expr)
    }
    df <- df %>% dplyr::group_by(!!!group_by_cols)
  }

  # info
  if (!quiet) {
    grouped_by <- ""
    if(length(grps <- dplyr::group_vars(df)) > 0) 
      grouped_by <- sprintf(" (grouped by '%s')", paste(grps, collapse = "', '"))
    sprintf("Info: marking death phase for %d growth curves%s... ", 
            dplyr::n_groups(df), grouped_by) %>% 
      message()
  }
  
  # mark death phase
  df <- df %>% 
    # find time at N max
    mutate(
      ..peak_time.. = 
        if (all(is.na(!!N_col))) NA_real_
        else (!!time_col)[which(!!N_col == max(!!N_col, na.rm = TRUE))[1]],
      death_phase = ifelse(is.na(..peak_time..), FALSE, !!time_col > ..peak_time..)
    ) %>% 
    # cleanup
    select(-..peak_time..)
  
  # ungroup if group_by was provided
  if (!rlang::is_null(group_by_expr)) {
    df <- df %>% ungroup()
  }
  
  return(df)
} 

#' Calculates growth curve parameters.
#' 
#' Uses fit_logistic_curve and extract_logistic_fit_parameters to simplify fitting logistic curves in a safe manner across multiple experiments in a data frame (pre-group before calling this function or provide the \code{group_by} parameter).
#' 
#' @param df pre-grouped data frame
#' @param time the time column
#' @param N the population size column (background corrected if appropriate)
#' @param group_by what identifies an individual growth curve? use c(x, y, z) to group by multiple colums
#' @param extract_parameters whether to extract the fit parameters automatically from the fit (if TRUE, use \code{coefficient_select} and \code{summary_select} to specify extraction)
#' @inheritParams extract_logistic_fit_parameters
#' @param keep_fit whether to keep the nls fit as a column (NO by default if parameters are extracted, otherwise YES)
#' @param quiet whether to provide summary information on the fits
#' @return data frame with summarized growth curve parameters
estimate_growth_curve_parameters <- function(
  df, time, N, group_by = NULL,
  extract_parameters = TRUE, 
  summary_select = c(n_used = nobs, rmsd = deviance),
  coefficient_select = c(N0 = estimate_N0, N0_se = std.error_N0,
                         K = estimate_K, K_se = std.error_K, 
                         r = estimate_r, r_se = std.error_r),
  keep_fit = !extract_parameters, quiet = FALSE) {
  
  # safety checks
  if (missing(time)) stop("specify a 'time' column", call. = FALSE)
  if (missing(N)) stop("specify a population size 'N' column", call. = FALSE)
  
  # get expressions
  time_col <- rlang::enexpr(time)
  N_col <- rlang::enexpr(N)
  
  # data frame
  df <- df %>% filter(!is.na(!!time_col), !is.na(!!N_col))
  
  # grouping
  group_by_expr <- rlang::enexpr(group_by)
  if (!rlang::is_null(group_by_expr)) {
    if (rlang::is_call(group_by_expr) && rlang::call_name(group_by_expr) == "c") {
      group_by_cols <- rlang::call_args(group_by_expr)
    } else {
      group_by_cols <- list(group_by_expr)
    }
    df <- df %>% dplyr::group_by(!!!group_by_cols)
  }
  
  # info
  if (!quiet) {
    grouped_by <- ""
    if(length(grps <- dplyr::group_vars(df)) > 0) 
      grouped_by <- sprintf(" (grouped by '%s')", paste(grps, collapse = "', '"))
    sprintf("Info: estimating growth parameters for %d growth curves%s... ", 
            dplyr::n_groups(df), grouped_by) %>% 
      message(appendLF = FALSE)
  }
  
  # fit logistic curve
  safely_fit_logistic_curve <- purrr::safely(fit_logistic_curve)
  df_fits <-
    df %>% 
    summarize(
      time_min = min(!!time_col), # data range min
      time_max = max(!!time_col), # data range max
      n_datapoints = length(!!N_col), # number of total data points
      safe_fit = list(safely_fit_logistic_curve(!!time_col, !!N_col)),
      fit = map(safe_fit, ~.x$result),
      error = map_chr(safe_fit, ~{
        if (!is.null(.x$error)) .x$error$message
        else NA_character_
      }),
      .groups = "drop"
    ) %>% 
    select(-safe_fit)
  
  # info
  if (!quiet) {
    sprintf("%d curves fitted successfully.", sum(is.na(df_fits$error))) %>% 
      message()
  }
  
  # extract fit parameters
  if (extract_parameters) 
    df_fits <- df_fits %>% 
      mutate(params = map(
        fit, extract_logistic_fit_parameters, 
        summary_select = !!enquo(summary_select),
        coefficient_select = !!enquo(coefficient_select)
      )) %>% 
      unnest(cols = c(params), keep_empty = TRUE)
  
  # keep fit
  if (!keep_fit) 
    df_fits <- df_fits %>% select(-fit)
  
  # return result
  return(df_fits)
}

#' Fits a logistic curve to data.
#'
#' This function fits a logistic curve to the supplied data, namely
#' N(t) = K / (1 + ( (K - N0) / N0) * exp(-r * t), where
#' N(t) is the number of cells (or density) at time t,
#' K is the carrying capacity,
#' N0 is the initial cell count or density, and
#' r is the "growth rate".
#' @param data_t    A vector of timepoints (data_n must also
#'                  be provided and be the same length).
#' @param data_n    A vector of cell counts or absorbance readings.
#' @return          An object of class nls.
fit_logistic_curve <- function(data_t, data_n) {
  
  # safety checks
  if (!is.numeric(data_t) |!is.numeric(data_n))
    stop("the input data (data_t and data_n) must be numeric vectors", call. = FALSE)
  if (length(data_t) != length(data_n))
    stop("the input data (data_t and data_n) must have the same length", call. = FALSE)
  
  # create dataset
  df <- tibble(time = data_t, N = data_n) %>% 
    # remove entries that will not be useful for fitting (exclude all negative and 0 N)
    filter(!is.na(N), !is.na(time), N > 0)
  
  # check for enough data
  if (!nrow(df) >= 3)
    stop("need at least 3 positive data points, background correct if appropriate.", call. = FALSE)
  
  # estimate initial parameters
  K_init <- 2 * max(df$N) # carrying capacity is near the max but possibly above
  N0_init <- min(df$N) # initial population size is near the min
  
  # make an initial estimate for growth rate r
  glm_fit <- stats::glm(
    N / K_init ~ time,
    family = stats::quasibinomial("logit"),
    data = df
  )
  r_init <- stats::coef(glm_fit)[[2]]
  
  # make sure slope estimate is positive
  if (is.na(r_init) || r_init <= 0) r_init <- 0.001
  
  # nls fit
  fit <- 
    tryCatch(
      stats::nls(
        formula = N ~ K/(1 + ((K - N0)/N0) * exp(-r * time)), 
        data = df, 
        start = list(K = K_init, N0 = N0_init, r = r_init),
        control = list(maxiter = 500),
        algorithm = "port",
        lower = c(stats::median(df$N), 0, 0),
        upper = c(Inf, max(df$N), Inf),
      ),
      error = function(e) {
        stop("could not fit logistic equation")
      }
    )
  
  # checks whether results make sense
  if(coefficients(fit)[['K']] < coefficients(fit)[['N0']]) {
    stop("carrying capacity K smaller than N0")
  }
  
  # return value
  return(fit)
}

#' Extract parameter estimates and fit statistics from the logistic fit
#' @param nls_fit the nls fit object
#' @param coefficient_select a dplyr::select statement for the coefficient columns
#' @param summary_select a dplyr::select statement for the fit summary columns
#' @return a tibble with the fit parameters
extract_logistic_fit_parameters <- function(nls_fit, summary_select = everything(), coefficient_select = everything()) {
  
  # safety check
  if (is.null(nls_fit)) return(tibble())
  
  # extract coefficients
  coefs <- 
    broom::tidy(nls_fit) %>% 
    select(term, estimate, std.error) %>% 
    tidyr::pivot_wider(names_from = term, values_from = c(estimate, std.error)) %>% 
    select(!!enquo(coefficient_select))
  
  # extract summary
  sums <-
    broom::glance(nls_fit) %>% 
    select(!!enquo(summary_select))
  
  # safety checks
  stopifnot(nrow(coefs) == 1L)
  stopifnot(nrow(sums) == 1L)
  
  # return
  return(dplyr::bind_cols(sums, coefs))
}


#' Generate logistic curve
#' 
#' @param time colum name for time column (newly generated)
#' @param N colum name for population size column (newly generated)
#' @param time_min column name for start time of curve
#' @param time_max column name for end time of curve
#' @param N0 column name for initial populatio size (N0) parameter
#' @param K column name for carrying capacity (K) parameter
#' @param r column name for growth rate (r) parameter
generate_logistic_curve <- function(df, time, N, time_min = time_min, time_max = time_max, N0 = N0, K = K, r = r) {
  
  # safety checks
  if (missing(time)) stop("specify a 'time' column", call. = FALSE)
  if (missing(N)) stop("specify a population size 'N' column", call. = FALSE)
  
  time_col <- rlang::enexpr(time)
  N_col <- rlang::enexpr(N)
  time_min_col <- rlang::enexpr(time_min)
  time_max_col <- rlang::enexpr(time_max)
  N0_col <- rlang::enexpr(N0)
  K_col <- rlang::enexpr(K)
  r_col <- rlang::enexpr(r)
  
  df %>%
    mutate(!!time_col := purrr::map2(!!time_min_col, !!time_max_col, seq, length.out = 100)) %>%
    unnest(!!time_col) %>%
    mutate(!!N_col := !!K_col/(1 + ((!!K_col - !!N0_col)/!!N0_col) * exp(-!!r_col * !!time_col))) %>% 
    filter(!is.na(!!N_col))
}
