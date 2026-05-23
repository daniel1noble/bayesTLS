# Documentation for the case-study datasets shipped with bayesTLS.
# The datasets are built from the raw CSVs in inst/extdata/ by
# data-raw/make_datasets.R. Each is analysis-ready for the workflow
# standardize_data() -> fit_4pl() -> extract_tdt().

#' Brown shrimp lethal thermal-death-time data
#'
#' Replicate lethal-TDT trials for brown shrimp (\emph{Crangon crangon}). Each
#' row is one tank of individuals exposed to a fixed assay temperature for a
#' fixed duration; the response is the number that died. The model-ready frame
#' for Case Study 1 (lethal endpoint).
#'
#' @format A data frame with 148 rows and 6 variables:
#' \describe{
#'   \item{Date}{Experiment date (use as a grouping factor).}
#'   \item{Tank}{Holding-tank identifier (use as a grouping factor).}
#'   \item{Temperature_assay}{Assay temperature (degrees C).}
#'   \item{Duration_exposure_hours}{Exposure duration (hours).}
#'   \item{N_individuals_after_trial}{Number of individuals in the trial.}
#'   \item{Mortality_after_trial}{Number that died during the trial, out of
#'         \code{N_individuals_after_trial}.}
#' }
#' @source Brown shrimp lethal-TDT assay (Case Study 1). Raw file:
#'   \code{system.file("extdata", "data_lethal_TDT_brown_shrimp.csv", package = "bayesTLS")}.
#' @examples
#' std <- standardize_data(shrimp_lethal, temp = "Temperature_assay",
#'                         duration = "Duration_exposure_hours",
#'                         n_total = "N_individuals_after_trial",
#'                         mortality = "Mortality_after_trial",
#'                         random_effects = c("Date", "Tank"),
#'                         duration_unit = "hours")
"shrimp_lethal"

#' Brown shrimp sublethal time-to-knockdown data
#'
#' Sublethal TDT trials for brown shrimp (\emph{Crangon crangon}): each cup of
#' individuals contributes the elapsed time to loss of response to touch
#' (knockdown) at a fixed assay temperature. Cleaned from the raw clock-time
#' records (excluded rows dropped; start/stop times parsed to elapsed minutes).
#'
#' @format A data frame with 299 rows and 5 variables:
#' \describe{
#'   \item{assay_temp}{Assay temperature (degrees C).}
#'   \item{time_to_event}{Time to knockdown (minutes).}
#'   \item{date_experiment}{Experiment date (grouping factor).}
#'   \item{tank_ID}{Holding-tank identifier (grouping factor).}
#'   \item{cup_ID}{Cup identifier, \code{Trial_ID_Sample} (grouping factor).}
#' }
#' @source Brown shrimp sublethal time-to-knockdown assay (Case Study 1,
#'   sublethal endpoint). Raw file:
#'   \code{system.file("extdata", "data_sublethal_TDT_brown_shrimp.csv", package = "bayesTLS")}.
"shrimp_sublethal"

#' Zebrafish lethal thermal-death-time data across life stages
#'
#' Lethal-TDT trials for zebrafish (\emph{Danio rerio}) at three life stages.
#' Built from the raw daily survival sheet by summing the per-day
#' morning/afternoon mortality counts into one death count per trial and dropping
#' excluded rows. One row per assay trial. The model-ready frame for Case Study 2.
#'
#' @format A data frame with 323 rows and 7 variables:
#' \describe{
#'   \item{assay_temp}{Assay temperature (degrees C).}
#'   \item{duration_h}{Exposure duration (hours).}
#'   \item{n_total}{Number of individuals in the trial.}
#'   \item{n_surv}{Number that survived.}
#'   \item{n_dead}{Number that died (\code{n_total - n_surv}).}
#'   \item{life_stage}{Life stage, a factor with levels \code{young_embryos},
#'         \code{old_embryos}, \code{larvae}.}
#'   \item{Date_experiment}{Experiment date (grouping factor).}
#' }
#' @source Zebrafish lethal-TDT assay (Case Study 2). Raw file:
#'   \code{system.file("extdata", "data_lethal_TDT_zebrafish.csv", package = "bayesTLS")}.
"zebrafish_lethal"

#' Snow gum leaf PSII functional-impairment thermal-tolerance data
#'
#' Chlorophyll-fluorescence (\eqn{F_v/F_m}) measurements on excised snow gum
#' (\emph{Eucalyptus pauciflora} var. \emph{pauciflora}) leaf segments before and
#' 24 h after heat exposure. The response is the continuous proportion
#' \code{fvfm_prop} (post/pre ratio), modelled with a Beta likelihood. The
#' model-ready frame for Case Study 3 (sublethal, continuous-proportion endpoint).
#'
#' @format A data frame with 319 rows and 8 variables:
#' \describe{
#'   \item{Temp}{Assay temperature (degrees C).}
#'   \item{Time}{Exposure duration (minutes).}
#'   \item{initial_fvfm}{\eqn{F_v/F_m} measured before heat exposure.}
#'   \item{final_fvfm}{\eqn{F_v/F_m} measured 24 h after heat exposure.}
#'   \item{fvfm_prop}{Retained PSII function, \code{final_fvfm / initial_fvfm}
#'         (a proportion in [0, 1]).}
#'   \item{Unique_ID}{Leaf-segment / tree identifier.}
#'   \item{G_Room}{Glasshouse room (grouping factor).}
#'   \item{Day}{Measurement day (grouping factor).}
#' }
#' @source Snow gum leaf PSII TDT assay (Case Study 3). Raw file:
#'   \code{system.file("extdata", "data_function_PSII_TDT_snowgum.csv", package = "bayesTLS")}.
#' @examples
#' std <- standardize_data(snowgum_psii, temp = "Temp", duration = "Time",
#'                         proportion = "fvfm_prop",
#'                         random_effects = c("Day", "G_Room"),
#'                         duration_unit = "minutes")
"snowgum_psii"

#' Acacia seed lethal thermal-death-time data
#'
#' Lethal-TDT trials for imbibed seeds of mountain hickory (\emph{Acacia
#' penninervis} var. \emph{penninervis}): the number of seeds respiring (alive)
#' after heat exposure, out of five per temperature x duration combination,
#' measured by closed-system respirometry. The model-ready frame for Case Study
#' 4. A deliberately sparse design (5 temperatures x 4 durations, n = 5).
#'
#' @format A data frame with 20 rows and 6 variables:
#' \describe{
#'   \item{Species}{Species name.}
#'   \item{Temperature}{Assay temperature (degrees C).}
#'   \item{Duration}{Exposure duration (minutes).}
#'   \item{Seeds_Total}{Number of seeds in the trial (5).}
#'   \item{Seeds_Resp}{Number respiring (alive) after exposure.}
#'   \item{Proportion}{\code{Seeds_Resp / Seeds_Total}.}
#' }
#' @source Acacia seed lethal-TDT assay (Case Study 4). Raw file:
#'   \code{system.file("extdata", "data_lethal_TDT_acacia_seeds.csv", package = "bayesTLS")}.
#' @examples
#' std <- standardize_data(acacia_seeds, temp = "Temperature",
#'                         duration = "Duration", n_total = "Seeds_Total",
#'                         n_surv = "Seeds_Resp", duration_unit = "minutes")
"acacia_seeds"
