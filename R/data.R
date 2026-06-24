# Documentation for the case-study datasets shipped with bayesTLS.
# The datasets are built from the raw CSVs in inst/extdata/ by
# data-raw/make_datasets.R. Each is analysis-ready for the workflow
# standardize_data() -> fit_4pl() -> extract_tdt().

#' Brown shrimp lethal thermal-death-time data
#'
#' Replicate lethal-TDT trials for brown shrimp (\emph{Crangon crangon}). Each
#' row is one tank of individuals exposed to a fixed assay temperature for a
#' fixed duration; the response is the proportion that died. The model-ready
#' frame for Case Study 1 (lethal endpoint).
#'
#' @format A data frame with 148 rows and 6 variables:
#' \describe{
#'   \item{Date}{Experiment date (use as a grouping factor).}
#'   \item{Tank}{Holding-tank identifier (use as a grouping factor).}
#'   \item{Temperature_assay}{Assay temperature (degrees C).}
#'   \item{Duration_exposure_hours}{Exposure duration (hours).}
#'   \item{N_individuals_after_trial}{Number of individuals in the trial.}
#'   \item{Mortality_after_trial}{Proportion that died during the trial
#'         (deaths / \code{N_individuals_after_trial}), in the unit interval.
#'         Consumed by
#'         \code{standardize_data(mortality = "Mortality_after_trial")}.}
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
#' (\emph{Eucalyptus pauciflora}) leaf sections before and 16--24 h after heat
#' exposure, from Experiment 1 (post-heat light vs dark recovery) of Arnold et al.
#' (2026). Short branches were cut from six mature trees grown outdoors in
#' Canberra, ACT; \eqn{\sim}1 cm\eqn{^2} leaf sections were dark-adapted, given an
#' initial \eqn{F_v/F_m}, then submerged in a temperature-controlled water bath
#' under sub-saturating light across a grid of assay temperatures (30--56 degrees
#' C) and exposure durations (5--120 min). After heat, paired arrays were held for
#' 90 min in moderate light (\code{recovery = "Light"}) or in darkness
#' (\code{recovery = "Dark"}); a final \eqn{F_v/F_m} was taken 16--24 h later. The
#' response is the continuous proportion \code{fvfm_prop} (post/pre ratio),
#' modelled with a Beta likelihood. The model-ready frame for the leaf PSII case
#' study (sublethal, continuous-proportion endpoint). The light/dark contrast is a
#' two-group categorical moderator; in the source experiment post-heat light
#' lowered apparent heat tolerance.
#'
#' @format A data frame with 396 rows and 8 variables:
#' \describe{
#'   \item{Temp}{Assay temperature (degrees C); 30, 35, 40, 44, 48, 52, 56.}
#'   \item{Time}{Exposure duration (minutes); 5, 15, 30, 60, 120.}
#'   \item{recovery}{Post-heat recovery light condition: \code{"Dark"} (darkness
#'         immediately after heat) or \code{"Light"} (90 min moderate light
#'         post-heat). A two-level moderator.}
#'   \item{plant}{Replicate mature tree (factor, 6 levels); the natural random-effect
#'         grouping.}
#'   \item{meas_day}{Assay day (factor, 2 levels). Two levels only, so a poor
#'         random-effect grouping; better treated as fixed or omitted.}
#'   \item{initial_fvfm}{\eqn{F_v/F_m} measured before heat exposure.}
#'   \item{final_fvfm}{\eqn{F_v/F_m} measured 16--24 h after heat exposure.}
#'   \item{fvfm_prop}{Retained PSII function, \code{final_fvfm / initial_fvfm}
#'         (a proportion in the unit interval; 0 indicates complete loss of
#'         measurable PSII function).}
#' }
#' @source Arnold PA, Harris RJ, Aitken SM, Hoek MM, Cook AM, Leigh A, Nicotra AB
#'   (2026) Towards a standard approach to investigating the thermal load
#'   sensitivity of photosystem II via chlorophyll fluorescence. bioRxiv
#'   \doi{10.64898/2026.04.09.717599} (CC BY-NC 4.0), Experiment 1, snow gum
#'   slice. Raw file:
#'   \code{system.file("extdata", "data_function_PSII_TDT_snowgum.csv", package = "bayesTLS")}.
#' @examples
#' std <- standardize_data(snowgum_psii, temp = "Temp", duration = "Time",
#'                         proportion = "fvfm_prop",
#'                         random_effects = "plant",
#'                         duration_unit = "minutes")
"snowgum_psii"

#' Drosophila suzukii multi-trait thermal-tolerance data
#'
#' Per-individual thermal-tolerance assays for the spotted-wing fly
#' (\emph{Drosophila suzukii}), one row per fly, carrying three thermal-tolerance
#' endpoints measured under static heat exposures at 34--38 degrees C:
#' a lethal endpoint (\code{dead}), a sublethal knockdown time-to-event
#' (\code{t_coma}), and a sublethal reproductive endpoint (\code{prod}). The
#' model-ready frame for Case Study 4; aggregate \code{dead} to counts for the
#' beta-binomial lethal fit, or use \code{t_coma} / \code{prod} directly for the
#' sublethal endpoints. \code{lvl} indexes the exposure-duration grid as a
#' percentage of the estimated median time-to-coma from the authors' initial TDT
#' curves; \code{time} is the realised duration in minutes.
#'
#' @format A data frame with 1407 rows and 9 variables:
#' \describe{
#'   \item{id}{Unique individual identifier (\code{temp-lvl-sex-rep}).}
#'   \item{temp}{Assay temperature (degrees C).}
#'   \item{lvl}{Exposure duration as a percentage of the estimated median
#'         time-to-coma from the authors' initial TDT curves.}
#'   \item{time}{Exposure duration (minutes).}
#'   \item{sex}{Sex, a factor with levels \code{F}, \code{M}.}
#'   \item{rep}{Replicate vial within a temperature x lvl x sex cell.}
#'   \item{prod}{Reproductive productivity (offspring per female per day).}
#'   \item{dead}{Mortality indicator: \code{1} = died, \code{0} = survived.}
#'   \item{t_coma}{Time to heat coma (minutes); \code{NA} where no coma was
#'         recorded for that individual.}
#' }
#' @source \enc{Ørsted}{Orsted} M, Willot Q, Olsen AK, Kongsgaard V, Overgaard J
#'   (2024). Data for: Thermal limits of survival and reproduction depend on
#'   stress duration: a case study of \emph{Drosophila suzukii}. Zenodo,
#'   \doi{10.5281/zenodo.10602268} (distributed under CC BY 4.0). Associated
#'   article: \doi{10.1111/ele.14421}. Raw file:
#'   \code{system.file("extdata", "data_multitrait_TDT_drosophila_suzukii.csv", package = "bayesTLS")}.
#' @examples
#' # Lethal endpoint: aggregate per-individual deaths to cell counts, then
#' # standardise for the beta-binomial 4PL.
#' cells <- dplyr::summarise(
#'   dplyr::group_by(dsuzukii, temp, time, sex),
#'   n_total = dplyr::n(), n_dead = sum(dead), .groups = "drop")
#' std <- standardize_data(cells, temp = "temp", duration = "time",
#'                         n_total = "n_total", n_dead = "n_dead",
#'                         random_effects = "sex", duration_unit = "minutes")
"dsuzukii"

#' Zebrafish lethal-TDT data across an oxygen gradient
#'
#' Survival of zebrafish (\emph{Danio rerio}) larvae assayed for upper thermal
#' tolerance under three oxygen treatments, the model-ready frame for the
#' oxygen-gradient case study. Diploid and triploid larvae were held at assay
#' temperatures of 26 (control), 38, 39 and 40 degrees C for 3.8--240 minutes
#' under hypoxia, normoxia or hyperoxia, and scored alive/dead. One row per
#' assay group; \code{oxygen} is the categorical moderator. Fit \code{CTmax} and
#' \code{z} as functions of \code{oxygen} (optionally \code{ploidy}) in one joint
#' 4PL to compare thermal tolerance across the gradient with full posterior
#' uncertainty.
#'
#' @format A data frame with 905 rows and 10 variables:
#' \describe{
#'   \item{cohort}{Larval cohort identifier.}
#'   \item{ploidy}{Ploidy, a factor with levels \code{diploid}, \code{triploid}.}
#'   \item{oxygen}{Oxygen treatment, a factor with levels \code{hypoxia},
#'         \code{normoxia}, \code{hyperoxia} (the modelling moderator).}
#'   \item{o2_nominal}{Nominal oxygen target as percent air saturation
#'         (25 / 100 / 225).}
#'   \item{o2_measured}{Measured oxygen level (percent air saturation).}
#'   \item{temp}{Target assay temperature (degrees C).}
#'   \item{temp_measured}{Measured assay temperature (degrees C).}
#'   \item{duration_min}{Exposure duration (minutes).}
#'   \item{n_total}{Number of larvae in the assay group.}
#'   \item{n_surv}{Number surviving after exposure and recovery.}
#' }
#' @source Saruhashi S, Boerrigter JGJ, Hooymans MHL, Sinclair BJ, Verberk WCEP
#'   (2026). Data and code for: Oxygen availability and oxygen delivery but not
#'   oxidative stress shape heat tolerance in diploid and triploid zebrafish
#'   larvae. Zenodo, \doi{10.5281/zenodo.20075355} (distributed under CC BY 4.0).
#'   Associated article: \emph{Journal of Experimental Biology} 229(10): jeb251548,
#'   \doi{10.1242/jeb.251548}. Raw file:
#'   \code{system.file("extdata", "data_lethal_TDT_zebrafish_oxygen.csv", package = "bayesTLS")}.
#' @examples
#' \dontrun{
#' std <- standardize_data(zebrafish_o2, temp = "temp", duration = "duration_min",
#'                         n_total = "n_total", n_surv = "n_surv",
#'                         duration_unit = "minutes")
#' wf <- fit_4pl(std, ctmax = ~ 0 + oxygen, z = ~ 0 + oxygen, t_ref = 60)
#' tls(wf, by = "oxygen", lethal = TRUE)   # z, CTmax, T_crit per oxygen treatment
#' }
"zebrafish_o2"

#' Cereal-aphid lethal-TDT data, three species across three ages
#'
#' Survival of three cereal-aphid species across a broad range of stressful high
#' and low temperatures, the model-ready frame for the multi-species case study.
#' Aphids of three ages (2, 6, 12 days old) were exposed to a heat branch
#' (34--40 degrees C) or a cold branch (-11 to -3 degrees C) for a range of
#' durations and scored alive/dead after recovery. One row per assay group;
#' \code{branch} flags the heat vs cold series. Subset to one \code{branch} (and
#' typically one \code{age}) and fit \code{CTmax} and \code{z} as functions of
#' \code{species} in one joint 4PL to compare species with full posterior
#' uncertainty.
#'
#' @format A data frame with 3041 rows and 7 variables:
#' \describe{
#'   \item{species}{Species, a factor with levels \code{M_dirhodum}
#'         (\emph{Metopolophium dirhodum}), \code{S_avenae}
#'         (\emph{Sitobion avenae}), \code{R_padi}
#'         (\emph{Rhopalosiphum padi}).}
#'   \item{age}{Age in days, a factor with levels \code{2}, \code{6}, \code{12}.}
#'   \item{branch}{Stress branch, a factor with levels \code{heat}
#'         (34--40 degrees C), \code{cold} (-11 to -3 degrees C).}
#'   \item{temp}{Assay temperature (degrees C).}
#'   \item{duration_min}{Exposure duration (minutes).}
#'   \item{n_total}{Number of aphids treated.}
#'   \item{n_surv}{Number surviving after treatment and recovery.}
#' }
#' @source Li Y-J, Chen S-Y, \enc{Jørgensen}{Jorgensen} LB, Overgaard J, Renault D,
#'   Colinet H, Ma C-S (2023). Data for: Interspecific differences in thermal
#'   tolerance landscape explain aphid community abundance under climate change.
#'   Dryad, \doi{10.5061/dryad.mcvdnck4j} (Dryad CC0). Associated article:
#'   \emph{Journal of Thermal Biology} 114: 103583, \doi{10.1016/j.jtherbio.2023.103583}.
#'   Raw file:
#'   \code{system.file("extdata", "data_lethal_TDT_aphid.csv", package = "bayesTLS")}.
#' @examples
#' \dontrun{
#' a <- subset(aphid_tdt, branch == "heat" & age == "6")
#' std <- standardize_data(a, temp = "temp", duration = "duration_min",
#'                         n_total = "n_total", n_surv = "n_surv",
#'                         duration_unit = "minutes")
#' wf <- fit_4pl(std, ctmax = ~ 0 + species, z = ~ 0 + species, t_ref = 60)
#' tls(wf, by = "species", lethal = TRUE)   # z, CTmax, T_crit per species
#' }
"aphid_tdt"
