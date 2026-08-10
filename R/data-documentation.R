#' Story Charcoal Matrix
#'
#' A matrix of charcoal accumulation rates from Story Lake, Indiana.
#'    This matrix has been binned at a 50-year resolution. Values are on the
#'    original measurement scale; scale them (e.g. with [scale()]) before
#'    passing them to [mnGLMM()] or [mnTS()].
#'
#' @format A `1` column matrix with `160` rows (observations).
#' @source For details see vignette
#'
#' @examples
#' data(story_char_matrix)
#' head(story_char_matrix)
"story_char_matrix"


#' Story Charcoal Wide
#'
#' A ~10,000-year data set from Story Lake, Indiana, of modelled ages, sample
#'    depths, and charcoal accumulation rate from a sediment core sample.
#'
#' @format A tibble with `736` rows and `3` columns.
#' @source For details see https://doi.org/10.1111/1365-2745.14289
#'
#' @examples
#' data(story_char_wide)
#' str(story_char_wide)
"story_char_wide"


#' Story Pollen Matrix
#'
#' A matrix of pollen taxa counts from Story Lake, Indiana. Data are organised
#'    into a reference group of "other" taxa, a functional group of hardwood
#'    taxa, and target taxa of Fagus grandifolia, Ulmus spp., and Quercus spp.
#'    Data have been binned at a 50-year resolution with empty bins filled with
#'    0-values (i.e., rows of 0 are empty bins, see vignette for details).
#'
#' @format A matrix with `160` rows (samples) and `5` columns (taxa).
#' @source For details see vignette
#'
#' @examples
#' data(story_pollen_matrix)
#' summary(story_pollen_matrix)
"story_pollen_matrix"



#' Story Pollen Wide
#'
#' A ~10,000-year data set from Story Lake, Indiana, of modelled ages, sample
#'    depths, pollen taxa. Data are organised into a reference group of
#'    "other" taxa, a functional group of hardwood taxa, and target taxa
#'    of Fagus grandifolia, Ulmus spp., and Quercus spp.
#'
#' @format A tibble with `95` rows (observations) and `7` columns.
#' @source For details see https://doi.org/10.1111/1365-2745.14289
#'
#' @examples
#' data(story_pollen_wide)
#' names(story_pollen_wide)
"story_pollen_wide"



#' Story 100 bootstraps
#'
#' 100 bootstraps of the story data as modelled in the vignette. This is the
#'    \code{coef.table.boot} element of a \code{\link{boot.mnTS}} result, not
#'    the full bootstrap object.
#'
#' @format A matrix with `29` rows (coefficients) and `8` columns:
#'    `boot_mean`, `se`, `t`, `P`, `upper68`, `lower68`, `n_fail_converged`
#'    and `n_total`.
#' @source Produced by bootstrapping the Story Lake model; see vignette.
#'
#' @examples
#' data(story_100_bootstraps)
#' head(story_100_bootstraps)
"story_100_bootstraps"


#' Sunfish Pond pollen
#'
#' A ~13,800 yr pollen record organised into
#' 5 focal taxa, a reference group, and an estimated age per sample
#'
#' @format A tibble with `69` rows (samples) and `7` columns: `age`, the
#'    reference group `other`, and the focal taxa `Betula`, `Fagus`,
#'    `P.strobus`, `Quercus` and `Tsuga`.
#' @source Johnson et al., in review. See the Sunfish Pond vignette.
#'
#' @examples
#' data(sunfish_pollen_wide)
"sunfish_pollen_wide"


#' Sunfish Pond level
#'
#' A ~14,950 yr lake level record from Sunfish Pond, used as a proxy for
#' hydroclimatic variability and plant water availability.
#'
#' @format A tibble with `299` rows and `4` columns: `age`, `lake_lvl`
#'    (cm relative to modern), and the reconstruction bounds `UpperBound`
#'    and `LowerBound`.
#' @source Johnson et al., in review. See the Sunfish Pond vignette.
#'
#' @examples
#' data(sunfish_ll)
"sunfish_ll"

#' Sunfish Pond bootstrapped sensitivity
#'
#' Bootstrapped sensitivity to binning width
#'
#' @format A list with one element per binning resolution (75- and 150-year
#'    bins), each the output of \code{\link{boot.mnTS}}.
#' @source Produced by the code in the Sunfish Pond vignette.
#'
#' @examples
#' data(bootstrap_sensitivity)
#' names(bootstrap_sensitivity)
"bootstrap_sensitivity"