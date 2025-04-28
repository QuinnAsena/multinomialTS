#' Story Charcoal Matrix
#'
#' A matrix of charcoal accumulation rates from Story Lake, Indiana.
#'    This matrix has been binned at a 50-year resolution and scaled.
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
