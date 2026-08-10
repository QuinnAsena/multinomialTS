# Maintenance script for the bundled datasets in data/.
#
# Run from the package root with:
#   Rscript data-raw/fix-datasets.R
#
# This records the corrections applied for 1.1.0 so that any future edit to the
# shipped data is reproducible rather than a one-off manual step. It reads and
# rewrites data/*.rda in place, so re-running it on already-corrected data is a
# no-op apart from recompression.

library(dplyr)

stopifnot(dir.exists("data"))

e <- new.env()
for (f in list.files("data", full.names = TRUE)) load(f, envir = e)

# 1. Two tibbles shipped as grouped_df, carrying a stray dplyr grouping that
#    surprised downstream group_by()/summarise() calls.
story_pollen_wide <- ungroup(get("story_pollen_wide", e))
sunfish_pollen_wide <- ungroup(get("sunfish_pollen_wide", e))

# 2. The Pinus strobus column name was truncated to "P.strobu" during data
#    preparation, and the Sunfish vignette renamed it on the fly every time.
if ("P.strobu" %in% names(sunfish_pollen_wide)) {
  sunfish_pollen_wide <- rename(sunfish_pollen_wide, P.strobus = "P.strobu")
}

save(story_pollen_wide, file = "data/story_pollen_wide.rda", compress = "xz")
save(sunfish_pollen_wide, file = "data/sunfish_pollen_wide.rda", compress = "xz")

# 3. Recompress everything with xz. R CMD check warns otherwise.
tools::resaveRdaFiles("data", compress = "xz")

cat("data/ is now", round(sum(file.size(list.files("data", full.names = TRUE))) / 1024),
    "KB\n")
