library(dplyr)
library(stringr)

#' Deconstruct plant IDs into their component parts
#'
#' Handles all of these forms (and combinations of them):
#'   A1.STIPUL.1        basic: plot A1, species STIPUL, individual 1
#'   23.A1.STIPUL.1     cohort prefix (1 yr old at collection)
#'   A1.STIPUL.1-G6     substitute individual taken from plot G6
#'   A1.STIPUL.1-2      second sample from the same individual
#'   A1.STIPUL.1a / 1b  two subsamples of the same leaf
#'   23.A1.STIPUL.1-G6-2  any combination of the above
#'
#' @param df     data frame containing the ID column
#' @param id_col name of the ID column (default "Sample ID")
#' @return df with parsed columns inserted right after the ID column,
#'         all original columns retained
parse_plant_id <- function(df, id_col = "Sample ID") {

  # find the ID column, tolerating odd spacing/case in the header
  # (e.g. "Sample ID " or a non-breaking space from Excel)
  norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
  hit  <- names(df)[norm(names(df)) == norm(id_col)]
  if (length(hit) != 1) {
    stop("Couldn't find a column called '", id_col, "'. Columns are:\n  ",
         paste(names(df), collapse = "\n  "),
         "\nPass the right one with id_col = \"...\"", call. = FALSE)
  }
  id_col <- hit
  ids <- str_trim(as.character(df[[id_col]]))

  # 1 cohort | 2 structure | 3 plot no. | 4 genus | 5 epithet | 6 individual | 7 suffixes
  pattern <- paste0(
    "^(?:(\\d{2})\\.)?",              # optional 2-digit cohort prefix, e.g. "23."
    "([A-Za-z]+)(\\d+)\\.",           # plot: structure letter(s) + plot number
    "([A-Za-z]{3})([A-Za-z]{3})\\.",  # 6-letter species code (genus + epithet)
    "(\\d+)",                         # individual number
    "([A-Za-z]?)",                    # optional a/b = subsample of the same leaf
    "((?:-[A-Za-z0-9]+)*)$"           # zero or more "-xx" suffixes
  )
  m <- str_match(ids, pattern)

  suffix <- m[, 9]
  # "-G6" style suffix (letters + digits) = plot the substitute came from
  donor  <- str_match(suffix, "-([A-Za-z]+\\d+)(?=-|$)")[, 2]
  # "-2" style suffix (digits only) = repeat sample from same individual
  rep    <- str_match(suffix, "-(\\d+)(?=-|$)")[, 2]

  parsed <- tibble(
    parse_ok       = !is.na(m[, 1]),
    cohort         = m[, 2],
    one_year_old   = !is.na(m[, 2]),
    structure      = toupper(m[, 3]),
    plot_num       = as.integer(m[, 4]),
    plot           = paste0(toupper(m[, 3]), m[, 4]),
    species_code   = toupper(paste0(m[, 5], m[, 6])),
    genus_code     = toupper(m[, 5]),
    epithet_code   = toupper(m[, 6]),
    individual     = as.integer(m[, 7]),
    leaf_subsample = na_if(tolower(m[, 8]), ""),  # "a"/"b", NA if none
    is_subsample   = !is.na(leaf_subsample),
    substitute     = !is.na(donor),
    source_plot    = toupper(donor),     # NA unless "-G6" style suffix
    sample_rep     = as.integer(coalesce(rep, "1")),
    is_resample    = !is.na(rep)
  ) %>%
    mutate(
      # clean up rows that didn't parse (standards, blanks, typos)
      across(-parse_ok, ~ if_else(parse_ok, .x, NA)),
      plot = if_else(parse_ok, plot, NA_character_),
      # plot the tissue physically came from
      collected_from = coalesce(source_plot, plot),
      # one key per physical plant (drops a/b and the "-2" replicate),
      # handy for averaging all samples from a plant
      plant_key = if_else(
        parse_ok,
        paste0(
          if_else(is.na(cohort), "", paste0(cohort, ".")),
          plot, ".", species_code, ".", individual,
          if_else(is.na(source_plot), "", paste0("-", source_plot))
        ),
        NA_character_
      ),
      # one key per leaf sample: keeps the "-2" replicate, drops a/b,
      # so a/b subsamples of the same leaf share a key
      leaf_key = if_else(parse_ok, paste0(plant_key, "_rep", sample_rep),
                         NA_character_)
    )

  if (any(!parsed$parse_ok)) {
    bad <- unique(ids[!parsed$parse_ok])
    warning(length(bad), " ID(s) did not match any known format: ",
            paste(head(bad, 10), collapse = ", "),
            if (length(bad) > 10) ", ..." else "", call. = FALSE)
  }

  # keep every original column, parsed columns sit right after the ID
  bind_cols(df, parsed) %>%
    relocate(all_of(names(parsed)), .after = all_of(id_col))
}

# ---- usage ----
# cn_parsed <- parse_plant_id(cn)                       # uses "Sample ID"
# cn_parsed <- parse_plant_id(cn, id_col = "Plant_ID")  # if you renamed it
#
# filter(cn_parsed, !parse_ok)   # inspect anything that didn't parse
#
# # average a/b subsamples to one row per leaf sample:
# cn_parsed %>%
#   group_by(leaf_key, plant_key, plot, species_code, cohort) %>%
#   summarise(across(c(`%N`, `δ13CVPDB (‰)`, `δ15NAir (‰)`), mean), .groups = "drop")
#
# # average all samples to one row per plant:
# cn_parsed %>%
#   group_by(plant_key, plot, species_code, cohort) %>%
#   summarise(across(c(`%N`, `δ13CVPDB (‰)`, `δ15NAir (‰)`), mean), .groups = "drop")
