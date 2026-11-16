# =============================================================================
# data-raw/make_datasets.R
# Builds the package datasets in data/ from the public source files in
# data-raw/sources/. Run from the package root:  Rscript data-raw/make_datasets.R
#
# Modelling units: LAD May 2024 geography with two merges forced by the data:
#   Cornwall + Isles of Scilly -> E06000052 (births published combined)
#   City of London + Hackney   -> E09000012 (Fingertips publishes combined)
# =============================================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(sf); library(spdep)
})
src <- "data-raw/sources"

merge_map <- c(E06000053 = "E06000052", E09000001 = "E09000012")
unit_of   <- function(code) ifelse(code %in% names(merge_map), merge_map[code], code)
recode_2025 <- c(E08000038 = "E08000016", E08000039 = "E08000019")   # Barnsley, Sheffield
pre2023 <- c(   # Census 2021 LTLAs abolished April 2023
  E07000026 = "E06000063", E07000028 = "E06000063", E07000029 = "E06000063",
  E07000027 = "E06000064", E07000030 = "E06000064", E07000031 = "E06000064",
  E07000163 = "E06000065", E07000164 = "E06000065", E07000165 = "E06000065",
  E07000166 = "E06000065", E07000167 = "E06000065", E07000168 = "E06000065", E07000169 = "E06000065",
  E07000187 = "E06000066", E07000188 = "E06000066", E07000189 = "E06000066", E07000246 = "E06000066")

# ---- geometry and adjacency ----------------------------------------------------
geo <- st_read(file.path(src, "LAD_May2024_UK_BGC.gpkg"), quiet = TRUE) |>
  filter(startsWith(LAD24CD, "E")) |> st_transform(27700) |> st_make_valid()
geo$unit <- unit_of(geo$LAD24CD)
lad_sf <- geo |> group_by(unit) |>
  summarise(name = paste(LAD24NM, collapse = " & "), .groups = "drop") |> st_make_valid()
lad_sf$name[lad_sf$unit == "E06000052"] <- "Cornwall & Isles of Scilly"
lad_sf$name[lad_sf$unit == "E09000012"] <- "Hackney & City of London"
lad_sf <- lad_sf[order(lad_sf$unit), ]
lad_nb <- poly2nb(lad_sf, queen = TRUE, snap = 50)   # adjacency from full-resolution BGC
# Isle of Wight: link to the mainland LADs served by its ferry routes
iow <- which(lad_sf$unit == "E06000046")
stopifnot(length(lad_nb[[iow]]) == 1 && lad_nb[[iow]][1] == 0L)
for (j in match(c("E06000044", "E06000045", "E07000091"), lad_sf$unit)) {
  lad_nb[[iow]] <- sort(unique(c(setdiff(lad_nb[[iow]], 0L), j)))
  lad_nb[[j]]   <- sort(unique(c(lad_nb[[j]], iow)))
}
attr(lad_nb, "region.id") <- lad_sf$unit
stopifnot(n.comp.nb(lad_nb)$nc == 1)
# simplify only the stored geometry (for mapping), after adjacency has been computed
lad_sf <- st_make_valid(st_simplify(lad_sf, dTolerance = 20, preserveTopology = TRUE))
lad_sf <- st_as_sf(as.data.frame(lad_sf))   # plain data.frame-based sf

# ---- upper-tier membership ---------------------------------------------------
lad_utla_lookup <- read.csv(file.path(src, "lad_to_utla_lookup.csv"), check.names = FALSE) |>
  filter(`Parent Code` != "E92000001") |>
  transmute(code = `Area Code`, utla = `Parent Code`, utla_name = `Parent Name`) |>
  mutate(unit = unit_of(code)) |>
  distinct(unit, .keep_all = TRUE) |> select(unit, utla, utla_name)
lad_utla_lookup$utla[lad_utla_lookup$unit == "E09000012"] <- "E09000012"
lad_utla_lookup$utla[lad_utla_lookup$unit == "E06000052"] <- "E06000052"

# ---- covariates that do not depend on the births year -------------------------
imd_raw <- read_excel(file.path(src, "IoD2025_IMD_LSOA.xlsx"), sheet = "IMD25")
names(imd_raw) <- c("lsoa", "lsoa_name", "lad", "lad_name", "imd_rank", "imd_decile")
N_lsoa <- nrow(imd_raw)
imd <- imd_raw |> mutate(unit = unit_of(lad)) |> group_by(unit) |>
  summarise(imd_dep_pctile = mean(100 * (1 - (imd_rank - 0.5) / N_lsoa)),
            pct_lsoa_most_dep_quintile = 100 * mean(imd_decile <= 2), .groups = "drop")

e <- read_excel(file.path(src, "ethnicities.xlsx"), sheet = "Dataset")
names(e) <- c("code", "name", "grp_code", "grp", "obs")
south_asian_rx <- paste0("^Asian, Asian British or Asian Welsh: (Indian|Pakistani|Bangladeshi|",
                         "Sri Lankan|Kashmiri|Punjabi|Sikh|Tamil|Nepali|Gujarati|Goan|Bengali)")
eth <- e |> filter(grp_code >= 0) |>
  mutate(code = ifelse(code %in% names(pre2023), pre2023[code], code), unit = unit_of(code)) |>
  group_by(unit) |>
  summarise(pop = sum(obs),
            wb = sum(obs[grepl("^White: English, Welsh, Scottish, Northern Irish or British$", grp)]),
            black = sum(obs[grepl("^Black", grp)]),
            sa = sum(obs[grepl(south_asian_rx, grp)]), .groups = "drop") |>
  transmute(unit, pct_black = 100 * black / pop, pct_south_asian = 100 * sa / pop,
            pct_non_white_british = 100 * (1 - wb / pop))

fp <- read.csv(file.path(src, "fingertips_maternity_indicators.csv"), check.names = FALSE) |>
  filter(`Area Type` == "Districts & UAs (from Apr 2023)")
fp_get <- function(id, period) fp |> filter(`Indicator ID` == id, `Time period` == period) |>
  transmute(code = `Area Code`, value = as.numeric(Value), count = as.numeric(Count),
            denom = as.numeric(Denominator), note = `Value note`) |>
  filter(!code %in% c("E06000053", "E09000001"))      # combined values carried by parent

fingertips_covariates <- function(period) {
  # maternal age 35+: counts; suppressed LADs (Torbay) use the latest published year
  age_all <- fp |> filter(`Indicator ID` == 92250, !is.na(Count)) |>
    transmute(code = `Area Code`, period_fb = `Time period`,
              count_fb = as.numeric(Count), denom_fb = as.numeric(Denominator)) |>
    filter(period_fb <= period) |> group_by(code) |> slice_max(period_fb, n = 1) |> ungroup()
  age <- fp_get(92250, period) |> left_join(age_all, by = "code") |>
    mutate(age_period = ifelse(is.na(count), period_fb, period),
           count = coalesce(count, count_fb), denom = coalesce(denom, denom_fb),
           age_dq_flag = grepl("data quality", note, ignore.case = TRUE), unit = unit_of(code)) |>
    group_by(unit) |>
    summarise(pct_age35plus = 100 * sum(count) / sum(denom), age_dq_flag = any(age_dq_flag),
              age_period = first(age_period), .groups = "drop")
  smk <- fp_get(93085, period) |>
    transmute(unit = unit_of(code), pct_smoking = value,
              smoking_dq_flag = grepl("data quality", note, ignore.case = TRUE))
  obs <- fp_get(93881, period) |> transmute(unit = unit_of(code), pct_obese_adult = value)
  age |> left_join(smk, by = "unit") |> left_join(obs, by = "unit")
}

# ---- births & stillbirths by year ---------------------------------------------
read_births <- function(file) {
  b <- read_excel(file, sheet = "Table_2", col_names = FALSE, skip = 6)
  names(b) <- c("code", "name", "geog", "age", "live_births", "stillbirths", "rate", "unreliable")
  b |> filter(age == "All ages",
              geog %in% c("Unitary Authority", "Non-metropolitan District",
                          "Metropolitan District", "London Borough"),
              startsWith(code, "E")) |>
    mutate(code = ifelse(code == "E06000052, E06000053", "E06000052", code),
           code = ifelse(code %in% names(recode_2025), recode_2025[code], code),
           unit = unit_of(code)) |>
    group_by(unit) |>
    summarise(live_births = sum(as.numeric(live_births)), stillbirths = sum(as.numeric(stillbirths)),
              .groups = "drop") |>
    mutate(total_births = live_births + stillbirths)
}

build_year <- function(year, fp_period) {
  out <- st_drop_geometry(lad_sf) |>
    left_join(read_births(file.path(src, paste0(year, "birthregistrations.xlsx"))), by = "unit") |>
    left_join(imd, by = "unit") |> left_join(eth, by = "unit") |>
    left_join(fingertips_covariates(fp_period), by = "unit") |>
    left_join(lad_utla_lookup |> select(unit, utla), by = "unit") |>
    mutate(year = year, fingertips_period = fp_period) |>
    relocate(year, .after = name) |> as.data.frame()
  core <- c("live_births", "stillbirths", "imd_dep_pctile", "pct_black", "pct_south_asian",
            "pct_age35plus", "pct_smoking", "pct_obese_adult", "utla")
  stopifnot(!anyNA(out[core]), identical(out$unit, lad_sf$unit))
  out
}
stillbirths_lad_2024 <- build_year(2024, "2024/25")
stillbirths_lad_2025 <- build_year(2025, "2024/25")

# ---- obesity in early pregnancy (upper-tier, as published) ----------------------
obesity_pregnancy_utla <- read.csv(file.path(src, "fingertips_obesity_pregnancy_utla.csv"),
                                   check.names = FALSE) |>
  filter(`Area Type` != "England", !is.na(Value)) |>
  transmute(utla = `Area Code`, utla_name = `Area Name`, period = `Time period`,
            pct_obese_preg = as.numeric(Value), n_obese = as.numeric(Count),
            n_bookings = as.numeric(Denominator)) |>
  filter(utla %in% lad_utla_lookup$utla) |> as.data.frame()

save(stillbirths_lad_2024, file = "data/stillbirths_lad_2024.rda", compress = "xz")
save(stillbirths_lad_2025, file = "data/stillbirths_lad_2025.rda", compress = "xz")
save(lad_sf, file = "data/lad_sf.rda", compress = "xz")
save(lad_nb, file = "data/lad_nb.rda", compress = "xz")
save(lad_utla_lookup, file = "data/lad_utla_lookup.rda", compress = "xz")
save(obesity_pregnancy_utla, file = "data/obesity_pregnancy_utla.rda", compress = "xz")

for (y in list(stillbirths_lad_2024, stillbirths_lad_2025))
  cat(y$year[1], ": units", nrow(y), " live births", sum(y$live_births), " stillbirths",
      sum(y$stillbirths), " rate/1000 total births",
      round(1000 * sum(y$stillbirths) / sum(y$total_births), 3), "\n")
cat("obesity_pregnancy_utla rows:", nrow(obesity_pregnancy_utla), " UTLAs:",
    length(unique(obesity_pregnancy_utla$utla)), "\n")
