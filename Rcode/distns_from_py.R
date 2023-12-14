library(reticulate)
library(tidyverse)

sys <- import("sys")
bi <- import_builtins()
pickle <- import("pickle")
states <- c("CA", "GA")


# Linelist delays ---------------------------------------------------------

linelist <- lapply(states, function(state) {
  state_file <- paste0(
    "LI_Dec_pres_plot_gen/Pycharm/data/", state, 
    "_Empshrink_delay_distribution_d60c_nov5_1.25.p")
  st_list <- pickle$load(bi$open(state_file, "rb"))
  df <- map(st_list, function(st) {
    tibble(
      state = state,
      dist = unname(unlist(st)),
      delay = 0:(length(st) - 1)
    )})
  names(df) <- names(st_list)
  bind_rows(df, .id = "Date") |>
    mutate(Date = as.Date(Date))
})

report_delay <- bind_rows(linelist)
write_rds(report_delay, "data/report_delay_ca_ga.rds")


# variant mix -------------------------------------------------------------

vmix <- read_csv("LI_Dec_pres_plot_gen/R/data/seq_df_post_decon_nov5.csv")
write_rds(vmix, "data/variant_mix.rds")

# Incubation period by state ----------------------------------------------

incubation <- lapply(states, function(state) {
  state_file <- paste0(
    "LI_Dec_pres_plot_gen/Pycharm/data/", state, 
    "_inc_distribution_nov5.p")
  st_list <- pickle$load(bi$open(state_file, "rb"))
  df <- map(st_list, function(st) {
    tibble(
      state = state,
      dist = unname(unlist(st)),
      delay = 0:(length(st) - 1)
    )})
  names(df) <- names(st_list)
  bind_rows(df, .id = "Date") |>
    mutate(Date = as.Date(Date))
})

incubation <- bind_rows(incubation)
write_rds(incubation, "data/incubation_ca_ga.rds")


# Convolved distributions -------------------------------------------------

convolved <- lapply(states, function(state) {
  state_file <- paste0(
    "LI_Dec_pres_plot_gen/Pycharm/data/", state, 
    "_conv_distribution_d60c_nov5.p")
  st_list <- pickle$load(bi$open(state_file, "rb"))
  df <- map(st_list, function(st) {
    tibble(
      state = state,
      dist = unname(unlist(st)),
      delay = 0:(length(st) - 1)
    )})
  names(df) <- names(st_list)
  bind_rows(df, .id = "Date") |>
    mutate(Date = as.Date(Date))
})

convolved <- bind_rows(convolved)
write_rds(convolved, "data/convolved_ca_ga.rds")


# deconvolved cases -------------------------------------------------------

np <- import("numpy")

deconvolved <- lapply(states, function(state) {
  state_file <- paste0(
    "LI_Dec_pres_plot_gen/Pycharm/data/", state, 
    "_deconvolution_res_d60c_nov5.npy")
  infects <- np$load(state_file)
  tibble(
    geo_value = tolower(state),
    infections = as.vector(infects),
    time_value = seq(
      as.Date("2020-03-01"), 
      length.out = length(infects), by = "1 day")
  )
}) |>
  bind_rows()

jhu <- epidatr::pub_covidcast(
  source = "jhu-csse", 
  signals = "confirmed_7dav_incidence_num", 
  time_type = "day",
  time_value = epidatr::epirange(20200301,20230301),
  geo_type = "state",
  geo_value = "ca,ga",
  asof = as.Date("2023-07-06")
) |>
  select(geo_value, time_value, cases = value) |>
  mutate()

deconvolved |>
  left_join(jhu) |>
  relocate(geo_value, time_value) |>
  write_rds("data/deconvolved_ca_ga.rds")
