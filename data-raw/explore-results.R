library(springRunDSM)
library(tidyverse)
library(plotly)
r2r_seeds <- springRunDSM::spring_run_model(scenario = NULL, mode = "seed",
                                            seeds = NULL, ..params = springRunDSM::r_to_r_baseline_params,
                                            delta_surv_inflation = FALSE)

r2r_model_results <- springRunDSM::spring_run_model(mode = "simulate", 
                                                    ..params = springRunDSM::r_to_r_kitchen_sink_params,
                                                    seeds = r2r_seeds,
                                                    delta_surv_inflation = TRUE)

r2r_model_results$spawners
r2r_model_results$phos

spawn <- dplyr::as_tibble(r2r_model_results$spawners) |>
  dplyr::mutate(location = fallRunDSM::watershed_labels) |>
  pivot_longer(cols = c(`1`:`20`), values_to = 'spawners', names_to = "year") |>
  group_by(year, location) |>
  summarize(total_spawners = sum(spawners)) |>
  # filter(!location %in% c("Feather River", "Butte Creek")) |>
  mutate(year = as.numeric(year)) |>
  ggplot(aes(year, total_spawners, color = location)) +
  geom_line() +
  theme_minimal() +
  labs(y = "Spawners",
       x = "Year") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = 1:20) +
  theme(text = element_text(size = 20))

plotly::ggplotly(spawn)


# CHECK against grandtab
grandtab_totals <- dplyr::as_tibble(DSMCalibrationData::grandtab_observed$spring)|> #change which results to look at diff plots
  dplyr::mutate(location = fallRunDSM::watershed_labels) |>
  pivot_longer(cols = c(`1998`:`2017`), values_to = 'spawners', names_to = "year") |>
  # filter(!location %in% non_spawn_regions) |>
  group_by(year,
           location
  ) |>
  summarize(total_spawners = sum(spawners, na.rm = TRUE)) |>
  mutate(year = as.numeric(year)) |>
  ggplot(aes(year, total_spawners,
             color = location
  )) +
  geom_line() +
  theme_minimal() +
  labs(y = "Spawners",
       x = "Year") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = 1:20) +
  theme(text = element_text(size = 20))

ggplotly(grandtab_totals)


plot_total_spawners <- function(model_results,
                                result_type = c("Total Spawners", "Total Natural Spawners")) {
  if (result_type == "Total Natural Spawners") {
    spawn <- dplyr::as_tibble(model_results$spawners * model_results$proportion_natural) |>
      dplyr::mutate(location = fallRunDSM::watershed_labels)
  }
  else {
    spawn <- dplyr::as_tibble(model_results$spawners) |>
      dplyr::mutate(location = fallRunDSM::watershed_labels)
  }
  
  spawn |>
    pivot_longer(cols = c(`1`:`20`), values_to = 'spawners', names_to = "year") |>
    group_by(year) |>
    summarize(total_spawners = sum(spawners)) |>
    mutate(year = as.numeric(year))  |> 
    ggplot(aes(year, total_spawners)) +
    geom_line() +
    theme_minimal() +
    labs(y = result_type,
         x = "Year") +
    scale_y_continuous(labels = scales::comma) +
    scale_x_continuous(breaks = 1:20) +
    theme(text = element_text(size = 20))
}

plot_total_spawners(r2r_model_results,"Total Spawners")
