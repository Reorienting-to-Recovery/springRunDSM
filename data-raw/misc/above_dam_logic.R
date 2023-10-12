remotes::install_github("Reorienting-to-Recovery/springRunDSM@add-r2r-params", force = TRUE)
library(springRunDSM)

spring_run_seeds <- spring_run_model(mode = "seed")
before_mods <- spring_run_model(scenario = DSMscenario::scenarios$ONE,
                                mode = "simulate",
                                seeds = spring_run_seeds, 
                                ..params = springRunDSM::r_to_r_baseline_params) 

# BUILD UPDATED PACKAGE
library(springRunDSM)
new_spring_run_seeds <- spring_run_model(mode = "seed")
# run the 20 year simulation
new_results <- spring_run_model(scenario = DSMscenario::scenarios$ONE,
                            mode = "simulate",
                            seeds = new_spring_run_seeds, 
                            ..params = springRunDSM::r_to_r_baseline_params)


new_results$spawners == before_mods$spawners

x <- surv_adult_prespawn(c(-5000:5000), r_to_r_baseline_params$..surv_adult_prespawn_int, 
                    .deg_day = r_to_r_baseline_params$.adult_prespawn_deg_day)

qplot(x = c(-5000:5000), y = x)

# plots -------------------------------------------------------------------
library(tidyverse)
spawning_watersheds <- c(DSMhabitat::watershed_species_present %>%
  filter(!(watershed_name %in% c("Upper Sacramento River", "Upper Mid Sac Region")),
         spawn) %>%
  pull(watershed_name), "San Joaquin River")

before_mods_tidy <- dplyr::as_tibble(before_mods$spawners) |>
  dplyr::mutate(location = fallRunDSM::watershed_labels) |>
  pivot_longer(cols = c(`1`:`20`), values_to = "value", names_to = "year") |>
  mutate(year = as.numeric(year)) |> 
  filter(location %in% spawning_watersheds) |> 
  mutate(scenario = "r2r-add-params")

ggplot(before_mods_tidy) + 
  geom_line(aes(x = year, y = value, color = location))

after_mods_tidy <- dplyr::as_tibble(new_results$spawners) |>
  dplyr::mutate(location = fallRunDSM::watershed_labels) |>
  pivot_longer(cols = c(`1`:`20`), values_to = "value", names_to = "year") |>
  mutate(year = as.numeric(year)) |> 
  filter(location %in% spawning_watersheds) |> 
  mutate(scenario = "with above dam logic")

ggplot(after_mods_tidy) + 
  geom_line(aes(x = year, y = value, color = location))

# both plotted together to see differences
after_mods_tidy |> 
  bind_rows(before_mods_tidy) |> 
ggplot() + 
  geom_line(aes(x = year, y = value, color = scenario)) +
  facet_wrap(~location, scales = "free_y") +
  ggtitle("above dam proportion applied; temperature above dam never exceeds 13C") +
  ylab("spawners") 

ggsave(filename = "data-raw/misc/spawners_with_above_dam_logic_applied.png", plot = last_plot(), 
       width = 10, height = 8)

# juveniles: 
before_mods_tidy_juv <- dplyr::as_tibble(before_mods$juveniles) |>
  mutate(scenario = "r2r-add-params")

after_mods_tidy_juv <- dplyr::as_tibble(new_results$juveniles) |>
  mutate(scenario = "with above dam logic")

before_mods_tidy_juv |> 
  bind_rows(after_mods_tidy_juv) |> 
  filter(watershed == "Feather River" & size == "s") |> 
ggplot() + 
  geom_line(aes(x = year, y = juveniles, color = scenario, linetype = size)) +
  facet_wrap(~watershed,  scales = "free_y") +
  ggtitle("above dam proportion applied; temperature above dam never exceeds 13C") +
  ylab("juveniles") 

ggsave(filename = "data-raw/misc/juvs_with_above_dam_logic_applied.png", plot = last_plot(), 
       width = 10, height = 8)


