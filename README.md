# Spring Chinook Salmon Reorienting to Recovery Model

This model is adapted from the [CVPIA SIT Spring Run Model](https://github.com/CVPIA-OSC/springRunDSM) for use in the Reorienting to Recovery Structured Decision Making process.

## License

[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

IP-117068

## Usage

### Package Installation

The `springRunDSM` package depends on a number of packages developed by the [Reorienting to Recovery Organization](https://github.com/reorienting-to-recovery). To install `springRunDSM` and additional `CVPIA-OSC` packages use the `remotes::install_github()` function.

``` r
# install.packages("remotes")
remotes::install_github("reorienting-to-recovery/springRunDSM")
remotes::install_github("reorienting-to-recovery/DSMscenario")

# optional - need if wanting to explore or modify flow, habitat, and temperature inputs
remotes::install_github("reorienting-to-recovery/DSMflow")
remotes::install_github("reorienting-to-recovery/DSMhabitat")
remotes::install_github("reorienting-to-recovery/DSMtemperature")
```

### Run Model

The `spring_run_model()` is a Spring Run Chinook life cycle model used for Reorienting to Recovery Structured Decision Making Process. Running the model simulates Spring Run Chinook population dynamics across 31 watersheds in California over a 20 year period.

The following code runs the spring run model:

``` r
# seed the model
spring_run_seeds <- spring_run_model(mode = "seed")

# run the 20 year simulation
results <- spring_run_model(mode = "simulate",
                            seeds = spring_run_seeds)
```

## Details on Supporting Data

### Dependencies

The `springRunDSM` package uses data from several other packages within the [Reorienting to Recovery Organization](https://github.com/reorienting-to-recovery). These relationships are visualized in the dependency graph below.

<img src="man/figures/dependencyChain.svg" width="100%"/>

### Flow, Habitat, and Temperature Data

All data used in the `springRunDSM` is passed in as a argument to `spring_run_model()` from a `springRunDSM::r_to_r_baseline_params` data list that is composed of data objects from the following packages:

-   **Flow Data**: View detailed documentation of flow data inputs at [DSMflow](https://reorienting-to-recovery.github.io/DSMflow/). Flow inputs to the `springRunDSM` are generated using CalSim 2 data.
-   **Habitat Data**: View detailed documentation of habitat data inputs at [DSMhabitat](https://reorienting-to-recovery.github.io/DSMhabitat/). Modeling details for each stream can be viewed [here](https://reorienting-to-recovery.github.io/DSMhabitat/reference/habitat_data.html#modeling-details-for-streams).
-   **Temperature Data**: View detailed documentation of temperature data inputs at [DSMtemperature](https://reorienting-to-recovery.github.io/DSMtemperature/). Modeling details for each stream can be viewed [here](https://reorienting-to-recovery.github.io/DSMtemperature/reference/stream_temperature.html#watershed-modeling-details).

### Calibration Data

This model is calibrated using methodology and calibration data prepared for the CVPIA SIT model in the `DSMCalibration` package:

1.  [GrandTab](https://wildlife.ca.gov/Conservation/Fishes/Chinook-Salmon/Anadromous-Assessment) estimated escapement data for the years 1998-2017. The GrandTab data is prepared as `DSMCalibrationData::grandtab_observed` and is used to measure the difference between model predictions and observed escapements. Grandtab data is additionally prepared as `DSMCalibrationData::grandtab_imputed` and is used to calculate the number of juveniles during the 20 year simulation.

2.  Proxy years are used to select Habitat, Flow, and Temperature data for 1998-2017 to correspond with the years of GrandTab escapement data. The data inputs to the DSM are for years 1980-1999. We selected proxy years for 1998-2017 from the 1980-1999 model inputs by [comparing the DWR water year indices](https://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST).

For a detailed overview of the calibration process see the [calibration markdown.](https://cvpia-osc.github.io/springRunDSM/articles/calibration-2021.html)