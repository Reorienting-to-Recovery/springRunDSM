#' Juvenile Month Dynamics
#' @description Function to run through monthly juvenile routing and rearing logic.
#' @param hypothesis movement hypothesis to test
#' @param fish list object containing tracking matrices for juveniles, yearlings, juveniles at chipps, and adults in ocean
#' @param year year
#' @param month month
#' @param mode mode
#' @param rearing_survival rearing survival for year and month
#' @param migratory_survival migratory survival for year and month
#' @param habitat habitat for year and month
#' @param ..params full params list from model parms
#' @examples
#' fish = list(juveniles = juveniles,
#'             yearlings = yearlings,
#'             north_delta_fish = north_delta_fish,
#'             south_delta_fish = south_delta_fish,
#'             juveniles_at_chipps = juveniles_at_chipps,
#'             adults_in_ocean = adults_in_ocean)
#'@export
juvenile_month_dynamic <- function(fish, year = year, month = month,
                                   mode,
                                   rearing_survival = rearing_survival,
                                   migratory_survival = migratory_survival,
                                   habitat = habitat, ..params = ..params,
                                   avg_ocean_transition_month = avg_ocean_transition_month,
                                   stochastic = stochastic,
                                   ic_growth, fp_growth, delta_growth,
                                   filling_fn = fallRunDSM::fill_natal,
                                   filling_args = NULL,
                                   filling_regional_fn = fallRunDSM::fill_regional,
                                   filling_regional_args = NULL,
                                   movement_fn = NULL,
                                   movement_args = NULL,
                                   movement_months = NULL) {
  
  juveniles <- fish$juveniles
  # TODO confirm yearlings are empty every time
  yearlings <- matrix(0, ncol = 4, nrow = 31, dimnames = list(springRunDSM::watershed_labels, springRunDSM::size_class_labels))
  lower_mid_sac_fish <- fish$lower_mid_sac_fish
  lower_sac_fish <- fish$lower_sac_fish
  upper_mid_sac_fish <- fish$upper_mid_sac_fish
  sutter_fish <- fish$sutter_fish
  yolo_fish <- fish$yolo_fish
  san_joaquin_fish <- fish$san_joaquin_fish
  north_delta_fish <- fish$north_delta_fish
  south_delta_fish <- fish$south_delta_fish
  juveniles_at_chipps <- fish$juveniles_at_chipps
  adults_in_ocean <- fish$adults_in_ocean
  
  migrants <- matrix(0, nrow = 31, ncol = 4, dimnames = list(fallRunDSM::watershed_labels, fallRunDSM::size_class_labels))
  
  # yearling logic
  if (month == 5) {
    # yearling logic here
    # 1 - 15, 18-20, 23, 25:30
    yearlings[c(1:15, 18:20, 23, 25:30), 1:2] <- juveniles[c(1:15, 18:20, 23, 25:30), 1:2]
    juveniles[c(1:15, 18:20, 23, 25:30), 1:2] <- 0 # set all to zero since they are yearlings now
    
    # all remaining fish outmigrate
    
    migrants <- juveniles
    
    sutter_fish <- migrate(sutter_fish, migratory_survival$sutter, stochastic = stochastic)
    upper_mid_sac_fish <- migrate(upper_mid_sac_fish + migrants[1:15, ], migratory_survival$uppermid_sac, stochastic = stochastic)
    migrants[1:15, ] <- upper_mid_sac_fish + sutter_fish
    
    lower_mid_sac_fish <- migrate(lower_mid_sac_fish + migrants[1:20, ], migratory_survival$lowermid_sac, stochastic = stochastic)
    yolo_fish <- migrate(yolo_fish, migratory_survival$yolo, stochastic = stochastic)
    migrants[1:20, ] <- lower_mid_sac_fish + yolo_fish
    
    lower_sac_fish <- migrate(lower_sac_fish + migrants[1:27, ], migratory_survival$lower_sac, stochastic = stochastic)
    
    san_joaquin_fish <- migrate(migrants[28:31, ] + san_joaquin_fish, migratory_survival$san_joaquin, stochastic = stochastic)
    migrants[28:31, ] <- san_joaquin_fish
    
    delta_fish <- route_and_rear_deltas(year = year, month = month,
                                        migrants = round(migrants),
                                        north_delta_fish = north_delta_fish,
                                        south_delta_fish = south_delta_fish,
                                        north_delta_habitat = habitat$north_delta,
                                        south_delta_habitat = habitat$south_delta,
                                        freeport_flows = ..params$freeport_flows,
                                        cc_gates_days_closed = ..params$cc_gates_days_closed,
                                        rearing_survival_delta = rearing_survival$delta,
                                        migratory_survival_delta = migratory_survival$delta,
                                        migratory_survival_bay_delta = migratory_survival$bay_delta,
                                        juveniles_at_chipps = juveniles_at_chipps,
                                        growth_rates = delta_growth,
                                        territory_size = ..params$territory_size,
                                        stochastic = stochastic)
    
    juveniles_at_chipps <- delta_fish$juveniles_at_chipps
    migrants_at_golden_gate <- delta_fish$migrants_at_golden_gate
  } else {
    
    if (month == 11 & year > 1) {
      # applying summer year to the yearlings and send them out to the ocean
      for (summer_months in 5:10) {
        # we only care for floodplain and inchannel
        yearling_habitat <- get_habitat(year, summer_months,
                                        inchannel_habitat_fry = ..params$inchannel_habitat_fry,
                                        inchannel_habitat_juvenile = ..params$inchannel_habitat_juvenile,
                                        floodplain_habitat = ..params$floodplain_habitat,
                                        sutter_habitat = ..params$sutter_habitat,
                                        yolo_habitat = ..params$yolo_habitat,
                                        delta_habitat = ..params$delta_habitat)
        
        # we only care for floodplain and inchannel
        yearlings_survival_rates <- get_rearing_survival(year, summer_months,
                                                         survival_adjustment = ..params$survival_adjustment,
                                                         mode = mode,
                                                         avg_temp = ..params$avg_temp,
                                                         avg_temp_delta = ..params$avg_temp_delta,
                                                         prob_strand_early = ..params$prob_strand_early,
                                                         prob_strand_late = ..params$prob_strand_late,
                                                         proportion_diverted = ..params$proportion_diverted,
                                                         total_diverted = ..params$total_diverted,
                                                         delta_proportion_diverted = ..params$delta_proportion_diverted,
                                                         delta_total_diverted = ..params$delta_total_diverted,
                                                         weeks_flooded = ..params$weeks_flooded,
                                                         prop_high_predation = ..params$prop_high_predation,
                                                         contact_points = ..params$contact_points,
                                                         delta_contact_points = ..params$delta_contact_points,
                                                         delta_prop_high_predation = ..params$delta_prop_high_predation,
                                                         ..surv_juv_rear_int = ..params$..surv_juv_rear_int,
                                                         .surv_juv_rear_contact_points = ..params$.surv_juv_rear_contact_points,
                                                         ..surv_juv_rear_contact_points = ..params$..surv_juv_rear_contact_points,
                                                         .surv_juv_rear_prop_diversions = ..params$.surv_juv_rear_prop_diversions,
                                                         ..surv_juv_rear_prop_diversions = ..params$..surv_juv_rear_prop_diversions,
                                                         .surv_juv_rear_total_diversions = ..params$.surv_juv_rear_total_diversions,
                                                         ..surv_juv_rear_total_diversions = ..params$..surv_juv_rear_total_diversions,
                                                         ..surv_juv_bypass_int = ..params$..surv_juv_bypass_int,
                                                         ..surv_juv_delta_int = ..params$..surv_juv_delta_int,
                                                         .surv_juv_delta_contact_points = ..params$.surv_juv_delta_contact_points,
                                                         ..surv_juv_delta_contact_points = ..params$..surv_juv_delta_contact_points,
                                                         .surv_juv_delta_total_diverted = ..params$.surv_juv_delta_total_diverted,
                                                         ..surv_juv_delta_total_diverted = ..params$..surv_juv_delta_total_diverted,
                                                         .surv_juv_rear_avg_temp_thresh = ..params$.surv_juv_rear_avg_temp_thresh,
                                                         .surv_juv_rear_high_predation = ..params$.surv_juv_rear_high_predation,
                                                         .surv_juv_rear_stranded = ..params$.surv_juv_rear_stranded,
                                                         .surv_juv_rear_medium = ..params$.surv_juv_rear_medium,
                                                         .surv_juv_rear_large = ..params$.surv_juv_rear_large,
                                                         .surv_juv_rear_floodplain = ..params$.surv_juv_rear_floodplain,
                                                         .surv_juv_bypass_avg_temp_thresh = ..params$.surv_juv_bypass_avg_temp_thresh,
                                                         .surv_juv_bypass_high_predation = ..params$.surv_juv_bypass_high_predation,
                                                         .surv_juv_bypass_medium = ..params$.surv_juv_bypass_medium,
                                                         .surv_juv_bypass_large = ..params$.surv_juv_bypass_large,
                                                         .surv_juv_bypass_floodplain = ..params$.surv_juv_bypass_floodplain,
                                                         .surv_juv_delta_avg_temp_thresh = ..params$.surv_juv_delta_avg_temp_thresh,
                                                         .surv_juv_delta_high_predation = ..params$.surv_juv_delta_high_predation,
                                                         .surv_juv_delta_prop_diverted = ..params$.surv_juv_delta_prop_diverted,
                                                         .surv_juv_delta_medium = ..params$.surv_juv_delta_medium,
                                                         .surv_juv_delta_large = ..params$.surv_juv_delta_large,
                                                         min_survival_rate = ..params$min_survival_rate,
                                                         stochastic = stochastic)
        
        yearlings <- fill_natal(juveniles = yearlings, inchannel_habitat = yearling_habitat$inchannel, 
                                floodplain_habitat = yearling_habitat$floodplain,
                                territory_size = ..params$yearling_territory_size, 
                                up_to_size_class = 4, yearlings = TRUE)
        
        if (summer_months %in% 9:10) {
          growth_ic <- ic_growth
          growth_fp <- fp_growth
        } else {
          growth_ic <- diag(1, 4, 4)
          growth_fp <- replicate(4, diag(1, 4, 4))
        }
        
        yearlings <- rear(juveniles = yearlings$inchannel, survival_rate = yearlings_survival_rates$inchannel, 
                          growth = ic_growth,
                          floodplain_juveniles = yearlings$floodplain,
                          floodplain_survival_rate = yearlings_survival_rates$floodplain, 
                          floodplain_growth = fp_growth,
                          weeks_flooded = ..params$weeks_flooded, 
                          stochastic)
        
        yearlings <- round(yearlings$inchannel + yearlings$floodplain)
        
      }
      
      sutter_detoured <- t(sapply(1:nrow(yearlings[1:15, ]), function(i) {
        if (stochastic) {
          rbinom(n = 4,
                 size = round(yearlings[i, ]),
                 prob = ..params$proportion_flow_bypass[month, year, 1])              
        } else {
          round(yearlings[i, ] * ..params$proportion_flow_bypass[month, year, 1])
        }
        
      }))
      
      yearlings_at_uppermid <- rbind(
        migrate(yearlings[1:15, ] - sutter_detoured, migratory_survival$uppermid_sac,
                stochastic = stochastic),
        matrix(0, ncol = 4, nrow = 2)
      )
      
      yearlings_at_sutter <- rbind(
        migrate(sutter_detoured, migratory_survival$sutter,
                stochastic = stochastic),
        matrix(0, ncol = 4, nrow = 2)
      )
      
      yearlings_at_uppermid <- yearlings_at_sutter + yearlings_at_uppermid
      
      yearlings_at_lowermid <- rbind(yearlings_at_uppermid, yearlings[18:20, ])
      
      yolo_detoured <- t(sapply(1:nrow(yearlings_at_lowermid), function(i) {
        if (stochastic) {
          rbinom(n = 4,
                 size = round(yearlings_at_lowermid[i, ]),
                 prob = ..params$proportion_flow_bypass[month, year, 2])  
        } else {
          round(yearlings_at_lowermid[i, ] * ..params$proportion_flow_bypass[month, year, 2])
        }
        
      }))
      
      yearlings_at_lowersac <- rbind(
        migrate(yearlings_at_lowermid - yolo_detoured, migratory_survival$lowermid_sac,
                stochastic = stochastic),
        matrix(0, ncol = 4, nrow = 3)
      )
      
      yearlings_at_lowersac[23, ] <- yearlings[23, ]
      
      yearlings_at_lowersac <- migrate(yearlings_at_lowersac, migratory_survival$lower_sac,
                                       stochastic = stochastic)
      
      prop_delta_fish_entrained <- route_to_south_delta(freeport_flow = ..params$freeport_flows[[month, year]] * 35.3147,
                                                        dcc_closed = ..params$cc_gates_days_closed[month],
                                                        month = month)
      
      sac_not_entrained <- t(sapply(1:nrow(yearlings_at_lowersac), function(i) {
        if (stochastic) {
          rbinom(n = 4, yearlings_at_lowersac[i, ], prob = 1 - prop_delta_fish_entrained)
        } else {
          yearlings_at_lowersac[i, ] * (1 - prop_delta_fish_entrained)
        }
      }))
      
      yearlings_at_north_delta <- sac_not_entrained +
        rbind(migrate(yolo_detoured, migratory_survival$yolo, stochastic = stochastic),
              matrix(0, ncol = 4, nrow = 3))
      
      yearlings_at_north_delta <- rbind(yearlings_at_north_delta,
                                        matrix(0, ncol = 4, nrow = 8))
      
      yearlings_at_south_delta <- rbind(
        matrix(0, ncol = 4, nrow = 24), # 24 rows for north delta/sac origin fish
        yearlings[25:27, ], # delta tribs
        migrate(yearlings[28:30,], migratory_survival$san_joaquin,
                stochastic = stochastic),
        matrix(0, ncol = 4, nrow = 1) # SJR
      ) +
        rbind(
          yearlings_at_lowersac - sac_not_entrained,
          matrix(0, ncol = 4, nrow = 8)
        )
      
      # estimate fish at Golden Gate Bridge and Chipps Island
      yearling_holding_south_delta <- matrix(0, nrow = 31, ncol = 4, dimnames = list(springRunDSM::watershed_labels, springRunDSM::size_class_labels))
      
      yearling_holding_south_delta[1:24, ] <- t(sapply(1:24, function(i) {
        if (stochastic) {
          rbinom(n = 4, size = round(yearlings_at_south_delta[i, ]), prob = migratory_survival$delta[1, ])              
        } else {
          round(yearlings_at_south_delta[i, ] * migratory_survival$delta[1, ])
        }
        
      }))
      
      yearling_holding_south_delta[26:27, ] <- t(sapply(26:27, function(i) {
        if (stochastic) {
          rbinom(n = 4, size = round(yearlings_at_south_delta[i, ]), prob = migratory_survival$delta[2, ])
        } else {
          round(yearlings_at_south_delta[i, ] * migratory_survival$delta[2, ])
        }
      }))
      
      yearling_holding_south_delta[25, ] <- if (stochastic) {
        rbinom(n = 4,
               yearlings_at_south_delta[25, , drop = F],
               prob = migratory_survival$delta[3, ]) 
      } else {
        round(yearlings_at_south_delta[25, , drop = F] * migratory_survival$delta[3, ])
      }
      
      yearling_holding_south_delta[28:31, ] <- t(sapply(28:31, function(i) {
        if (stochastic) {
          rbinom(n = 4, size = round(yearlings_at_south_delta[i, ]), prob = migratory_survival$delta[4, ])
        } else {
          round(yearlings_at_south_delta[i, ] * migratory_survival$delta[4, ])
        }
      }))
      
      
      survived_yearlings_out <- t(sapply(1:nrow(yearlings_at_north_delta), function(i) {
        if (stochastic) {
          rbinom(n = 4,
                 size = round(yearlings_at_north_delta[i, ]),
                 prob = migratory_survival$bay_delta)  
        } else {
          round(yearlings_at_north_delta[i, ] * migratory_survival$bay_delta)
        }
        
      }))
      
      survived_yearling_holding_south_delta <- t(sapply(1:nrow(yearling_holding_south_delta), function(i) {
        if (stochastic) {fyeralings_
          rbinom(n = 4,
                 size = round(yearling_holding_south_delta[i, ]),
                 prob = migratory_survival$bay_delta)              
        } else {
          round(yearling_holding_south_delta[i, ] * migratory_survival$bay_delta)
        }
        
      }))
      
      yearlings_at_golden_gate <- survived_yearlings_out + survived_yearling_holding_south_delta
      
      # TODO we needed a delta_fish object here. Is this appropriate? 
      delta_fish <- route_and_rear_deltas(year = year, month = month,
                                          migrants = round(migrants),
                                          north_delta_fish = north_delta_fish,
                                          south_delta_fish = south_delta_fish,
                                          north_delta_habitat = habitat$north_delta,
                                          south_delta_habitat = habitat$south_delta,
                                          freeport_flows = ..params$freeport_flows,
                                          cc_gates_days_closed = ..params$cc_gates_days_closed,
                                          rearing_survival_delta = rearing_survival$delta,
                                          migratory_survival_delta = migratory_survival$delta,
                                          migratory_survival_bay_delta = migratory_survival$bay_delta,
                                          juveniles_at_chipps = juveniles_at_chipps,
                                          growth_rates = delta_growth,
                                          territory_size = ..params$territory_size,
                                          stochastic = stochastic)
      
      juveniles_at_chipps <- juveniles_at_chipps + yearlings_at_north_delta + yearling_holding_south_delta
      juveniles_at_chipps <- delta_fish$juveniles_at_chipps
      adults_in_ocean <- adults_in_ocean + ocean_entry_success(migrants = yearlings_at_golden_gate,
                                                               month = 11,
                                                               avg_ocean_transition_month = avg_ocean_transition_month,
                                                               stochastic = stochastic)
      
      yearlings <- matrix(0, ncol = 4, nrow = 31, dimnames = list(springRunDSM::watershed_labels, springRunDSM::size_class_labels))
    }
    # if month < 8
    # route northern natal fish stay and rear or migrate downstream ------
    upper_sac_trib_fish <-  route(year = year,
                                  month = month,
                                  juveniles = juveniles[1:15, ],
                                  inchannel_habitat = habitat$inchannel[1:15],
                                  floodplain_habitat = habitat$floodplain[1:15],
                                  prop_pulse_flows = ..params$prop_pulse_flows[1:15, ],
                                  .pulse_movement_intercept = ..params$.pulse_movement_intercept,
                                  .pulse_movement_proportion_pulse = ..params$.pulse_movement_proportion_pulse,
                                  .pulse_movement_medium = ..params$.pulse_movement_medium,
                                  .pulse_movement_large = ..params$.pulse_movement_large,
                                  .pulse_movement_vlarge = ..params$.pulse_movement_vlarge,
                                  .pulse_movement_medium_pulse = ..params$.pulse_movement_medium_pulse,
                                  .pulse_movement_large_pulse = ..params$.pulse_movement_large_pulse,
                                  .pulse_movement_very_large_pulse = ..params$.pulse_movement_very_large_pulse,
                                  territory_size = ..params$territory_size,
                                  stochastic = stochastic)
    
    upper_sac_trib_rear <- rear(juveniles = upper_sac_trib_fish$inchannel,
                                survival_rate = rearing_survival$inchannel[1:15, ],
                                growth = ic_growth[,,1:15],
                                floodplain_juveniles = upper_sac_trib_fish$floodplain,
                                floodplain_survival_rate = rearing_survival$floodplain[1:15, ],
                                floodplain_growth = fp_growth[,,1:15],
                                weeks_flooded = ..params$weeks_flooded[1:15, month, year], 
                                stochastic = stochastic)
    
    juveniles[1:15, ] <- upper_sac_trib_rear$inchannel + upper_sac_trib_rear$floodplain
    
    # route migrant fish into Upper-mid Sac Region (fish from watersheds 1:15)
    # regional fish stay and rear
    # or migrate further downstream or in sutter bypass
    
    upper_mid_sac_fish <- route_regional(month = month,
                                         year = year,
                                         migrants = upper_mid_sac_fish + upper_sac_trib_fish$migrants,
                                         inchannel_habitat = habitat$inchannel[16],
                                         floodplain_habitat = habitat$floodplain[16],
                                         prop_pulse_flows = ..params$prop_pulse_flows[16, , drop = FALSE],
                                         migration_survival_rate = migratory_survival$uppermid_sac,
                                         proportion_flow_bypass = ..params$proportion_flow_bypass,
                                         detour = 'sutter',
                                         territory_size = ..params$territory_size,
                                         stochastic = stochastic)
    
    
    sutter_fish <- route_bypass(bypass_fish = sutter_fish + upper_mid_sac_fish$detoured,
                                bypass_habitat = habitat$sutter,
                                migration_survival_rate = migratory_survival$sutter,
                                territory_size = ..params$territory_size,
                                stochastic = stochastic)
    
    migrants[1:15, ] <- upper_mid_sac_fish$migrants + sutter_fish$migrants
    
    upper_mid_sac_fish <- rear(juveniles = upper_mid_sac_fish$inchannel,
                               survival_rate = rearing_survival$inchannel[16, ],
                               growth = ic_growth[,,16],
                               floodplain_juveniles = upper_mid_sac_fish$floodplain,
                               floodplain_survival_rate = rearing_survival$floodplain[16, ],
                               floodplain_growth = fp_growth[,,16],
                               weeks_flooded = rep(..params$weeks_flooded[16, month, year], nrow(upper_mid_sac_fish$inchannel)),
                               stochastic = stochastic)
    
    upper_mid_sac_fish <- upper_mid_sac_fish$inchannel + upper_mid_sac_fish$floodplain
    
    sutter_fish <- rear(juveniles = sutter_fish$inchannel,
                        survival_rate = rearing_survival$sutter[1,],
                        growth = ic_growth[,,17],
                        stochastic = stochastic)
    
    
    # route migrant fish into Lower-mid Sac Region (fish from watersheds 18:20, and migrants from Upper-mid Sac Region)
    # regional fish stay and rear
    # or migrate further downstream  or in yolo bypass
    lower_mid_sac_trib_fish <- route(year = year,
                                     month = month,
                                     juveniles = juveniles[18:20, ],
                                     inchannel_habitat = habitat$inchannel[18:20],
                                     floodplain_habitat = habitat$floodplain[18:20],
                                     prop_pulse_flows =  ..params$prop_pulse_flows[18:20, ],
                                     .pulse_movement_intercept = ..params$.pulse_movement_intercept,
                                     .pulse_movement_proportion_pulse = ..params$.pulse_movement_proportion_pulse,
                                     .pulse_movement_medium = ..params$.pulse_movement_medium,
                                     .pulse_movement_large = ..params$.pulse_movement_large,
                                     .pulse_movement_vlarge = ..params$.pulse_movement_vlarge,
                                     .pulse_movement_medium_pulse = ..params$.pulse_movement_medium_pulse,
                                     .pulse_movement_large_pulse = ..params$.pulse_movement_large_pulse,
                                     .pulse_movement_very_large_pulse = ..params$.pulse_movement_very_large_pulse,
                                     territory_size = ..params$territory_size,
                                     stochastic = stochastic)
    
    lower_mid_sac_trib_rear <- rear(juveniles = lower_mid_sac_trib_fish$inchannel,
                                    survival_rate = rearing_survival$inchannel[18:20, ],
                                    growth = ic_growth[,,18:20],
                                    floodplain_juveniles = lower_mid_sac_trib_fish$floodplain,
                                    floodplain_survival_rate = rearing_survival$floodplain[18:20, ],
                                    floodplain_growth = fp_growth[,,18:20],
                                    weeks_flooded = ..params$weeks_flooded[18:20, month, year],
                                    stochastic = stochastic)
    
    juveniles[18:20, ] <- lower_mid_sac_trib_rear$inchannel + lower_mid_sac_trib_rear$floodplain
    migrants[18:20, ] <- lower_mid_sac_trib_fish$migrants
    
    lower_mid_sac_fish <- route_regional(month = month,
                                         year = year,
                                         migrants = lower_mid_sac_fish + migrants[1:20, ],
                                         inchannel_habitat = habitat$inchannel[21],
                                         floodplain_habitat = habitat$floodplain[21],
                                         prop_pulse_flows = ..params$prop_pulse_flows[21, , drop = FALSE],
                                         migration_survival_rate = migratory_survival$lowermid_sac,
                                         proportion_flow_bypass = ..params$proportion_flow_bypass,
                                         detour = 'yolo',
                                         territory_size = ..params$territory_size,
                                         stochastic = stochastic)
    
    yolo_fish <- route_bypass(bypass_fish = yolo_fish + lower_mid_sac_fish$detoured,
                              bypass_habitat = habitat$yolo,
                              migration_survival_rate = migratory_survival$yolo,
                              territory_size = ..params$territory_size,
                              stochastic = stochastic)
    
    migrants[1:20, ] <- lower_mid_sac_fish$migrants + yolo_fish$migrants
    
    lower_mid_sac_fish <- rear(juveniles = lower_mid_sac_fish$inchannel,
                               survival_rate = rearing_survival$inchannel[21, ],
                               growth = ic_growth[,,21],
                               floodplain_juveniles = lower_mid_sac_fish$floodplain,
                               floodplain_survival_rate = rearing_survival$floodplain[21, ],
                               floodplain_growth = fp_growth[,,21],
                               weeks_flooded = rep(..params$weeks_flooded[21, month, year], nrow(lower_mid_sac_fish$inchannel)),
                               stochastic = stochastic)
    
    lower_mid_sac_fish <- lower_mid_sac_fish$inchannel + lower_mid_sac_fish$floodplain
    
    yolo_fish <- rear(juveniles = yolo_fish$inchannel,
                      survival_rate = rearing_survival$yolo[1,],
                      growth = ic_growth[,,22],
                      stochastic = stochastic)
    
    
    # route migrant fish into Lower Sac Region (fish from watershed 23, and migrants from Lower-mid Sac Region)
    # regional fish stay and rear
    # or migrate north delta
    lower_sac_trib_fish <- route(year = year,
                                 month = month,
                                 juveniles = juveniles[23, , drop = FALSE],
                                 inchannel_habitat = habitat$inchannel[23],
                                 floodplain_habitat = habitat$floodplain[23],
                                 prop_pulse_flows =  ..params$prop_pulse_flows[23, , drop = FALSE],
                                 .pulse_movement_intercept = ..params$.pulse_movement_intercept,
                                 .pulse_movement_proportion_pulse = ..params$.pulse_movement_proportion_pulse,
                                 .pulse_movement_medium = ..params$.pulse_movement_medium,
                                 .pulse_movement_large = ..params$.pulse_movement_large,
                                 .pulse_movement_vlarge = ..params$.pulse_movement_vlarge,
                                 .pulse_movement_medium_pulse = ..params$.pulse_movement_medium_pulse,
                                 .pulse_movement_large_pulse = ..params$.pulse_movement_large_pulse,
                                 .pulse_movement_very_large_pulse = ..params$.pulse_movement_very_large_pulse,
                                 territory_size = ..params$territory_size,
                                 stochastic = stochastic)
    
    lower_sac_trib_rear <- rear(juveniles = lower_sac_trib_fish$inchannel,
                                survival_rate = rearing_survival$inchannel[23, , drop = FALSE],
                                growth = ic_growth[,,23],
                                floodplain_juveniles = lower_sac_trib_fish$floodplain,
                                floodplain_survival_rate = rearing_survival$floodplain[23, , drop = FALSE],
                                floodplain_growth = fp_growth[,,23],
                                weeks_flooded = ..params$weeks_flooded[23, month, year],
                                stochastic = stochastic)
    
    juveniles[23, ] <- lower_sac_trib_rear$inchannel + lower_sac_trib_rear$floodplain
    
    migrants[23, ] <- lower_sac_trib_fish$migrants
    
    lower_sac_fish <- route_regional(month = month,
                                     year = year,
                                     migrants = lower_sac_fish + migrants[1:27, ],
                                     inchannel_habitat = habitat$inchannel[24],
                                     floodplain_habitat = habitat$floodplain[24],
                                     prop_pulse_flows = ..params$prop_pulse_flows[24, , drop = FALSE],
                                     migration_survival_rate = migratory_survival$lower_sac,
                                     territory_size = ..params$territory_size,
                                     stochastic = stochastic)
    
    migrants[1:27, ] <- lower_sac_fish$migrants
    
    lower_sac_fish <- rear(juveniles = lower_sac_fish$inchannel,
                           survival_rate = rearing_survival$inchannel[24, ],
                           growth = ic_growth[,,24],
                           floodplain_juveniles = lower_sac_fish$floodplain,
                           floodplain_survival_rate = rearing_survival$floodplain[24, ],
                           floodplain_growth = fp_growth[,,24],
                           weeks_flooded = rep(..params$weeks_flooded[24, month, year], nrow(lower_sac_fish$inchannel)),
                           stochastic = stochastic)
    
    lower_sac_fish <- lower_sac_fish$inchannel + lower_sac_fish$floodplain
    
    # route southern natal fish stay and rear or migrate downstream ------
    
    # route migrant fish into South Delta Region (fish from watersheds 25:27)
    # regional fish stay and rear
    # or migrate to south delta
    south_delta_trib_fish <- route(year = year,
                                   month = month,
                                   juveniles = juveniles[25:27, ],
                                   inchannel_habitat = habitat$inchannel[25:27],
                                   floodplain_habitat = habitat$floodplain[25:27],
                                   prop_pulse_flows =  ..params$prop_pulse_flows[25:27, ],
                                   .pulse_movement_intercept = ..params$.pulse_movement_intercept,
                                   .pulse_movement_proportion_pulse = ..params$.pulse_movement_proportion_pulse,
                                   .pulse_movement_medium = ..params$.pulse_movement_medium,
                                   .pulse_movement_large = ..params$.pulse_movement_large,
                                   .pulse_movement_vlarge = ..params$.pulse_movement_vlarge,
                                   .pulse_movement_medium_pulse = ..params$.pulse_movement_medium_pulse,
                                   .pulse_movement_large_pulse = ..params$.pulse_movement_large_pulse,
                                   .pulse_movement_very_large_pulse = ..params$.pulse_movement_very_large_pulse,
                                   territory_size = ..params$territory_size,
                                   stochastic = stochastic)
    
    south_delta_trib_rear <- rear(juveniles = south_delta_trib_fish$inchannel,
                                  survival_rate = rearing_survival$inchannel[25:27, ],
                                  growth = ic_growth[,,25:27],
                                  floodplain_juveniles = south_delta_trib_fish$floodplain,
                                  floodplain_survival_rate = rearing_survival$floodplain[25:27, ],
                                  floodplain_growth = fp_growth[,,25:27],
                                  weeks_flooded = ..params$weeks_flooded[25:27, month, year],
                                  stochastic = stochastic)
    
    juveniles[25:27, ] <- south_delta_trib_rear$inchannel + south_delta_trib_rear$floodplain
    
    migrants[25:27, ] <- south_delta_trib_fish$migrants
    
    # route migrant fish into San Joquin River (fish from watersheds 28:30)
    # regional fish stay and rear
    # or migrate to south delta
    
    san_joaquin_trib_fish <- route(year = year,
                                   month = month,
                                   juveniles = juveniles[28:31, ],
                                   inchannel_habitat = habitat$inchannel[28:31],
                                   floodplain_habitat = habitat$floodplain[28:31],
                                   prop_pulse_flows =  ..params$prop_pulse_flows[28:31, ],
                                   .pulse_movement_intercept = ..params$.pulse_movement_intercept,
                                   .pulse_movement_proportion_pulse = ..params$.pulse_movement_proportion_pulse,
                                   .pulse_movement_medium = ..params$.pulse_movement_medium,
                                   .pulse_movement_large = ..params$.pulse_movement_large,
                                   .pulse_movement_vlarge = ..params$.pulse_movement_vlarge,
                                   .pulse_movement_medium_pulse = ..params$.pulse_movement_medium_pulse,
                                   .pulse_movement_large_pulse = ..params$.pulse_movement_large_pulse,
                                   .pulse_movement_very_large_pulse = ..params$.pulse_movement_very_large_pulse,
                                   territory_size = ..params$territory_size,
                                   stochastic = stochastic)
    
    san_joaquin_trib_rear <- rear(juveniles = san_joaquin_trib_fish$inchannel,
                                  survival_rate = rearing_survival$inchannel[28:31, ],
                                  growth = ic_growth[,,28:31],
                                  floodplain_juveniles = san_joaquin_trib_fish$floodplain,
                                  floodplain_survival_rate = rearing_survival$floodplain[28:31, ],
                                  floodplain_growth = fp_growth[,,28:31],
                                  weeks_flooded = ..params$weeks_flooded[28:31, month, year],
                                  stochastic = stochastic)
    
    juveniles[28:31, ] <- san_joaquin_trib_rear$inchannel + san_joaquin_trib_rear$floodplain
    
    san_joaquin_fish <- route_regional(month = month,
                                       year = year,
                                       migrants = san_joaquin_fish + san_joaquin_trib_fish$migrants,
                                       inchannel_habitat = habitat$inchannel[31],
                                       floodplain_habitat = habitat$floodplain[31],
                                       prop_pulse_flows = ..params$prop_pulse_flows[31, , drop = FALSE],
                                       migration_survival_rate = migratory_survival$san_joaquin,
                                       territory_size = ..params$territory_size,
                                       stochastic = stochastic)
    
    migrants[28:31, ] <- san_joaquin_fish$migrants
    
    san_joaquin_fish <- rear(juveniles = san_joaquin_fish$inchannel,
                             survival_rate = rearing_survival$inchannel[31, ],
                             growth = ic_growth[,,31],
                             floodplain_juveniles = san_joaquin_fish$floodplain,
                             floodplain_survival_rate = rearing_survival$floodplain[31, ],
                             floodplain_growth = fp_growth[,,31],
                             weeks_flooded = rep(..params$weeks_flooded[31, month, year], nrow(san_joaquin_fish$inchannel)),
                             stochastic = stochastic)
    
    san_joaquin_fish <- san_joaquin_fish$inchannel + san_joaquin_fish$floodplain
    
    delta_fish <- route_and_rear_deltas(year = year, month = month,
                                        migrants = round(migrants),
                                        north_delta_fish = north_delta_fish,
                                        south_delta_fish = south_delta_fish,
                                        north_delta_habitat = habitat$north_delta,
                                        south_delta_habitat = habitat$south_delta,
                                        freeport_flows = ..params$freeport_flows,
                                        cc_gates_days_closed = ..params$cc_gates_days_closed,
                                        rearing_survival_delta = rearing_survival$delta,
                                        migratory_survival_delta = migratory_survival$delta,
                                        migratory_survival_bay_delta = migratory_survival$bay_delta,
                                        juveniles_at_chipps = juveniles_at_chipps,
                                        growth_rates = delta_growth,
                                        territory_size = ..params$territory_size,
                                        stochastic = stochastic)
    
    migrants_at_golden_gate <- delta_fish$migrants_at_golden_gate
    north_delta_fish <- delta_fish$north_delta_fish
    south_delta_fish <- delta_fish$south_delta_fish
    juveniles_at_chipps <- delta_fish$juveniles_at_chipps
  }
  
  adults_in_ocean <- adults_in_ocean + ocean_entry_success(migrants = migrants_at_golden_gate,
                                                           month = month,
                                                           avg_ocean_transition_month = avg_ocean_transition_month,
                                                           .ocean_entry_success_length = ..params$.ocean_entry_success_length,
                                                           ..ocean_entry_success_int = ..params$..ocean_entry_success_int,
                                                           .ocean_entry_success_months = ..params$.ocean_entry_success_months,
                                                           stochastic = stochastic)
  return(list(juveniles = juveniles,
              lower_mid_sac_fish = lower_mid_sac_fish,
              lower_sac_fish = lower_sac_fish,
              upper_mid_sac_fish = upper_mid_sac_fish,
              sutter_fish = sutter_fish,
              yolo_fish = yolo_fish,
              san_joaquin_fish = san_joaquin_fish,
              north_delta_fish = north_delta_fish,
              south_delta_fish = south_delta_fish,
              juveniles_at_chipps = juveniles_at_chipps,
              adults_in_ocean = adults_in_ocean,
              migrants_at_golden_gate = migrants_at_golden_gate)
  )
  
}

convert_number_to_word <- function(number) {
  number_to_word <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
  word_form <- number_to_word[number]
  return(word_form)
}

create_fish_df <- function(fish_df, month, year) {
  
  hypothesis <- convert_number_to_word(fish_df$hypothesis)
  
  tmp <- data.frame("s" = rep(0, 31),
                    "m" = rep(0, 31),
                    "l" = rep(0, 31),
                    'vl' = rep(0, 31))
  
  north_delta_fish <-  data.frame(fish_df$north_delta_fish)
  south_delta_fish <- data.frame(fish_df$south_delta_fish)
  
  tmp[1:nrow(north_delta_fish), ] <- north_delta_fish
  fish_df <- data.frame(tmp + south_delta_fish) |>
    dplyr::mutate(watershed = fallRunDSM::watershed_labels[1:31],
                  month = month,
                  year = year,
                  hypothesis = hypothesis)
  rownames(fish_df) <- NULL
  
  return(fish_df)
  
}
