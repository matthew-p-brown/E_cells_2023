# Function for analyzing activity and sleep between ZT0-6, ZT6-12, and ZT12-24

dusk_behavior_suite <- function(dt){
  
  # Add phase information to dt.
  # Separately labels morning (ZT0-6), afternoon (ZT6-9), evening (ZT9-12), twilight (ZT12-15), and night (ZT15-24) for baseline/LD and activation/DD days.
  dt <- dt[, phase := ifelse(t %between% c(hours(0), hours(6)), "base morn", 
                             ifelse(t %between% c(hours(6), hours(9)), "base after", 
                                    ifelse(t %between% c(hours(9), hours(12)), "base eve", 
                                           ifelse(t %between% c(hours(12), hours(15)), "base twilight", 
                                                  ifelse(t %between% c(hours(15), hours(24)), "base night", 
                                                         ifelse(t %between% c(hours(24), hours(30)), "activate morn", 
                                                                ifelse(t %between% c(hours(30), hours(33)), "activate after", 
                                                                       ifelse(t %between% c(hours(33), hours(36)), "activate eve", 
                                                                              ifelse(t %between% c(hours(36), hours(39)), "activate twilight",  "activate night")))))))))]
  
  # Summarize activity and sleep per phase per fly.
  # Adds up activity and sleep amount time (in mins) for phases listed above.
  # Afternoon below includes both afternoon and evening phases from above (ZT6-12).
  dt_summ <- rejoin(dt[,
                       .(
                         # this is where the computation happens
                         sleep_sum_base_morn = sum(asleep[phase == "base morn"]),
                         sleep_sum_base_after = sum(asleep[phase %in% c("base after", "base eve")]),
                         sleep_sum_base_eve = sum(asleep[phase == "base eve"]), 
                         sleep_sum_base_twi = sum(asleep[phase == "base twilight"]), 
                         sleep_sum_base_night = sum(asleep[phase %in% c("base twilight", "base night")]), 
                         sleep_sum_activ_morn = sum(asleep[phase == "activate morn"]), 
                         sleep_sum_activ_after = sum(asleep[phase %in% c("activate after", "activate eve")]), 
                         sleep_sum_activ_eve = sum(asleep[phase == "activate eve"]), 
                         sleep_sum_activ_twi = sum(asleep[phase == "activate twilight"]), 
                         sleep_sum_activ_night = sum(asleep[phase %in% c("activate twilight", "activate night")]), 
                         activity_sum_base_morn = sum(activity[phase == "base morn"]),
                         activity_sum_base_after = sum(activity[phase %in% c("base after", "base eve")]),
                         activity_sum_base_eve = sum(activity[phase == "base eve"]), 
                         activity_sum_base_twi = sum(activity[phase == "base twilight"]), 
                         activity_sum_base_night = sum(activity[phase %in% c("base twilight", "base night")]), 
                         activity_sum_activ_morn = sum(activity[phase == "activate morn"]), 
                         activity_sum_activ_after = sum(activity[phase %in% c("activate after", "activate eve")]), 
                         activity_sum_activ_eve = sum(activity[phase == "activate eve"]), 
                         activity_sum_activ_twi = sum(activity[phase == "activate twilight"]), 
                         activity_sum_activ_night = sum(activity[phase %in% c("activate twilight", "activate night")]), 
                         activity_slope_base_after = lm(activity[phase %in% c("base after", "base eve")] ~ t[phase %in% c("base after", "base eve")])$coefficients[2] * 1800, 
                         activity_slope_activ_after = lm(activity[phase %in% c("activate after", "activate eve")] ~ t[phase %in% c("activate after", "activate eve")])$coefficients[2] * 1800, 
                         activity_intercept_base_after = (-lm(activity[phase %in% c("base after", "base eve")] ~ t[phase %in% c("base after", "base eve")])$coefficients[1] / lm(activity[phase %in% c("base after", "base eve")] ~ t[phase %in% c("base after", "base eve")])$coefficients[2]) / 3600, 
                         activity_intercept_activ_after = ((-lm(activity[phase %in% c("activate after", "activate eve")] ~ t[phase %in% c("activate after", "activate eve")])$coefficients[1] / lm(activity[phase %in% c("activate after", "activate eve")] ~ t[phase %in% c("activate after", "activate eve")])$coefficients[2]) / 3600) - 24
                         ),
                       ,by=id])
  
  
  # Trim original data table to only include the nighttime of both nights for bout analysis
  dt_trimmed <- dt[phase %in% c("base twilight", "base night", "activate twilight", "activate night")]
  
  # Perform bout analysis from Rethomics
  dt_bouts <- bout_analysis(asleep, dt_trimmed)
  
  # Remove redundant columns
  dt_bouts <- dt_bouts[asleep == TRUE, -"asleep"]
  
  # Summarize bouts for the baseline/LD day
  dt_bouts_baseline_day <- dt_bouts[t %between%  c(hours(12), days(1))]
  baseline_bout_summary <- dt_bouts_baseline_day[,.(
    base_latency = t[1] / 3600, # the first bout is at t[1] 
    base_first_bout_length = duration[1] / 60,
    base_latency_to_longest_bout = t[which.max(duration)] / 3600,
    base_length_longest_bout = max(duration) / 60,
    base_n_bouts = .N,
    base_mean_bout_length = mean(duration) / 60
  ),
  by=id]
  
  # Summarize bouts for the activation/DD day
  dt_bouts_second_day <- dt_bouts[t %between%  c(hours(36), days(2))]
  second_day_bout_summary <- dt_bouts_second_day[,.(
    second_latency = (t[1] - hours(24)) / 3600, # the first bout is at t[1] 
    second_first_bout_length = duration[1] / 60,
    second_latency_to_longest_bout = (t[which.max(duration)] - hours(24)) / 3600,
    second_length_longest_bout = max(duration) / 60,
    second_n_bouts = .N,
    second_mean_bout_length = mean(duration) / 60
  ),
  by=id]
  
  dt_output <- merge(dt_summ, baseline_bout_summary, by = "id")
  dt_output <- merge(dt_output, second_day_bout_summary, by = "id")
  dt_output
}