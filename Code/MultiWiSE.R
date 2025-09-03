# MultiWiSE.R
# This code is associated with the manuscript 'Multiyear Wildfire Smoke Exposure (MultiWiSE) metrics: 
# A data-driven approach to characterizing episodic PM2.5 exposures for epidemiologic research on chronic health effects'
# The script allows for replication of the results included in the manuscript. It includes code for the data-driven approach  
# used to process the 2010-2023 BC dissemination area-level PM2.5 data and generate the 12 MultiWiSE Metrics for all census 
# subdivisions in the province. It also includes code to generate the figures and tables included in the manuscript. 

####################################################################################################################################
############################### PACKAGES, WORKING DIRECTORY, & SHAPE FILES #########################################################
####################################################################################################################################

#### Install/Load packages ####
packages <- c(
  "cancensus", "dplyr", "sf", "data.table",
  "lubridate", "ggplot2", "tidyr", "ggtext", 
  "purrr", "readr", "stringr", "gt",
  "patchwork", "ggcorrplot",
  "ggspatial", "gridExtra", "ggrepel",
  "cowplot"
)

installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install)) {
  install.packages(to_install)
}

lapply(packages, library, character.only = TRUE)
rm(packages, installed, to_install)

#### Set a working directory ####
setwd("<working_directory>")

#### Get DA and CSD shapefiles from cancensus ####
da_census <- get_census(
  dataset    = "CA16",
  regions    = list(PR = "59"),
  level      = "DA",
  geo_format = "sf"
) %>% select(GeoUID, CSD_UID, Population)

csd_census <- get_census(
  dataset    = "CA16",
  regions    = list(PR = "59"),
  level      = "CSD",
  geo_format = "sf"
) %>% mutate(GeoUID = as.character(GeoUID))

####################################################################################################################################
###################################### READ IN, COMBINE, & CLEAN DA-LEVEL DATA #####################################################
####################################################################################################################################

#### Read in and combine DA-level PM2.5 data ####
# Set up input directory
input_dir <- file.path(getwd(), "Data")

# Get list of file names - one file per year
da_files  <- list.files(input_dir, pattern = "DA_CanOSSEM_\\d{4}_WFS.*\\.RData$", 
                        full.names = TRUE)

all_da_list <- vector("list", length(da_files))
names(all_da_list) <- basename(da_files)

# Loop through each year-file and store it in a list
for (file in da_files) {
  year_extracted <- sub(".*DA_CanOSSEM_(\\d{4})_WFS.*\\.RData$", "\\1", basename(file))
  load(file)
  pm25_da <- DA_CanOSSEM_WFS
  pm25_da <- merge(pm25_da, st_drop_geometry(da_census),
                   by.x = "DAUID", by.y = "GeoUID")
  pm25_da <- pm25_da %>% filter(!is.na(Population) & Population > 0)
  pm25_da$year <- as.integer(year_extracted)
  all_da_list[[basename(file)]] <- pm25_da
}

# Combine list into a data.table
all_da <- rbindlist(all_da_list, use.names = TRUE, fill = TRUE)

#### Convert date variable into continuous epiweeks ####

all_da <- all_da %>% mutate(DATE = as.Date(DATE))

# Create a complete sequence of dates from start to end of study period
dates <- seq(
  from = as.Date("2010-01-01"),    
  to   = as.Date("2023-12-31"),      
  by   = "days"                     
)

# Build a helper data.frame mapping each date to its epiweek number
weeks_df <- data.frame(
  date = dates,
  epiweek_num = epiweek(as.POSIXct(dates))
)

# Compute a continuous epiweek variable
weeks_df$cont_epiweek_num <- rleid(weeks_df$epiweek_num) - 1

# Remove incomplete weeks at beginning and end of study period
weeks_df <- weeks_df %>% 
  group_by(cont_epiweek_num) %>% 
  mutate(n = n()) %>% 
  ungroup() %>%
  mutate(cont_epiweek_num = ifelse(n < 7, NA, cont_epiweek_num))

# Join the continuous epiweek back onto the main data
all_da$epiweek <- weeks_df$cont_epiweek_num[
  match(all_da$DATE, weeks_df$date) 
]

# Clean up: drop rows where epiweek is NA (incomplete weeks at start/end)
all_da <- all_da %>%
  filter(!is.na(epiweek)) %>%    
  mutate(epiweek = as.integer(epiweek)) # convert epiweek to integer

# Save combined DA PM2.5 data as an .RDS file
saveRDS(all_da, file = "./Processed_Data/combined_DA_daily_PM25_2010-2023.RDS")

####################################################################################################################################
################################### AGGREGATE DATA TO CSD AND WEEKLY LEVEL #########################################################
####################################################################################################################################

#### Aggregate daily DA-level data to weekly CSD-level data ####

# Read in combined DA data
all_da <- readRDS("./Processed_Data/combined_DA_daily_PM25_2010-2023.RDS")

# Calculate population‑weighted daily PM2.5 for each census subdivision
csd_daily <- all_da %>%
  group_by(CSD_UID, DATE, epiweek) %>%
  summarize(
    popw_pm25 = sum(DA_PM25 * Population, na.rm = TRUE) /
      sum(Population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(csd = as.character(CSD_UID))

# Aggregate to weekly level
csd_weekly <- csd_daily %>%
  group_by(csd, epiweek) %>%
  summarize(
    weekly_n_days    = n(),                 # How many days had data for this week
    weekly_avg_pm25  = mean(popw_pm25),     # Average of daily PM2.5 values
    weekly_sum_pm25  = sum(popw_pm25),      # Total sum of daily PM2.5 values
    year_start       = min(year(DATE)),     # Year of the first day in the epiweek
    start_date       = min(DATE),           # First calendar date in the epiweek
    end_date         = max(DATE),           # Last calendar date in the epiweek
    .groups          = "drop"              
  )

saveRDS(csd_weekly, "./Processed_Data/CSD_weekly_PM25_2010-2023.RDS")

####################################################################################################################################
############################ DEFINE FUNCTIONS FOR PROCESSING DATA & GENERATING METRICS #############################################
####################################################################################################################################

# Function to calculate modified Z-Score
calc_modified_z <- function(values) {
  med <- median(values, na.rm = TRUE)
  mad_val <- mad(values, constant = 1.4826, na.rm = TRUE)
  if (mad_val == 0) return(rep(NA_real_, length(values)))
  (values - med) / mad_val
}

# Function to identify smoke episodes
identify_smoke_episodes <- function(
    is_smoke_impacted,                 # Vector identifying if an epiweek is WFS-impacted or not (1 or 0) 
    cont_week,                         # Vector of continuous epiweeks
    weekly_sum_wfs_pm25_values,        # Vector of the sum of WFS PM2.5 for each epiweek
    episode_threshold,                 # Value that must be exceeded to count as an episode (0 for WFS episode, 250 for severe episode)
    num_weeks_in_episode,              # Number of WFS-impacted weeks that need to be present to count as an episode (2 for WFS episode, 1 for severe episode)
    max_no_smoke_gap,                  # Maximum number of continuous non-WFS-impacted weeks that can be present and still count as an episode (3 weeks)
    weeks_needed_above_threshold       # Number of weeks above the threshold needed to count as an episode (1 week)
) {
  # Total number of epiweeks
  n <- length(is_smoke_impacted)
  
  # Initialize output vector (episode IDs), current episode counter, and flag
  epi_id <- integer(n)
  cur_id <- 0
  in_epi <- FALSE
  
  # Trackers within each candidate episode
  smoke_wks <- 0      # total smoke-impacted weeks
  wks_above <- 0      # weeks where PM2.5 > threshold
  gap       <- 0      # current non-smoke gap length
  bridge    <- integer(0)  # indices of non-smoke weeks within the gap
  
  # Helper to abandon a failed episode (reset its epi_id entries to 0)
  fail <- function(id) epi_id[epi_id == id] <<- 0
  
  # Iterate week by week
  for (i in seq_len(n)) {
    
    # Break episode if non-consecutive week index
    if (i > 1 && cont_week[i] - cont_week[i - 1] > 1) {
      # If in an episode that doesn’t meet criteria, drop it
      if (in_epi &&
          (smoke_wks < num_weeks_in_episode ||
           wks_above < weeks_needed_above_threshold)) {
        fail(cur_id)
      }
      # Reset all episode trackers
      in_epi    <- FALSE
      smoke_wks <- 0
      wks_above <- 0
      gap       <- 0
      bridge    <- integer(0)
    }
    
    # Handle a smoke-impacted week
    if (is_smoke_impacted[i] == 1) {
      
      # Start a new episode if not already in one
      if (!in_epi) {
        cur_id <- cur_id + 1
        in_epi <- TRUE
        smoke_wks <- 0
        wks_above <- 0
        gap <- 0
      }
      
      # Any bridged non-smoke weeks now belong to this episode
      if (length(bridge) > 0) {
        epi_id[bridge] <- cur_id
        bridge <- integer(0)
      }
      
      # Update counters
      smoke_wks <- smoke_wks + 1
      if (weekly_sum_wfs_pm25_values[i] > episode_threshold) {
        wks_above <- wks_above + 1
      }
      
      epi_id[i] <- cur_id
      gap <- 0  # reset gap
      
    # Handle a non-smoke week within an ongoing episode
    } else if (in_epi) {
      gap <- gap + 1
      bridge <- c(bridge, i)  # tentatively include this week
      
      # If gap is too long, close out the episode
      if (gap > max_no_smoke_gap) {
        
        # Abandon if criteria not met
        if (smoke_wks < num_weeks_in_episode ||
            wks_above < weeks_needed_above_threshold) {
          fail(cur_id)
        }
        in_epi <- FALSE
        bridge <- integer(0)
      }
    }
  }
  
  # Final check at end of series
  if (in_epi &&
      (smoke_wks < num_weeks_in_episode ||
       wks_above < weeks_needed_above_threshold)) {
    fail(cur_id)
  }
  
  # Trim trailing non-smoke weeks after last smoke week in each episode
  for (id in setdiff(unique(epi_id), 0)) {
    idx <- which(epi_id == id)
    smoke_idx <- idx[is_smoke_impacted[idx] == 1]
    epi_id[idx[idx > max(smoke_idx)]] <- 0
  }
  
  # Reset episode IDs to 1,2,3, etc. in order of discovery
  uniq <- setdiff(unique(epi_id), 0)
  for (k in seq_along(uniq)) {
    epi_id[epi_id == uniq[k]] <- k
  }
  epi_id
}

# Function to number metrics from #1-12
metric_numbering <- c(
  "Cumulative WFS PM2.5 (mg/m³)"           = "1. Cumulative WFS PM₂.₅ (mg/m³)",
  "WFS Fraction (%)"                       = "2. WFS Fraction (%)",
  "Average WFS PM2.5 (µg/m³)"              = "3. Average WFS PM₂.₅ (µg/m³)",
  "Any WFS (# weeks)"                      = "4. Any WFS (# weeks)",
  "WFS PM2.5 > 5 µg/m³ (# weeks)"          = "5. WFS PM₂.₅ > 5 µg/m³ (# weeks)",
  "Total PM2.5 > 25 µg/m³ (# weeks)"       = "6. Total PM₂.₅ > 25 µg/m³ (# weeks)",
  "WFS Episodes (# episodes)"              = "7. WFS Episodes (# episodes)",
  "Severe Episodes (# episodes)"           = "8. Severe Episodes (# episodes)",
  "Longest Episode (# weeks)"              = "9. Longest Episode (# weeks)",
  "Worst Episode (µg/m³)"                  = "10. Worst Episode (µg/m³)",
  "WFS from Severe Episodes (%)"           = "11. WFS from Severe Episodes (%)",
  "Average Recovery (# weeks)"             = "12. Average Recovery (# weeks)"
)
apply_metric_numbers <- function(df) {
  rename_with(df, ~ metric_numbering[.x] %||% .x)
}

# Function to process data for a given CSD to generate WFS and non-WFS PM2.5 and identify WFS episodes
process_csd <- function(
    csd_df,                                     # data.frame with weekly CSD-level PM2.5 data
    threshold_episode               = 0,        # Value that must be exceeded to count as a WFS episode (0 ug/m3)
    num_wfs_weeks_in_episode        = 2,        # Number of WFS-impacted weeks that need to be present to count as a WFS episode (2 weeks)
    threshold_severe_episode        = 250,      # Value that must be exceeded to count as a severe episode (250 ug/m3)
    num_wfs_weeks_in_severe_episode = 1,        # Number of WFS-impacted weeks that need to be present to count as a severe episode (1 weeks)
    max_no_smoke_gap                = 3,        # Maximum number of continuous non-WFS-impacted weeks that can be present and still count as an episode (3 weeks)
    weeks_needed_above_threshold    = 1,        # Maximum number of continuous non-WFS-impacted weeks that can be present and still count as an episode (3 weeks)
    min_smoke_thresh                = 0.005     # Minimum value for WFS PM2.5 - used to remove near-zero estimates
) {

  # Identify wildfire season and calculate modified z-score
  csd_df <- csd_df %>%
    arrange(epiweek) %>%
    mutate(
      total_cumulative = cumsum(weekly_sum_pm25),
      is_wildfire_season = (end_date >= as.Date(paste0(year_start, "-05-01")) &
                              start_date < as.Date(paste0(year_start, "-11-01"))),
      modified_z = calc_modified_z(weekly_sum_pm25)
    )
  
  # Calculate counterfactual value
  cf_slope <- csd_df %>%
    filter(!is.na(modified_z),
           dplyr::between(modified_z, -2, 2),
           weekly_n_days == 7) %>%
    pull(weekly_sum_pm25) %>% median(na.rm = TRUE) %>% replace_na(0)
  
  # Separate WFS and non-WFS components
  csd_df <- csd_df %>%
    mutate(
      non_wfs_weekly_sum_pm25 = ifelse(modified_z > 2 & is_wildfire_season,
                                       cf_slope, weekly_sum_pm25),
      non_wfs_weekly_avg_pm25 = ifelse(modified_z > 2 & is_wildfire_season,
                                       cf_slope / 7, weekly_avg_pm25),
      wfs_weekly_sum_pm25 = weekly_sum_pm25 - non_wfs_weekly_sum_pm25,
      wfs_weekly_avg_pm25 = weekly_avg_pm25 - non_wfs_weekly_avg_pm25,
      non_wfs_cumulative = cumsum(non_wfs_weekly_sum_pm25),
      wfs_cumulative = cumsum(wfs_weekly_sum_pm25)
    ) %>%
    mutate(
      wfs_weekly_avg_pm25 = if_else(wfs_weekly_avg_pm25 < min_smoke_thresh, 0, wfs_weekly_avg_pm25),
      wfs_weekly_sum_pm25 = if_else(wfs_weekly_sum_pm25 < min_smoke_thresh, 0, wfs_weekly_sum_pm25),
      is_smoke_impacted = as.integer(wfs_weekly_avg_pm25 > 0)
    )
  
  # Identify WFS episodes & severe episodes
  norm_ids <- identify_smoke_episodes(
    csd_df$is_smoke_impacted, csd_df$epiweek, csd_df$wfs_weekly_sum_pm25,
    episode_threshold = threshold_episode,
    num_weeks_in_episode = num_wfs_weeks_in_episode,
    max_no_smoke_gap = max_no_smoke_gap,
    weeks_needed_above_threshold = weeks_needed_above_threshold
  )
  
  sev_ids <- identify_smoke_episodes(
    csd_df$is_smoke_impacted, csd_df$epiweek, csd_df$wfs_weekly_sum_pm25,
    episode_threshold = threshold_severe_episode,
    num_weeks_in_episode = num_wfs_weeks_in_severe_episode,
    max_no_smoke_gap = max_no_smoke_gap,
    weeks_needed_above_threshold = weeks_needed_above_threshold
  )
  
  epi_id <- integer(nrow(csd_df))
  cur    <- 0; in_epi <- FALSE
  for (i in seq_len(nrow(csd_df))) {
    if (norm_ids[i] > 0 || sev_ids[i] > 0) {
      if (!in_epi) { cur <- cur + 1; in_epi <- TRUE }
      epi_id[i] <- cur
    } else in_epi <- FALSE
  }
  
  csd_df$episode_id <- epi_id
  csd_df$severe_episode_id <- sev_ids
  
  # Return processed data
  csd_df
}

# Function to compute the 12 MultiWiSE metrics for a given CSD
compute_metrics_csd <- function(
    csd_df,                   # data.frame with weekly CSD-level PM2.5 and episode data
    threshold_low    = 5,     # Threshold used for metric #5
    threshold_high   = 25,    # Threshold used for metric #6
    min_smoke_thresh = 0.005
) {
  
  # Conversion for ug to mg
  MICRO_TO_MILLI <- 1/1000
  
  csd_df <- csd_df %>%
    mutate(
      wfs_weekly_avg_pm25 = if_else(wfs_weekly_avg_pm25 < min_smoke_thresh, 0, wfs_weekly_avg_pm25),
      wfs_weekly_sum_pm25 = if_else(wfs_weekly_sum_pm25 < min_smoke_thresh, 0, wfs_weekly_sum_pm25)
    )
  
  # 1. Cumulative WFS PM2.5
  wfs_cumulative <- sum(csd_df$wfs_weekly_sum_pm25, na.rm = TRUE) * MICRO_TO_MILLI
  
  # 2. WFS Fraction
  wfs_fraction <- 100 * wfs_cumulative / (sum(csd_df$weekly_sum_pm25, na.rm = TRUE) * MICRO_TO_MILLI)
  
  # 3. Average WFS PM2.5
  mean_smoke_week <- if (sum(csd_df$wfs_weekly_avg_pm25 > 0, na.rm = TRUE) > 0)
    mean(csd_df$wfs_weekly_avg_pm25[csd_df$wfs_weekly_avg_pm25 > 0])
  else 0
  
  # 4. Any WFS
  total_smoke_weeks <- sum(csd_df$wfs_weekly_avg_pm25 > 0, na.rm = TRUE)
  
  # 5. WFS PM2.5 > 5
  weeks_wfs_over_5 <- sum(csd_df$wfs_weekly_avg_pm25 > threshold_low, na.rm = TRUE)
  
  # 6. Total PM2.5 > 25
  weeks_over_25 <- sum(csd_df$wfs_weekly_avg_pm25 > 0 & csd_df$weekly_avg_pm25 > threshold_high, na.rm = TRUE)
  
  # 7. WFS Episodes
  n_normal <- max(csd_df$episode_id, na.rm = TRUE)
  
  # 8. Severe Episodes
  n_severe <- max(csd_df$severe_episode_id, na.rm = TRUE)
  
  # 9. Longest Episode
  ep_summary_all <- csd_df %>%
    filter(episode_id > 0) %>%
    group_by(episode_id) %>%
    summarise(
      active_weeks = n(),
      .groups = "drop"
    )
  longest_episode_len <- ifelse(nrow(ep_summary_all) > 0, max(ep_summary_all$active_weeks), 0)
  
  # 10. Worst Episode
  ep_summary_smk <- csd_df %>%
    filter(episode_id > 0, is_smoke_impacted == 1) %>%
    group_by(episode_id) %>%
    summarise(
      avg_pm25 = mean(wfs_weekly_avg_pm25),
      .groups = "drop"
    )
  worst_exposure <- ifelse(nrow(ep_summary_smk) > 0, max(ep_summary_smk$avg_pm25), 0)
  
  # 11. WFS from Severe Episodes
  severe_pm25_cum <- csd_df %>%
    filter(severe_episode_id > 0, is_smoke_impacted == 1) %>%
    summarise(sum(wfs_weekly_sum_pm25, na.rm = TRUE)) %>% pull()
  severe_pm25_cum <- severe_pm25_cum * MICRO_TO_MILLI
  pct_wfs_from_severe <- ifelse(wfs_cumulative > 0, 100 * severe_pm25_cum / wfs_cumulative, 0)
  
  # 12. Average Recovery 
  period_start <- min(csd_df$start_date, na.rm = TRUE)
  period_end   <- max(csd_df$end_date,   na.rm = TRUE)
  if (n_normal == 0) {    # No episodes -> recovery is the full window
    avg_time_between <- round(as.numeric(period_end - period_start + 1) / 7, 2)
  } else {
    episode_info <- csd_df %>%
      filter(episode_id > 0, is_smoke_impacted == 1) %>%
      group_by(episode_id) %>%
      summarise(
        start_date = min(start_date),
        end_date = max(end_date),
        .groups = "drop"
      ) %>%
      arrange(start_date)
    gap_starts <- c(period_start, episode_info$end_date + 1) # first day after each episode
    gap_ends <- c(episode_info$start_date - 1, period_end) # last day before next episode
    gaps_days <- pmax(0, as.numeric(gap_ends - gap_starts + 1)) 
    avg_time_between <- round(mean(gaps_days) / 7, 2)
  }
  
  # Create a table (tibble) that assigns the metrics to the official metric name
  numbered <- tibble(
    `Cumulative WFS PM2.5 (mg/m³)` = round(wfs_cumulative, 2),
    `WFS Fraction (%)` = round(wfs_fraction, 2),
    `Average WFS PM2.5 (µg/m³)`= round(mean_smoke_week, 2),
    `Any WFS (# weeks)` = total_smoke_weeks,
    `WFS PM2.5 > 5 µg/m³ (# weeks)` = weeks_wfs_over_5,
    `Total PM2.5 > 25 µg/m³ (# weeks)` = weeks_over_25,
    `WFS Episodes (# episodes)` = n_normal,
    `Severe Episodes (# episodes)` = n_severe,
    `Longest Episode (# weeks)` = longest_episode_len,
    `Worst Episode (µg/m³)` = round(worst_exposure, 2),
    `WFS from Severe Episodes (%)` = round(pct_wfs_from_severe, 1),
    `Average Recovery (# weeks)` = avg_time_between
  ) |> apply_metric_numbers()
  
  # Additional variables (total, non-WFS, and counterfactual) to be included in results
  extras <- tibble::tibble(
    `Cumulative_Total_PM25`     = round(sum(csd_df$non_wfs_weekly_sum_pm25,  na.rm = TRUE) * MICRO_TO_MILLI, 2),
    `Average_Total_PM25`        = round(mean(csd_df$weekly_avg_pm25, na.rm = TRUE), 2),
    `Cumulative_NonWFS_PM25`    = round(sum(csd_df$weekly_sum_pm25, na.rm = TRUE) * MICRO_TO_MILLI, 2),
    `Average_NonWFS_PM25`       = round(mean(csd_df$non_wfs_weekly_avg_pm25, na.rm = TRUE), 2),
    `Counterfactual_Value`      = round(unique(csd_df$non_wfs_weekly_avg_pm25[csd_df$is_smoke_impacted == 1]), 2)
  )
  
  dplyr::bind_cols(numbered, extras)
}

####################################################################################################################################
############################ PROCESS WEEKLY CSD DATA & GENERATE METRICS ############################################################
####################################################################################################################################

csd_weekly <- readRDS("./Processed_Data/CSD_weekly_PM25_2010-2023.RDS")

# Exclude CSDs with less than 75% of data for the study period
csd_weekly <- csd_weekly %>%
  group_by(csd) %>%
  filter(n() >= 0.75 * 730) %>%
  ungroup()

# Process CSD data to get estimates of WFS and non-WFS PM2.5 and episodes
csd_processed <- csd_weekly %>%
  group_by(csd) %>%
  group_split() %>%
  map_dfr(process_csd)

# Save processed CSD data to a .RDS and .csv file
saveRDS(csd_processed, "./Processed_Data/CSD_weekly_PM25_2010-2023_processed.RDS")
write_csv(csd_processed[,c('csd',
                          'epiweek',
                          'weekly_n_days',
                          'start_date',
                          'end_date',
                          'weekly_avg_pm25',
                          'wfs_weekly_avg_pm25',
                          'non_wfs_weekly_avg_pm25',
                          'weekly_sum_pm25',
                          'wfs_weekly_sum_pm25',
                          'non_wfs_weekly_sum_pm25',
                          'modified_z',
                          'is_wildfire_season',
                          'is_smoke_impacted',
                          'episode_id',
                          'severe_episode_id')],
          './Results/BC_CSD_Weekly_PM25_Episode_Estimates_2010-2023.csv')

# Calculate the 12 MultiWiSE metrics
all_csd_metrics <- csd_processed %>%
  group_by(csd) %>%
  group_split() %>%
  map_dfr(~ tibble(csd = unique(.x$csd)) %>% bind_cols(compute_metrics_csd(.x)))

# Save CSD metrics to .csv file and .RDS file
saveRDS(all_csd_metrics[,1:13], "./Processed_Data/CSD_MultiWiSE_metrics_2010-2023.RDS")
names(all_csd_metrics) <- c('CSD', '1_Cumulative_WFS_PM25', 
                            '2_WFS_Fraction', 
                            '3_Average_WFS_PM25',
                            '4_Any_WFS' ,
                            '5_WFS_PM25_exceeds_5',
                            '6_Total_PM25_exceeds_25',
                            '7_WFS_Episodes',
                            '8_Severe_Episodes',
                            '9_Longest_Episode',
                            '10_Worst_Episode',
                            '11_WFS_from_Severe_Episodes',
                            '12_Average_Recovery',
                            'Cumulative_Total_PM25',
                            'Average_Total_PM25',
                            'Cumulative_NonWFS_PM25',
                            'Average_NonWFS_PM25',
                            'Counterfactual_Value'
)
write_csv(all_csd_metrics, "./Results/BC_CSD_MultiWiSE_Metrics_2010-2023.csv")


####################################################################################################################################
############################ GENERATE FIGURES & TABLES #############################################################################
####################################################################################################################################

# Define custom ggplot theme for maps
theme_map <- function(base_size = 26) {
  theme_void(base_size = base_size) +
    theme(
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(1.2, "cm"),
      legend.position = c(0.84, 0.70),
      legend.background = element_blank(),
      plot.title = element_markdown(hjust = 0.5, face = "bold"),
      plot.margin = margin(5,5,5,5)
    )
}

# Function to generate cumulative PM2.5 plot
plot_cumulative_pm25_by_year <- function(csd_df, y_limits, line_thin = 0.25) {
  
  csd_df <- csd_df %>% 
    mutate(total_cumulative = total_cumulative / 1000,
           non_wfs_cumulative = non_wfs_cumulative / 1000, 
           wfs_cumulative = wfs_cumulative / 1000)
  
  ggplot(csd_df, aes(x = start_date)) +
    geom_line(aes(y = total_cumulative,       color = "Total"),    linewidth = 1) +
    geom_line(aes(y = non_wfs_cumulative, color = "Non-WFS"), linewidth = 1) +
    geom_line(aes(y = wfs_cumulative,  color = "WFS"),     linewidth = 1) +
    scale_color_manual(
      values = c("Total" = "black", "Non-WFS" = "#B0C6DF", "WFS" = "#FF8572"),
      name   = NULL
    ) +
    scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand      = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = y_limits
    ) +
    labs( title = "A",
          x = "Year",
          y = expression("Cumulative PM"[2.5]*" (mg/m"^3*")")
    ) +
    theme_classic(base_size = 26) +
    theme(
      legend.position = "top",
      axis.line = element_line(linewidth = line_thin),
      plot.title = element_text(hjust = 0),
      axis.ticks  = element_line(linewidth = line_thin)
    )
}

# Function to generate PM2.5 time-series plot
plot_weekly_avg_pm25_by_year <- function(csd_df, y_limits, line_thin = 0.25) {
  
  ggplot(csd_df, aes(x = start_date)) +
    geom_line(aes(y = weekly_avg_pm25, color = "Total"),    linewidth = 1) +
    geom_line(aes(y = non_wfs_weekly_avg_pm25, color = "Non-WFS"), linewidth = 1) +
    geom_line(aes(y = wfs_weekly_avg_pm25, color = "WFS"),     linewidth = 1) +
    scale_color_manual(
      values = c("Total" = "black", "Non-WFS" = "#B0C6DF", "WFS" = "#FF8572"),
      name   = NULL
    ) +
    scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = y_limits
    ) +
    labs(title = "B",
         x = "Year",
         y = expression("Mean PM"[2.5]*" ("*mu*"g/m"^3*")")
    ) +
    theme_classic(base_size = 26) +
    theme(
      legend.position = "top",
      axis.line = element_line(linewidth = line_thin),
      plot.title = element_text(hjust = 0),
      axis.ticks  = element_line(linewidth = line_thin)
    )
}

# Function to generate histogram of weekly sums (used to identify counterfactual value)
plot_slope_histogram <- function(csd_df, y_limits, x_limits, binwidth = 8, line_thin = 0.25) {
  
  median_slope <- median(csd_df$weekly_sum_pm25, na.rm = TRUE)
  
  z_vals <- csd_df$modified_z
  z_vals_norm <- which(!is.na(z_vals) & z_vals >= -2 & z_vals <= 2)
  if (length(z_vals_norm) > 0) {
    min_norm_slope <- min(csd_df[z_vals_norm,]$weekly_sum_pm25, na.rm = TRUE)
    max_norm_slope <- max(csd_df[z_vals_norm,]$weekly_sum_pm25, na.rm = TRUE)
  } else {
    min_norm_slope <- max_norm_slope <- NA_real_
  }
  
  norm_slopes_rect <- data.frame(
    xmin = min_norm_slope, xmax = max_norm_slope,
    ymin = 0, ymax = max(y_limits)
  )
  
  ggplot() +
    # Z‐range rectangle
    geom_rect(
      data = norm_slopes_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "-2 < Z-Score < 2"),
      alpha = 0.3, inherit.aes = FALSE
    ) +
    # Histogram bars
    geom_histogram(
      data        = csd_df,
      aes(x = weekly_sum_pm25),
      binwidth    = binwidth,
      boundary    = 0,
      fill        = "grey90",
      color       = "black",
      linewidth   = line_thin,
      show.legend = FALSE
    ) +
    # Cutoff lines (no legend)
    { if (!is.na(min_norm_slope)) geom_vline(xintercept = min_norm_slope, linetype = "dashed", show.legend = FALSE) } +
    { if (!is.na(max_norm_slope)) geom_vline(xintercept = max_norm_slope, linetype = "dashed", show.legend = FALSE) } +
    # Median line
    geom_vline(
      aes(xintercept = median_slope, color = "Median"),
      linewidth = 1.2
    ) +
    # Manual scales
    scale_fill_manual(
      name   = NULL,
      values = c("-2 < Z-Score < 2" = "yellow")
    ) +
    scale_color_manual(
      name   = NULL,
      values = c("Median" = "#0085E4")
    ) +
    guides(
      fill  = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0), oob = scales::oob_keep) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0), oob = scales::oob_keep) +
    labs(
      title = "C",
      x     = expression("Weekly Sum of PM"[2.5]*" ("*mu*"g/m"^3*" per epiweek)"),
      y     = "Count"
    ) +
    theme_classic(base_size = 26) +
    theme(
      plot.title      = element_text(hjust = 0),
      legend.position = "top",
      axis.line       = element_line(linewidth = line_thin),
      panel.grid      = element_blank(),
      axis.ticks      = element_line(linewidth = line_thin)
    )
}

#### Table 1 - Summary of Metrics across all CSDs ####
all_metrics <- readRDS("./Processed_Data/CSD_MultiWiSE_metrics_2010-2023.RDS")

metrics_order <- c(
  "1. Cumulative WFS PM₂.₅ (mg/m³)",
  "2. WFS Fraction (%)",
  "3. Average WFS PM₂.₅ (µg/m³)",
  "4. Any WFS (# weeks)",
  "5. WFS PM₂.₅ > 5 µg/m³ (# weeks)",
  "6. Total PM₂.₅ > 25 µg/m³ (# weeks)",
  "7. WFS Episodes (# episodes)",
  "8. Severe Episodes (# episodes)",
  "9. Longest Episode (# weeks)",
  "10. Worst Episode (µg/m³)",
  "11. WFS from Severe Episodes (%)",
  "12. Average Recovery (# weeks)"
)

get_summary <- function(x) {
  c(sprintf("%.4f (%.4f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)),
    sprintf("%.4f (%.4f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE)),
    sprintf("(%.4f, %.4f)", min(x, na.rm = TRUE),  max(x, na.rm = TRUE)))
}

table_data <- metrics_order %>%
  setNames(nm = .) %>%
  lapply(function(metric) {
    if (!metric %in% names(all_metrics)) return(NULL)
    tibble(
      Metric            = metric,
      `Mean (SD)`       = get_summary(all_metrics[[metric]])[1],
      `Median (IQR)`    = get_summary(all_metrics[[metric]])[2],
      `Range (Min, Max)`= get_summary(all_metrics[[metric]])[3]
    )
  }) %>% bind_rows() %>%
  mutate(Category = case_when(
    Metric %in% metrics_order[1:2]   ~ "Cumulative Exposure",
    Metric %in% metrics_order[3:6]   ~ "Weekly Exposure",
    Metric %in% metrics_order[7:12]  ~ "Episode Exposure"
  ))

table_one <- table_data %>%
  gt(rowname_col = "Metric", groupname_col = "Category") %>%
  tab_stubhead(label = "Metric") %>%
  tab_style(style = cell_text(indent = px(20)), locations = cells_stub())

print(table_one)

#### Figure 2 - Kitimat Stikine-C VS Vernon plots ####
csd_processed <- readRDS("./Processed_Data/CSD_weekly_PM25_2010-2023_processed.RDS")

csd_ks <- filter(csd_processed, csd == "5949013") # Kitimat‑StikineC
csd_vn <- filter(csd_processed, csd == "5937014") # Vernon

# Determine relevant ranges for cumulative data across both CSDs so plots can share y axis
combo_cumul <- bind_rows(
  csd_ks %>% transmute(
    epiweek,
    total_cumulative = total_cumulative / 1000
  ),
  csd_vn %>% transmute(
    epiweek,
    total_cumulative = total_cumulative / 1000
  )
)

y_max_cumul   <- max(combo_cumul$total_cumulative, na.rm = TRUE)
y_limits_cumul <- c(0, y_max_cumul * 1.05)

combo_mean <- bind_rows(
  csd_ks %>% select(epiweek, weekly_avg_pm25),
  csd_vn %>% select(epiweek, weekly_avg_pm25)
)

y_max_mean <- max(combo_mean$weekly_avg_pm25, na.rm = TRUE)
y_limits_mean <- c(0, y_max_mean * 1.05)

combo_slope <- bind_rows(
  csd_ks %>% select(weekly_sum_pm25, modified_z),
  csd_vn %>% select(weekly_sum_pm25, modified_z)
)
x_max_slope <- max(combo_slope$weekly_sum_pm25, na.rm = TRUE)
x_limits_slope <- c(0, x_max_slope)

get_max_bin <- function(data, binwidth = 8) {
  tmp <- ggplot_build(ggplot(data, aes(x = weekly_sum_pm25)) + geom_histogram(binwidth = binwidth))
  max(tmp$data[[1]]$count, na.rm = TRUE)
}
max_bin_ks <- get_max_bin(csd_ks)
max_bin_vn <- get_max_bin(csd_vn)
y_limits_slope <- c(0, max(max_bin_ks, max_bin_vn) * 1.1)

# Build the 6 plots using the plotting functions
p1_vn <- plot_cumulative_pm25_by_year(csd_vn, y_limits = y_limits_cumul) +
  labs(title = NULL)
p2_vn <- plot_weekly_avg_pm25_by_year(csd_vn, y_limits = y_limits_mean) +
  labs(title = NULL)
p3_vn <- plot_slope_histogram(csd_vn, y_limits = y_limits_slope, x_limits = x_limits_slope) +
  labs(title = NULL)

p1_ks <- plot_cumulative_pm25_by_year(csd_ks, y_limits = y_limits_cumul) 
p2_ks <- plot_weekly_avg_pm25_by_year(csd_ks, y_limits = y_limits_mean) 
p3_ks <- plot_slope_histogram(csd_ks, y_limits = y_limits_slope, x_limits = x_limits_slope) 

# Combine with patchwork and add titles
title_ks <- ggplot() +
  geom_text(aes(x = 0.5, y = 0.5, label = "Kitimat-Stikine C"), 
            size = 12, fontface = "bold") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

title_vn <- ggplot() +
  geom_text(aes(x = 0.5, y = 0.5, label = "Vernon"), 
            size = 12, fontface = "bold") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

combined_plot <- (
  (title_ks | title_vn) /
    (p1_ks | p1_vn) /
    (p2_ks | p2_vn) /
    (p3_ks | p3_vn)
) +
  plot_layout(
    heights = c(0.15, 1, 1, 1)
  )

# Preview plot
print(combined_plot)

# Save as a .png
ggsave("./Results/Figure2.png",
       combined_plot, width = 28, height = 28, dpi = 300)

#### Figure 4 - Correlation Matrix ####
all_csd_metrics <- readRDS("./Processed_Data/CSD_MultiWiSE_metrics_2010-2023.RDS")

# Make HTML-safe names with subscripts and escaped periods
html_names <- names(all_csd_metrics) %>%
  str_replace_all("PM₂\\.₅", "PM<sub>2.5</sub>") %>%
  str_replace("^([0-9]+)\\.", "\\1\\\\.") %>%
  str_replace_all(" > ", " > ")
names(all_csd_metrics) <- html_names

# Create order of metric display
corr_order <- c(
  "12\\. Average Recovery (# weeks)",
  "11\\. WFS from Severe Episodes (%)",
  "10\\. Worst Episode (µg/m³)",
  "9\\. Longest Episode (# weeks)",
  "8\\. Severe Episodes (# episodes)",
  "7\\. WFS Episodes (# episodes)",
  "6\\. Total PM<sub>2.5</sub> > 25 µg/m³ (# weeks)",
  "5\\. WFS PM<sub>2.5</sub> > 5 µg/m³ (# weeks)",
  "4\\. Any WFS (# weeks)",
  "3\\. Average WFS PM<sub>2.5</sub> (µg/m³)",
  "2\\. WFS Fraction (%)",
  "1\\. Cumulative WFS PM<sub>2.5</sub> (mg/m³)"
)

# Build the Pearson correlation matrix
metrics_numeric <- all_csd_metrics %>%
  select(all_of(corr_order)) %>%
  as.data.frame()

cor_matrix <- cor(metrics_numeric,
                  use    = "complete.obs",
                  method = "pearson")

# Mask the upper triangle
masked_mat <- cor_matrix
masked_mat[upper.tri(masked_mat, diag = FALSE)] <- NA

corr_plot <- ggcorrplot(
  masked_mat,
  method        = "circle",
  type          = "full",
  lab           = TRUE,
  lab_size      = 3,
  ggtheme       = theme_minimal(),
  outline.color = "black",
  legend.title  = "Correlation"
) +
  scale_fill_gradient2(
    low      = "#B0C6DF",
    mid      = "white",
    high     = "#FF8572",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Correlation",
    guide    = guide_colorbar(frame.colour = "black",
                              ticks.colour = "black")
  ) +
  scale_size_continuous(range = c(12, 12), guide = FALSE) +
  theme(
    plot.title       = element_markdown(size = 18, hjust = 0.5),
    axis.text.x      = element_markdown(size = 14,
                                        angle = 45,
                                        hjust = 1),
    axis.text.y      = element_markdown(size = 14),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black",
                                    linewidth = 0.5)
  )

# Print plot
print(corr_plot)

# Save as .jpeg
ggsave(
  filename = "./Results/Figure4.jpeg",
  plot     = corr_plot,
  width    = 11,
  height   = 9,
  dpi      = 300
)

#### Figure 5 - MultiWiSE Metric Maps ####
# Get census shapefiles
bc_csd <- csd_census %>% 
  mutate(csd = GeoUID) %>%
  st_transform(crs = 26910)

# Load and convert metrics to long format
all_csd_metrics <- readRDS("./Processed_Data/CSD_MultiWiSE_metrics_2010-2023.RDS")

all_csd_metrics_long <- all_csd_metrics %>%
  pivot_longer(
    cols      = -csd,
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  select(csd, Metric, Value) %>%
  mutate(csd = str_extract(csd, "\\d+$"))

joined <- left_join(bc_csd, all_csd_metrics_long, by = "csd")

# 1. Cumulative WFS PM2.5 (mg/m3)
map1 <- joined %>%
  filter(Metric == "1. Cumulative WFS PM₂.₅ (mg/m³)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "mg/m³",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    ) 
  ) +
  labs(title = "<span>1. Cumulative WFS PM<sub>2.5</sub></span>") +
  theme_void() +
  theme_map(base_size = 24) +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.2,              
    height        = unit(0.5, "cm"),   
    line_width    = 1,                 
    text_cex      = 1.2,               
    pad_x         = unit(0.5, "cm"),  
    pad_y         = unit(0.5, "cm"),   
    bar_cols      = c("black", "white"),
    unit_category = "metric"
  )

# 2. WFS Fraction (%)
map2 <- joined %>%
  filter(Metric == "2. WFS Fraction (%)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "%",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>2. WFS Fraction</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 3. Average WFS PM2.5 (ug/m3)
map3 <- joined %>%
  filter(Metric == "3. Average WFS PM₂.₅ (µg/m³)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "µg/m³",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>3. Average WFS PM<sub>2.5</sub></span>") +
  theme_void() +
  theme_map(base_size = 24)  

# 4. Any WFS (# weeks)
map4 <- joined %>%
  filter(Metric == "4. Any WFS (# weeks)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "# weeks",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>4. Any WFS</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 5. WFS PM2.5 > 5 ug/m3 (# weeks)
map5 <- joined %>%
  filter(Metric == "5. WFS PM₂.₅ > 5 µg/m³ (# weeks)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "# weeks",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>5. WFS PM<sub>2.5</sub> > 5 µg/m³</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 6. Total PM2.5 > 25 ug/m3 (# weeks)
map6 <- joined %>%
  filter(Metric == "6. Total PM₂.₅ > 25 µg/m³ (# weeks)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "# weeks",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>6. Total PM<sub>2.5</sub> > 25 µg/m³</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 7. WFS Episodes (# episodes)
map7 <- joined %>%
  filter(Metric == "7. WFS Episodes (# episodes)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "# episodes",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>7. WFS Episodes</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 8. Severe Episodes (# episodes)
map8 <- joined %>%
  filter(Metric == "8. Severe Episodes (# episodes)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "# episodes",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>8. Severe Episodes</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 9. Longest Episode (# weeks)
map9 <- joined %>%
  filter(Metric == "9. Longest Episode (# weeks)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "# weeks",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>9. Longest Episode</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# 10. Worst Episode (ug/m3)
map10 <- joined %>%
  filter(Metric == "10. Worst Episode (µg/m³)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "µg/m³",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>10. Worst Episode</span>") +
  theme_void() +
  theme_map(base_size = 24) 


# 11. WFS from Severe Episodes (%)
map11 <- joined %>%
  filter(Metric == "11. WFS from Severe Episodes (%)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80", name = "%",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>11. WFS from Severe Episodes</span>") +
  theme_void() +
  theme_map(base_size = 24) 


# 12. Average Recovery (# weeks)
map12 <- joined %>%
  filter(Metric == "12. Average Recovery (# weeks)") %>%
  ggplot(aes(geometry = geometry, fill = Value)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = -1,
    na.value = "grey80", name = "# weeks",
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2.5,
      barheight      = 16
    )   ) +
  labs(title = "<span>12. Average Recovery</span>") +
  theme_void() +
  theme_map(base_size = 24) 

# Save combined plot as a .jpg
jpeg(
  filename = "./Results/Figure5.jpg",
  width    = 30,
  height   = 32,
  units    = "in",
  res      = 300
)

# Arrange the 12 maps on a grid (4 by 3)
grid.arrange(
  map1,  map2,  map3,  map4,
  map5,  map6,  map7,  map8,
  map9,  map10, map11, map12,
  ncol = 3
)

dev.off()

# Supplemental Figures & Tables
#### Table S1 - Kitimat Stikine-C VS Vernon table ####
csd_processed <- readRDS("./Processed_Data/CSD_weekly_PM25_2010-2023_processed.RDS")
csd_ks <- filter(csd_processed, csd == "5949013") # Kitimat‑StikineC
csd_vn <- filter(csd_processed, csd == "5937014") # Vernon
kitimat <- compute_metrics_csd(csd_ks)[,1:12]
vernon  <-  compute_metrics_csd(csd_vn)[,1:12]

subset_metrics <- left_join(
  pivot_longer(kitimat, everything(), names_to = "Metric", values_to = "Kitimat"),
  pivot_longer(vernon,  everything(), names_to = "Metric", values_to = "Vernon"),
  by = "Metric"
)

metric_order <- c(
  "1. Cumulative WFS PM₂.₅ (mg/m³)",
  "2. WFS Fraction (%)",
  "3. Average WFS PM₂.₅ (µg/m³)",
  "4. Any WFS (# weeks)",
  "5. WFS PM₂.₅ > 5 µg/m³ (# weeks)",
  "6. Total PM₂.₅ > 25 µg/m³ (# weeks)",
  "7. WFS Episodes (# episodes)",
  "8. Severe Episodes (# episodes)",
  "9. Longest Episode (# weeks)",
  "10. Worst Episode (µg/m³)",
  "11. WFS from Severe Episodes (%)",
  "12. Average Recovery (# weeks)"
)

subset_metrics_cat <- subset_metrics %>%
  mutate(
    Metric = str_replace_all(Metric, "PM2\\.5", "PM₂.₅"),
    Category = case_when(
      Metric %in% metric_order[1:2]  ~ "Cumulative Exposure",
      Metric %in% metric_order[3:6]  ~ "Weekly Exposure",
      Metric %in% metric_order[7:12] ~ "Episode Exposure",
      TRUE                           ~ "Other"
    ),
    Metric   = factor(Metric, levels = metric_order),
    Category = factor(
      Category,
      levels = c("Cumulative Exposure", "Weekly Exposure", "Episode Exposure", "Other")
    )
  ) %>%
  arrange(Category, Metric)

subset_metrics_gt <- subset_metrics_cat %>%
  gt(rowname_col = "Metric", groupname_col = "Category") %>%
  cols_hide("Category") %>%
  cols_label(Kitimat = "Kitimat‑Stikine C", Vernon = "Vernon") %>%
  fmt_number(columns = c(Kitimat, Vernon), decimals = 2, use_seps = FALSE) %>%
  tab_style(
    style     = cell_text(indent = px(20)),
    locations = cells_stub(rows = everything())
  ) %>%
  cols_align(
    align   = "left",
    columns = everything()
  )

print(subset_metrics_gt)

#### Figure S2 - Location of Kitimat Stikine-C VS Vernon on BC Map ####
bc_csd <- csd_census %>% 
  mutate(csd = GeoUID) %>%
  st_transform(crs = 26910)

target_csds <- c("5949013", "5937014")
bc_csd <- bc_csd %>%
  mutate(
    CSD_Type = case_when(
      GeoUID == "5949013" ~ "Kitimat-Stikine C",
      GeoUID == "5937014" ~ "Vernon",
      TRUE ~ "Non-Target"
    ),
    LabelName = case_when(
      GeoUID == "5949013" ~ "Kitimat-Stikine C",
      GeoUID == "5937014" ~ "Vernon",
      TRUE ~ NA_character_
    )
  )

label_centroids <- bc_csd %>%
  filter(CSD_Type != "Non-Target") %>%
  st_centroid()

vernon_bbox <- bc_csd %>% 
  filter(CSD_Type == "Vernon") %>% 
  st_bbox()

x_range <- vernon_bbox["xmax"] - vernon_bbox["xmin"]
y_range <- vernon_bbox["ymax"] - vernon_bbox["ymin"]

inset_xlim <- c(vernon_bbox["xmin"] - 0.1 * x_range, vernon_bbox["xmax"] + 0.1 * x_range)
inset_ylim <- c(vernon_bbox["ymin"] - 0.1 * y_range, vernon_bbox["ymax"] + 0.1 * y_range)

main_map <- ggplot() +
  geom_sf(data = filter(bc_csd, CSD_Type == "Non-Target"),
          fill = "grey88", color = NA) +
  geom_sf(data = filter(bc_csd, CSD_Type == "Kitimat-Stikine C"),
          fill = "#A79296",
          color = "#57494C",
          size = 0.8) +
  geom_sf(data = filter(bc_csd, CSD_Type == "Vernon"),
          fill = "#F49461",
          color = "#D06B35",
          size = 0.8) +
  geom_label_repel(
    data = label_centroids,
    aes(label = LabelName, geometry = geometry),
    stat = "sf_coordinates",
    size = 3,
    fill = "white",
    color = "black",
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.15, "lines"),
    label.size = 0.4,
    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
    box.padding = 1.5,
    min.segment.length = 0
  ) +
  annotation_scale(
    location      = "br",
    width_hint    = 0.2,
    bar_cols      = c("black", "white"),
    unit_category = "metric"
  ) +
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    style       = north_arrow_fancy_orienteering()
  ) +
  coord_sf(datum = NA) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

vernon_inset <- ggplot() +
  geom_sf(data = filter(bc_csd, CSD_Type == "Non-Target"),
          fill = "white", color = NA) +
  geom_sf(data = filter(bc_csd, CSD_Type == "Vernon"),
          fill = "#F49461", color = "#D06B35", size = 0.8) +
  coord_sf(xlim = inset_xlim, ylim = inset_ylim, datum = NA) +
  theme_minimal() +
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    panel.grid  = element_blank(),
    plot.background = element_rect(color = "black", linewidth = 0.5)
  ) +
  labs(title = "Vernon")

final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(vernon_inset, x = 0.61, y = 0.05, width = 0.25, height = 0.25)

print(final_map)

ggsave(
  filename = "./Results/FigureS2.jpg",
  plot     = final_map,
  width    = 16,
  height   = 9,
  dpi      = 300
)

#### Figure S3 - Map of the counterfactual weekly value ####
csd_data <- readRDS("./Processed_Data/CSD_weekly_PM25_2010-2023_processed.RDS")

bc_csd <- csd_census %>% 
  mutate(csd = GeoUID) %>%
  st_transform(crs = 26910)

# calculate counterfactual
counterfactual <- csd_data %>%
  filter(
    !is.na(modified_z),
    modified_z >= -2, modified_z <= 2,
    weekly_n_days == 7
  ) %>%
  group_by(csd) %>%
  summarize(
    median_slope  = median(weekly_sum_pm25, na.rm = TRUE),
    counterfactual_pm25 = median_slope / 7
  ) %>%
  ungroup()

s3 <- bc_csd %>%
  left_join(counterfactual, by = "csd") %>%
  ggplot(aes(geometry = geometry, fill = counterfactual_pm25)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette   = "YlOrRd", direction = 1, na.value = "grey80",
    name = expression(mu * g / m^3),
    guide  = guide_colorbar(barwidth = 2.5, barheight = 16)
  ) +
  theme_void() +
  theme(
    plot.title            = element_markdown(size = 30, hjust = 0.5),
    legend.title          = element_text(size = 26),
    legend.position       = "inside",
    legend.position.inside = c(0.78, 0.7),
    legend.text = element_text(size = 24)
  ) +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.2,              
    height        = unit(0.5, "cm"),  
    line_width    = 1,                 
    text_cex      = 1.2,             
    pad_x         = unit(0.5, "cm"),  
    pad_y         = unit(0.5, "cm"),   
    bar_cols      = c("black", "white"),
    unit_category = "metric"
  )

# save as JPEG
ggsave(
  filename = "./Results/FigureS3.jpg",
  plot     = s3,
  width    = 10,
  height   = 9,
  units    = "in",
  dpi      = 300
)

print(s3)

#### Figure S4 - Map of the weekly average non-wildfire smoke PM2.5 ####
csd_data <- readRDS("./Processed_Data/CSD_weekly_PM25_2010-2023_processed.RDS")

bc_csd <- csd_census %>% 
  mutate(csd = GeoUID) %>%
  st_transform(crs = 26910)

# calculate average non-wfs PM2.5
mean_nonwfs <- csd_data %>%
  group_by(csd) %>%
  summarize(mean_nonwfs_pm25 = mean(non_wfs_weekly_avg_pm25, na.rm = TRUE)) 

s4 <- bc_csd %>%
  left_join(mean_nonwfs, by = "csd") %>%
  ggplot(aes(geometry = geometry, fill = mean_nonwfs_pm25)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette   = "YlOrRd", direction = 1, na.value = "grey80",
    name = expression(mu * g / m^3),
    guide  = guide_colorbar(barwidth = 2.5, barheight = 16)
  ) +
  theme_void() +
  theme(
    plot.title            = element_markdown(size = 30, hjust = 0.5),
    legend.title          = element_text(size = 26),
    legend.position       = "inside",
    legend.position.inside = c(0.78, 0.7),
    legend.text = element_text(size = 24)
  ) +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.2,             
    height        = unit(0.5, "cm"),  
    line_width    = 1,                 
    text_cex      = 1.2,              
    pad_x         = unit(0.5, "cm"),   
    pad_y         = unit(0.5, "cm"),  
    bar_cols      = c("black", "white"),
    unit_category = "metric"
  )

# save as JPEG
ggsave(
  filename = "./Results/FigureS4.jpg",
  plot     = s4,
  width    = 10,      
  height   = 9,       
  units    = "in",
  dpi      = 300
)

print(s4)


#### Figure S5 - Map of the weekly average total PM2.5 ####
csd_processed <- readRDS("./Processed_Data/CSD_weekly_PM25_2010-2023_processed.RDS")

bc_csd <- csd_census %>% 
  mutate(csd = GeoUID) %>%
  st_transform(crs = 26910)

# calculate average total PM2.5
mean_total_pm25 <- csd_processed %>%
  group_by(csd) %>%
  summarise(total_pm25 = round(mean(weekly_avg_pm25),2), .groups = "drop")

map_data <- left_join(bc_csd, mean_total_pm25, by = "csd")

s5 <- ggplot(map_data, aes(geometry = geometry, fill = total_pm25)) +
  geom_sf(color = NA) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    na.value = "grey80",
    name = "mg/m³",
    guide = guide_colorbar(barwidth = 2.5, barheight = 16)
  ) +
  theme_void(base_size = 24) +
  theme(
    plot.title            = element_markdown(size = 30, hjust = 0.5),
    legend.title          = element_text(size = 26),
    legend.position       = "inside",
    legend.position.inside = c(0.78, 0.7),
    legend.text = element_text(size = 24)
  ) +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.2,
    height        = unit(0.5, "cm"),
    line_width    = 1,
    text_cex      = 1.2,
    pad_x         = unit(0.5, "cm"), 
    pad_y         = unit(0.5, "cm"),
    bar_cols      = c("black", "white"),
    unit_category = "metric"
  )

print(s5)

# Save as .jpg
ggsave(
  filename = "./Results/FigureS5.jpg",
  plot     = s5,
  width    = 10,
  height   = 9,
  dpi      = 300
)

