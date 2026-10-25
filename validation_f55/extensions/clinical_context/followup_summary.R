invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
source_dir <- args[1]
output_dir <- args[2]
eligible <- readRDS(file.path(source_dir, 'eligible_primary.rds'))
datasets <- readRDS(file.path(source_dir, 'analysis_m40.rds'))
stopifnot(nrow(eligible) == 7289, sum(eligible$event_lm) == 948,
          length(datasets) == 40, !anyDuplicated(eligible$subject_id),
          !anyNA(eligible$time_lm_days), all(eligible$time_lm_days > 0),
          all(eligible$time_lm_days <= 29))
for (d in datasets) {
    ix <- match(eligible$stay_id, d$stay_id)
    stopifnot(!anyNA(ix), nrow(d) == nrow(eligible),
              identical(eligible$event_lm, d$event_lm[ix]),
              isTRUE(all.equal(eligible$time_lm_days, d$time_lm_days[ix])))
}
# Reconstruct the existing endpoint from cached death timing, without fitting models.
expected_event <- as.integer(!is.na(eligible$death_hours) & eligible$death_hours > 24 & eligible$death_hours <= 720)
expected_time <- ifelse(expected_event == 1, (eligible$death_hours - 24) / 24, 29)
stopifnot(identical(as.integer(eligible$event_lm), expected_event),
          isTRUE(all.equal(eligible$time_lm_days, expected_time)),
          all(eligible$time_lm_days[eligible$event_lm == 0] == 29))
t <- eligible$time_lm_days
summary <- data.frame(n = length(t), deaths = sum(eligible$event_lm),
                      censored = sum(eligible$event_lm == 0),
                      mean_days = mean(t), median_days = median(t),
                      q1_days = unname(quantile(t, .25)), q3_days = unname(quantile(t, .75)),
                      total_person_days = sum(t), min_days = min(t), max_days = max(t))
write.csv(summary, file.path(output_dir, 'followup_summary.csv'), row.names = FALSE)
writeLines(c('All 40 completed analysis datasets match the original eligible cohort on identity, event and follow-up time.',
             'The saved event and follow-up fields match reconstruction from saved death_hours.',
             'Summary uses elapsed time from hour 24 to recorded death or administrative censoring at ICU-admission day 30.',
             'No model or imputation was rerun. Only aggregate output is written.'), file.path(output_dir, 'followup_verification.txt'))
print(summary, digits = 12)
