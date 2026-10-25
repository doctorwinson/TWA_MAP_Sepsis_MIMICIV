# Aggregate-only extension of the retained F55 cohort; no exposure imputation is rerun.
# Usage: Rscript context_revision.R <F55 analysis folder> <aggregate output folder>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages({library(data.table); library(survival); library(mice)})
src <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- as.data.table(readRDS(file.path(src, "restricted_cache/eligible_primary.rds")))
dl <- readRDS(file.path(src, "restricted_cache/analysis_m40.rds"))
sp <- readRDS(file.path(src, "restricted_cache/analysis_specs.rds"))
stopifnot(nrow(d) == 7289L, sum(d$event_lm) == 948L, length(dl) == 40L,
          identical(d$stay_id, dl[[1]]$stay_id), !anyNA(d$onset_hours))
d[, map_q := as.character(dl[[1]]$map_q)]
d[, recording := ifelse(n_hour_obs == 24L, "Complete 24 hours", "18-23 observed hours")]
num <- function(x, digits = 1) {
  q <- quantile(x, c(.5,.25,.75), na.rm = TRUE)
  sprintf(paste0("%.",digits,"f [%.",digits,"f, %.",digits,"f]"),q[1],q[2],q[3])
}
bin <- function(x) sprintf("%d (%.1f%%)", sum(x), 100*mean(x))
group_summary <- function(z, label) {
  data.table(group = label, n = nrow(z), deaths = sum(z$event_lm),
             onset_hours = num(z$onset_hours, 2), onset_at_or_before_0 = bin(z$onset_hours <= 0),
             onset_after_0 = bin(z$onset_hours > 0), onset_after_6 = bin(z$onset_hours > 6),
             onset_after_12 = bin(z$onset_hours > 12), onset_after_18 = bin(z$onset_hours > 18),
             first_observed_hour = num(apply(!is.na(as.matrix(z[,sprintf("h%02d",0:23),with=FALSE])),1,function(x) which(x)[1]-1L)),
             age = num(z$age,0), sofa = num(z$sofa,0), mortality = bin(z$event_lm == 1),
             ventilation = bin(z$itvtn_24h_vent_tag == 1), vasopressor = bin(z$ne_equiv_mean_0_24h > 0),
             neuro_service = bin(z$service_group == "Neuro service"),
             cardiac_icu = bin(z$careunit_group == "Cardiac ICU"),
             surgical_trauma_icu = bin(z$careunit_group == "Surgical/Trauma ICU"),
             medical_icu = bin(z$careunit_group == "Medical ICU"),
             early_emergency = bin(z$admission_group == "Emergency/Direct emergent"))
}
summ <- rbindlist(c(list(group_summary(d,"Overall")),
  lapply(c("Q1","Q2","Q3","Q4"),function(g) group_summary(d[map_q == g],g)),
  lapply(c("Complete 24 hours","18-23 observed hours"),function(g) group_summary(d[recording == g],g))))
fwrite(summ,file.path(out,"onset_recording_summary.csv"))

# A precise hospital death timestamp can locate death relative to the index ICU stay,
# but neither this timestamp nor diagnosis codes establish a cause of death.
d[, death_location := fifelse(is.na(deathtime), "Date-only death; ICU timing unresolved",
  fifelse(deathtime >= intime & deathtime <= outtime, "During index ICU stay", "Outside index ICU stay"))]
loc <- d[event_lm == 1L, .N, by = .(map_q,death_location)]
loc[, denominator := sum(N), by = map_q]
loc[, percent := 100*N/denominator]
fwrite(loc,file.path(out,"death_location.csv"))

rows <- list()
for(v in c("careunit_group","admission_group")) {
  for(l in sort(unique(as.character(d[[v]])))) {
    vals <- c(bin(d[[v]] == l),sapply(c("Q1","Q2","Q3","Q4"),function(q) bin(d[map_q==q][[v]] == l)))
    rows[[length(rows)+1L]] <- data.table(variable=v,level=l,Overall=vals[1],Q1=vals[2],Q2=vals[3],Q3=vals[4],Q4=vals[5])
  }
}
v <- "advanced_monitor"
vals <- c(bin(d[[v]] == 1),sapply(c("Q1","Q2","Q3","Q4"),function(q) bin(d[map_q==q][[v]] == 1)))
rows[[length(rows)+1L]] <- data.table(variable=v,level="PiCCO/NICOM, n (%)",Overall=vals[1],Q1=vals[2],Q2=vals[3],Q3=vals[4],Q4=vals[5])
fwrite(rbindlist(rows),file.path(out,"table1_context_rows.csv"))
stopifnot(sum(d$admission_group == "Emergency/Direct emergent") == 2763L,
          sum(d$careunit_group == "Cardiac ICU") == 4199L)

# Refit the already reported all-covariate piecewise specification, extracting the
# Q4-by-late-window contrast directly. Do not compare overlap of marginal CIs.
est <- rbindlist(lapply(seq_along(dl),function(i) {
  x <- dl[[i]]
  mm <- model.matrix(as.formula(paste("~map_q +",paste(sp$primary,collapse=" + "))),data=x)[,-1,drop=FALSE]
  z <- data.frame(time=x$time_lm_days,event=x$event_lm,mm,check.names=FALSE)
  z <- survSplit(Surv(time,event)~.,data=z,cut=6,episode="window")
  cn <- colnames(mm)
  for(v in cn) z[[paste0(v,"_late")]] <- z[[v]]*as.integer(z$window==2)
  fit <- coxph(as.formula(paste("Surv(tstart,time,event) ~",paste(c(cn,paste0(cn,"_late"),"strata(window)"),collapse=" + "))),data=z)
  b <- coef(fit); V <- vcov(fit); q <- "map_qQ4"; t <- "map_qQ4_late"
  stopifnot(all(is.finite(b)))
  data.table(imputation=i,contrast=c("Early HR","Late HR","Late/early HR ratio"),
             b=c(b[q],b[q]+b[t],b[t]),u=c(V[q,q],V[q,q]+V[t,t]+2*V[q,t],V[t,t]),
             events=fit$nevent,k=length(b))
}))
pooled <- est[, {
  p <- pool.scalar(b,u,n=events[1],k=k[1])
  se <- sqrt(p$t); crit <- qt(.975,p$df)
  list(estimate=exp(p$qbar),lo=exp(p$qbar-crit*se),hi=exp(p$qbar+crit*se),
       p=2*pt(-abs(p$qbar/se),p$df),df=p$df,imputations=.N)
},by=contrast]
stopifnot(round(pooled[contrast=="Early HR",estimate],2)==2.25,
          round(pooled[contrast=="Late HR",estimate],2)==1.94)
fwrite(pooled,file.path(out,"piecewise_comparison.csv"))
fwrite(est,file.path(out,"piecewise_estimates_by_imputation.csv"))
capture.output(sessionInfo(), file=file.path(out,"R_sessionInfo.txt"))
cat("AGGREGATE SUMMARY\n"); print(summ); print(pooled); print(loc)
