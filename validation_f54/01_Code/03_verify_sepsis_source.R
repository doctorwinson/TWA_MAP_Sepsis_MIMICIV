invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(DBI);library(odbc);library(data.table)})
con<-dbConnect(odbc(),dsn="mimic4_v31",timeout=15)
dbExecute(con,"SET default_transaction_read_only=on")
x<-as.data.table(dbGetQuery(con,paste(readLines("01_Code/03_verify_sepsis_source.sql"),collapse="\n")))
dbDisconnect(con)
d<-as.data.table(readRDS("restricted_cache/source_enriched.rds"))
x<-x[stay_id %in% d$stay_id]
x[,sepsis3:=tolower(as.character(sepsis3)) %in% c("1","t","true")]
stopifnot(nrow(x)==nrow(d),!anyDuplicated(x$stay_id))
x[,both_pair_by24:=antibiotic_time<=intime+86400 & culture_time<=intime+86400]
x[,eligible_by24:=source_onset<=intime+86400]
x[,source_match:=source_infection_time==suspected_infection_time & source_sofa_time==sofa_time & source_onset==pmax(suspected_infection_time,sofa_time)]
fwrite(data.table(metric=c("matched_source_rows","official_source_mismatch","official_sepsis_false_or_missing","SOFA_below2","early_with_pair_completed_after24","pair_timing_missing"),n=c(nrow(x),sum(!x$source_match|is.na(x$source_match)),sum(is.na(x$sepsis3)|!x$sepsis3),sum(x$sofa_score<2,na.rm=TRUE),sum(x$eligible_by24 & !x$both_pair_by24,na.rm=TRUE),sum(is.na(x$both_pair_by24)))),"04_QC/official_sepsis_audit.csv")
stopifnot(all(x$source_match),all(x$sepsis3),all(x$sofa_score>=2))
saveRDS(x,"restricted_cache/sepsis_source_verification.rds")
print(fread("04_QC/official_sepsis_audit.csv"))
