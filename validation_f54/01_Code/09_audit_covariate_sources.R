invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(DBI);library(odbc);library(data.table)})
con<-dbConnect(odbc(),dsn="mimic4_v31",timeout=15)
dbExecute(con,"SET default_transaction_read_only=on")
c<-as.data.table(dbGetQuery(con,"SELECT table_name,column_name FROM information_schema.columns WHERE table_schema='mimiciv_derived' AND table_name IN ('age','charlson','first_day_weight','weight_durations','first_day_sofa','first_day_rrt','first_day_ventilation') ORDER BY table_name,ordinal_position"))
print(c);fwrite(c,"04_QC/official_covariate_columns.csv")
x<-as.data.table(dbGetQuery(con,"SELECT p.stay_id,p.age,ag.age official_age,p.weight,w.weight official_weight,w.weight_admit,p.charlson,c.charlson_comorbidity_index official_charlson,p.itvtn_24h_rrt_tag,r.dialysis_active,r.dialysis_present FROM bdmcc.bdmcc_population p LEFT JOIN mimiciv_derived.age ag USING(hadm_id) LEFT JOIN mimiciv_derived.first_day_weight w USING(stay_id) LEFT JOIN mimiciv_derived.charlson c ON c.hadm_id=p.hadm_id LEFT JOIN mimiciv_derived.first_day_rrt r USING(stay_id) WHERE p.crtr_sepsis3=1"))
ids<-readRDS("restricted_cache/eligible_primary.rds")$stay_id
x<-x[stay_id %in% ids]
stopifnot(nrow(x)==length(ids))
results<-data.table(comparison=c("age vs official age","weight vs first_day_weight.weight","weight vs first_day_weight.weight_admit","Charlson vs official total","RRT vs dialysis_active","RRT vs dialysis_present"),
 differences=c(sum(abs(x$age-x$official_age)>1e-8,na.rm=TRUE),sum(abs(x$weight-x$official_weight)>1e-8,na.rm=TRUE),sum(abs(x$weight-x$weight_admit)>1e-8,na.rm=TRUE),sum(x$charlson!=x$official_charlson,na.rm=TRUE),sum(x$itvtn_24h_rrt_tag!=x$dialysis_active,na.rm=TRUE),sum(x$itvtn_24h_rrt_tag!=x$dialysis_present,na.rm=TRUE)))
print(results);fwrite(results,"04_QC/covariate_provenance_comparisons.csv")
dbDisconnect(con)
