invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice)})
setDTthreads(2)
scope <- Sys.getenv("F54_MI_SCOPE","primary")
target_iterations <- as.integer(Sys.getenv("F54_MAXIT","30"))
stopifnot(scope %in% c("primary","enriched"))
logcon <- file(paste0("04_QC/imputation_",scope,".log"),"wt")
sink(logcon)
cat("Start",format(Sys.time()),"scope",scope,"\n")
d <- as.data.table(readRDS("restricted_cache/source_enriched.rds"))
spec <- readRDS("restricted_cache/pre_spec.rds")
primary <- spec$primary; hours <- spec$hours
d <- d[!is.na(onset_hours) & onset_hours<=24]
d <- d[complete.cases(d[,..primary])]
d[,gender:=factor(gender)]
stopifnot(nrow(d)==7001L,all(d$n_hour_obs>=18),all(d$time_lm_days>0))
pred_cov <- primary
if(scope=="enriched") {
  domains <- paste0("sofa_",c("respiration","coagulation","liver","cardiovascular","cns","renal"))
  for(v in domains) {
    mv <- paste0(v,"_missing")
    d[,(mv):=as.integer(is.na(get(v)))]
    set(d,which(is.na(d[[v]])),v,0)
  }
  stopifnot(all(rowSums(d[,..domains])==d$sofa))
  pred_cov <- unique(c(setdiff(primary,"sofa"),domains,paste0(domains,"_missing"),"icd_hypertension","icd_stroke","careunit_group","admission_group","cardio_proc_device","advanced_monitor","ne_log1p","service_group"))
  pred_cov <- pred_cov[vapply(d[,..pred_cov],function(x)uniqueN(x)>1L,logical(1))]
}
stopifnot(all(complete.cases(d[,..pred_cov])))
cols <- unique(c("stay_id","time_lm_days","event_lm",pred_cov,hours))
dat <- droplevels(as.data.frame(d[,..cols]))
ini <- mice(dat,maxit=0,printFlag=FALSE)
meth <- ini$method; meth[] <- ""; meth[hours] <- "pmm"
pm <- ini$predictorMatrix;pm[,] <- 0
for(h in 0:23) {
  target <- sprintf("h%02d",h)
  neighbors <- sprintf("h%02d",setdiff(max(0,h-2):min(23,h+2),h))
  pm[target,c(neighbors,pred_cov,"time_lm_days","event_lm")] <- 1
}
set.seed(42)
resume_path <- Sys.getenv("F54_RESUME_MI","")
if(nzchar(resume_path)) {
  imp <- readRDS(resume_path)
  stopifnot(identical(imp$data,dat),imp$m==40L,imp$iteration<target_iterations)
  imp <- mice.mids(imp,maxit=target_iterations-imp$iteration,printFlag=TRUE)
}else{
  imp <- mice(dat,m=40,maxit=target_iterations,method=meth,predictorMatrix=pm,ridge=.01,remove.collinear=TRUE,remove.constant=TRUE,printFlag=TRUE)
}
for(i in 1:40) {
  z <- complete(imp,i)
  stopifnot(nrow(z)==nrow(dat),!anyNA(z[,hours]))
  for(h in hours) stopifnot(identical(z[[h]][!is.na(dat[[h]])],dat[[h]][!is.na(dat[[h]])]))
}
saveRDS(imp,paste0("restricted_cache/mice40_",scope,".rds"))
if(scope=="primary")saveRDS(d,"restricted_cache/eligible_primary.rds")
if(!is.null(imp$loggedEvents))fwrite(as.data.table(imp$loggedEvents),paste0("04_QC/MI_",scope,"_logged_events.csv"))
fwrite(data.table(scope=scope,n=nrow(dat),m=imp$m,iterations=imp$iteration,events=sum(dat$event_lm),residual_missing_hours=0,observed_values_unchanged=TRUE),paste0("04_QC/MI_",scope,"_QC.csv"))
if(scope=="primary") {
  z <- complete(imp,1)
  vals <- rowMeans(z[,hours])
  cuts <- as.numeric(quantile(vals,c(0,.25,.5,.75,1)))
  knots <- as.numeric(quantile(vals,c(.10,.50,.90)))
  spec$cutpoints <- cuts;spec$knots <- knots;spec$boundary <- range(vals);spec$reference <- median(vals)
  saveRDS(spec,"restricted_cache/model_spec.rds")
  fwrite(data.table(boundary=c("min","q25","q50","q75","max"),value=cuts),"02_Results/cutpoints.csv")
  fwrite(data.table(hour=hours,missing=vapply(dat[,hours],function(x)sum(is.na(x)),integer(1)),n=nrow(dat)),"02_Results/hourly_missingness.csv")
}
capture.output(sessionInfo(),file=paste0("04_QC/session_MI_",scope,".txt"))
cat("End",format(Sys.time()),"\n");sink();close(logcon)
cat("Completed",scope,"m=40\n")
