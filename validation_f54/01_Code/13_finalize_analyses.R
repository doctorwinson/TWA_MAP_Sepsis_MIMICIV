invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
scripts<-c("05_splines_diagnostics.R","06_enriched_and_tables.R","08_time_varying_covariates.R","12_source_sensitivities.R","07_figures.R")
for(s in scripts){
 cat("Starting",s,"\n")
 status<-system2(file.path(R.home("bin"),"Rscript.exe"),file.path("01_Code",s),stdout=file.path("04_QC",paste0(s,".log")),stderr=file.path("04_QC",paste0(s,".stderr.log")))
 if(status!=0)stop("Stage failed: ",s)
}
cat("Canonical 30-iteration analyses and plots complete\n")
