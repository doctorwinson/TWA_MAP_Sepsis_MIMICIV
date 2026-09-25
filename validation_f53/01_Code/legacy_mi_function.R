run_mi_hours_only <- function(df, covars, hour_cols, m = 40, maxit = 10, seed = 42, max_retry_drop_na = 2) {
    keep_cols <- unique(c("stay_id", "time_lm_days", "event_lm", covars, hour_cols))
    dat0 <- as.data.frame(df[, ..keep_cols])
    ini <- mice(dat0, maxit = 0, printFlag = FALSE)
    meth <- ini$method
    pred <- ini$predictorMatrix
    meth[] <- ""
    meth[hour_cols] <- "pmm"
    for (v in c("stay_id", "time_lm_days", "event_lm")) {
        if (v %in% names(meth))
            meth[v] <- ""
    }
    pred[, ] <- 0
    cov_for_pred <- intersect(covars, colnames(pred))
    for (h in 0:23) {
        target <- sprintf("h%02d", h)
        if (!target %in% rownames(pred))
            next
        nb <- (h - 2):(h + 2)
        nb <- nb[nb >= 0 & nb <= 23]
        nb_cols <- sprintf("h%02d", nb)
        nb_cols <- nb_cols[nb_cols %in% colnames(pred)]
        pred[target, nb_cols] <- 1
        pred[target, target] <- 0
        if (length(cov_for_pred) > 0)
            pred[target, cov_for_pred] <- 1
        if ("time_lm_days" %in% colnames(pred))
            pred[target, "time_lm_days"] <- 1
        if ("event_lm" %in% colnames(pred))
            pred[target, "event_lm"] <- 1
    }
    set.seed(seed)
    imp <- mice(dat0, m = m, maxit = maxit, method = meth, predictorMatrix = pred, ridge = 0.01, remove.collinear = TRUE,
        remove.constant = TRUE, printFlag = TRUE)
    bad_ids <- c()
    for (i in seq_len(m)) {
        di <- mice::complete(imp, i)
        na_row <- which(rowSums(is.na(di[, hour_cols, drop = FALSE])) > 0)
        if (length(na_row) > 0)
            bad_ids <- union(bad_ids, di$stay_id[na_row])
    }
    bad_ids <- unique(bad_ids)
    if (length(bad_ids) == 0)
        return(list(imp = imp, dropped_stay_id = character(0)))
    if (max_retry_drop_na <= 0)
        return(list(imp = imp, dropped_stay_id = bad_ids))
    df2 <- as.data.table(df)[!stay_id %in% bad_ids]
    res2 <- run_mi_hours_only(df = df2, covars = covars, hour_cols = hour_cols, m = m, maxit = maxit, seed = seed, max_retry_drop_na = max_retry_drop_na -
        1)
    res2$dropped_stay_id <- unique(c(bad_ids, res2$dropped_stay_id))
    res2
}
