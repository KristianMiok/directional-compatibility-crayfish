#!/usr/bin/env Rscript
# Stage 3 production SDM runner.
#
# Differs from R/01_run_stage2_esm_single_core.R in three ways:
#   1. Extracts maxSSS thresholds (biomod2 cutoff column) per CV run.
#   2. Sets the biomod2 working directory to a persistent stage3 workspace,
#      so trained models accumulate in one known location for Stage 3b to
#      load via BIOMOD_LoadModels().
#   3. Writes a per-species per-(predictor_set, core) metadata JSON
#      summarizing AUC/TSS/Boyce/maxSSS aggregates and QC flags.
#
# All Boyce sentinel handling matches Stage 2 (filter ecospat -1.0 with AUC > 0.5).

suppressPackageStartupMessages({
  library(biomod2)
  library(ecospat)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i[length(i)] + 1]
}

species       <- get_arg("--species")
predictor_set <- get_arg("--predictor-set")
core          <- toupper(get_arg("--core", "GBM"))
input_dir     <- get_arg("--input-dir",     "data/processed/stage3/biomod_inputs")
workspace     <- get_arg("--workspace",     "data/processed/stage3/biomod_workspace")
artifacts_dir <- get_arg("--artifacts-dir", "data/processed/stage3/sdm_artifacts")
cv_reps       <- as.integer(get_arg("--cv-reps", "5"))
seed          <- as.integer(get_arg("--seed", "42"))

if (is.null(species) || is.null(predictor_set)) stop("--species and --predictor-set are required")
if (!core %in% c("GLM", "GBM")) stop("--core must be GLM or GBM")

slug <- sub("_$", "", sub("^_", "", gsub("[^a-z0-9]+", "_", tolower(species))))
input_csv <- file.path(input_dir,
                       sprintf("%s__%s_biomod_pilot.csv", slug, predictor_set))
if (!file.exists(input_csv)) stop(sprintf("Missing biomod input: %s", input_csv))

dir.create(workspace, showWarnings = FALSE, recursive = TRUE)
dir.create(artifacts_dir, showWarnings = FALSE, recursive = TRUE)
species_artifact_dir <- file.path(artifacts_dir, slug)
dir.create(species_artifact_dir, showWarnings = FALSE, recursive = TRUE)

# Critical: make biomod2 drop its workspace folders into our known location
orig_wd <- getwd()
setwd(workspace)
on.exit(setwd(orig_wd), add = TRUE)

# --- helper: write empty/failed metadata JSON and exit ----------------------
write_failed_metadata <- function(reason) {
  meta <- list(
    species = species, species_slug = slug,
    predictor_set = predictor_set, core = core,
    fit_status = reason,
    n_records = NA, n_cv_runs = 0,
    auc_mean = NA, auc_std = NA,
    tss_mean = NA, tss_std = NA,
    boyce_mean = NA, boyce_std = NA, boyce_n_valid = 0,
    maxsss_mean = NA, maxsss_std = NA,
    qc_pass = FALSE, unreliable_boyce = TRUE
  )
  out_json <- file.path(species_artifact_dir,
                        sprintf("%s__%s_metadata.json", predictor_set, tolower(core)))
  write_json(meta, out_json, auto_unbox = TRUE, pretty = TRUE, na = "null")
  cat(sprintf("Wrote failed metadata: %s (%s)\n", out_json, reason))
}

# --- load input -------------------------------------------------------------
# input_csv may be absolute (when --input-dir was passed absolute) or relative
# (when defaulted). file.path() doesn't auto-detect, so check explicitly.
input_csv_abs <- if (substr(input_csv, 1, 1) == "/") input_csv else file.path(orig_wd, input_csv)
df <- read.csv(input_csv_abs, stringsAsFactors = FALSE)
resp_col <- "resp"; lon_col <- "long_or"; lat_col <- "lat_or"

if (!resp_col %in% names(df)) stop(sprintf("Missing response column: %s", resp_col))
if (!all(c(lon_col, lat_col) %in% names(df))) stop("Missing long_or / lat_or columns")

predictors <- setdiff(names(df), c("Crayfish_scientific_name", "subc_id", "basin_id",
                                    lon_col, lat_col, resp_col))
n_pres <- sum(df[[resp_col]] == 1, na.rm = TRUE)
n_abs  <- sum(df[[resp_col]] == 0, na.rm = TRUE)
cat(sprintf("Records: presences=%d absences=%d predictors=%d\n", n_pres, n_abs, length(predictors)))
cat(sprintf("Predictors: %s\n", paste(predictors, collapse = ", ")))

# --- format biomod data -----------------------------------------------------
biomod_data <- BIOMOD_FormatingData(
  resp.var = df[[resp_col]],
  expl.var = df[, predictors, drop = FALSE],
  resp.xy = df[, c(lon_col, lat_col)],
  resp.name = slug,
  filter.raster = FALSE
)

# --- build modeling options (Stage 2 patch carried over) --------------------
opt_user_val <- NULL
if (core == "GBM") {
  gbm_params <- list(
    n.trees = 500, interaction.depth = 2, shrinkage = 0.01,
    n.minobsinnode = 3, bag.fraction = 0.75, train.fraction = 1.0,
    cv.folds = 0, keep.data = TRUE, verbose = FALSE
  )
  opt_user_val <- list(GBM.binary.gbm.gbm = list(for_all_datasets = gbm_params))
}

modeling_id <- sprintf("%s__%s__esm_%s", slug, predictor_set, tolower(core))

mods <- tryCatch(
  BIOMOD_Modeling(
    bm.format = biomod_data,
    modeling.id = modeling_id,
    models = core,
    CV.strategy = "random",
    CV.nb.rep = cv_reps,
    CV.perc = 0.8,
    OPT.strategy = if (is.null(opt_user_val)) "bigboss" else "user.defined",
    OPT.user.val = opt_user_val,
    metric.eval = c("AUCroc", "TSS", "KAPPA"),
    var.import = 0,
    seed.val = seed
  ),
  error = function(e) {
    cat(sprintf("BIOMOD_Modeling failed: %s\n", conditionMessage(e)))
    NULL
  }
)

if (is.null(mods)) {
  write_failed_metadata("fit_error_biomod_modeling")
  quit(status = 0)
}

# --- evaluations: KEEP the cutoff column this time --------------------------
eval_df <- tryCatch(get_evaluations(mods), error = function(e) {
  cat(sprintf("get_evaluations failed: %s\n", conditionMessage(e)))
  NULL
})

if (is.null(eval_df) || !is.data.frame(eval_df) || nrow(eval_df) == 0) {
  write_failed_metadata("no_evaluations_returned")
  quit(status = 0)
}

# Map biomod2 metric names to ours
metric_map <- c("AUCroc" = "AUC", "TSS" = "TSS", "KAPPA" = "Kappa")
eval_df$metric <- ifelse(
  as.character(eval_df$metric.eval) %in% names(metric_map),
  metric_map[as.character(eval_df$metric.eval)],
  as.character(eval_df$metric.eval)
)

# Reshape KEEPING cutoff (this is the Stage 3 addition over Stage 2)
wide_val <- reshape(
  eval_df[, c("full.name", "PA", "run", "algo", "metric", "validation")],
  idvar = c("full.name", "PA", "run", "algo"),
  timevar = "metric", direction = "wide"
)
names(wide_val) <- sub("^validation\\.", "", names(wide_val))

wide_cut <- reshape(
  eval_df[, c("full.name", "PA", "run", "algo", "metric", "cutoff")],
  idvar = c("full.name", "PA", "run", "algo"),
  timevar = "metric", direction = "wide"
)
names(wide_cut) <- sub("^cutoff\\.", "cutoff_", names(wide_cut))

wide <- merge(wide_val, wide_cut,
              by = c("full.name", "PA", "run", "algo"), all = TRUE)

wide$species <- species
wide$species_slug <- slug
wide$model_core <- core
wide$predictor_set <- predictor_set

if (!("AUC"   %in% names(wide))) wide$AUC   <- NA_real_
if (!("TSS"   %in% names(wide))) wide$TSS   <- NA_real_
if (!("Kappa" %in% names(wide))) wide$Kappa <- NA_real_

# biomod2's cutoff column is on the 0-1000 scale; convert to 0-1
# We use the TSS-optimal cutoff as our maxSSS proxy
# (biomod2's "TSS" cutoff is the one maximizing Sens + Spec - 1 = Youden's J,
# which IS the max-Sens+Spec threshold = maxSSS)
wide$maxSSS <- if ("cutoff_TSS" %in% names(wide)) wide$cutoff_TSS / 1000 else NA_real_

wide$Boyce <- NA_real_
wide$MPA <- NA_real_

# --- Boyce computation with sentinel filter (carried from Stage 2) ----------
pred_df <- try(get_predictions(mods), silent = TRUE)
if (!inherits(pred_df, "try-error") && is.data.frame(pred_df) && "pred" %in% names(pred_df)) {
  vals_all <- as.numeric(pred_df$pred)
  if (max(vals_all, na.rm = TRUE) > 1.01) vals_all <- vals_all / 1000

  for (i in seq_len(nrow(wide))) {
    run_name <- as.character(wide$full.name[i])
    sub_pred <- pred_df[as.character(pred_df$full.name) == run_name, ]
    if (nrow(sub_pred) == 0) next
    vals <- as.numeric(sub_pred$pred)
    if (max(vals, na.rm = TRUE) > 1.01) vals <- vals / 1000
    obs <- as.numeric(sub_pred$points)
    obs_present <- df[[resp_col]][as.integer(sub_pred$points)]

    ok <- !is.na(vals) & !is.na(obs_present)
    if (sum(ok) < 5 || sum(obs_present[ok] == 1) < 3) next

    eb <- try(ecospat.boyce(fit = vals[ok],
                            obs = vals[ok][obs_present[ok] == 1],
                            PEplot = FALSE), silent = TRUE)
    if (!inherits(eb, "try-error")) {
      boyce_val <- NA_real_
      if (!is.null(eb$cor)) boyce_val <- as.numeric(eb$cor)
      else if (!is.null(eb$Spearman.cor)) boyce_val <- as.numeric(eb$Spearman.cor)
      auc_this <- if (!is.na(wide$AUC[i])) wide$AUC[i] else 0
      if (!is.na(boyce_val) && boyce_val == -1.0 && auc_this > 0.5) boyce_val <- NA_real_
      wide$Boyce[i] <- boyce_val
    }
    wide$MPA[i] <- mean(vals[ok] >= 0.5)
  }
}

names(wide)[names(wide) == "full.name"] <- "dataset"
names(wide)[names(wide) == "PA"] <- "pa"
wide$fit_status <- "ok"

# --- write per-run evaluations CSV (Stage 3 location, not Stage 2) ----------
eval_dir <- file.path(orig_wd, "data/processed/stage3/evaluations")
dir.create(eval_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(eval_dir, sprintf("%s__%s_esm_%s_evaluations.csv",
                                         slug, predictor_set, tolower(core)))
keep_cols <- intersect(c("species","species_slug","predictor_set","model_core",
                          "dataset","pa","run","algo","AUC","TSS","Boyce","Kappa",
                          "MPA","maxSSS","fit_status"), names(wide))
write.csv(wide[, keep_cols], out_csv, row.names = FALSE)
cat(sprintf("Wrote %s\n", out_csv))

# --- compute per-species aggregates and write metadata JSON -----------------
per_run <- wide[grepl("^RUN", as.character(wide$run)), ]
auc_vals <- per_run$AUC[!is.na(per_run$AUC)]
tss_vals <- per_run$TSS[!is.na(per_run$TSS)]
boyce_vals <- per_run$Boyce[!is.na(per_run$Boyce)]
maxsss_vals <- per_run$maxSSS[!is.na(per_run$maxSSS)]

auc_mean <- if (length(auc_vals)) mean(auc_vals) else NA
tss_mean <- if (length(tss_vals)) mean(tss_vals) else NA
qc_pass <- isTRUE(!is.na(auc_mean) && auc_mean > 0.7 &&
                  !is.na(tss_mean) && tss_mean > 0.4)
unreliable_boyce <- length(boyce_vals) < 3

meta <- list(
  species = species, species_slug = slug,
  predictor_set = predictor_set, core = core,
  modeling_id = modeling_id,
  fit_status = "ok",
  n_records = n_pres + n_abs,
  n_presences = n_pres, n_pseudo_absences = n_abs,
  n_predictors = length(predictors),
  predictors = predictors,
  n_cv_runs = nrow(per_run),
  auc_mean = auc_mean,
  auc_std = if (length(auc_vals) > 1) sd(auc_vals) else NA,
  tss_mean = tss_mean,
  tss_std = if (length(tss_vals) > 1) sd(tss_vals) else NA,
  boyce_mean = if (length(boyce_vals)) mean(boyce_vals) else NA,
  boyce_std = if (length(boyce_vals) > 1) sd(boyce_vals) else NA,
  boyce_n_valid = length(boyce_vals),
  maxsss_mean = if (length(maxsss_vals)) mean(maxsss_vals) else NA,
  maxsss_std = if (length(maxsss_vals) > 1) sd(maxsss_vals) else NA,
  qc_pass = qc_pass,
  unreliable_boyce = unreliable_boyce,
  biomod_models_path = file.path(workspace, slug)
)

out_json <- file.path(species_artifact_dir,
                      sprintf("%s__%s_metadata.json", predictor_set, tolower(core)))
write_json(meta, out_json, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat(sprintf("Wrote metadata: %s\n", out_json))
cat(sprintf("  AUC=%.3f  TSS=%.3f  Boyce=%s (n_valid=%d)  maxSSS=%.3f  qc_pass=%s\n",
            auc_mean, tss_mean,
            if (length(boyce_vals)) sprintf("%.3f", mean(boyce_vals)) else "NA",
            length(boyce_vals),
            if (length(maxsss_vals)) mean(maxsss_vals) else NA,
            qc_pass))