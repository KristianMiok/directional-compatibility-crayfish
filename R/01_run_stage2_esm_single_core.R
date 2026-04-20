suppressPackageStartupMessages({
  library(ecospat)
  library(biomod2)
  library(terra)
  library(dismo)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i[length(i)] + 1]
}

species <- get_arg("--species")
predictor_set <- get_arg("--predictor-set")
core <- toupper(get_arg("--core", "GLM"))
input_csv <- get_arg("--input-csv")
out_dir <- get_arg("--output-dir", "data/processed/stage2_pilot")
pa_reps <- as.integer(get_arg("--pa-reps", "5"))
cv_reps <- as.integer(get_arg("--cv-reps", "5"))
seed <- as.integer(get_arg("--seed", "42"))

if (is.null(species) || is.null(predictor_set)) {
  stop("Required: --species and --predictor-set")
}
if (!(core %in% c("GLM", "GBM"))) {
  stop("--core must be GLM or GBM")
}

slug <- gsub("[^a-z0-9]+", "_", tolower(species))
if (is.null(input_csv)) {
  input_csv <- file.path("data/processed/stage2_pilot", sprintf("%s__%s_biomod_pilot.csv", slug, predictor_set))
}
if (!file.exists(input_csv)) {
  stop(sprintf("Input CSV not found: %s", input_csv))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_csv, check.names = FALSE)

presence_candidates <- c("presence", "resp", "occ", "occurrence", "pa", "Presence", "pres_pa", "response")
presence_col <- presence_candidates[presence_candidates %in% names(dat)][1]
if (is.na(presence_col)) stop(sprintf("No presence column found. Available columns: %s", paste(names(dat), collapse = ", ")))

lon_candidates <- c("lon", "x", "X_WGS84_DD", "longitude", "Longitude", "decimalLongitude", "long_or", "LONG_OR")
lat_candidates <- c("lat", "y", "Y_WGS84_DD", "latitude", "Latitude", "decimalLatitude", "lat_or", "LAT_OR")
lon_col <- lon_candidates[lon_candidates %in% names(dat)][1]
lat_col <- lat_candidates[lat_candidates %in% names(dat)][1]
if (is.na(lon_col) || is.na(lat_col)) stop(sprintf("Could not detect coordinate columns. Available columns: %s", paste(names(dat), collapse = ", ")))

pred_cols <- setdiff(names(dat), c(
  "species","species_slug","predictor_set",
  presence_candidates,
  lon_candidates, lat_candidates,
  "CellID","FID","cell_id","fid","subc_id","basin_id","Crayfish_scientific_name"
))
if (length(pred_cols) < 2) stop("Need at least 2 predictor columns")

resp <- dat[[presence_col]]
xy <- as.matrix(dat[, c(lon_col, lat_col)])
env <- dat[, pred_cols, drop = FALSE]

cat(sprintf("Using presence column: %s\n", presence_col))
cat(sprintf("Using longitude column: %s\n", lon_col))
cat(sprintf("Using latitude column: %s\n", lat_col))
cat(sprintf("Predictor count: %d\n", length(pred_cols)))

has_true_absences <- any(resp == 0, na.rm = TRUE)
has_presences <- any(resp == 1, na.rm = TRUE)

if (!has_presences) stop("No presences found in response column")
if (!has_true_absences) {
  biomod_data <- BIOMOD_FormatingData(
    resp.var = resp,
    expl.var = env,
    resp.xy = xy,
    resp.name = slug,
    PA.nb.rep = pa_reps,
    PA.nb.absences = sum(resp == 1, na.rm = TRUE),
    PA.strategy = "random",
    na.rm = TRUE
  )
} else {
  cat("Detected true absences in input; skipping pseudo-absence generation\n")
  biomod_data <- BIOMOD_FormatingData(
    resp.var = resp,
    expl.var = env,
    resp.xy = xy,
    resp.name = slug,
    na.rm = TRUE
  )
}

modeling_id <- sprintf("%s__%s__esm_%s", slug, predictor_set, tolower(core))

write_empty_result <- function(reason) {
  empty <- data.frame(
    species = species,
    species_slug = slug,
    predictor_set = predictor_set,
    model_core = core,
    dataset = NA,
    pa = NA,
    run = NA,
    algo = core,
    AUC = NA_real_,
    TSS = NA_real_,
    Boyce = NA_real_,
    Kappa = NA_real_,
    MPA = NA_real_,
    fit_status = reason,
    stringsAsFactors = FALSE
  )
  outfile_local <- file.path(out_dir, sprintf("%s__%s_esm_%s_evaluations.csv", slug, predictor_set, tolower(core)))
  write.csv(empty, outfile_local, row.names = FALSE)
  cat(sprintf("Wrote %s (status: %s)\n", outfile_local, reason))
}

opt_user_val <- NULL
if (core == "GBM") {
  gbm_params <- list(
    n.trees = 500,
    interaction.depth = 2,
    shrinkage = 0.01,
    n.minobsinnode = 3,
    bag.fraction = 0.75,
    train.fraction = 1.0,
    cv.folds = 0,
    keep.data = TRUE,
    verbose = FALSE
  )
  opt_user_val <- list(GBM.binary.gbm.gbm = list(for_all_datasets = gbm_params))
}

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
  write_empty_result("fit_error_biomod_modeling")
  quit(status = 0)
}

eval_df <- tryCatch(get_evaluations(mods), error = function(e) {
  cat(sprintf("get_evaluations failed: %s\n", conditionMessage(e)))
  NULL
})

if (is.null(eval_df) || !is.data.frame(eval_df) || nrow(eval_df) == 0) {
  write_empty_result("no_evaluations_returned")
  quit(status = 0)
}

cat("\nEvaluation data frame columns:\n"); print(names(eval_df))
cat("First few rows:\n"); print(head(eval_df, 3))

metric_map <- c(AUCroc = "AUC", ROC = "AUC", TSS = "TSS", KAPPA = "Kappa")
eval_df$metric <- ifelse(
  as.character(eval_df$metric.eval) %in% names(metric_map),
  metric_map[as.character(eval_df$metric.eval)],
  as.character(eval_df$metric.eval)
)

wide <- reshape(
  eval_df[, c("full.name", "PA", "run", "algo", "metric", "validation")],
  idvar = c("full.name", "PA", "run", "algo"),
  timevar = "metric",
  direction = "wide"
)
names(wide) <- sub("^validation\\.", "", names(wide))
wide$species <- species
wide$species_slug <- slug
wide$model_core <- core
wide$predictor_set <- predictor_set

if (!("AUC" %in% names(wide))) wide$AUC <- NA_real_
if (!("TSS" %in% names(wide))) wide$TSS <- NA_real_
if (!("Kappa" %in% names(wide))) wide$Kappa <- NA_real_

wide$Boyce <- NA_real_
wide$MPA <- NA_real_

pred_df <- try(get_predictions(mods), silent = TRUE)
if (!inherits(pred_df, "try-error") && is.data.frame(pred_df) && "pred" %in% names(pred_df)) {
  cat("Predictions data frame columns:\n"); print(names(pred_df))
  vals_all <- as.numeric(pred_df$pred)
  if (max(vals_all, na.rm = TRUE) > 1.01) vals_all <- vals_all / 1000

  for (i in seq_len(nrow(wide))) {
    run_name <- as.character(wide$full.name[i])
    sub_pred <- pred_df[as.character(pred_df$full.name) == run_name, , drop = FALSE]
    if (nrow(sub_pred) == 0) next
    vals <- as.numeric(sub_pred$pred)
    if (max(vals, na.rm = TRUE) > 1.01) vals <- vals / 1000

    if ("points" %in% names(sub_pred)) {
      pts <- as.integer(sub_pred$points)
      obs <- resp[pts]
    } else if (length(vals) == length(resp)) {
      obs <- resp
    } else {
      next
    }
    ok <- is.finite(vals) & is.finite(obs)
    if (sum(ok) < 10 || length(unique(obs[ok])) < 2) next

    eb <- try(ecospat.boyce(fit = vals[ok],
                            obs = vals[ok][obs[ok] == 1],
                            PEplot = FALSE), silent = TRUE)
    if (!inherits(eb, "try-error")) {
      boyce_val <- NA_real_
      if (!is.null(eb$cor)) boyce_val <- as.numeric(eb$cor)
      else if (!is.null(eb$Spearman.cor)) boyce_val <- as.numeric(eb$Spearman.cor)
      # Treat ecospat.boyce sentinel returns as NA:
      #   - NaN: classification cell had zero variance
      #   - exactly -1.0 co-occurring with AUC > 0.5: known ecospat degeneracy
      #     (all obs in one pred-suitability class, monotone F-ratio by chance)
      auc_this <- if (!is.na(wide$AUC[i])) wide$AUC[i] else 0
      if (!is.na(boyce_val) && boyce_val == -1.0 && auc_this > 0.5) boyce_val <- NA_real_
      wide$Boyce[i] <- boyce_val
    }
    wide$MPA[i] <- mean(vals[ok] >= 0.5)
  }
}

names(wide)[names(wide) == "full.name"] <- "dataset"
names(wide)[names(wide) == "PA"] <- "pa"

outfile <- file.path(out_dir, sprintf("%s__%s_esm_%s_evaluations.csv", slug, predictor_set, tolower(core)))
wide$fit_status <- "ok"
write.csv(wide[, intersect(c("species","species_slug","predictor_set","model_core","dataset","pa","run","algo","AUC","TSS","Boyce","Kappa","MPA","fit_status"), names(wide))], outfile, row.names = FALSE)
cat(sprintf("Wrote %s\n", outfile))
