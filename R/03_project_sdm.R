# Stage 3b — Directional SDM projection
#
# Loads one trained biomod2 SDM for species A, applies it to species B's
# native occurrence env vectors for every target B in the target_species_file,
# writes a long-form parquet with columns:
#   target_species, target_cell_idx, subc_id, predicted_suitability
#
# Predictions are per-CV-run averaged (biomod2 returns 5-6 values per cell,
# one per RUN + the allRun summary; we average across RUNs, excluding allRun).
# Predictions come out on biomod2's 0-1000 scale and are rescaled to 0-1.
#
# Invoked by scripts/30_run_stage3b_projections.py.

suppressPackageStartupMessages({
  library(biomod2)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i[length(i)] + 1]
}

source_species  <- get_arg("--source-species")
predictor_set   <- get_arg("--predictor-set")
core            <- toupper(get_arg("--core", "GBM"))
workspace       <- get_arg("--workspace",         "data/processed/stage3/biomod_workspace")
target_csv      <- get_arg("--target-env-csv")     # long-form CSV with all target species' env vectors
output_parquet  <- get_arg("--output-parquet")
source_slug     <- get_arg("--source-slug")        # e.g. "astacus_astacus" (underscore form)

if (is.null(source_species) || is.null(predictor_set) || is.null(target_csv) ||
    is.null(output_parquet) || is.null(source_slug)) {
  stop("Required: --source-species --predictor-set --target-env-csv --output-parquet --source-slug")
}

# biomod2 uses dots, not underscores, in species folder names
dot_slug <- gsub("_", ".", source_slug)

# Locate the trained model .out file
modeling_id <- sprintf("%s__%s__esm_%s", source_slug, predictor_set, tolower(core))
models_out_path <- file.path(workspace, dot_slug,
                             sprintf("%s.%s.models.out", dot_slug, modeling_id))

if (!file.exists(models_out_path)) {
  stop(sprintf("Model .out file not found: %s", models_out_path))
}

cat(sprintf("Loading model: %s\n", basename(models_out_path)))

# biomod2 serializes via R's native save(). The .out file contains an object
# whose name is derived from modeling_id. Load into a fresh env and grab it.
load_env <- new.env()
loaded_names <- load(models_out_path, envir = load_env)
if (length(loaded_names) == 0) stop("Nothing loaded from .out file")
mods <- load_env[[loaded_names[1]]]

cat(sprintf("Model loaded. Class: %s\n", class(mods)[1]))
cat(sprintf("Expected predictors: %s\n",
            paste(mods@expl.var.names, collapse = ", ")))

# Read the target env CSV (all target species' presence env vectors, long form)
cat(sprintf("Reading target env vectors from: %s\n", target_csv))
targets <- read.csv(target_csv, stringsAsFactors = FALSE)
cat(sprintf("Target rows: %d across %d species\n",
            nrow(targets), length(unique(targets$target_species))))

# Extract the env predictor matrix in the order the model expects
predictor_cols <- mods@expl.var.names
missing_cols <- setdiff(predictor_cols, names(targets))
if (length(missing_cols) > 0) {
  stop(sprintf("Target CSV missing predictor columns: %s",
               paste(missing_cols, collapse = ", ")))
}

env_matrix <- as.data.frame(targets[, predictor_cols, drop = FALSE])

# Drop rows with any NA in predictors (can't project where env is undefined)
complete_rows <- complete.cases(env_matrix)
n_dropped <- sum(!complete_rows)
if (n_dropped > 0) {
  cat(sprintf("Dropping %d rows with NA in predictors (%d remaining)\n",
              n_dropped, sum(complete_rows)))
}
env_matrix_ok <- env_matrix[complete_rows, , drop = FALSE]
targets_ok    <- targets[complete_rows, , drop = FALSE]

# Predict via BIOMOD_Projection — the cleanest path that handles all CV runs
cat("Running BIOMOD_Projection...\n")
proj_id <- sprintf("proj_%s_to_targets", source_slug)

# Work inside the species' workspace folder so biomod2 writes its intermediate
# projection files in a contained location
orig_wd <- getwd()
setwd(file.path(workspace, dot_slug))
on.exit(setwd(orig_wd), add = TRUE)

proj <- tryCatch(
  BIOMOD_Projection(
    bm.mod          = mods,
    proj.name       = proj_id,
    new.env         = env_matrix_ok,
    models.chosen   = "all",
    metric.binary   = NULL,
    metric.filter   = NULL,
    build.clamping.mask = FALSE,
    nb.cpu          = 1,
    do.stack        = FALSE,
    keep.in.memory  = TRUE,
    output.format   = ".RData"
  ),
  error = function(e) {
    cat(sprintf("BIOMOD_Projection failed: %s\n", conditionMessage(e)))
    NULL
  }
)

setwd(orig_wd)

if (is.null(proj)) {
  stop("Projection failed, see error above")
}

# get_predictions returns a data.frame for tabular input
pred_df <- get_predictions(proj)
cat(sprintf("get_predictions returned: class=%s, nrow=%d\n",
            class(pred_df)[1], nrow(pred_df)))

# pred_df structure: one row per (input_point, model_run), columns include
# points/pred/full.name. We average pred across CV runs (excluding allRun).
pred_df$full.name <- as.character(pred_df$full.name)
is_cv_run <- grepl("_RUN[0-9]+_", pred_df$full.name)
pred_df_cv <- pred_df[is_cv_run, , drop = FALSE]

if (nrow(pred_df_cv) == 0) {
  stop("No CV-run predictions found (only allRun?). Cannot average.")
}

# Aggregate: mean pred per point across CV runs
# The 'points' column in biomod2's get_predictions corresponds to row index in env_matrix_ok
agg <- aggregate(pred ~ points, data = pred_df_cv, FUN = mean, na.rm = TRUE)
names(agg)[names(agg) == "pred"] <- "predicted_suitability"

# Rescale biomod2 0-1000 to 0-1
if (max(agg$predicted_suitability, na.rm = TRUE) > 1.01) {
  agg$predicted_suitability <- agg$predicted_suitability / 1000
}

# Join with target metadata (species, cell idx, subc_id)
targets_ok$points <- seq_len(nrow(targets_ok))
out <- merge(targets_ok, agg, by = "points", all.x = TRUE)

# Assemble final long-form output
result <- data.frame(
  source_species       = source_species,
  target_species       = out$target_species,
  target_cell_idx      = out$target_cell_idx,
  subc_id              = out$subc_id,
  predicted_suitability = out$predicted_suitability,
  predictor_set        = predictor_set,
  core                 = core,
  stringsAsFactors     = FALSE
)

# Write parquet
dir.create(dirname(output_parquet), showWarnings = FALSE, recursive = TRUE)
write_parquet(result, output_parquet)
cat(sprintf("Wrote %s (%d rows, %d target species, %d cells projected)\n",
            output_parquet, nrow(result),
            length(unique(result$target_species)),
            sum(!is.na(result$predicted_suitability))))