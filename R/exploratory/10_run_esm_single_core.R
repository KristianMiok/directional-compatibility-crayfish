args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript scripts_r/10_run_esm_single_core.R <csv_path> <GLM|GBM>")
}

csv_path <- args[1]
model_core <- toupper(args[2])

if (!model_core %in% c("GLM", "GBM")) {
  stop("model_core must be GLM or GBM")
}

library(ecospat)
library(biomod2)

df <- read.csv(csv_path)

meta_cols <- c("resp", "long_or", "lat_or", "subc_id", "basin_id", "Crayfish_scientific_name")
preds <- setdiff(colnames(df), meta_cols)

resp_name <- gsub("_biomod_pilot\\.csv$", "", basename(csv_path))
resp_name <- gsub("[^A-Za-z0-9_]", "_", resp_name)

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = as.numeric(df$resp),
  expl.var  = df[, preds],
  resp.xy   = df[, c("long_or", "lat_or")],
  resp.name = resp_name
)

cat("Starting ecospat.ESM.Modeling for", resp_name, "with", model_core, "...\n")
cat("Predictor columns:\n")
print(preds)

myESM <- ecospat.ESM.Modeling(
  data = myBiomodData,
  models = c(model_core),
  NbRunEval = 2,
  DataSplit = 70,
  weighting.score = c("AUC"),
  parallel = FALSE
)

esm_out <- file.path(dirname(csv_path), paste0(resp_name, "_esm_", tolower(model_core), ".rds"))
saveRDS(myESM, esm_out)

esm_ens <- ecospat.ESM.EnsembleModeling(
  ESM.modeling.output = myESM,
  weighting.score = c("AUC"),
  threshold = NULL,
  models = c(model_core)
)

ens_out <- file.path(dirname(csv_path), paste0(resp_name, "_esm_", tolower(model_core), "_ensemble.rds"))
saveRDS(esm_ens, ens_out)

eval_out <- file.path(dirname(csv_path), paste0(resp_name, "_esm_", tolower(model_core), "_evaluations.csv"))
write.csv(esm_ens$ESM.evaluations, eval_out, row.names = FALSE)

cat("\nSaved modeling object:", esm_out, "\n")
cat("Saved ensemble object:", ens_out, "\n")
cat("Saved evaluation table:", eval_out, "\n")
cat("\nESM evaluations:\n")
print(esm_ens$ESM.evaluations)
