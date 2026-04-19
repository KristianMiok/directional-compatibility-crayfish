library(ecospat)

esm <- readRDS("data/processed/stage2_pilot/astacus_astacus_esm_gbm.rds")

esm_ens <- ecospat.ESM.EnsembleModeling(
  ESM.modeling.output = esm,
  weighting.score = c("AUC"),
  threshold = NULL,
  models = c("GBM")
)

saveRDS(esm_ens, "data/processed/stage2_pilot/astacus_astacus_esm_gbm_ensemble.rds")

cat("Saved: data/processed/stage2_pilot/astacus_astacus_esm_gbm_ensemble.rds\n")
cat("\nTop-level names:\n")
print(names(esm_ens))

cat("\nESM evaluations:\n")
print(esm_ens$ESM.evaluations)
