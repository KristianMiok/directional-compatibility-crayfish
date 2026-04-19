library(ecospat)

esm <- readRDS("data/processed/stage2_pilot/astacus_astacus_esm_glm.rds")
esm_ens <- readRDS("data/processed/stage2_pilot/astacus_astacus_esm_glm_ensemble.rds")

cat("=== Modeling object ===\n")
print(names(esm))

cat("\n=== Ensemble object ===\n")
print(names(esm_ens))

cat("\n=== Ensemble evaluations ===\n")
print(esm_ens$ESM.evaluations)

cat("\n=== Weight summary ===\n")
if (!is.null(esm_ens$weights)) {
  print(summary(esm_ens$weights))
} else {
  cat("No 'weights' slot found\n")
}

cat("\n=== Failed bivariate models ===\n")
if (!is.null(esm_ens$failed)) {
  print(esm_ens$failed)
} else {
  cat("No failed slot found in ensemble object\n")
}
