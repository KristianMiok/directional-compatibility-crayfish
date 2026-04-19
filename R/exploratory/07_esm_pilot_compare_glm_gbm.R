library(ecospat)

glm_ens <- readRDS("data/processed/stage2_pilot/astacus_astacus_esm_glm_ensemble.rds")
gbm_ens <- readRDS("data/processed/stage2_pilot/astacus_astacus_esm_gbm_ensemble.rds")

glm_eval <- glm_ens$ESM.evaluations
glm_eval$model_family <- "GLM"

gbm_eval <- gbm_ens$ESM.evaluations
gbm_eval$model_family <- "GBM"

cmp <- rbind(glm_eval, gbm_eval)

write.csv(cmp, "data/processed/stage2_pilot/astacus_astacus_esm_glm_vs_gbm.csv", row.names = FALSE)

cat("Saved: data/processed/stage2_pilot/astacus_astacus_esm_glm_vs_gbm.csv\n\n")
print(cmp)

cat("\nMean metrics by family:\n")
print(aggregate(cbind(AUC, TSS, Boyce, Kappa) ~ model_family, data = cmp, FUN = mean))
