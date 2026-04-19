library(ecospat)
library(biomod2)

infile <- "data/processed/stage2_pilot/astacus_astacus_biomod_pilot.csv"
df <- read.csv(infile)

preds <- c("l_CLI3","l_CLI23","l_CLI19","l_CLI15","l_CLI47","l_CLI59","l_TOP15","l_TOP122")

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = as.numeric(df$resp),
  expl.var  = df[, preds],
  resp.xy   = df[, c("long_or", "lat_or")],
  resp.name = "Astacus_astacus"
)

cat("Starting ecospat.ESM.Modeling with GBM only...\n")

myESM <- ecospat.ESM.Modeling(
  data = myBiomodData,
  models = c("GBM"),
  NbRunEval = 2,
  DataSplit = 70,
  weighting.score = c("AUC"),
  parallel = FALSE
)

saveRDS(myESM, "data/processed/stage2_pilot/astacus_astacus_esm_gbm.rds")

cat("\nESM GBM pilot finished\n")
cat("Saved: data/processed/stage2_pilot/astacus_astacus_esm_gbm.rds\n")

cat("\nTop-level names in object:\n")
print(names(myESM))

if (!is.null(myESM$failed)) {
  cat("\nFailed bivariate models:\n")
  print(myESM$failed)
}
