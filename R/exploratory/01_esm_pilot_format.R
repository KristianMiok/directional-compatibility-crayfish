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

cat("BIOMOD data object created successfully\n")
cat("Rows in source table:", nrow(df), "\n")
cat("Presences:", sum(df$resp == 1), "\n")
cat("Background:", sum(df$resp == 0), "\n")
print(myBiomodData)
