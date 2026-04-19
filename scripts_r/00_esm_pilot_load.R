library(ecospat)

infile <- "data/processed/stage2_pilot/astacus_astacus_native_training.csv"
df <- read.csv(infile)

cat("Loaded file:", infile, "\n")
cat("Rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n\n")

cat("Column names:\n")
print(colnames(df))

cat("\nSummary of predictors:\n")
preds <- c("l_CLI3","l_CLI23","l_CLI19","l_CLI15","l_CLI47","l_CLI59","l_TOP15","l_TOP122")
print(summary(df[, preds]))

cat("\nUnique basins:", length(unique(df$basin_id)), "\n")
cat("Unique segments:", length(unique(df$subc_id)), "\n")
