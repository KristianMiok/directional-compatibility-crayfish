# Stage 2 pilot summary — Astacus astacus

A pilot ESM workflow was tested successfully on native-only, filtered, deduplicated records for *Astacus astacus*.

## Mean evaluation metrics

| model_family | AUC | TSS | Boyce | Kappa |
|---|---:|---:|---:|---:|
| GBM | 0.9850 | 0.9030 | 0.9925 | 0.8590 |
| GLM | 0.9820 | 0.8905 | 0.9930 | 0.8440 |

## Interpretation

Both ESM cores ran successfully:
- GLM-based ESM
- GBM-based ESM

GBM performed slightly better on this pilot, but both approaches gave very strong results.
This supports the feasibility of the dual-core ESM strategy proposed for Stage 2.