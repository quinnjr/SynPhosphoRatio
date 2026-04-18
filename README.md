# SynPhosphoRatio

Combine per-sample *total* α-synuclein and *pS129-phosphorylated* α-synuclein
quantifications (typically produced by two runs of `SynProteinFilter` with
`form=total` and `form=ps129`) into the pS129/total ratio — a biomarker of
Lewy-body pathology in Parkinson's disease and related synucleinopathies.

**Research tool, not a clinical diagnostic.**

## Input

Two outputs from `SynProteinFilter` (same sample set, same order not required):

- `total` — single-feature sample×value CSV of total α-syn intensity
- `phospho` — single-feature sample×value CSV of pS129 α-syn intensity

## Parameters (`parameters.synphospho.txt`, tab-delimited)

| Key | Required | Default | Meaning |
|---|---|---|---|
| `total` | yes | — | Path to total α-syn CSV |
| `phospho` | yes | — | Path to pS129 α-syn CSV |
| `log_inputs` | no | `true` | Reverse the default `SynProteinFilter` log2 transform before dividing |
| `zero_guard` | no | `1e-6` | Minimum denominator to avoid divide-by-zero |

## Output (TSV)

```
sample    ps129_ratio    total    phospho
PD_001    0.278          1842000  512000
CTRL_001  0.050          820000   41000
```

Chain the output into ML classifiers (`RandomForest`, `LogReg`, `SVM`) or
`PlotROC` for case/control discrimination.
