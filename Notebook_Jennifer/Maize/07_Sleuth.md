# 07 Sleuth
Will need a `metadata.tsv` file similar to the following. Give the path to kallisto output folders.

| sample | condition | path |
| :--|:--|:--|
|sample1 | control | kallisto/samp1_quant/ |
|sample2 | control | kallisto/samp2_quant/ |
|sample3 | treatment | kallisto/samp3_quant/ |
|sample4 | treatment | kallisto/samp4_quant/ |

```
library(sleuth)
library(tidyverse)

# Input: sample to condition file
s2c <- readr::read_delim("metadata.tsv", delim="\t")

# (1/4) Initialize the Sleuth Object
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# (2/4) Fitting a full model
so_full <- sleuth_fit(so, ~condition, 'full')

# (3/4) Fitting a reduced model
so_reduced <- sleuth_fit(so_full, ~1, 'reduced')

# (4/4) lrt test
so_lrt <- sleuth_lrt(so_reduced, 'reduced', 'full')

# Output: significant gene table
sleuth_table <- sleuth_results(so_lrt, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write_delim(sleuth_significant, "bee_Sleuth_results.tsv", delim="\t")
```