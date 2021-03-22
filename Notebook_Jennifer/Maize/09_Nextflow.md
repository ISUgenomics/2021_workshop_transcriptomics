# 09 Nextflow

Wrapped [01_Quality-Control](01_Quality-Control), [03_GSNAP](03_GSNAP.md), [06_Kallisto](06_Kallisto.md), and [08_Salmon](08_Salmon.md) in a Nextflow script.

```
N E X T F L O W  ~  version 20.10.0
Launching `maize_run.nf` [naughty_keller] - revision: 6370a27351
executor >  slurm (94)
[0b/95fa5c] process > fastqc (batched)               [100%] 6 of 6 ✔
[58/608969] process > multiqc                        [100%] 1 of 1 ✔
[2e/cfe950] process > kallisto_index (GCF_9021671... [100%] 1 of 1 ✔
[f3/7ca926] process > kallisto_quant (SRR1573504)    [100%] 24 of 24 ✔
[6c/12e9b0] process > salmon_index (GCF_902167145... [100%] 1 of 1 ✔
[05/9679fc] process > salmon_quant (SRR1573504)      [100%] 24 of 24 ✔
[7d/8f8217] process > gsnap_index (GCF_902167145)    [100%] 1 of 1 ✔
[46/12c26d] process > gsnap_quant (SRR1573510)       [100%] 5 of 5, failed: 1
[84/37bae6] process > featureCounts_gene (SRR1573... [100%] 3 of 3
[7c/abf2e1] process > featureCounts_mRNA (SRR1573... [100%] 3 of 3
[50/1af38a] process > featureCounts_geneMult (SRR... [100%] 3 of 3
Error executing process > 'gsnap_quant (SRR1573511)'

```

Ran pretty quickly over a weekend on Atlas HPC, ~X minutes.

![](results/assets/screenshot_maize_timeline.png)

View entire [bee_timeline.html](results/maize_timeline.html) and scroll down to see GSNAP and salmon run times. Kallisto was drastically faster than the other two.

To learn Nextflow, visit [https://www.nextflow.io/](https://www.nextflow.io/)
