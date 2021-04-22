# 09 Nextflow

Wrapped [01_Quality-Control](01_Quality-Control), [03_GSNAP](03_GSNAP.md), [06_Kallisto](06_Kallisto.md), and [08_Salmon](08_Salmon.md) in a Nextflow script.

```
N E X T F L O W  ~  version 20.10.0
Launching `maize_run.nf` [naughty_keller] - revision: 6370a27351
executor >  slurm (154)
[65/9f5b8d] process > fastqc (batched)               [100%] 6 of 6 ✔
[a6/f76504] process > multiqc                        [100%] 1 of 1 ✔
[3f/cb5455] process > kallisto_index (GCF_9021671... [100%] 1 of 1 ✔
[96/940514] process > kallisto_quant (SRR1573504)    [100%] 24 of 24 ✔
[2d/ab4e35] process > salmon_index (GCF_902167145... [100%] 1 of 1 ✔
[b5/260d0a] process > salmon_quant (SRR1573509)      [100%] 24 of 24 ✔
[42/d3600c] process > gsnap_index (GCF_902167145)    [100%] 1 of 1 ✔
[17/b68c59] process > gsnap_quant (SRR1573511)       [100%] 24 of 24 ✔
[51/de2e96] process > featureCounts_gene (SRR1573... [100%] 24 of 24 ✔
[5b/3c04f4] process > featureCounts_mRNA (SRR1573... [100%] 24 of 24 ✔
[7c/7ccc62] process > featureCounts_geneMult (SRR... [100%] 24 of 24 ✔
Completed at: 22-Mar-2021 04:51:51
Duration    : 4h 20m 7s
CPU hours   : 71.3
Succeeded   : 154
```

Total pipeline ran in 4 hours and 20 minutes on Atlas HPC.

![](results/assets/screenshot_maize_timeline.png)

View entire [maize_timeline.html](results/maize_timeline.html) and scroll down to see GSNAP and salmon run times. Kallisto and salmon was drastically faster than GSNAP.

To learn Nextflow, visit [https://www.nextflow.io/](https://www.nextflow.io/)
