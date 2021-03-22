# 01 Quality Control

(1) Run fastqc on each read file

```
module load fastqc parallel
PROC=16
parallel -j${PROC} "fastqc -o 01_Quality-Control {1}" ::: 00_Raw-Data/*.fastq.gz
```

(2) Merge into one html report

```
cd 01_Quality-Control
multiqc .
```

(3) Look at [multiqc_report.html](results/multiqc_report.html)

![](results/assets/screenshot_multiqc.png)

