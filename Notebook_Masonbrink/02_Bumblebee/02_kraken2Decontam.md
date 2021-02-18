#  Run Kraken on reads to identify contamination


## Install Kraken2
```
#/work/gif/remkv6/Baum

ml minconda3; conda create -n kraken2
source activate kraken2
conda install -c bioconda kraken2
```

## Create Kraken2 database
```
#/work/gif/remkv6/Baum

ml miniconda3;source activate kraken2
kraken2-build --download-taxonomy --db BacteriaVirusArchae
kraken2-build --download-library viral --db BacteriaVirusArchaea/
kraken2-build --download-library bacteria --db BacteriaVirusArchaea/
kraken2-build --download-library archaea --db BacteriaVirusArchaea/
kraken2-build --download-library fungi --db BacteriaVirusArchaea/
kraken2-build --download-library protozoa --db BacteriaVirusArchaea/
kraken2-build --add-to-library BterrestrisTax.fa --db BacteriaVirusArchaea/
kraken2-build --build --db BacteriaVirusArchaea/ --threads 36

```

## Run reads through the kraken database
```
#/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/01_Align
paste <(ls -1 *R1*gz) |while read line; do echo "kraken2 -db /work/gif/remkv6/Baum/BacteriaVirusArchaea --threads 36 --report "$line".report --gzip-compressed  --unclassified-out "${line%.*}"unclassified.fq --classified-out "${line%.*}"classified.fq "$line" > "${line%.*}"Kraken.out" ;done  >kraken.sh
```



### Summarize the classification results
```

#taxons with 0.01% of reads or greater, and at least 100 reads atttributed to each taxon.
for f in *report; do echo "awk '\$1>0 && \$3>100 {print \$0,FILENAME}' "$f" >"$f".summary" ;done >summarizer.sh


#This prints the number of reads, the taxonomic name, and the sample name.
cat *report.summary |sed 's/_L002_R1_001.fastq.gz.report//g' |awk '{print $3,$6,$7}' |sed 's/1-/\t1-/1' |grep -v "unclassified" |grep -v "root" |grep -v "cellular" |grep -v "Eukaryota" |grep -v "Bombus" |grep -v "Opisthokonta" |grep -v "Bacteria" |sed 's/^[ ]*//g' >FinalContaminationNetwork.tab
```

### General Results

There are multiple taxons that are present in the bee data, though it does not appear to be related to a treatment.  

##### General abundance within the 60 samples.
```
 less FinalContaminationNetwork.tab |awk '{print $2}' |sort |uniq -c |sort -k1,1nr |less
```
<details>
  <summary>Abundance Summary</summary>
  <pre>
| 60 | Sarcocystidae |
| 60 | Sporisoriumgraminicola |
| 57 | Apicomplexa |
| 54 | Sar |
| 51 | leotiomyceta |
| 48 | AbyssogenaphaseoliformissymbiontOG214 |
| 34 | DictyosteliumdiscoideumAX4 |
| 22 | Dikarya |
| 14 | Plasmodiumvivax |
| 11 | Choristoneurafumiferanagranulovirus |
| 7 | Plasmodiumrelictum |
| 6 | Besnoitiabesnoiti |
| 6 | StaphylococcusphageAndhra |
| 6 | ToxoplasmagondiiME49 |
| 5 | Encephalitozoon |
| 5 | Plasmodium |
| 4 | ZymoseptoriatriticiIPO323 |
| 3 | Bacilli |
| 3 | Brevibacillus |
| 3 | Pseudomonastolaasii |
| 3 | Staphylococcus |
| 3 | Staphylococcusaureus |
| 2 | Enterobacteriaceae |
| 2 | Mycobacterium |
| 2 | Terrabacteriagroup |
| 1 | Actinomycetia |
| 1 | BabesiabovisT2Bo |
| 1 | Corynebacteriales |
| 1 | Enterobacterales |
| 1 | Lactobacillales |
| 1 | Latilactobacillus |
| 1 | Latilactobacilluscurvatus |
| 1 | Mycobacteriaceae |
| 1 | Mycobacteriumavium |
| 1 | Proteobacteria |
| 1 | Staphylococcaceae |
| 1 | Weissellaparamesenteroides |
| 1 | Zygosaccharomycesrouxii |
</details>
</pre>
