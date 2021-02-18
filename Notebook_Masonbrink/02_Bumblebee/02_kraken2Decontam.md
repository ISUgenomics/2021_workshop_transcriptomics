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


##### Read content identified as Bombus terrestris

```
grep "terrestris" *.report|awk '$5=="S1"' |awk '$2>0' |sed 's/_/\t/2'  |cut -f 1,3- |sed 's/[ ]*//g' |sed 's/terrestris/ terrestris/g' |less
```
<details>
<p>

| Sample   | Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|----------|----------|----|---------|-------------------------------|
| 1-A01-A1_S7   | 4538578  | 4538578  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A02-A2_S8   | 2910840  | 2910840  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A03-A3_S9   | 2952334  | 2952334  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A04-A4_S10  | 3186146  | 3186146  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A05-A5_S11  | 2645088  | 2645088  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A06-A6_S12  | 3646012  | 3646012  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A07-A7_S13  | 4533142  | 4533142  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A08-A8_S14  | 4631837  | 4631837  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A09-A9_S15  | 4798282  | 4798282  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A10-A10_S16 | 4080344  | 4080344  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A11-A11_S17 | 4131962  | 4131962  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-A12-A12_S18 | 3995787  | 3995787  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B01-A13_S19 | 4703231  | 4703231  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B02-A14_S20 | 2702453  | 2702453  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B03-A16_S21 | 15673979 | 15673979 | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B04-A17_S22 | 2453892  | 2453892  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B05-A18_S23 | 2032239  | 2032239  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B06-A19_S24 | 4071425  | 4071425  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B07-B2_S25  | 4556365  | 4556365  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B08-B3_S26  | 4186330  | 4186330  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B09-B4_S27  | 3163793  | 3163793  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B10-B5_S28  | 3780311  | 3780311  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B11-B6_S29  | 1421950  | 1421950  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-B12-B7_S30  | 5409627  | 5409627  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C01-B8_S31  | 4288291  | 4288291  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C02-B10_S32 | 3136857  | 3136857  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C03-C1_S33  | 2804020  | 2804020  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C04-C2_S34  | 2124400  | 2124400  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C05-C3_S35  | 3376622  | 3376622  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C06-C4_S36  | 3840422  | 3840422  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C07-C6_S37  | 2802287  | 2802287  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C08-C7_S38  | 3546076  | 3546076  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C09-C8_S39  | 3439938  | 3439938  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C10-C9_S40  | 4165227  | 4165227  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C11-D1_S41  | 4166471  | 4166471  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-C12-D2_S42  | 3316185  | 3316185  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D01-D5_S43  | 3862527  | 3862527  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D02-D6_S44  | 2600531  | 2600531  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D03-D7_S45  | 2846932  | 2846932  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D04-D8_S46  | 1990119  | 1990119  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D05-D9_S47  | 2865459  | 2865459  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D06-D10_S48 | 2763776  | 2763776  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D07-D11_S49 | 5864803  | 5864803  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D08-E1_S50  | 3020576  | 3020576  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D09-E2_S51  | 3288241  | 3288241  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D10-E5_S52  | 2963740  | 2963740  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D11-E7_S53  | 2077059  | 2077059  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-D12-E8_S54  | 4606376  | 4606376  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E01-E9_S55  | 2886014  | 2886014  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E02-E10_S56 | 4003479  | 4003479  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E03-F1_S57  | 4512959  | 4512959  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E04-F2_S58  | 6787519  | 6787519  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E05-F3_S59  | 2057803  | 2057803  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E06-F4_S60  | 2237088  | 2237088  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E07-F5_S61  | 273212   | 273212   | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E08-F6_S62  | 3612650  | 3612650  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E09-F7_S63  | 3289933  | 3289933  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E10-F8_S64  | 3017514  | 3017514  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E11-F9_S65  | 3344462  | 3344462  | S1 | 1255232 | Bombus terrestris terrestris |
| 1-E12-F10_S66 | 3567365  | 3567365  | S1 | 1255232 | Bombus terrestris terrestris |

</p>
</details>


###### Read content unknown
```
grep "unclassified" *.report|awk '$5=="U"' |awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |less
```
<details>
<p>

| Sample   | Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|----------|----------|----|---------|-------------------------------|
| 1-A01-A1_S7   | 24.84 | 1516516 | 1516516 | U | 0 | unclassified |
| 1-A02-A2_S8   | 22.45 | 853087  | 853087  | U | 0 | unclassified |
| 1-A03-A3_S9   | 20.02 | 745608  | 745608  | U | 0 | unclassified |
| 1-A04-A4_S10  | 29.52 | 1359534 | 1359534 | U | 0 | unclassified |
| 1-A05-A5_S11  | 23.69 | 831482  | 831482  | U | 0 | unclassified |
| 1-A06-A6_S12  | 25.02 | 1230674 | 1230674 | U | 0 | unclassified |
| 1-A07-A7_S13  | 22.73 | 1345876 | 1345876 | U | 0 | unclassified |
| 1-A08-A8_S14  | 26.09 | 1657168 | 1657168 | U | 0 | unclassified |
| 1-A09-A9_S15  | 20.29 | 1233397 | 1233397 | U | 0 | unclassified |
| 1-A10-A10_S16 | 22.54 | 1203840 | 1203840 | U | 0 | unclassified |
| 1-A11-A11_S17 | 24.49 | 1353557 | 1353557 | U | 0 | unclassified |
| 1-A12-A12_S18 | 26.80 | 1478075 | 1478075 | U | 0 | unclassified |
| 1-B01-A13_S19 | 26.56 | 1728738 | 1728738 | U | 0 | unclassified |
| 1-B02-A14_S20 | 29.18 | 1137401 | 1137401 | U | 0 | unclassified |
| 1-B03-A16_S21 | 30.92 | 7211713 | 7211713 | U | 0 | unclassified |
| 1-B04-A17_S22 | 30.03 | 1080381 | 1080381 | U | 0 | unclassified |
| 1-B05-A18_S23 | 29.51 | 871026  | 871026  | U | 0 | unclassified |
| 1-B06-A19_S24 | 25.68 | 1430203 | 1430203 | U | 0 | unclassified |
| 1-B07-B2_S25  | 25.35 | 1573941 | 1573941 | U | 0 | unclassified |
| 1-B08-B3_S26  | 27.82 | 1644335 | 1644335 | U | 0 | unclassified |
| 1-B09-B4_S27  | 22.66 | 940452  | 940452  | U | 0 | unclassified |
| 1-B10-B5_S28  | 27.41 | 1454053 | 1454053 | U | 0 | unclassified |
| 1-B11-B6_S29  | 42.15 | 1072356 | 1072356 | U | 0 | unclassified |
| 1-B12-B7_S30  | 22.86 | 1617682 | 1617682 | U | 0 | unclassified |
| 1-C01-B8_S31  | 18.25 | 962481  | 962481  | U | 0 | unclassified |
| 1-C02-B10_S32 | 25.43 | 1081827 | 1081827 | U | 0 | unclassified |
| 1-C03-C1_S33  | 31.48 | 1325328 | 1325328 | U | 0 | unclassified |
| 1-C04-C2_S34  | 26.83 | 793783  | 793783  | U | 0 | unclassified |
| 1-C05-C3_S35  | 18.94 | 795132  | 795132  | U | 0 | unclassified |
| 1-C06-C4_S36  | 18.96 | 906997  | 906997  | U | 0 | unclassified |
| 1-C07-C6_S37  | 31.46 | 1318571 | 1318571 | U | 0 | unclassified |
| 1-C08-C7_S38  | 27.25 | 1346398 | 1346398 | U | 0 | unclassified |
| 1-C09-C8_S39  | 26.00 | 1226647 | 1226647 | U | 0 | unclassified |
| 1-C10-C9_S40  | 27.96 | 1639424 | 1639424 | U | 0 | unclassified |
| 1-C11-D1_S41  | 22.58 | 1229723 | 1229723 | U | 0 | unclassified |
| 1-C12-D2_S42  | 24.67 | 1096310 | 1096310 | U | 0 | unclassified |
| 1-D01-D5_S43  | 30.86 | 1755718 | 1755718 | U | 0 | unclassified |
| 1-D02-D6_S44  | 26.30 | 942340  | 942340  | U | 0 | unclassified |
| 1-D03-D7_S45  | 37.81 | 1816044 | 1816044 | U | 0 | unclassified |
| 1-D04-D8_S46  | 37.11 | 1212096 | 1212096 | U | 0 | unclassified |
| 1-D05-D9_S47  | 29.46 | 1227819 | 1227819 | U | 0 | unclassified |
| 1-D06-D10_S48 | 34.70 | 1509722 | 1509722 | U | 0 | unclassified |
| 1-D07-D11_S49 | 26.60 | 2166354 | 2166354 | U | 0 | unclassified |
| 1-D08-E1_S50  | 35.46 | 1697129 | 1697129 | U | 0 | unclassified |
| 1-D09-E2_S51  | 27.77 | 1285615 | 1285615 | U | 0 | unclassified |
| 1-D10-E5_S52  | 35.50 | 1682169 | 1682169 | U | 0 | unclassified |
| 1-D11-E7_S53  | 22.45 | 607659  | 607659  | U | 0 | unclassified |
| 1-D12-E8_S54  | 19.49 | 1123230 | 1123230 | U | 0 | unclassified |
| 1-E01-E9_S55  | 21.46 | 797977  | 797977  | U | 0 | unclassified |
| 1-E02-E10_S56 | 30.78 | 1807438 | 1807438 | U | 0 | unclassified |
| 1-E03-F1_S57  | 21.83 | 1274251 | 1274251 | U | 0 | unclassified |
| 1-E04-F2_S58  | 30.93 | 3101800 | 3101800 | U | 0 | unclassified |
| 1-E05-F3_S59  | 23.89 | 655411  | 655411  | U | 0 | unclassified |
| 1-E06-F4_S60  | 32.32 | 1091287 | 1091287 | U | 0 | unclassified |
| 1-E07-F5_S61  | 43.26 | 224512  | 224512  | U | 0 | unclassified |
| 1-E08-F6_S62  | 26.64 | 1332765 | 1332765 | U | 0 | unclassified |
| 1-E09-F7_S63  | 25.58 | 1145841 | 1145841 | U | 0 | unclassified |
| 1-E10-F8_S64  | 28.02 | 1201251 | 1201251 | U | 0 | unclassified |
| 1-E11-F9_S65  | 25.83 | 1187496 | 1187496 | U | 0 | unclassified |
| 1-E12-F10_S66 | 29.91 | 1556936 | 1556936 | U | 0 | unclassified |

</p>
</details>

##### General contaminant abundance within the 60 samples.
```
 less FinalContaminationNetwork.tab |awk '{print $2}' |sort |uniq -c |sort -k1,1nr |less
```
<details>
<p>

| Samples | Taxonomic Group                         |
|----|---------------------------------------|
| 60 | Sarcocystidae                         |
| 60 | Sporisorium graminicola                |
| 57 | Apicomplexa                           |
| 54 | Sar                                   |
| 51 | leotiomyceta                          |
| 48 | Abyssogena phaseoliformis symbiont OG214 |
| 34 | Dictyostelium discoideum AX4            |
| 22 | Dikarya                               |
| 14 | Plasmodium vivax                       |
| 11 | Choristoneura fumiferana granulovirus   |
| 7  | Plasmodium relictum                    |
| 6  | Besnoitiabes noiti                     |
| 6  | Staphylococcus phage Andhra             |
| 6  | Toxoplasma gondii ME49                  |
| 5  | Encephalitozoon                       |
| 5  | Plasmodium                            |
| 4  | ZymoseptoriatriticiIPO323             |
| 3  | Bacilli                               |
| 3  | Brevibacillus                         |
| 3  | Pseudomonastolaasii                   |
| 3  | Staphylococcus                        |
| 3  | Staphylococcusaureus                  |
| 2  | Enterobacteriaceae                    |
| 2  | Mycobacterium                         |
| 2  | Terrabacteriagroup                    |
| 1  | Actinomycetia                         |
| 1  | BabesiabovisT2Bo                      |
| 1  | Corynebacteriales                     |
| 1  | Enterobacterales                      |
| 1  | Lactobacillales                       |
| 1  | Latilactobacillus                     |
| 1  | Latilactobacilluscurvatus             |
| 1  | Mycobacteriaceae                      |
| 1  | Mycobacteriumavium                    |
| 1  | Proteobacteria                        |
| 1  | Staphylococcaceae                     |
| 1  | Weissellaparamesenteroides            |
| 1  | Zygosaccharomycesrouxii               |

</p>
</details>


###### Sarcoystidae Read Counts
```
 grep "Sarcocystidae" *.report|awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |less  
```

This is a family of cyst-forming apicomplexan protozoan parasites. Kraken could not pinpoint a species.
<details>
<p>

| Sample   | Proportion of reads |Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|------|--------|--------|---|------|---------------|
| 1-A01-A1_S7   | 0.23 | 13849  | 13493  | F | 5809 | Sarcocystidae |
| 1-A02-A2_S8   | 0.37 | 13951  | 13718  | F | 5809 | Sarcocystidae |
| 1-A03-A3_S9   | 0.24 | 8960   | 8759   | F | 5809 | Sarcocystidae |
| 1-A04-A4_S10  | 0.44 | 20374  | 20033  | F | 5809 | Sarcocystidae |
| 1-A05-A5_S11  | 0.39 | 13611  | 13427  | F | 5809 | Sarcocystidae |
| 1-A06-A6_S12  | 0.28 | 13809  | 13559  | F | 5809 | Sarcocystidae |
| 1-A07-A7_S13  | 0.23 | 13331  | 13081  | F | 5809 | Sarcocystidae |
| 1-A08-A8_S14  | 0.34 | 21361  | 21035  | F | 5809 | Sarcocystidae |
| 1-A09-A9_S15  | 0.27 | 16270  | 15955  | F | 5809 | Sarcocystidae |
| 1-A10-A10_S16 | 0.40 | 21608  | 21221  | F | 5809 | Sarcocystidae |
| 1-A11-A11_S17 | 0.31 | 17290  | 17099  | F | 5809 | Sarcocystidae |
| 1-A12-A12_S18 | 0.16 | 8641   | 8386   | F | 5809 | Sarcocystidae |
| 1-B01-A13_S19 | 0.28 | 18252  | 17908  | F | 5809 | Sarcocystidae |
| 1-B02-A14_S20 | 0.56 | 21783  | 21513  | F | 5809 | Sarcocystidae |
| 1-B03-A16_S21 | 0.74 | 172800 | 170695 | F | 5809 | Sarcocystidae |
| 1-B04-A17_S22 | 0.71 | 25563  | 25238  | F | 5809 | Sarcocystidae |
| 1-B05-A18_S23 | 0.70 | 20774  | 20589  | F | 5809 | Sarcocystidae |
| 1-B06-A19_S24 | 0.50 | 27679  | 27354  | F | 5809 | Sarcocystidae |
| 1-B07-B2_S25  | 0.52 | 32577  | 32322  | F | 5809 | Sarcocystidae |
| 1-B08-B3_S26  | 0.58 | 34108  | 33739  | F | 5809 | Sarcocystidae |
| 1-B09-B4_S27  | 0.43 | 17780  | 17519  | F | 5809 | Sarcocystidae |
| 1-B10-B5_S28  | 0.57 | 30463  | 30228  | F | 5809 | Sarcocystidae |
| 1-B11-B6_S29  | 0.03 | 732    | 525    | F | 5809 | Sarcocystidae |
| 1-B12-B7_S30  | 0.26 | 18475  | 18062  | F | 5809 | Sarcocystidae |
| 1-C01-B8_S31  | 0.19 | 9784   | 9651   | F | 5809 | Sarcocystidae |
| 1-C02-B10_S32 | 0.29 | 12516  | 12332  | F | 5809 | Sarcocystidae |
| 1-C03-C1_S33  | 0.76 | 31824  | 31382  | F | 5809 | Sarcocystidae |
| 1-C04-C2_S34  | 0.56 | 16570  | 16338  | F | 5809 | Sarcocystidae |
| 1-C05-C3_S35  | 0.30 | 12480  | 12374  | F | 5809 | Sarcocystidae |
| 1-C06-C4_S36  | 0.33 | 15641  | 15497  | F | 5809 | Sarcocystidae |
| 1-C07-C6_S37  | 0.74 | 31066  | 30854  | F | 5809 | Sarcocystidae |
| 1-C08-C7_S38  | 0.46 | 22605  | 22408  | F | 5809 | Sarcocystidae |
| 1-C09-C8_S39  | 0.56 | 26409  | 26272  | F | 5809 | Sarcocystidae |
| 1-C10-C9_S40  | 0.37 | 21653  | 21450  | F | 5809 | Sarcocystidae |
| 1-C11-D1_S41  | 0.41 | 22175  | 21975  | F | 5809 | Sarcocystidae |
| 1-C12-D2_S42  | 0.28 | 12319  | 12184  | F | 5809 | Sarcocystidae |
| 1-D01-D5_S43  | 0.41 | 23350  | 22776  | F | 5809 | Sarcocystidae |
| 1-D02-D6_S44  | 0.58 | 20728  | 20591  | F | 5809 | Sarcocystidae |
| 1-D03-D7_S45  | 0.92 | 44304  | 43577  | F | 5809 | Sarcocystidae |
| 1-D04-D8_S46  | 0.52 | 16964  | 16682  | F | 5809 | Sarcocystidae |
| 1-D05-D9_S47  | 0.78 | 32542  | 32200  | F | 5809 | Sarcocystidae |
| 1-D06-D10_S48 | 0.62 | 26940  | 26620  | F | 5809 | Sarcocystidae |
| 1-D07-D11_S49 | 0.63 | 50925  | 50494  | F | 5809 | Sarcocystidae |
| 1-D08-E1_S50  | 0.39 | 18830  | 18548  | F | 5809 | Sarcocystidae |
| 1-D09-E2_S51  | 0.59 | 27221  | 27049  | F | 5809 | Sarcocystidae |
| 1-D10-E5_S52  | 0.93 | 44089  | 43799  | F | 5809 | Sarcocystidae |
| 1-D11-E7_S53  | 0.39 | 10625  | 10542  | F | 5809 | Sarcocystidae |
| 1-D12-E8_S54  | 0.18 | 10475  | 10234  | F | 5809 | Sarcocystidae |
| 1-E01-E9_S55  | 0.39 | 14549  | 14347  | F | 5809 | Sarcocystidae |
| 1-E02-E10_S56 | 0.33 | 19594  | 19261  | F | 5809 | Sarcocystidae |
| 1-E03-F1_S57  | 0.39 | 22930  | 22626  | F | 5809 | Sarcocystidae |
| 1-E04-F2_S58  | 0.62 | 61776  | 61021  | F | 5809 | Sarcocystidae |
| 1-E05-F3_S59  | 0.51 | 14118  | 14024  | F | 5809 | Sarcocystidae |
| 1-E06-F4_S60  | 0.53 | 17959  | 17758  | F | 5809 | Sarcocystidae |
| 1-E07-F5_S61  | 0.51 | 2656   | 2598   | F | 5809 | Sarcocystidae |
| 1-E08-F6_S62  | 0.50 | 24870  | 24628  | F | 5809 | Sarcocystidae |
| 1-E09-F7_S63  | 0.39 | 17367  | 17204  | F | 5809 | Sarcocystidae |
| 1-E10-F8_S64  | 0.62 | 26464  | 26076  | F | 5809 | Sarcocystidae |
| 1-E11-F9_S65  | 0.62 | 28318  | 28008  | F | 5809 | Sarcocystidae |
| 1-E12-F10_S66 | 0.51 | 26303  | 25935  | F | 5809 | Sarcocystidae |

</p>
</details>


###### Sporisorium graminicola Read Counts

```
grep "Sporisorium graminicola" *.report|awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |sed 's/graminicola/_graminicola/g' |less
```
Not much information out there. The basidiomycete Sporisorium graminicola (formally Pseudozyma graminicola) strain CBS10092 was originally isolated from an herbaceous. Highly likely this is truly present in every sample, and not just a contaminant.
<details>
<p>

| Sample   | Proportion of reads |Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|------|-------|-------|---|--------|-------------------------|
| 1-A01-A1_S7   | 0.07 | 4115  | 4115  | S | 280036 | Sporisorium_graminicola |
| 1-A02-A2_S8   | 0.10 | 3989  | 3989  | S | 280036 | Sporisorium_graminicola |
| 1-A03-A3_S9   | 0.07 | 2520  | 2520  | S | 280036 | Sporisorium_graminicola |
| 1-A04-A4_S10  | 0.12 | 5448  | 5448  | S | 280036 | Sporisorium_graminicola |
| 1-A05-A5_S11  | 0.10 | 3646  | 3646  | S | 280036 | Sporisorium_graminicola |
| 1-A06-A6_S12  | 0.08 | 4068  | 4068  | S | 280036 | Sporisorium_graminicola |
| 1-A07-A7_S13  | 0.06 | 3673  | 3673  | S | 280036 | Sporisorium_graminicola |
| 1-A08-A8_S14  | 0.09 | 5498  | 5498  | S | 280036 | Sporisorium_graminicola |
| 1-A09-A9_S15  | 0.07 | 4477  | 4477  | S | 280036 | Sporisorium_graminicola |
| 1-A10-A10_S16 | 0.10 | 5339  | 5339  | S | 280036 | Sporisorium_graminicola |
| 1-A11-A11_S17 | 0.09 | 5221  | 5221  | S | 280036 | Sporisorium_graminicola |
| 1-A12-A12_S18 | 0.06 | 3268  | 3268  | S | 280036 | Sporisorium_graminicola |
| 1-B01-A13_S19 | 0.08 | 5192  | 5192  | S | 280036 | Sporisorium_graminicola |
| 1-B02-A14_S20 | 0.14 | 5641  | 5641  | S | 280036 | Sporisorium_graminicola |
| 1-B03-A16_S21 | 0.17 | 39231 | 39231 | S | 280036 | Sporisorium_graminicola |
| 1-B04-A17_S22 | 0.16 | 5861  | 5861  | S | 280036 | Sporisorium_graminicola |
| 1-B05-A18_S23 | 0.18 | 5215  | 5215  | S | 280036 | Sporisorium_graminicola |
| 1-B06-A19_S24 | 0.11 | 6380  | 6380  | S | 280036 | Sporisorium_graminicola |
| 1-B07-B2_S25  | 0.13 | 8138  | 8138  | S | 280036 | Sporisorium_graminicola |
| 1-B08-B3_S26  | 0.14 | 8120  | 8120  | S | 280036 | Sporisorium_graminicola |
| 1-B09-B4_S27  | 0.11 | 4399  | 4399  | S | 280036 | Sporisorium_graminicola |
| 1-B10-B5_S28  | 0.15 | 7830  | 7830  | S | 280036 | Sporisorium_graminicola |
| 1-B11-B6_S29  | 0.01 | 205   | 205   | S | 280036 | Sporisorium_graminicola |
| 1-B12-B7_S30  | 0.07 | 4916  | 4916  | S | 280036 | Sporisorium_graminicola |
| 1-C01-B8_S31  | 0.05 | 2804  | 2804  | S | 280036 | Sporisorium_graminicola |
| 1-C02-B10_S32 | 0.09 | 3631  | 3631  | S | 280036 | Sporisorium_graminicola |
| 1-C03-C1_S33  | 0.21 | 8692  | 8692  | S | 280036 | Sporisorium_graminicola |
| 1-C04-C2_S34  | 0.14 | 4069  | 4069  | S | 280036 | Sporisorium_graminicola |
| 1-C05-C3_S35  | 0.07 | 3085  | 3085  | S | 280036 | Sporisorium_graminicola |
| 1-C06-C4_S36  | 0.08 | 4059  | 4059  | S | 280036 | Sporisorium_graminicola |
| 1-C07-C6_S37  | 0.19 | 7823  | 7823  | S | 280036 | Sporisorium_graminicola |
| 1-C08-C7_S38  | 0.11 | 5518  | 5518  | S | 280036 | Sporisorium_graminicola |
| 1-C09-C8_S39  | 0.14 | 6761  | 6761  | S | 280036 | Sporisorium_graminicola |
| 1-C10-C9_S40  | 0.11 | 6303  | 6303  | S | 280036 | Sporisorium_graminicola |
| 1-C11-D1_S41  | 0.10 | 5694  | 5694  | S | 280036 | Sporisorium_graminicola |
| 1-C12-D2_S42  | 0.08 | 3389  | 3389  | S | 280036 | Sporisorium_graminicola |
| 1-D01-D5_S43  | 0.11 | 6024  | 6024  | S | 280036 | Sporisorium_graminicola |
| 1-D02-D6_S44  | 0.15 | 5211  | 5211  | S | 280036 | Sporisorium_graminicola |
| 1-D03-D7_S45  | 0.24 | 11701 | 11701 | S | 280036 | Sporisorium_graminicola |
| 1-D04-D8_S46  | 0.15 | 4773  | 4773  | S | 280036 | Sporisorium_graminicola |
| 1-D05-D9_S47  | 0.19 | 8035  | 8035  | S | 280036 | Sporisorium_graminicola |
| 1-D06-D10_S48 | 0.16 | 6930  | 6930  | S | 280036 | Sporisorium_graminicola |
| 1-D07-D11_S49 | 0.15 | 11848 | 11848 | S | 280036 | Sporisorium_graminicola |
| 1-D08-E1_S50  | 0.11 | 5080  | 5080  | S | 280036 | Sporisorium_graminicola |
| 1-D09-E2_S51  | 0.14 | 6651  | 6651  | S | 280036 | Sporisorium_graminicola |
| 1-D10-E5_S52  | 0.20 | 9598  | 9598  | S | 280036 | Sporisorium_graminicola |
| 1-D11-E7_S53  | 0.11 | 2892  | 2892  | S | 280036 | Sporisorium_graminicola |
| 1-D12-E8_S54  | 0.05 | 3106  | 3106  | S | 280036 | Sporisorium_graminicola |
| 1-E01-E9_S55  | 0.10 | 3857  | 3857  | S | 280036 | Sporisorium_graminicola |
| 1-E02-E10_S56 | 0.10 | 6118  | 6118  | S | 280036 | Sporisorium_graminicola |
| 1-E03-F1_S57  | 0.10 | 5909  | 5909  | S | 280036 | Sporisorium_graminicola |
| 1-E04-F2_S58  | 0.15 | 15490 | 15490 | S | 280036 | Sporisorium_graminicola |
| 1-E05-F3_S59  | 0.12 | 3201  | 3201  | S | 280036 | Sporisorium_graminicola |
| 1-E06-F4_S60  | 0.14 | 4770  | 4770  | S | 280036 | Sporisorium_graminicola |
| 1-E07-F5_S61  | 0.13 | 680   | 680   | S | 280036 | Sporisorium_graminicola |
| 1-E08-F6_S62  | 0.12 | 6235  | 6235  | S | 280036 | Sporisorium_graminicola |
| 1-E09-F7_S63  | 0.11 | 4799  | 4799  | S | 280036 | Sporisorium_graminicola |
| 1-E10-F8_S64  | 0.15 | 6417  | 6417  | S | 280036 | Sporisorium_graminicola |
| 1-E11-F9_S65  | 0.15 | 6793  | 6793  | S | 280036 | Sporisorium_graminicola |
| 1-E12-F10_S66 | 0.14 | 7316  | 7316  | S | 280036 | Sporisorium_graminicola |

</p>
</details>


###### leotiomyceta Read Counts

```
grep "Sporisorium graminicola" *.report|awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |sed 's/graminicola/_graminicola/g' |less
```
Half of these reads go to the species Zymoseptoria tritici IPO323
<details>
<p>

| Sample   | Proportion of reads |Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|------|-------|-------|---|--------|-------------------------|
1-A01-A1_S7     0.01    429     157     P3      716546  leotiomyceta
1-A02-A2_S8     0.01    307     172     P3      716546  leotiomyceta
1-A04-A4_S10    0.01    508     213     P3      716546  leotiomyceta
1-A05-A5_S11    0.01    222     92      P3      716546  leotiomyceta
1-A06-A6_S12    0.01    386     185     P3      716546  leotiomyceta
1-A07-A7_S13    0.01    424     192     P3      716546  leotiomyceta
1-A08-A8_S14    0.01    637     261     P3      716546  leotiomyceta
1-A09-A9_S15    0.01    327     167     P3      716546  leotiomyceta
1-A10-A10_S16   0.01    403     201     P3      716546  leotiomyceta
1-A11-A11_S17   0.01    370     176     P3      716546  leotiomyceta
1-A12-A12_S18   0.01    449     95      P3      716546  leotiomyceta
1-B01-A13_S19   0.01    582     313     P3      716546  leotiomyceta
1-B02-A14_S20   0.01    446     223     P3      716546  leotiomyceta
1-B03-A16_S21   0.01    2440    1031    P3      716546  leotiomyceta
1-B04-A17_S22   0.01    513     215     P3      716546  leotiomyceta
1-B05-A18_S23   0.01    320     153     P3      716546  leotiomyceta
1-B06-A19_S24   0.01    564     263     P3      716546  leotiomyceta
1-B07-B2_S25    0.01    522     216     P3      716546  leotiomyceta
1-B08-B3_S26    0.01    479     161     P3      716546  leotiomyceta
1-B09-B4_S27    0.01    284     120     P3      716546  leotiomyceta
1-B10-B5_S28    0.01    532     235     P3      716546  leotiomyceta
1-B11-B6_S29    0.02    504     9       P3      716546  leotiomyceta
1-B12-B7_S30    0.01    608     335     P3      716546  leotiomyceta
1-C01-B8_S31    0.01    293     172     P3      716546  leotiomyceta
1-C02-B10_S32   0.01    373     204     P3      716546  leotiomyceta
1-C03-C1_S33    0.01    629     260     P3      716546  leotiomyceta
1-C04-C2_S34    0.01    313     157     P3      716546  leotiomyceta
1-C05-C3_S35    0.01    210     114     P3      716546  leotiomyceta
1-C06-C4_S36    0.01    369     168     P3      716546  leotiomyceta
1-C07-C6_S37    0.01    601     307     P3      716546  leotiomyceta
1-C08-C7_S38    0.01    410     208     P3      716546  leotiomyceta
1-C09-C8_S39    0.01    375     211     P3      716546  leotiomyceta
1-C10-C9_S40    0.01    562     265     P3      716546  leotiomyceta
1-C11-D1_S41    0.01    326     173     P3      716546  leotiomyceta
1-C12-D2_S42    0.01    381     231     P3      716546  leotiomyceta
1-D01-D5_S43    0.01    570     173     P3      716546  leotiomyceta
1-D02-D6_S44    0.01    313     171     P3      716546  leotiomyceta
1-D03-D7_S45    0.02    1033    385     P3      716546  leotiomyceta
1-D04-D8_S46    0.02    696     257     P3      716546  leotiomyceta
1-D05-D9_S47    0.01    489     232     P3      716546  leotiomyceta
1-D06-D10_S48   0.01    631     282     P3      716546  leotiomyceta
1-D07-D11_S49   0.01    911     500     P3      716546  leotiomyceta
1-D08-E1_S50    0.02    730     249     P3      716546  leotiomyceta
1-D09-E2_S51    0.01    500     272     P3      716546  leotiomyceta
1-D10-E5_S52    0.01    416     52      P3      716546  leotiomyceta
1-D11-E7_S53    0.01    158     86      P3      716546  leotiomyceta
1-E01-E9_S55    0.01    262     116     P3      716546  leotiomyceta
1-E02-E10_S56   0.01    643     296     P3      716546  leotiomyceta
1-E03-F1_S57    0.01    468     237     P3      716546  leotiomyceta
1-E04-F2_S58    0.01    1137    559     P3      716546  leotiomyceta
1-E05-F3_S59    0.01    168     52      P3      716546  leotiomyceta
1-E06-F4_S60    0.01    477     231     P3      716546  leotiomyceta
1-E07-F5_S61    0.02    128     26      P3      716546  leotiomyceta
1-E08-F6_S62    0.01    449     144     P3      716546  leotiomyceta
1-E09-F7_S63    0.01    333     128     P3      716546  leotiomyceta
1-E10-F8_S64    0.01    444     156     P3      716546  leotiomyceta
1-E11-F9_S65    0.01    447     204     P3      716546  leotiomyceta
1-E12-F10_S66   0.01    561     161     P3      716546  leotiomyceta

</p>
</details>


###### Abyssogena phaseoliformis symbiont OG214

```
grep "Abyssogena" *.report|awk '$5=="S1"'|awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |sed 's/graminicola/_graminicola/g' |less
```
This is and endosymbiont of a vesicomyid bivalve very little known.

<details>
<p>

| Sample   | Proportion of reads |Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|------|-------|-------|---|--------|-------------------------|
| 1-A01-A1_S7   | 0.01 | 356  | 356  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A02-A2_S8   | 0.01 | 253  | 253  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A03-A3_S9   | 0.01 | 399  | 399  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A04-A4_S10  | 0.01 | 307  | 307  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A05-A5_S11  | 0.01 | 295  | 295  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A06-A6_S12  | 0.01 | 488  | 488  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A07-A7_S13  | 0.01 | 368  | 368  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A09-A9_S15  | 0.01 | 550  | 550  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A10-A10_S16 | 0.01 | 369  | 369  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-A11-A11_S17 | 0.01 | 308  | 308  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B01-A13_S19 | 0.01 | 401  | 401  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B02-A14_S20 | 0.01 | 336  | 336  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B03-A16_S21 | 0.01 | 2979 | 2979 | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B04-A17_S22 | 0.01 | 384  | 384  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B05-A18_S23 | 0.01 | 329  | 329  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B06-A19_S24 | 0.01 | 545  | 545  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B07-B2_S25  | 0.01 | 592  | 592  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B08-B3_S26  | 0.01 | 441  | 441  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B09-B4_S27  | 0.01 | 396  | 396  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-B10-B5_S28  | 0.01 | 351  | 351  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C03-C1_S33  | 0.01 | 337  | 337  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C04-C2_S34  | 0.01 | 272  | 272  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C05-C3_S35  | 0.01 | 261  | 261  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C06-C4_S36  | 0.01 | 279  | 279  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C07-C6_S37  | 0.01 | 429  | 429  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C08-C7_S38  | 0.01 | 263  | 263  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C09-C8_S39  | 0.01 | 262  | 262  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C10-C9_S40  | 0.01 | 420  | 420  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-C11-D1_S41  | 0.01 | 325  | 325  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D01-D5_S43  | 0.01 | 547  | 547  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D02-D6_S44  | 0.01 | 199  | 199  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D03-D7_S45  | 0.01 | 477  | 477  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D04-D8_S46  | 0.01 | 190  | 190  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D05-D9_S47  | 0.01 | 406  | 406  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D06-D10_S48 | 0.01 | 334  | 334  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D07-D11_S49 | 0.01 | 559  | 559  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D09-E2_S51  | 0.01 | 399  | 399  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D10-E5_S52  | 0.01 | 356  | 356  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D11-E7_S53  | 0.01 | 204  | 204  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-D12-E8_S54  | 0.01 | 329  | 329  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E01-E9_S55  | 0.01 | 245  | 245  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E04-F2_S58  | 0.01 | 576  | 576  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E05-F3_S59  | 0.01 | 175  | 175  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E06-F4_S60  | 0.01 | 181  | 181  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E08-F6_S62  | 0.01 | 572  | 572  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E09-F7_S63  | 0.01 | 255  | 255  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E11-F9_S65  | 0.01 | 456  | 456  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |
| 1-E12-F10_S66 | 0.01 | 449  | 449  | S1 | 1235283 | AbyssogenaphaseoliformissymbiontOG214 |

</p>
</details>


###### Dictyostelium discoideum AX4
```
 grep "Dictyostelium" *.report|awk '$5=="S1"'|awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |sed 's/graminicola/_graminicola/g' |less
```
The amoebozoa are a richly diverse group of organisms whose genomes remain largely unexplored. The soil-dwelling social amoeba Dictyostelium discoideum has been actively studied for the past 50 years and has contributed greatly to our understanding of cellular motility, signalling and interaction1. Dictyostelium amoebae inhabit forest soil and consume bacteria and yeast, which they track by chemotaxis. Starvation, however, prompts the solitary cells to aggregate and develop as a true multicellular organism, producing a fruiting body comprised of a cellular, cellulosic stalk supporting a bolus of spores. Thus, Dictyostelium has evolved mechanisms that direct the differentiation of a homogeneous population of cells into distinct cell types, regulate the proportions between tissues and orchestrate the construction of an effective structure for the dispersal of spores. Many of the genes necessary for these processes in Dictyostelium were also inherited by Metazoa and fashioned through evolution for use within many different modes of development.

The amoebozoa are also noteworthy as representing one of the earliest branches from the last common ancestor of all eukaryotes. Each of the surviving branches of the crown group of eukaryotes provides an example of the ways in which the ancestral genome has been sculpted and adapted by lineage-specific gene duplication, divergence and deletion. Comparison between representatives of these branches promises to shed light not only on the nature and content of the ancestral eukaryotic genome, but on the diversity of ways in which its components have been adapted to meet the needs of complex organisms.

[Reference](https://www.hgsc.bcm.edu/microbiome/dictyostelium-discoideum-ax4)

<details>
<p>

| Sample   | Proportion of reads |Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|------|-------|-------|---|--------|-------------------------|
| 1-A01-A1_S7   | 0.01 | 309  | 309  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-A04-A4_S10  | 0.01 | 387  | 387  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-A06-A6_S12  | 0.01 | 270  | 270  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-A08-A8_S14  | 0.01 | 401  | 401  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-A11-A11_S17 | 0.01 | 307  | 307  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-A12-A12_S18 | 0.01 | 424  | 424  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-B01-A13_S19 | 0.01 | 333  | 333  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-B02-A14_S20 | 0.01 | 251  | 251  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-B03-A16_S21 | 0.01 | 1356 | 1356 | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-B04-A17_S22 | 0.01 | 190  | 190  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-B10-B5_S28  | 0.01 | 276  | 276  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-B11-B6_S29  | 0.02 | 413  | 413  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-C02-B10_S32 | 0.01 | 273  | 273  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-C03-C1_S33  | 0.01 | 285  | 285  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-C07-C6_S37  | 0.01 | 220  | 220  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-C08-C7_S38  | 0.01 | 252  | 252  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-C10-C9_S40  | 0.01 | 367  | 367  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-C12-D2_S42  | 0.01 | 244  | 244  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D01-D5_S43  | 0.01 | 300  | 300  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D02-D6_S44  | 0.01 | 200  | 200  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D03-D7_S45  | 0.01 | 409  | 409  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D04-D8_S46  | 0.01 | 347  | 347  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D05-D9_S47  | 0.01 | 222  | 222  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D06-D10_S48 | 0.01 | 338  | 338  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D08-E1_S50  | 0.01 | 557  | 557  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D10-E5_S52  | 0.01 | 292  | 292  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-D11-E7_S53  | 0.01 | 144  | 144  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E02-E10_S56 | 0.01 | 566  | 566  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E04-F2_S58  | 0.01 | 692  | 692  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E06-F4_S60  | 0.01 | 292  | 292  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E07-F5_S61  | 0.01 | 72   | 72   | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E08-F6_S62  | 0.01 | 287  | 287  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E09-F7_S63  | 0.01 | 245  | 245  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E10-F8_S64  | 0.01 | 282  | 282  | S1 | 352472 | DictyosteliumdiscoideumAX4 |
| 1-E12-F10_S66 | 0.01 | 331  | 331  | S1 | 352472 | DictyosteliumdiscoideumAX4 |


</p>
</details>


###### Plasmodium Genus Content
```
grep "Plasmodium" *.report|awk '$5=="G"'|awk '$2>0' |sed 's/_/\t/2'  |sed 's/: /\t/1'  |cut -f 1,3-|sed 's/[ ]*//g' |sed 's/graminicola/_graminicola/g' |less
```

This guy is generally known for malaria, but kraken had a hard time pinpointing a species, half reads between Plasmodium vivax and most of the other half to Plasmodium relictum. Bee malaria?
<details>
<p>

| Sample   | Proportion of reads |Taxon group reads assigned  | taxon specific reads assigned | Taxon code | Taxonomy id | classification |
|---------------|------|-------|-------|---|--------|-------------------------|
| 1-A01-A1_S7   | 0.01 | 570  | 49  | G | 5820 | Plasmodium |
| 1-A02-A2_S8   | 0.01 | 308  | 28  | G | 5820 | Plasmodium |
| 1-A03-A3_S9   | 0.01 | 198  | 19  | G | 5820 | Plasmodium |
| 1-A04-A4_S10  | 0.01 | 574  | 46  | G | 5820 | Plasmodium |
| 1-A05-A5_S11  | 0.01 | 278  | 18  | G | 5820 | Plasmodium |
| 1-A06-A6_S12  | 0.01 | 415  | 36  | G | 5820 | Plasmodium |
| 1-A07-A7_S13  | 0.01 | 387  | 29  | G | 5820 | Plasmodium |
| 1-A08-A8_S14  | 0.01 | 630  | 65  | G | 5820 | Plasmodium |
| 1-A09-A9_S15  | 0.01 | 379  | 27  | G | 5820 | Plasmodium |
| 1-A10-A10_S16 | 0.01 | 459  | 59  | G | 5820 | Plasmodium |
| 1-A11-A11_S17 | 0.01 | 323  | 33  | G | 5820 | Plasmodium |
| 1-A12-A12_S18 | 0.01 | 480  | 62  | G | 5820 | Plasmodium |
| 1-B01-A13_S19 | 0.01 | 587  | 44  | G | 5820 | Plasmodium |
| 1-B02-A14_S20 | 0.01 | 413  | 42  | G | 5820 | Plasmodium |
| 1-B03-A16_S21 | 0.01 | 3409 | 389 | G | 5820 | Plasmodium |
| 1-B04-A17_S22 | 0.01 | 456  | 44  | G | 5820 | Plasmodium |
| 1-B05-A18_S23 | 0.01 | 306  | 29  | G | 5820 | Plasmodium |
| 1-B06-A19_S24 | 0.01 | 575  | 67  | G | 5820 | Plasmodium |
| 1-B07-B2_S25  | 0.01 | 548  | 63  | G | 5820 | Plasmodium |
| 1-B08-B3_S26  | 0.01 | 705  | 69  | G | 5820 | Plasmodium |
| 1-B09-B4_S27  | 0.01 | 293  | 35  | G | 5820 | Plasmodium |
| 1-B10-B5_S28  | 0.01 | 426  | 39  | G | 5820 | Plasmodium |
| 1-B11-B6_S29  | 0.01 | 301  | 109 | G | 5820 | Plasmodium |
| 1-B12-B7_S30  | 0.01 | 464  | 42  | G | 5820 | Plasmodium |
| 1-C02-B10_S32 | 0.01 | 333  | 28  | G | 5820 | Plasmodium |
| 1-C03-C1_S33  | 0.02 | 673  | 72  | G | 5820 | Plasmodium |
| 1-C04-C2_S34  | 0.01 | 334  | 31  | G | 5820 | Plasmodium |
| 1-C07-C6_S37  | 0.01 | 479  | 54  | G | 5820 | Plasmodium |
| 1-C08-C7_S38  | 0.01 | 270  | 26  | G | 5820 | Plasmodium |
| 1-C09-C8_S39  | 0.01 | 280  | 21  | G | 5820 | Plasmodium |
| 1-C10-C9_S40  | 0.01 | 488  | 61  | G | 5820 | Plasmodium |
| 1-C12-D2_S42  | 0.01 | 272  | 28  | G | 5820 | Plasmodium |
| 1-D01-D5_S43  | 0.01 | 837  | 61  | G | 5820 | Plasmodium |
| 1-D02-D6_S44  | 0.01 | 305  | 25  | G | 5820 | Plasmodium |
| 1-D03-D7_S45  | 0.02 | 1154 | 138 | G | 5820 | Plasmodium |
| 1-D04-D8_S46  | 0.02 | 685  | 112 | G | 5820 | Plasmodium |
| 1-D05-D9_S47  | 0.01 | 451  | 49  | G | 5820 | Plasmodium |
| 1-D06-D10_S48 | 0.01 | 601  | 64  | G | 5820 | Plasmodium |
| 1-D07-D11_S49 | 0.01 | 707  | 89  | G | 5820 | Plasmodium |
| 1-D08-E1_S50  | 0.01 | 605  | 78  | G | 5820 | Plasmodium |
| 1-D09-E2_S51  | 0.01 | 352  | 29  | G | 5820 | Plasmodium |
| 1-D10-E5_S52  | 0.01 | 488  | 81  | G | 5820 | Plasmodium |
| 1-D12-E8_S54  | 0.01 | 325  | 23  | G | 5820 | Plasmodium |
| 1-E01-E9_S55  | 0.01 | 344  | 40  | G | 5820 | Plasmodium |
| 1-E02-E10_S56 | 0.01 | 754  | 66  | G | 5820 | Plasmodium |
| 1-E03-F1_S57  | 0.01 | 386  | 27  | G | 5820 | Plasmodium |
| 1-E04-F2_S58  | 0.01 | 1191 | 103 | G | 5820 | Plasmodium |
| 1-E05-F3_S59  | 0.01 | 208  | 18  | G | 5820 | Plasmodium |
| 1-E06-F4_S60  | 0.01 | 437  | 40  | G | 5820 | Plasmodium |
| 1-E07-F5_S61  | 0.02 | 95   | 16  | G | 5820 | Plasmodium |
| 1-E08-F6_S62  | 0.01 | 425  | 45  | G | 5820 | Plasmodium |
| 1-E09-F7_S63  | 0.01 | 305  | 29  | G | 5820 | Plasmodium |
| 1-E10-F8_S64  | 0.01 | 402  | 37  | G | 5820 | Plasmodium |
| 1-E11-F9_S65  | 0.01 | 424  | 48  | G | 5820 | Plasmodium |
| 1-E12-F10_S66 | 0.01 | 687  | 72  | G | 5820 | Plasmodium |

</p>
</details>
