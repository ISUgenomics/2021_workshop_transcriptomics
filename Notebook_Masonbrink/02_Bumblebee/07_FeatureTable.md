# Create a table that has all data assessed

## create gene to transcript name mapping
```
less GCF_000188095.3_BIMP_2.2_genomic.gff |awk '$3=="mRNA"' |cut -f 9 |sed 's/;/\t/1'|sed 's/;/\t/1' |sed 's/ID=//g' |sed 's/Parent=//g' |cut -f 1,2 >GenetomRNA.tab

```
## Grab the functional annotations mapped to mRNA names
```
less GCF_000188095.3_BIMP_2.2_genomic.gff |awk '$3=="mRNA"' |sed 's/product=/\t/g' |sed 's/;/\t/1' |cut -f 1,4,5,9,11 |awk '{print $4"\t"$1":"$2"-"$3"\t"$5,$6,$7,
$14}' |sed 's/ID=//g' >Mrnafeatures.tab
```

### create vlookup for the deseq output files to create one conglomerate file
Note cannot paste from terminal, must be saved and opened in notepad or the dang thing wont paste in excel
```
x=1; while [ $x -le 90 ]; do echo "=VLOOKUP(\$b2,AllExposedvsAllControlGeneMult.csv\!\$a\$1:\$bq\$40934,$x,FALSE)";x=$(( $x +1));done |sed 's/\\//g' |tr "\n" "\t" >getter.tab
```

Merged all three deseq experiments into one table with the above information. 
