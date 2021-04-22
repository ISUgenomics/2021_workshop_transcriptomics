# Need some ontology germs associated with genes to determine enrichment

### Split and run interproscan on Hymenoptera proteins for Bombus impatiens
```
#/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/07_InterproScan


#get sequences
ml cufflinks; gffread -g Bombus_impatiens_GCF_000188095.3_BIMP_2.2_genomic.fa  Bombus_impatiens_bimp_OGSv1.0_liftover_BIMP_2.2.gff3  -x BombusImpatiens_Transcript.fa  -y  BombusImpatiens_Protein.fa

#remove period from fasta
sed -i 's/\.//g' BombusImpatiens_Protein.fa

#split fastas
fasta-splitter.pl --n-parts 28 BombusImpatiens_Protein.fa

#create temp folders
for f in *part*; do mkdir $f.TEMP;done

#run interproscan
for f in *part*.fa; do echo "ml interproscan/5.38-76.0-py3-7ahddsa; interproscan.sh -cpu 18 -f TSV,GFF3 -goterms -i "$f" -pa -T "$f".TEMP"; done >IPRScan.sh
```
##### Investigate interpro results
```

#How many annotations were found total?
cat *fa.tsv  |awk '$3!="polypeptide"' |awk '$4!="Coils"' |awk '$4!="MobiDBLite"' |awk '$4!="CDD"' |wc
    498    7199   68441

#Number of genes that were annotated?    
[remkv6@nova052 07_InterproScan]$ cat *fa.tsv  |awk '$3!="polypeptide"' |awk '$4!="Coils"' |awk '$4!="MobiDBLite"' |awk '$4!="CDD" {print $1}' |sort|uniq|wc
    226     226    2938

#How many GO term annotations
cat *fa.tsv  |awk '$3!="polypeptide"' |awk '$4!="Coils"' |awk '$4!="MobiDBLite"' |awk '$4!="CDD" ' |grep -c "GO"
147

#How many genes own these 147 GO terms?
cat *fa.tsv  |awk '$3!="polypeptide"' |awk '$4!="Coils"' |awk '$4!="MobiDBLite"' |awk '$4!="CDD" ' |grep  "GO" |awk '{print $1}' |sort|uniq|wc
     61      61     793
```


Perhaps I should try to see how the NCBI annotations fare with Interproscan.
