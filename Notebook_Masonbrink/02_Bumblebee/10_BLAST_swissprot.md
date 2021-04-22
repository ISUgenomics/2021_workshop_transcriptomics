# Need GO annotations, BLAST to swissprot -- hymenopterabase annotations


### Create SwissProt blast database
```
/work/gif/databases/Blast/SwissProt4-19-21/

#download and create blast database
wget https://ftp.ncbi.nlm.nih.gov/blast/db/v5/FASTA/swissprot.gz
gunzip swissprot.gz
makeblastdb -in swissprot -dbtype prot -out swissprot

```


### Run the blastp of our predicted peptides
```
/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/08_SwissProtAnnot

#Softlink annotation parts
for f in  ../07_InterproScan/*part*fa; do ln -s $f;done

#set up the blast
for f in *; do echo "ml blast-plus; blastp -db /work/gif/databases/Blast/SwissProt4-19-21/swissprot -query "$f" -outfmt 6 -num_threads 36 -out "${f%.*}"_blastp.out";done  >blast.sh
```


### Explore the blast and  create tab file
```
#How many genes were annotated with swissprot above an evalue of .001?
cat BombusImpatiens_Protein.part-*_blastp.out |awk '$11<.001' |awk '{print $1}' |sort|uniq|wc
  10687   10687  138931

cat BombusImpatiens_Protein.part-*_blastp.out |awk '$11<.001' |sort -u -k1,1 >One2OneBlastout.txt

How many of the annotations were unique? i.e. same annotation for different gene
cut -f 2  One2OneBlastout.txt |sort|uniq|wc
   8228    8228   74176


cut -f 2  One2OneBlastout.txt |sort|uniq >AnnotateMe.txt

#convert to swissprot ids

# also have all annotations for
uniprotOut.tab

#three annotations were not found by ID ??
less UniprotJustGOs |cut -f 1,7 |awk 'NR>1' >SPAnnotDB.list

get rid of those from my blast, as they mess up the alignment for my geneid "\t" annotations
 awk '{print $1 }' SPAnnotDB.list |cat - <(less One2OneBlastout.txt |awk '{print substr($2,1,length($2)-2)}') |sort |uniq -c |awk '$1==1 {print $2}' |grep -v -f - One2OneBlastout.txt  >FilteredOne2OneBlastout.txt

 less FilteredOne2OneBlastout.txt |awk '{print substr($2,1,length($2)-2)}' |while read line; do echo "awk '{if(\$1~\""$line"\"){print \$0} else {next;;}}' SPAnnotDB.list >>BlastAnnotTest"; done >grabAnnots.sh
 sh grabAnnots.sh

```


## Swissprot annotation for NCBI


###### Get the sequences
```
/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/09_SwissProtAnnotNCBI

 ln -s  ../01_Align/GCF_000188095.3_BIMP_2.2_genomic.gff
 ln -s ../01_Align/Bumblebee.fasta

 ml cufflinks; gffread -g Bumblebee.fasta GCF_000188095.3_BIMP_2.2_genomic.gff  -x BombusImpatiensNCBI_Transcript.fa  -y  BombusImpatiensNCBI_Protein.fa
```

### Run the blast to swissprot
```
/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/09_SwissProtAnnotNCBI

fasta-splitter.pl 8 BombusImpatiensNCBI_Protein.fa

#set up the blast
for f in *part*; do echo "ml blast-plus; blastp -db /work/gif/databases/Blast/SwissProt4-19-21/swissprot -query "$f" -outfmt 6 -num_threads 36 -out "${f%.*}"_blastp.out";done  >blast.sh
```


### Explore the NCBI annot blast and  create tab file
```
#How many genes were annotated with swissprot above an evalue of .001?
cat *part-*_blastp.out |awk '$11<.001' |awk '{print $1}' |sort|uniq|wc
  21737   21737  413003


cat *part-*_blastp.out  |awk '$11<.001' |sort -u -k1,1 >One2OneBlastout.txt

How many of the annotations were unique? i.e. same annotation for different gene
cut -f 2  One2OneBlastout.txt |sort|uniq|wc
   8540    8540   76980



cut -f 2  One2OneBlastout.txt |sort|uniq >AnnotateMe.txt

#convert to swissprot ids
NCBIAnnots_swissprot.tab


#three annotations were not found by ID ??
less NCBIAnnots_swissprot.tab |cut -f 1,7 |awk 'NR>1' >NCBIAnnotDB.list

#get rid of those from my blast, as they mess up the alignment for my geneid "\t" annotations
 awk '{print $1 }'  NCBIAnnotDB.list |cat - <(less One2OneBlastout.txt |awk '{print substr($2,1,length($2)-2)}') |sort |uniq -c |awk '$1==1 {print $2}' |grep -w -v -f - One2OneBlastout.txt  >FilteredOne2OneBlastout.txt
less BlastAnnotTest|awk '{print $1}' |grep -v -f -  FilteredOne2OneBlastout.txt > ReFilteredOne2OneBlastout.txt


 less FilteredOne2OneBlastout.txt |awk '{print substr($2,1,length($2)-2)}' |while read line; do echo "awk '{if(\$1==\""$line"\"){print \$0} else {next;;}}' NCBIAnnotDB.list >>BlastAnnotTest"; done >grabAnnots.sh
 sh grabAnnots.sh

 paste  <(cut -f 1 ReFilteredOne2OneBlastout.txt)  <(cut -f 2 BlastAnnotTest ) >mRNAAnnots.tab

#convert mrna names to gene names
awk '$3=="mRNA"' GCF_000188095.3_BIMP_2.2_genomic.gff |cut -f 9 |sed 's/;/\t/g' |cut -f 1,2 |sed 's/ID=//g' |sed 's/Parent=//g' >mRNA2Gene.list

less mRNAAnnots.tab |awk '{print $1}' |while read line; do grep -w "$line" mRNA2Gene.list ;done  |paste - <(cut -f 2 mRNAAnnots.tab )|sort -k2,2 -u |cut -f 2,3>GeneAnnots.tab
```
