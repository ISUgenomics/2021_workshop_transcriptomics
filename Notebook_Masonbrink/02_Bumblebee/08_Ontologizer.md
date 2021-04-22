#  See what kind of gene functions are up and downregulated upon heavy metal exposure




### Get GO terms

NCBI annotations only have the nucleotide entry for a definition.  Need to convert
```
/work/gif/TranscriptomicsWorkshop/remkv6/02_BumbleBee/06_Ontologizer


wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz

less GCF_000188095.3_BIMP_2.2_genomic.gff | awk '$3=="mRNA"' |sed 's/transcript_id=/\t/g' |cut -f 10 |while read line; do grep "$line" bombusOnly.accession >>Convert1.tab ;done &

```


### Ontologizer setup
```
wget http://ontologizer.de/cmdline/Ontologizer.jar
wget http://purl.obolibrary.org/obo/go.obo

ln -s ../09_SwissProtAnnotNCBI/GeneAnnots.tab
sed 's/; /,/g' GeneAnnots.tab >ReformatGeneAnnots.tab
awk '$2!=""' ReformatGeneAnnots.ids >NoEmptyGeneAnnots.ids
#make sure you change the file extension above or you will see no end to errors. Note there can be no terms in your association file to genes without a go term.

less ../03_Deseq_gene/AllExposedvsAllControlGene.csv |sed 's/,/\t/g' |awk '$7<.050000000000000000000000001' |cut -f 1 >GenesUniqSignificant.list

#genes upregulated in exposed bees/downregulated in ctrl bees
less ../03_Deseq_gene/AllExposedvsAllControlGene.csv |sed 's/,/\t/g' |awk '$7<.050000000000000000000000001 && $3>0' |cut -f 1 >UpGenesUniqSignificant.list

#genes downregulated in exposed bees/upregulated in ctrl bees
less ../03_Deseq_gene/AllExposedvsAllControlGene.csv |sed 's/,/\t/g' |awk '$7<.050000000000000000000000001 && $3<0' |cut -f 1 >DownGenesUniqSignificant.list
wc UpGenesUniqSignificant.list
 349  349 6282 UpGenesUniqSignificant.list
[remkv6@nova014 06_Ontologizer]$ wc DownGenesUniqSignificant.list
 214  214 3852 DownGenesUniqSignificant.list



#gene uniq expressed only
less ../03_Deseq_gene/AllExposedvsAllControlGene.csv |sed 's/,/\t/g' |awk '$2>10' |cut -f 1 |awk 'NR>1' >GenesExpressedPopulation.list

#All genes population
less ../04_Deseq_Mult/AllExposedvsAllControlGeneMult.csv |sed 's/,/\t/g' |cut -f 1 |awk 'NR>1' >AllGenesPopulation.list

```




### Regular significant genes to unique reads
##### compare all up and downregulated significant genes against only genes that are expressed with 10 or more normalized reads
```
java -jar Ontologizer.jar -d -g go.obo -a NoEmptyGeneAnnots.ids -p GenesExpressedPopulation.list -s GenesUniqSignificant.list -m Benjamini-Hochberg -d 0.05 -r 1000
#No significant terms, plenty that were not trivial
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0045595      7742    10      563     4       90      9       3       false   0.007737456091803563    0.8142987400530527      1.748054486444731E-13   "regulation of cell differentiation"
GO:0051253      7742    14      563     4       58      5       4       false   0.010049069032735159    0.8142987400530527      9.859073675355085E-14   "negative regulation of RNA metabolic process"
GO:0000041      7742    2       563     2       14      2       1       false   0.010989010989010973    0.8142987400530527      0.010989010989010973    "transition metal ion transport"
GO:0010628      7742    11      563     4       67      7       3       false   0.011366989443621262    0.8142987400530527      7.781717560880857E-13   "positive regulation of gene expression"
GO:0034708      7742    2       563     2       121     15      2       false   0.014462809917355188    0.8142987400530527      1.3774104683195873E-4   "methyltransferase complex"
GO:0008283      7742    16      563     5       120     13      1       false   0.015000006364228026    0.8142987400530527      3.2212283363904E-20     "cell population proliferation"
GO:0045934      7742    15      563     4       70      6       4       false   0.016757644005038548    0.8142987400530527      1.3860384767736189E-15  "negative regulation of nucleobase-containing compound metabolic process"
GO:0035097      7742    2       563     2       49      7       3       false   0.017857142857143102    0.8142987400530527      8.503401360544278E-4    "histone methyltransferase complex"
GO:0050661      7742    2       563     2       34      6       1       false   0.026737967914438263    0.8142987400530527      0.0017825311942958834   "NADP binding"
GO:0071559      7742    3       563     2       15      2       2       false   0.0285714285714285      0.8142987400530527      0.002197802197802196    "response to transforming growth factor beta"

```
##### Compare all up and downregulated significant genes against all genes, regardless of expression
```
 java -jar Ontologizer.jar -d -g go.obo -a NoEmptyGeneAnnots.ids -p AllGenesPopulation.list -s GenesUniqSignificant.list  -m Benjamini-Hochberg -d 0.05 -r 1000
#none significant
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0000041      12925   2       563     2       20      2       1       false   0.005263157894736846    0.7524882811049166      0.005263157894736846    "transition metal ion transport"
GO:0032270      12925   19      563     5       101     8       4       false   0.00556216590543168     0.7524882811049166      6.134740354731204E-21   "positive regulation of cellular protein metabolic process"
GO:0051247      12925   20      563     5       99      8       4       false   0.007895465035864614    0.7524882811049166      2.332161908415525E-21   "positive regulation of protein metabolic process"
GO:0034708      12925   2       563     2       174     16      2       false   0.007972892166633555    0.7524882811049166      6.64407680552764E-5     "methyltransferase complex"
GO:0008283      12925   19      563     5       180     14      1       false   0.008344646215190678    0.7524882811049166      4.5962676768007216E-26  "cell population proliferation"
GO:0045595      12925   13      563     4       127     10      3       false   0.0100785260088023      0.7524882811049166      5.2574134715946324E-18  "regulation of cell differentiation"
GO:0035097      12925   2       563     2       63      7       3       false   0.010752688172043159    0.7524882811049166      5.120327700972899E-4    "histone methyltransferase complex"
GO:0010628      12925   12      563     4       94      9       3       false   0.014079455927745814    0.7524882811049166      2.0940533695286887E-15  "positive regulation of gene expression"
GO:0050661      12925   2       563     2       46      6       1       false   0.014492753623188337    0.7524882811049166      9.66183574879222E-4     "NADP binding"
GO:1902494      12925   34      563     5       97      6       1       false   0.019100951019507907    0.7524882811049166      6.0849798289196194E-27  "catalytic complex"
```

##### Genes that are upregulated in gene uniq comparisons -- compare only to expressed genes -- up in bees exposed/down in ctrl
```
java -jar Ontologizer.jar -d -g go.obo -a NoEmptyGeneAnnots.ids -p GenesExpressedPopulation.list -s UpGenesUniqSignificant.list -m Benjamini-Hochberg -d 0.05 -r 1000
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0008283      12925   19      349     5       180     10      1       false   0.001391160536632669    0.113   4.5962676768007216E-26  "cell population proliferation"
GO:0034708      12925   2       349     2       174     11      2       false   0.0036542422430406603   0.261   6.64407680552764E-5     "methyltransferase complex"
GO:0045595      12925   13      349     4       127     8       3       false   0.003798893924469157    0.266   5.2574134715946324E-18  "regulation of cell differentiation"
GO:0035097      12925   2       349     2       63      5       3       false   0.00512032770097285     0.334   5.120327700972899E-4    "histone methyltransferase complex"
GO:0009890      12925   21      349     4       106     5       3       false   0.005220736398607474    0.338   1.255613223786561E-22   "negative regulation of biosynthetic process"
GO:0000041      12925   2       349     2       20      2       1       false   0.005263157894736846    0.34    0.005263157894736846    "transition metal ion transport"
GO:0045934      12925   21      349     4       102     5       4       false   0.006064640077452856    0.378   3.0805726315004322E-22  "negative regulation of nucleobase-containing compound metabolic process"
GO:0031327      12925   21      349     4       102     5       4       false   0.006064640077452856    0.378   3.0805726315004322E-22  "negative regulation of cellular biosynthetic process"
GO:0008134      12925   11      349     3       84      4       1       false   0.00641357532336086     0.388   5.383819446067372E-14   "transcription factor binding"
GO:0005730      12925   18      349     4       103     6       2       false   0.008162480603122229    0.458   1.8213000235557105E-20  "nucleolus"
GO:0051253      12925   19      349     4       83      5       4       false   0.008944288009297394    0.485   3.912071938929928E-19   "negative regulation of RNA metabolic process"
GO:0032270      12925   19      349     4       101     6       4       false   0.010930630197460552    0.571   6.134740354731204E-21   "positive regulation of cellular protein metabolic process"
GO:0004197      12925   3       349     2       22      2       2       false   0.01298701298701296     0.662   6.493506493506473E-4    "cysteine-type endopeptidase activity"
GO:0042127      12925   17      349     4       122     8       2       false   0.013092452806845185    0.664   3.8950633928039094E-21  "regulation of cell population proliferation"
GO:1902494      12925   34      349     4       97      4       1       false   0.013384745038731834    0.686   6.0849798289196194E-27  "catalytic complex"
GO:0010558      12925   21      349     4       83      5       4       false   0.013481217243161807    0.692   4.075074936385274E-20   "negative regulation of macromolecule biosynthetic process"
GO:0051247      12925   20      349     4       99      6       4       false   0.014449440666811954    0.721   2.332161908415525E-21   "positive regulation of protein metabolic process"
GO:0051301      12925   11      349     3       180     10      1       false   0.01624645158850002     0.771   8.482373013288667E-18   "cell division"
GO:1901681      12925   4       349     2       156     9       1       false   0.016799140138795816    0.785   4.212527931692659E-8    "sulfur compound binding"
GO:0051287      12925   3       349     2       46      4       1       false   0.016864295125164665    0.785   6.587615283267417E-5    "NAD binding"
GO:0045596      12925   3       349     2       84      7       4       false   0.01733764325595097     0.798   1.049494143822672E-5    "negative regulation of cell differentiation"
GO:0051093      12925   9       349     3       107     8       3       false   0.018957610819996404    0.838   2.78905312540007E-13    "negative regulation of developmental process"
GO:0032774      12925   49      349     5       104     5       3       false   0.02073544744097849     0.889   7.498766231381724E-31   "RNA biosynthetic process"

```
##### Genes that are upregulated in gene uniq comparisons -- compare to all genes  -- up in bees exposed/down in ctrl
```
java -jar Ontologizer.jar -d -g go.obo -a NoEmptyGeneAnnots.ids -p AllGenesPopulation.list -s UpGenesUniqSignificant.list -m Benjamini-Hochberg -d 0.05 -r 1000
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0008283      12925   19      349     5       180     10      1       false   0.001391160536632669    0.117   4.5962676768007216E-26  "cell population proliferation"
GO:0034708      12925   2       349     2       174     11      2       false   0.0036542422430406603   0.286   6.64407680552764E-5     "methyltransferase complex"
GO:0045595      12925   13      349     4       127     8       3       false   0.003798893924469157    0.293   5.2574134715946324E-18  "regulation of cell differentiation"
GO:0035097      12925   2       349     2       63      5       3       false   0.00512032770097285     0.357   5.120327700972899E-4    "histone methyltransferase complex"
GO:0009890      12925   21      349     4       106     5       3       false   0.005220736398607474    0.359   1.255613223786561E-22   "negative regulation of biosynthetic process"
GO:0000041      12925   2       349     2       20      2       1       false   0.005263157894736846    0.365   0.005263157894736846    "transition metal ion transport"
GO:0045934      12925   21      349     4       102     5       4       false   0.006064640077452856    0.403   3.0805726315004322E-22  "negative regulation of nucleobase-containing compound metabolic process"
GO:0031327      12925   21      349     4       102     5       4       false   0.006064640077452856    0.403   3.0805726315004322E-22  "negative regulation of cellular biosynthetic process"
GO:0008134      12925   11      349     3       84      4       1       false   0.00641357532336086     0.423   5.383819446067372E-14   "transcription factor binding"
GO:0005730      12925   18      349     4       103     6       2       false   0.008162480603122229    0.494   1.8213000235557105E-20  "nucleolus"
GO:0051253      12925   19      349     4       83      5       4       false   0.008944288009297394    0.523   3.912071938929928E-19   "negative regulation of RNA metabolic process"
GO:0032270      12925   19      349     4       101     6       4       false   0.010930630197460552    0.593   6.134740354731204E-21   "positive regulation of cellular protein metabolic process"
GO:0004197      12925   3       349     2       22      2       2       false   0.01298701298701296     0.666   6.493506493506473E-4    "cysteine-type endopeptidase activity"
GO:0042127      12925   17      349     4       122     8       2       false   0.013092452806845185    0.668   3.8950633928039094E-21  "regulation of cell population proliferation"
GO:1902494      12925   34      349     4       97      4       1       false   0.013384745038731834    0.686   6.0849798289196194E-27  "catalytic complex"
GO:0010558      12925   21      349     4       83      5       4       false   0.013481217243161807    0.688   4.075074936385274E-20   "negative regulation of macromolecule biosynthetic process"
GO:0051247      12925   20      349     4       99      6       4       false   0.014449440666811954    0.713   2.332161908415525E-21   "positive regulation of protein metabolic process"
GO:0051301      12925   11      349     3       180     10      1       false   0.01624645158850002     0.769   8.482373013288667E-18   "cell division"
GO:1901681      12925   4       349     2       156     9       1       false   0.016799140138795816    0.78    4.212527931692659E-8    "sulfur compound binding"
GO:0051287      12925   3       349     2       46      4       1       false   0.016864295125164665    0.78    6.587615283267417E-5    "NAD binding"
GO:0045596      12925   3       349     2       84      7       4       false   0.01733764325595097     0.795   1.049494143822672E-5    "negative regulation of cell differentiation"
GO:0051093      12925   9       349     3       107     8       3       false   0.018957610819996404    0.855   2.78905312540007E-13    "negative regulation of developmental process"
GO:0032774      12925   49      349     5       104     5       3       false   0.02073544744097849     0.887   7.498766231381724E-31   "RNA biosynthetic process"

```
##### Genes that are downregulated in gene uniq comparisons -- compare only to expressed genes  -- up in ctrl bees/down in exposed bees
```
java -jar Ontologizer.jar -d -g go.obo -a NoEmptyGeneAnnots.ids -p GenesExpressedPopulation.list -s DownGenesUniqSignificant.list -m Benjamini-Hochberg -d 0.05 -r 1000
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0045934      7730    15      349     4       70      4       4       false   0.0014887200824521987   0.742496998799522       1.3860384767736189E-15  "negative regulation of nucleobase-containing compound metabolic process"
GO:0008283      7730    16      349     5       120     9       1       false   0.0020662109007705716   0.742496998799522       3.2212283363904E-20     "cell population proliferation"
GO:0051253      7730    14      349     4       58      4       4       false   0.0023593466424682644   0.742496998799522       9.859073675355085E-14   "negative regulation of RNA metabolic process"
GO:0045595      7730    10      349     4       90      7       3       false   0.0024181410578162734   0.742496998799522       1.748054486444731E-13   "regulation of cell differentiation"
GO:0034708      7730    2       349     2       121     10      2       false   0.006198347107438396    0.742496998799522       1.3774104683195873E-4   "methyltransferase complex"
GO:0035097      7730    2       349     2       49      5       3       false   0.008503401360544198    0.742496998799522       8.503401360544278E-4    "histone methyltransferase complex"
GO:0000041      7730    2       349     2       14      2       1       false   0.010989010989010973    0.742496998799522       0.010989010989010973    "transition metal ion transport"
GO:0042127      7730    14      349     4       87      7       2       false   0.011593013244187763    0.742496998799522       1.8488704923520847E-16  "regulation of cell population proliferation"
GO:0008134      7730    9       349     3       53      4       1       false   0.013052164261931077    0.742496998799522       2.2565144472039573E-10  "transcription factor binding"
GO:0048523      7730    46      349     7       120     9       3       false   0.01577769265038141     0.742496998799522       2.7209859937467587E-34  "negative regulation of cellular process"


```
##### Genes that are downregulated in gene uniq comparisons -- compare to all genes -- up in ctrl bees/down in exposed bees
```
java -jar Ontologizer.jar -d -g go.obo -a NoEmptyGeneAnnots.ids -p AllGenesPopulation.list -s DownGenesUniqSignificant.list -m Benjamini-Hochberg -d 0.05 -r 1000
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:1903317      12925   1       214     1       78      1       3       false   0.01282051282051256     0.7950913242009178      0.01282051282051256     "regulation of protein maturation"
GO:0010181      12925   1       214     1       49      1       2       false   0.020408163265306332    0.7950913242009178      0.020408163265306332    "FMN binding"
GO:1903513      12925   1       214     1       44      1       1       false   0.02272727272727276     0.7950913242009178      0.02272727272727276     "endoplasmic reticulum to cytosol transport"
GO:0044322      12925   1       214     1       191     5       2       false   0.026178010471206725    0.7950913242009178      0.005235602094241107    "endoplasmic reticulum quality control compartment"
GO:0006809      12925   1       214     1       76      2       3       false   0.02631578947368393     0.7950913242009178      0.013157894736841941    "nitric oxide biosynthetic process"
GO:0046167      12925   1       214     1       36      1       3       false   0.027777777777777804    0.7950913242009178      0.027777777777777804    "glycerol-3-phosphate biosynthetic process"
GO:2001057      12925   1       214     1       134     4       1       false   0.02985074626865628     0.7950913242009178      0.007462686567164042    "reactive nitrogen species metabolic process"
GO:1904813      12925   2       214     1       67      1       2       false   0.029850746268657073    0.7950913242009178      4.522840343735954E-4    "ficolin-1-rich granule lumen"
GO:0052646      12925   2       214     1       65      1       3       false   0.03076923076923048     0.7950913242009178      4.8076923076922174E-4   "alditol phosphate metabolic process"
GO:1903319      12925   1       214     1       29      1       4       false   0.034482758620689634    0.7950913242009178      0.034482758620689634    "positive regulation of protein maturation"
GO:0046903      12925   15      214     2       76      2       1       false   0.036842105263157475    0.7950913242009178      3.520294354604658E-16   "secretion"

```



### Gene Mult comparisons
