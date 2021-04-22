# Review: Gene Ontology

Review gene ontology from first principles

Download ontology file

* http://geneontology.org/docs/download-ontology/

```
wget http://purl.obolibrary.org/obo/go.obo

less go.obo
```

Each term looks similar to the following:

```
[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process
def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
synonym: "mitochondrial inheritance" EXACT []
is_a: GO:0048308 ! organelle inheritance
is_a: GO:0048311 ! mitochondrion distribution
```

which can be represented by a tree, or edgelist

```
from_node   to_node
GO:0048308  GO:0000001
GO:0048311  GO:0000001
```

and nodelist

```
id  name    namespace
GO:0000001  mitochondrion inheritance   biological_process
GO:0048308  organelle inheritance   biological_process
GO:0048311  mitochondrion distribution  biological_process
```

Generating the entire edge and node list from each Gene Ontology term gives us the <b>Gene Ontology Topology</b> graph. 

Topology-based Gene Ontology Enrichment is usually mapping Up/Down DEGs (differentially expressed genes) onto the Gene Ontology Topology and identifying regions of the network that are significantly have more Up or Down votes. Some normalization is necessary since certain terms are over-represented in all genes.
