# BioLlama
bioinformatics tools

## Installation
```bash
git clone <biollama_url>
cd biollama
python -m venv venv
source venv/bin/activate
python setup.py develop
python biollama/test/test_annotation.py run
```

## Annotate a dataframe with Ensembl
```python
import pandas as pd
from biollama.core.annotation import LlamaEnsembl, UCSCapi
llama = LlamaEnsembl()  # hg19 by default
#llama = LlamaEnsembl(genome='hg38')

# get all genes at a given position
llama.db.genes_at_locus('chr15', 32990907, 33011217)
llama.get_genes('chr15', 32990907, 33011217)

# annotate dataframe
df = pd.read_csv('table.txt', sep='\t')
ndf = llama.annotate_dataframe(df, chrom_col='CHROM', start_col='START', end_col='END')
mdf = pd.concat([df, ndf], axis=1)

# get rsid locations from dbSNP
df = llama.annotate_variants(["rs11111", "rs222222"])

# annotate with UCSC
llama = UCSCapi() # hg19 by default
# llama = UCSCapi(genome='hg38')
xdf = llama.annotate_dataframe(df)
```
