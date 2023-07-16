# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_BA2'], na_filter=None)
variants = variants.append(pd.read_csv(config['codon_variant_table_file_BQ11'], na_filter=None))
variants = variants.append(pd.read_csv(config['codon_variant_table_file_XBB15'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>BA2</td>
      <td>pool1</td>
      <td>AAAAAAAAAAACGCGA</td>
      <td>3</td>
      <td>ATT88GTT</td>
      <td>I88V</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA2</td>
      <td>pool1</td>
      <td>AAAAAAAAACAGCGAG</td>
      <td>7</td>
      <td>AGA163ATT</td>
      <td>R163I</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA2</td>
      <td>pool1</td>
      <td>AAAAAAAAACCACGAA</td>
      <td>5</td>
      <td>CAA79ATG</td>
      <td>Q79M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA2</td>
      <td>pool1</td>
      <td>AAAAAAAAACCGAACT</td>
      <td>3</td>
      <td>GGT83AAA</td>
      <td>G83K</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA2</td>
      <td>pool1</td>
      <td>AAAAAAAAACGCTATG</td>
      <td>6</td>
      <td>GAA135GAT</td>
      <td>E135D</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAGTTATTTGG</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAATAAAATTGTTTT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACAACATATTATAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAGTTGGTAAAAAAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AATTTAACTATAAAAT</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 43 duplicated barcodes.Started with 649567 barcodes:
    After removing duplicates, there are 649481 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_BA2'],
                                     feature_parse_specs=config['feature_parse_specs_BA2'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>230517</td>
      <td>690999</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s1_b1_S1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>230517</td>
      <td>1253223</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s1_b2_S2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>230517</td>
      <td>1865810</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s1_b3_S3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>230517</td>
      <td>5791174</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s1_b4_S4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>230517</td>
      <td>1043072</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s2_b1_S5_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>230517</td>
      <td>1668641</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s2_b2_S6_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>230517</td>
      <td>2079964</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s2_b3_S7_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>230517</td>
      <td>5453883</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s2_b4_S8_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>230517</td>
      <td>2104022</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s3_b1_S9_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>230517</td>
      <td>1203237</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s3_b2_S10_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>230517</td>
      <td>2197256</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s3_b3_S11_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>230517</td>
      <td>4758681</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s3_b4_S12_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>230517</td>
      <td>2343998</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s4_b1_S13_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>230517</td>
      <td>2212026</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s4_b2_S14_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>230517</td>
      <td>4129711</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s4_b3_S15_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>230517</td>
      <td>1640237</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s4_b4_S16_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>230517</td>
      <td>4370574</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s5_b1_S17_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>230517</td>
      <td>4056727</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s5_b2_S18_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>230517</td>
      <td>1824750</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s5_b3_S19_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>230517</td>
      <td>7668</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s5_b4_S20_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>230517</td>
      <td>7333866</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s6_b1_S21_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>230517</td>
      <td>3057202</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s6_b2_S22_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>230517</td>
      <td>16789</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s6_b3_S23_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>230517</td>
      <td>195</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s6_b4_S24_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>230517</td>
      <td>9725601</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s7_b1_S25_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>230517</td>
      <td>312586</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s7_b2_S26_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>230517</td>
      <td>541</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s7_b3_S27_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>230517</td>
      <td>179</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s7_b4_S28_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>230517</td>
      <td>9811543</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s8_b1_S29_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>230517</td>
      <td>217714</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s8_b2_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>230517</td>
      <td>488</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s8_b3_S31_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>230517</td>
      <td>189</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s8_b4_S32_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>230517</td>
      <td>9824994</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>230517</td>
      <td>205386</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>230517</td>
      <td>554</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>230517</td>
      <td>190</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>230518</td>
      <td>523376</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s1_b1_S49_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>230518</td>
      <td>575326</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s1_b2_S50_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>230518</td>
      <td>1424453</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s1_b3_S51_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>230518</td>
      <td>6801514</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s1_b4_S52_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>230518</td>
      <td>327274</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s2_b1_S53_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>230518</td>
      <td>1258791</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s2_b2_S54_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>230518</td>
      <td>1959051</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s2_b3_S55_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>230518</td>
      <td>6568947</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s2_b4_S56_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>230518</td>
      <td>1316469</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s3_b1_S57_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>230518</td>
      <td>879345</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s3_b2_S58_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>230518</td>
      <td>2317196</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s3_b3_S59_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>230518</td>
      <td>5740592</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s3_b4_S60_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>230518</td>
      <td>1938463</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s4_b1_S61_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X87_230712_A00421_0570_BHFTYKDRX3_S87_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X87_230712_A00421_0570_BHFTYKDRX3_S87_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>230518</td>
      <td>1893991</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s4_b2_S62_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X88_230712_A00421_0570_BHFTYKDRX3_S88_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X88_230712_A00421_0570_BHFTYKDRX3_S88_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>230518</td>
      <td>4232319</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s4_b3_S63_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X89_230712_A00421_0570_BHFTYKDRX3_S89_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X89_230712_A00421_0570_BHFTYKDRX3_S89_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>230518</td>
      <td>1546552</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s4_b4_S64_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>230518</td>
      <td>4075797</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s5_b1_S65_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>230518</td>
      <td>4405145</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s5_b2_S66_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X90_230712_A00421_0570_BHFTYKDRX3_S90_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X90_230712_A00421_0570_BHFTYKDRX3_S90_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>230518</td>
      <td>1691129</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s5_b3_S67_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X91_230712_A00421_0570_BHFTYKDRX3_S91_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X91_230712_A00421_0570_BHFTYKDRX3_S91_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>230518</td>
      <td>10122</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s5_b4_S68_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>230518</td>
      <td>6859685</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X92_230712_A00421_0570_BHFTYKDRX3_S92_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X92_230712_A00421_0570_BHFTYKDRX3_S92_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>230518</td>
      <td>2433550</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s6_b2_S70_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X93_230712_A00421_0570_BHFTYKDRX3_S93_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X93_230712_A00421_0570_BHFTYKDRX3_S93_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>230518</td>
      <td>15917</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s6_b3_S71_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>230518</td>
      <td>198</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s6_b4_S72_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>230518</td>
      <td>9527244</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X94_230712_A00421_0570_BHFTYKDRX3_S94_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X94_230712_A00421_0570_BHFTYKDRX3_S94_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>230518</td>
      <td>353513</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X95_230712_A00421_0570_BHFTYKDRX3_S95_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X95_230712_A00421_0570_BHFTYKDRX3_S95_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>230518</td>
      <td>672</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s7_b3_S39_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>230518</td>
      <td>226</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s7_b4_S40_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>230518</td>
      <td>10027026</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s8_b1_S41_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>230518</td>
      <td>287556</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s8_b2_S42_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>230518</td>
      <td>552</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s8_b3_S43_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>230518</td>
      <td>176</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s8_b4_S44_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>230518</td>
      <td>9595454</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s9_b1_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>230518</td>
      <td>258786</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s9_b2_S46_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>230518</td>
      <td>611</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s9_b3_S47_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>230518</td>
      <td>199</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_s9_b4_S48_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>230517</td>
      <td>1639896</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin1_1_S37_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin1_2_S38_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin1_3_S39_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X85_230712_A00421_0570_BHFTYKDRX3_S85_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X85_230712_A00421_0570_BHFTYKDRX3_S85_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>230517</td>
      <td>3200001</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X86_230712_A00421_0570_BHFTYKDRX3_S86_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X86_230712_A00421_0570_BHFTYKDRX3_S86_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>230517</td>
      <td>3668604</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin3_1_S43_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin3_2_S44_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin3_3_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>230517</td>
      <td>6032001</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin4_1_S46_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin4_2_S47_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230517_exp1_RBD_bin4_3_S48_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>230518</td>
      <td>1344078</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin1_1_S49_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin1_2_S50_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin1_3_S51_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X96_230712_A00421_0570_BHFTYKDRX3_S96_L001_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2023/230628_Aerium-Vir-mAbs_XBB-polish/21073R/Fastq/21073X96_230712_A00421_0570_BHFTYKDRX3_S96_L002_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>230518</td>
      <td>3350001</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin2_1_S52_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin2_2_S53_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin2_3_S54_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>230518</td>
      <td>3606351</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin3_1_S55_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin3_2_S56_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin3_3_S57_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>230518</td>
      <td>4804152</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin4_1_S58_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin4_2_S59_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230531_Omicron-xbb-bq_bc/230518_exp2_RBD_bin4_3_S60_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool2</td>
      <td>334677</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>314804</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>AAGGAAAATATTAACT</td>
      <td>2041</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>AATGATCTCTTATCAT</td>
      <td>1820</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>TTTTATCGGAAACGAC</td>
      <td>1426</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>ACGAACGTTTCGCTCA</td>
      <td>1337</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>GCAATTCAACGGCAAC</td>
      <td>1088</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>1222961</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>263712</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>177811</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>18933</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">pool1</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>275542</td>
      <td>225114</td>
      <td>24595</td>
      <td>1249400</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>4610055</td>
      <td>339262</td>
      <td>97434</td>
      <td>156067</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>998411</td>
      <td>1620852</td>
      <td>111551</td>
      <td>7512130</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>2015344</td>
      <td>3174280</td>
      <td>216855</td>
      <td>14529105</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>177811</td>
      <td>263712</td>
      <td>18933</td>
      <td>1222961</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>300809</td>
      <td>464717</td>
      <td>35299</td>
      <td>2170325</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>541159</td>
      <td>906438</td>
      <td>64957</td>
      <td>4121493</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1529519</td>
      <td>2534073</td>
      <td>163513</td>
      <td>11116734</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>254405</td>
      <td>359571</td>
      <td>25278</td>
      <td>1671724</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>403448</td>
      <td>612122</td>
      <td>45543</td>
      <td>2845316</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>515018</td>
      <td>828881</td>
      <td>60165</td>
      <td>3974453</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>1759670</td>
      <td>2688797</td>
      <td>185216</td>
      <td>12709995</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>574780</td>
      <td>868878</td>
      <td>64039</td>
      <td>4039863</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>289257</td>
      <td>463702</td>
      <td>34474</td>
      <td>2130458</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>1203660</td>
      <td>2023146</td>
      <td>140488</td>
      <td>9296384</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1366085</td>
      <td>2067734</td>
      <td>136253</td>
      <td>9685610</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>689200</td>
      <td>1124812</td>
      <td>80228</td>
      <td>4970823</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>681749</td>
      <td>1141312</td>
      <td>81756</td>
      <td>5244232</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>1340577</td>
      <td>2158217</td>
      <td>145493</td>
      <td>9866287</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>480833</td>
      <td>734619</td>
      <td>46107</td>
      <td>3341744</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>1464937</td>
      <td>2372469</td>
      <td>167186</td>
      <td>10741742</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>1476767</td>
      <td>2397611</td>
      <td>164699</td>
      <td>10933152</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>603693</td>
      <td>952628</td>
      <td>58640</td>
      <td>4197546</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>1725</td>
      <td>3403</td>
      <td>160</td>
      <td>14240</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>2254160</td>
      <td>3788985</td>
      <td>260056</td>
      <td>16720851</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>921444</td>
      <td>1446868</td>
      <td>94106</td>
      <td>6500578</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>4926</td>
      <td>10393</td>
      <td>564</td>
      <td>42627</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>744</td>
      <td>1414</td>
      <td>71</td>
      <td>5048</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>4618058</td>
      <td>7477213</td>
      <td>507023</td>
      <td>33336460</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>127894</td>
      <td>200611</td>
      <td>12915</td>
      <td>912617</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>516</td>
      <td>3044</td>
      <td>74</td>
      <td>3042</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>133</td>
      <td>318</td>
      <td>15</td>
      <td>667</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>4217865</td>
      <td>6661000</td>
      <td>456497</td>
      <td>30463880</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>111984</td>
      <td>180664</td>
      <td>12107</td>
      <td>813037</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>668</td>
      <td>1332</td>
      <td>31</td>
      <td>3932</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>93</td>
      <td>322</td>
      <td>3</td>
      <td>274</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>4079203</td>
      <td>6589682</td>
      <td>446180</td>
      <td>29599882</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>112596</td>
      <td>176062</td>
      <td>12330</td>
      <td>808521</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>63</td>
      <td>933</td>
      <td>14</td>
      <td>376</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>23</td>
      <td>149</td>
      <td>3</td>
      <td>58</td>
    </tr>
    <tr>
      <th rowspan="40" valign="top">pool2</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>281596</td>
      <td>262324</td>
      <td>23084</td>
      <td>1300031</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>534201</td>
      <td>704863</td>
      <td>51101</td>
      <td>3344208</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>545101</td>
      <td>786456</td>
      <td>56829</td>
      <td>3751848</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>967626</td>
      <td>1336114</td>
      <td>94091</td>
      <td>6386601</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>128010</td>
      <td>196606</td>
      <td>14748</td>
      <td>896466</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>133981</td>
      <td>206771</td>
      <td>15860</td>
      <td>943101</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>300923</td>
      <td>467663</td>
      <td>31834</td>
      <td>2151082</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1953879</td>
      <td>2843984</td>
      <td>185404</td>
      <td>12863647</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>85795</td>
      <td>127913</td>
      <td>9103</td>
      <td>569897</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>334913</td>
      <td>520630</td>
      <td>37628</td>
      <td>2358046</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>425215</td>
      <td>677159</td>
      <td>46782</td>
      <td>2995353</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>2075395</td>
      <td>3059196</td>
      <td>193053</td>
      <td>13465243</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>311823</td>
      <td>486341</td>
      <td>34991</td>
      <td>2196343</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>213261</td>
      <td>336366</td>
      <td>25158</td>
      <td>1532938</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>513518</td>
      <td>794521</td>
      <td>54255</td>
      <td>3628148</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1808572</td>
      <td>2615694</td>
      <td>165119</td>
      <td>11686366</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>289896</td>
      <td>396783</td>
      <td>36709</td>
      <td>2015253</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>261532</td>
      <td>343549</td>
      <td>32662</td>
      <td>1852702</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>690448</td>
      <td>974481</td>
      <td>77236</td>
      <td>4549941</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>257532</td>
      <td>356417</td>
      <td>21821</td>
      <td>1602651</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>629183</td>
      <td>978054</td>
      <td>67938</td>
      <td>4439191</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>711720</td>
      <td>940684</td>
      <td>74526</td>
      <td>4692489</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>280143</td>
      <td>338405</td>
      <td>26740</td>
      <td>1733861</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>1742</td>
      <td>2274</td>
      <td>92</td>
      <td>8903</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>1570493</td>
      <td>766542</td>
      <td>222027</td>
      <td>9189833</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>358980</td>
      <td>366112</td>
      <td>38877</td>
      <td>2235493</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>6748</td>
      <td>6704</td>
      <td>446</td>
      <td>27398</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>115</td>
      <td>101</td>
      <td>11</td>
      <td>299</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>7713326</td>
      <td>602714</td>
      <td>171070</td>
      <td>301735</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>95457</td>
      <td>51875</td>
      <td>14500</td>
      <td>614041</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>388</td>
      <td>625</td>
      <td>31</td>
      <td>2030</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>45</td>
      <td>393</td>
      <td>8</td>
      <td>94</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>3053115</td>
      <td>4387955</td>
      <td>304440</td>
      <td>20387637</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>123473</td>
      <td>170434</td>
      <td>12113</td>
      <td>820127</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>660</td>
      <td>1142</td>
      <td>52</td>
      <td>4156</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>1107</td>
      <td>346</td>
      <td>18</td>
      <td>170</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>3497367</td>
      <td>4816158</td>
      <td>346977</td>
      <td>23173990</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>132860</td>
      <td>184550</td>
      <td>13693</td>
      <td>884810</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>454</td>
      <td>6330</td>
      <td>27</td>
      <td>279</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>110</td>
      <td>287</td>
      <td>8</td>
      <td>529</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
