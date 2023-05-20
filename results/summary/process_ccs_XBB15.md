<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Process-CCSs" data-toc-modified-id="Process-CCSs-1">Process CCSs</a></span><ul class="toc-item"><li><span><a href="#Setup" data-toc-modified-id="Setup-1.1">Setup</a></span></li><li><span><a href="#PacBio-amplicons" data-toc-modified-id="PacBio-amplicons-1.2">PacBio amplicons</a></span></li><li><span><a href="#CCS-stats-for-PacBio-runs" data-toc-modified-id="CCS-stats-for-PacBio-runs-1.3">CCS stats for PacBio runs</a></span></li><li><span><a href="#Align-CCSs-to-amplicons" data-toc-modified-id="Align-CCSs-to-amplicons-1.4">Align CCSs to amplicons</a></span></li><li><span><a href="#Write-valid-CCSs" data-toc-modified-id="Write-valid-CCSs-1.5">Write valid CCSs</a></span></li></ul></li></ul></div>

This Python Jupyter notebook processes the PacBio circular consensus sequences (CCSs) to extract barcodes and call mutations in the gene. It then builds consensus sequences for barcoded variants from the mutations, to build a barcode-variant lookup table.

# Process CCSs
First, process the PacBio CCSs to extract barcodes and call mutations in the gene.

Define the background for CCS mapping


```python
background = "XBB15"
```

## Setup

Import Python modules

Plotting is done with [plotnine](https://plotnine.readthedocs.io/en/stable/), which uses ggplot2-like syntax.

The analysis uses the Bloom lab's [alignparse](https://jbloomlab.github.io/alignparse) and [dms_variants](https://jbloomlab.github.io/dms_variants) packages.


```python
import collections
import math
import os
import re
import time
import warnings

import alignparse
import alignparse.ccs
from alignparse.constants import CBPALETTE
import alignparse.minimap2
import alignparse.targets
import alignparse.consensus
import alignparse.utils

import dms_variants
import dms_variants.codonvarianttable
import dms_variants.plotnine_themes
import dms_variants.utils

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
print(f"Using pandas version {pd.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.2
    Using pandas version 1.3.5


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory for figures:


```python
os.makedirs(config['process_ccs_dir'], exist_ok=True)
```

## PacBio reads
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
pacbio_runs = (
    pd.read_csv(config['pacbio_runs'], dtype=str)
    .drop(columns=['ccs'])
    .assign(name=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'],
            fastq=lambda x: config['ccs_dir'] + '/' + x['name'] + '_ccs.fastq.gz'
            )
    )
 
pacbio_runs.drop(pacbio_runs[pacbio_runs['bg'] != background].index,inplace=True)

display(HTML(pacbio_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>bg</th>
      <th>run</th>
      <th>name</th>
      <th>fastq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>XBB15</td>
      <td>230520</td>
      <td>pool1_XBB15_230520</td>
      <td>results/ccs/pool1_XBB15_230520_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>XBB15</td>
      <td>230520</td>
      <td>pool2_XBB15_230520</td>
      <td>results/ccs/pool2_XBB15_230520_ccs.fastq.gz</td>
    </tr>
  </tbody>
</table>


## PacBio amplicons
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
print(f"Reading amplicons from {config['amplicons_' + background]}")
print(f"Reading feature parse specs from {config['feature_parse_specs_' + background]}")

targets = alignparse.targets.Targets(
                seqsfile=config['amplicons_' + background],
                feature_parse_specs=config['feature_parse_specs_' + background])
```

    Reading amplicons from data/PacBio_amplicon_XBB15.gb
    Reading feature parse specs from data/feature_parse_specs_XBB15.yaml


Draw the target amplicons:


```python
fig = targets.plot(ax_width=7,
                   plots_indexing='biopython',  # numbering starts at 0
                   ax_height=2,  # height of each plot
                   hspace=1.2,  # vertical space between plots
                   )

plotfile = os.path.join(config['process_ccs_dir'], 'amplicons_'+background+'.pdf')
print(f"Saving plot to {plotfile}")
fig.savefig(plotfile, bbox_inches='tight')
```

    Saving plot to results/process_ccs/amplicons_XBB15.pdf



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_22_1.png)
    


Write out the specs used to parse the features (these are the same specs provided as `feature_parse_specs` when initializing `targets`, but with defaults filled in):


```python
print(targets.feature_parse_specs('yaml'))
```

    XBB15:
      query_clip5: 4
      query_clip3: 4
      termini5:
        filter:
          clip5: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip3: 0
        return: []
      gene:
        filter:
          mutation_nt_count: 90
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - mutations
        - accuracy
      spacer:
        filter:
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return: []
      barcode:
        filter:
          mutation_nt_count: 0
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - sequence
        - accuracy
      termini3:
        filter:
          clip3: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
        return: []
    


## Align CCSs to amplicons
We now align the CCSs to the amplicon and parse features from the resulting alignments using the specs above.

First, we initialize an `alignparse.minimap2.Mapper` to align the reads to SAM files:


```python
mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)

print(f"Using `minimap2` {mapper.version} with these options:\n" +
      ' '.join(mapper.options))
```

    Using `minimap2` 2.18-r1015 with these options:
    -A2 -B4 -O12 -E2 --end-bonus=13 --secondary=no --cs


Next, we use `Targets.align_and_parse` to create the alignments and parse them:


```python
readstats, aligned, filtered = targets.align_and_parse(
        df=pacbio_runs,
        mapper=mapper,
        outdir=config['process_ccs_dir'],
        name_col='run',
        group_cols=['name', 'library'],
        queryfile_col='fastq',
        overwrite=True,
        ncpus=config['max_cpus'],
        )
```

First, examine the read stats from the alignment / parsing, both extracting alignment target name and getting stats aggregated by target:


```python
readstats = (
    readstats
    .assign(category_all_targets=lambda x: x['category'].str.split().str[0],
            target=lambda x: x['category'].str.split(None, 1).str[1],
            valid=lambda x: x['category_all_targets'] == 'aligned')
    )
```

Now plot the read stats by run (combining all targets and libraries within a run):


```python
ncol = 2
p = (
    ggplot(readstats
           .groupby(['name', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.85 * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_32_0.png)
    


And the read stats by library (combining all targets and runs within a library):


```python
p = (
    ggplot(readstats
           .groupby(['library', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.5 * pacbio_runs['library'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_34_0.png)
    


And the number of reads by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('target', 'count')) +
    geom_point(stat='identity', size=3) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_y_log10(name='number of reads')
    )
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_36_0.png)
    


And read stats by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index()
           .assign(total=lambda x: x.groupby('target')['count'].transform('sum'),
                   frac=lambda x: x['count'] / x['total'],
                   ), 
           aes('target', 'frac', fill='valid')) +
    geom_bar(stat='identity') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.5 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_38_0.png)
    


Now let's see **why** we filtered the reads.
First, we do some transformations on the `filtered` dict returned by `Targets.align_and_parse`.
Then we count up the number of CCSs filtered for each reason, and group together "unusual" reasons that represent less than some fraction of all filtering.
For now, we group together all targets to the stats represent all targets combined:


```python
other_cutoff = 0.02  # group as "other" reasons with <= this frac

filtered_df = (
    pd.concat(df.assign(target=target) for target, df in filtered.items())
    .groupby(['library', 'name', 'run', 'filter_reason'])
    .size()
    .rename('count')
    .reset_index()
    .assign(tot_reason_frac=lambda x: (x.groupby('filter_reason')['count']
                                       .transform('sum')) / x['count'].sum(),
            filter_reason=lambda x: numpy.where(x['tot_reason_frac'] > other_cutoff,
                                                x['filter_reason'],
                                                'other')
            )
    )
```

Now plot the filtering reason for all runs:


```python
ncol = 7
nreasons = filtered_df['filter_reason'].nunique()

p = (
    ggplot(filtered_df, aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.25 * nreasons * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          )
    )
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_42_0.png)
    


Finally, we take the successfully parsed alignments and read them into a data frame, keeping track of the target that each CCS aligns to.
We also drop the pieces of information we won't use going forward, and rename a few columns:


```python
aligned_df = (
    pd.concat(df.assign(target=target) for target, df in aligned.items())
    .drop(columns=['query_clip5', 'query_clip3', 'run','name'])
    .rename(columns={'barcode_sequence': 'barcode'})
    )

print(f"First few lines of information on the parsed alignments:")
display(HTML(aligned_df.head().to_html(index=False)))
```

    First few lines of information on the parsed alignments:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/4/ccs</td>
      <td>A595T A596T A597T</td>
      <td>0.999692</td>
      <td>TTCGAAATCACTACAA</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/25/ccs</td>
      <td>T367G A368C C369T</td>
      <td>1.000000</td>
      <td>TCAAACACTATATGAT</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/40/ccs</td>
      <td>A163T C164T</td>
      <td>1.000000</td>
      <td>ATGGTCTCAAGCACTG</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/57/ccs</td>
      <td>A214T T215C</td>
      <td>1.000000</td>
      <td>ACGACCGCCATTGCTA</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/61/ccs</td>
      <td>G454A G455A T456A</td>
      <td>1.000000</td>
      <td>TATCACTGACAAATGA</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
  </tbody>
</table>


## Write valid CCSs

Write the processed CCSs to a file:


```python
aligned_df.to_csv(config['processed_ccs_file' + '_' + background], index=False)

print("Barcodes and mutations for valid processed CCSs "
      f"have been written to {config['processed_ccs_file' + '_' + background]}.")
```

    Barcodes and mutations for valid processed CCSs have been written to results/process_ccs/processed_ccs_XBB15.csv.


Next, we analyze these processed CCSs to build the variants.

# Build barcode variant table
Builds consensus sequences for barcoded variants from the mutations called in the processed PacBio CCSs.
Uses these consensus sequences to build a codon variant table.

Make output directories if needed:


```python
os.makedirs(config['variants_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

Read the CSV file with the processed CCSs into a data frame, display first few lines:


```python
processed_ccs = pd.read_csv(config['processed_ccs_file' + '_' + background], na_filter=None)

nlibs = processed_ccs['library'].nunique()  # number of unique libraries

ntargets = processed_ccs['target'].nunique()  # number of unique targets

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")

display(HTML(processed_ccs.head().to_html(index=False)))
```

    Read 1654449 CCSs from 2 libraries and 1 targets.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/4/ccs</td>
      <td>A595T A596T A597T</td>
      <td>0.999692</td>
      <td>TTCGAAATCACTACAA</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/25/ccs</td>
      <td>T367G A368C C369T</td>
      <td>1.000000</td>
      <td>TCAAACACTATATGAT</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/40/ccs</td>
      <td>A163T C164T</td>
      <td>1.000000</td>
      <td>ATGGTCTCAAGCACTG</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/57/ccs</td>
      <td>A214T T215C</td>
      <td>1.000000</td>
      <td>ACGACCGCCATTGCTA</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>m64272e_230518_163630/61/ccs</td>
      <td>G454A G455A T456A</td>
      <td>1.000000</td>
      <td>TATCACTGACAAATGA</td>
      <td>1.0</td>
      <td>XBB15</td>
    </tr>
  </tbody>
</table>


Optional: discard reads that have a `fwd` or `rev` indicator suggestive of heteroduplexes (or we suspect in these runs, some form of pseudo-duplex).


```python
#start = len(processed_ccs)
#processed_ccs = processed_ccs.loc[~processed_ccs['query_name'].str.contains('fwd')]
#processed_ccs = processed_ccs.loc[~processed_ccs['query_name'].str.contains('rev')]
#end = len(processed_ccs)

#print(f"Went from {start} CCSs to {end} after discarded split/duplexed fwd and rev reads.")

#display(HTML(processed_ccs.head().to_html(index=False)))
```

Overall statistics on number of total CCSs and number of unique barcodes:


```python
display(HTML(
    processed_ccs
    .groupby(['target', 'library'])
    .aggregate(total_CCSs=('barcode', 'size'),
               unique_barcodes=('barcode', 'nunique'))
    .assign(avg_CCSs_per_barcode=lambda x: x['total_CCSs'] / x['unique_barcodes'])
    .round(2)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>total_CCSs</th>
      <th>unique_barcodes</th>
      <th>avg_CCSs_per_barcode</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">XBB15</th>
      <th>pool1</th>
      <td>817300</td>
      <td>101421</td>
      <td>8.06</td>
    </tr>
    <tr>
      <th>pool2</th>
      <td>837149</td>
      <td>120358</td>
      <td>6.96</td>
    </tr>
  </tbody>
</table>


### Filter processed CCSs
We have the PacBio `ccs` program's estimated "accuracy" for both the barcode and the gene sequence for each processed CCS.
We will filter the CCSs to only keep ones of sufficiently high accuracy.

First, we want to plot the accuracies.
It is actually visually easier to look at the error rate, which is one minus the accuracy.
Because we want to plot on a log scale (which can't show error rates of zero), we define a *error_rate_floor*, and set all error rates less than this to that value:


```python
error_rate_floor = 1e-7  # error rates < this set to this
if error_rate_floor >= config['max_error_rate']:
    raise ValueError('error_rate_floor must be < max_error_rate')

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: numpy.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: numpy.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            )
    )
```

Now plot the error rates, drawing a dashed vertical line at the threshold separating the CCSs we retain for consensus building versus those that we discard:


```python
_ = (
 ggplot(processed_ccs
        .melt(value_vars=['barcode_error', 'gene_error'],
              var_name='feature_type', value_name='error rate'),
        aes('error rate')) +
 geom_histogram(bins=25) +
 geom_vline(xintercept=config['max_error_rate'],
            linetype='dashed',
            color=CBPALETTE[1]) +
 facet_wrap('~ feature_type') +
 theme(figure_size=(4.5, 2)) +
 ylab('number of CCSs') +
 scale_x_log10()
 ).draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_61_0.png)
    


Flag the CCSs to retain, and indicate how many we are retaining and purging due to the accuracy filter:


```python
processed_ccs = (
    processed_ccs
    .assign(retained=lambda x: ((x['gene_error'] < config['max_error_rate']) &
                                (x['barcode_error'] < config['max_error_rate'])))
    )
```

Here are number of retained CCSs:


```python
_ = (
 ggplot(processed_ccs.assign(xlabel=lambda x: x['target'] + ', ' + x['library'])
                     .groupby(['xlabel', 'retained'])
                     .size()
                     .rename('count')
                     .reset_index(),
        aes('xlabel', 'count', color='retained', label='count')) +
 geom_point(size=3) +
 geom_text(va='bottom', size=7, ha='center',format_string='{:.3g}', nudge_y=0.2) +
 theme(figure_size=(0.5 * nlibs * ntargets, 3),
       panel_grid_major_x=element_blank(),
       axis_text_x=element_text(angle=90),
       ) +
 scale_y_log10(name='number of CCSs') +
 xlab('') +
 scale_color_manual(values=CBPALETTE[1:])
 ).draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_65_0.png)
    


### Sequences per barcode
How many times is each barcode sequenced?
This is useful to know for thinking about building the barcode consensus.

Plot the distribution of the number of times each **barcode** is observed among the retained CCSs:


```python
max_count = 8 # in plot, group all barcodes with >= this many counts

p = (
 ggplot(
    processed_ccs
     .query('retained')
     .groupby(['library', 'barcode'])
     .size()
     .rename('nseqs')
     .reset_index()
     .assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_count)),
    aes('nseqs')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of barcodes') +
 xlab('CCSs for barcode')
 )

_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_67_0.png)
    


### Empirical accuracy of CCSs
We want to directly estimate the accuracy of the gene-barcode link rather than relying on the PacBio `ccs` accuracy, which doesn't include inaccuracies due to things like strand exchange or the same barcode on different sequences.

One way to do this is to examine instances when we have multiple sequences for the same barcode. 
We can calculate the empirical accuracy of the sequences by looking at all instances of multiple sequences of the same barcode and determining how often they are identical.
This calculation is performed by `alignparse.consensus.empirical_accuracy` using the equations described in the docs for that function.

We will do this four for sets of sequences:

 1. All of the CCSs retained above.
 2. CCSs retained by applying a PacBio `ccs` accuracy filter 10-fold more stringent than the one above.
    The rationale is that if this improves the concordance (real accuracy) of the CCSs substantially then maybe we should make the accuracy filter more stringent.
 3. Like (1) but excluding all CCSs with indels.
    the rationale is that we only really care about substitutions, and will exclude sequences with indels anyway.
 4. Like (2) but excluding all CCSs with indels.
 
First, we annotate the sequences with the number of indels and whether they have an indel to enable categorization into the aforementioned sets. We want to convert the way we call *only* in-frame deletions to be substitution style (that is, `del511to513` becomes `C511- C512- G513-`), if the deletion is a single-codon in-frame deletion. We do this using alignparse util InFrameDeletionsToSubs. Below, when we parse indels, these in-frame deletions will be treated as "substitutions" and not add to indel counts, so we can then parse on indel count which will be count of out-of-frame indels.:


```python
targets = alignparse.targets.Targets(seqsfile=config['amplicons_' + background],
                                     feature_parse_specs=config['feature_parse_specs_' + background])
geneseq = targets.get_target(background).get_feature('gene').seq

print(f"Read gene of {len(geneseq)} nts for {background} from {config['amplicons_' + background]}")
```

    Read gene of 603 nts for XBB15 from data/PacBio_amplicon_XBB15.gb



```python
processed_ccs = processed_ccs.reset_index(drop=True)

deltosubs = alignparse.utils.InFrameDeletionsToSubs(geneseq)

processed_ccs['gene_mutations_with_del'] = processed_ccs.apply(lambda x: deltosubs.dels_to_subs(x['gene_mutations']), axis=1)

processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations_with_del',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)

processed_ccs.head(n=12)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
      <th>barcode_error</th>
      <th>gene_error</th>
      <th>retained</th>
      <th>gene_mutations_with_del</th>
      <th>n_indels</th>
      <th>has_indel</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/4/ccs</td>
      <td>A595T A596T A597T</td>
      <td>0.999692</td>
      <td>TTCGAAATCACTACAA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>3.075795e-04</td>
      <td>False</td>
      <td>A595T A596T A597T</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/25/ccs</td>
      <td>T367G A368C C369T</td>
      <td>1.000000</td>
      <td>TCAAACACTATATGAT</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>T367G A368C C369T</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/40/ccs</td>
      <td>A163T C164T</td>
      <td>1.000000</td>
      <td>ATGGTCTCAAGCACTG</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.047740e-07</td>
      <td>True</td>
      <td>A163T C164T</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/57/ccs</td>
      <td>A214T T215C</td>
      <td>1.000000</td>
      <td>ACGACCGCCATTGCTA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>A214T T215C</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/61/ccs</td>
      <td>G454A G455A T456A</td>
      <td>1.000000</td>
      <td>TATCACTGACAAATGA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>G454A G455A T456A</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>5</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/71/ccs</td>
      <td>A166C G168A</td>
      <td>1.000000</td>
      <td>AAACTTTCGAGAGCTA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>A166C G168A</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>6</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/78/ccs</td>
      <td>T376A C378G</td>
      <td>0.999834</td>
      <td>CCACCAAACTACTTAA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.658589e-04</td>
      <td>False</td>
      <td>T376A C378G</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>7</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/115/ccs</td>
      <td>T355A A356G C357A</td>
      <td>1.000000</td>
      <td>TCAACCGAAAAAGATC</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>T355A A356G C357A</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>8</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/124/ccs</td>
      <td>A166G A167G G168T</td>
      <td>1.000000</td>
      <td>TAACTTGGAACGAAAA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>A166G A167G G168T</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>9</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/127/ccs</td>
      <td>A310G T311C</td>
      <td>1.000000</td>
      <td>TCCAGGTATACCTTGA</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>A310G T311C</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>10</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/140/ccs</td>
      <td>T85A C86G T87A</td>
      <td>1.000000</td>
      <td>AAAAAAAAACACGAAT</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>T85A C86G T87A</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11</th>
      <td>pool1</td>
      <td>m64272e_230518_163630/160/ccs</td>
      <td>T532G C534A</td>
      <td>1.000000</td>
      <td>GTGATAAGGTTAGAGT</td>
      <td>1.0</td>
      <td>XBB15</td>
      <td>1.000000e-07</td>
      <td>1.000000e-07</td>
      <td>True</td>
      <td>T532G C534A</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



Plot how many sequences have indels:


```python
_ = (
 ggplot(processed_ccs,
        aes('retained', fill='has_indel')) +
 geom_bar(position='dodge') +
 geom_text(aes(label='..count..'), stat='count', va='bottom', size=7,
           position=position_dodge(width=0.9), format_string='{:.2g}') +
 theme(figure_size=(2.5 * nlibs, 3),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of CCSs') +
 scale_fill_manual(values=CBPALETTE[1:]) +
 facet_wrap('~ library', nrow=1)
 ).draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_72_0.png)
    


Now get the empirical accuracy for each of the CCS groups mentioned above:


```python
high_acc = config['max_error_rate'] / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)
    
    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: numpy.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )
```

Display table of the empirical accuracies:


```python
display(HTML(empirical_acc.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>accuracy</th>
      <th>description</th>
      <th>number_CCSs</th>
      <th>actual</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>0.986468</td>
      <td>retained</td>
      <td>712159</td>
      <td>False</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>0.985223</td>
      <td>retained</td>
      <td>728329</td>
      <td>False</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>0.992330</td>
      <td>retained, no indel</td>
      <td>704757</td>
      <td>True</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>0.990594</td>
      <td>retained, no indel</td>
      <td>720566</td>
      <td>True</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>0.988251</td>
      <td>10X accuracy</td>
      <td>685682</td>
      <td>False</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>0.987000</td>
      <td>10X accuracy</td>
      <td>701223</td>
      <td>False</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>0.992367</td>
      <td>10X accuracy, no indel</td>
      <td>679760</td>
      <td>False</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>0.990678</td>
      <td>10X accuracy, no indel</td>
      <td>694940</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Plot the empirical accuracies, using a different color to show the category that we will actually use:


```python
p = (
    ggplot(empirical_acc,
           aes('description', 'accuracy', color='actual', label='accuracy')
           ) +
    geom_point(size=3) +
    geom_text(va='bottom', size=9, format_string='{:.3g}', nudge_y=0.003) +
    facet_wrap('~ library', ncol=8) +
    theme(figure_size=(1.75 * nlibs, 2.25),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank(),
          ) +
    xlab('') +
    scale_y_continuous(name='empirical accuracy', limits=(0.2, 1.005)) +
    scale_color_manual(values=CBPALETTE, guide=False)
    )

plotfile = os.path.join(config['figs_dir'], 'empirical_CCS_accuracy.pdf')
print(f"Saving plot to {plotfile}")
_ = p.draw()
```

    Saving plot to results/figures/empirical_CCS_accuracy.pdf



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_78_1.png)
    


The above analysis shows that if we exclude sequences with indels (which we plan to do among our consensus sequences), then the accuracy of each CCS is around 99%. 
We do **not** get notably higher empirical accuracy by imposing a more stringent filter from the PacBio `ccs` program, indicating that the major sources of error are due to processes that are not modeled in this program's accuracy filter (perhaps strand exchange or barcode sharing).

Note that this empirical accuracy is for a **single** CCS.
When we build the consensus sequences for each barcode below, we will take the consensus of CCSs within a barcode.
So for barcodes with multiple CCSs, the actual accuracy of the consensus sequences will be higher than the empirical accuracy above due to capturing information from multiple CCSs.

### Consensus sequences for barcodes
We call the consensus sequence for each barcode using the simple method implemented in [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus).
The documentation for that function explains the method in detail, but basically it works like this:
 1. When there is just one CCS per barcode, the consensus is just that sequence.
 2. When there are multiple CCSs per barcode, they are used to build a consensus--however, the entire barcode is discarded if there are many differences between CCSs with the barcode, or high-frequency non-consensus mutations. The reason that barcodes are discarded in such cases as many differences between CCSs or high-frequency non-consensus mutations suggest errors such as barcode collisions or strand exchange.
 
First, call the consensus for each barcode including **all** retained sequences, even those with undesirable indels.


```python
consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations_with_del',
                        )
```

Here are the first few lines of the data frame of consensus sequences for each barcode.
In addition to giving the library, barcode, target, and mutations, it also has a column indicating how many CCSs support the variant call:


```python
display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations_with_del</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAAATGATT</td>
      <td>XBB15</td>
      <td>G574A C575T T576G</td>
      <td>5</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAACAGCGA</td>
      <td>XBB15</td>
      <td>A340G A341C G342T</td>
      <td>30</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAATACGTG</td>
      <td>XBB15</td>
      <td>A325G C327T</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAACACGAAT</td>
      <td>XBB15</td>
      <td>T85A C86G T87A</td>
      <td>8</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAGACAACT</td>
      <td>XBB15</td>
      <td>A424G T425A C426A</td>
      <td>3</td>
    </tr>
  </tbody>
</table>


Since we retain variants with substitutions and single-codon deletions, add information about substitution mutations, and for troubleshooting, number of indels:


```python
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations_with_del',
                    sub_str_col='substitutions',
                    del_str_col='deletions',
                    ins_str_col='insertions',
                    n_indel_col='number_of_indels',
                    overwrite_cols=True)

display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations_with_del</th>
      <th>variant_call_support</th>
      <th>substitutions</th>
      <th>deletions</th>
      <th>insertions</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAAATGATT</td>
      <td>XBB15</td>
      <td>G574A C575T T576G</td>
      <td>5</td>
      <td>G574A C575T T576G</td>
      <td></td>
      <td></td>
      <td>0</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAACAGCGA</td>
      <td>XBB15</td>
      <td>A340G A341C G342T</td>
      <td>30</td>
      <td>A340G A341C G342T</td>
      <td></td>
      <td></td>
      <td>0</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAATACGTG</td>
      <td>XBB15</td>
      <td>A325G C327T</td>
      <td>3</td>
      <td>A325G C327T</td>
      <td></td>
      <td></td>
      <td>0</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAACACGAAT</td>
      <td>XBB15</td>
      <td>T85A C86G T87A</td>
      <td>8</td>
      <td>T85A C86G T87A</td>
      <td></td>
      <td></td>
      <td>0</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAAAAAGACAACT</td>
      <td>XBB15</td>
      <td>A424G T425A C426A</td>
      <td>3</td>
      <td>A424G T425A C426A</td>
      <td></td>
      <td></td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Plot distribution of number of CCSs supporting each variant call (consensus), indicating whether or not there is an indel:


```python
max_variant_call_support = 6  # group variants with >= this much support

_ = (
 ggplot(consensus
        .assign(variant_call_support=lambda x: numpy.clip(x['variant_call_support'],
                                                          None,
                                                          max_variant_call_support),
                indel_state=lambda x: numpy.where(x['number_of_indels'] > 0,
                                                  'has indel', 'no indel')
                ),
        aes('variant_call_support')) +
 geom_bar() +
 ylab('number of variants') +
 facet_grid('indel_state ~ library') +
 theme(figure_size=(1.75 * nlibs, 3.5),
       panel_grid_major_x=element_blank(),
       ) 
 ).draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_87_0.png)
    


We see that most variant consensus sequences do **not** have indels, especially if we limit to the more "accurate" ones that have multiple CCSs supporting them.

We will ignore all consensus sequences with indels in the variant-barcode lookup table. 
We do this for two reasons:
 1. When there is just one CCS supporting a consensus, it is less likely to be accurate as indels are the main mode of PacBio error.
 2. For the purposes of our studies, we are interested in point mutations rather than indels anyway.
 
Here are number of valid consensus sequence (no indels) for each library and target:


```python
consensus = consensus.query('number_of_indels < 1')

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )

display(HTML(lib_target_counts.to_html(index=False)))

p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )

_ = p.draw()
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>target</th>
      <th>consensus sequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>XBB15</td>
      <td>92846</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>XBB15</td>
      <td>109743</td>
    </tr>
  </tbody>
</table>



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_89_1.png)
    


Below we write the retained consensus sequences to a CSV file that links the nucleotide mutations to the barcodes.
(The next section analyzes this variant table in detail, and provides have more precise information on the number of variants and relevant statistics):


```python
print(f"Writing nucleotide variants to {config['nt_variant_table_file' + '_' + background]}")
      
(consensus
 [['target', 'library', 'barcode', 'gene_mutations_with_del', 'substitutions', 'insertions','deletions', 'variant_call_support']]
 .to_csv(config['nt_variant_table_file' + '_' + background], index=False)
 )
      
print('Here are the first few lines of this file:')
display(HTML(
    pd.read_csv(config['nt_variant_table_file' + '_' + background], na_filter=None)
    .head()
    .to_html(index=False)
    ))
```

    Writing nucleotide variants to results/variants/nucleotide_variant_table_XBB15.csv
    Here are the first few lines of this file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>gene_mutations_with_del</th>
      <th>substitutions</th>
      <th>insertions</th>
      <th>deletions</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAAATGATT</td>
      <td>G574A C575T T576G</td>
      <td>G574A C575T T576G</td>
      <td></td>
      <td></td>
      <td>5</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAACAGCGA</td>
      <td>A340G A341C G342T</td>
      <td>A340G A341C G342T</td>
      <td></td>
      <td></td>
      <td>30</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAATACGTG</td>
      <td>A325G C327T</td>
      <td>A325G C327T</td>
      <td></td>
      <td></td>
      <td>3</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAACACGAAT</td>
      <td>T85A C86G T87A</td>
      <td>T85A C86G T87A</td>
      <td></td>
      <td></td>
      <td>8</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAGACAACT</td>
      <td>A424G T425A C426A</td>
      <td>A424G T425A C426A</td>
      <td></td>
      <td></td>
      <td>3</td>
    </tr>
  </tbody>
</table>


What happened to the barcodes that we "dropped" because we could not construct a reliable consensus?
The `dropped` data frame from [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus) has this information:


```python
display(HTML(dropped.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>drop_reason</th>
      <th>nseqs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAAACAATTAACTA</td>
      <td>XBB15</td>
      <td>minor subs too frequent</td>
      <td>22</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAACGAGCAAGGA</td>
      <td>XBB15</td>
      <td>subs diff too large</td>
      <td>12</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAACGCGCCAGTA</td>
      <td>XBB15</td>
      <td>subs diff too large</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAACGGGACTAAA</td>
      <td>XBB15</td>
      <td>minor subs too frequent</td>
      <td>23</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAAACTAAAAGAGT</td>
      <td>XBB15</td>
      <td>subs diff too large</td>
      <td>5</td>
    </tr>
  </tbody>
</table>


Summarize the information in this data frame on dropped barcodes with the plot below.
This plot shows several things.
First, we see that the total number of barcodes dropped is modest (just a few thousand per library) relative to the total number of barcodes per library (seen above to be on the order of hundreds of thousands).
Second, the main reason that barcodes are dropped is that there are CCSs within the same barcode with suspiciously large numbers of mutations relative to the consensus---which we use as a filter to discard the entire barcode as it could indicate strand exchange or some other issue.
In any case, the modest number of dropped barcodes indicates that there probably isn't much of a need to worry: 


```python
max_nseqs = 8  # plot together all barcodes with >= this many sequences

_ = (
 ggplot(
    dropped.assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_nseqs)),
    aes('nseqs')) + 
 geom_bar() + 
 scale_x_continuous(limits=(1, None)) +
 xlab('number of sequences for barcode') +
 ylab('number of barcodes') +
 facet_grid('library ~ drop_reason') +
 theme(figure_size=(10, 1.5 * nlibs),
       panel_grid_major_x=element_blank(),
       )
 ).draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_95_0.png)
    


## Create barcode-variant table
We now create a [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) that stores and processes all the information about the variant consensus sequences.
Below we initialize such a table, and then analyze information about its composition. **Note, I believe this table eliminates our indel mutants, so we probably need to go back to the nt variants table to parse indels ourselves. Nonetheless, we'll keep this here so we can see what our (amino acid) mutant coverage looks like**

### Initialize codon variant table
In order to initialize the codon variant table, we need two pieces of information:
  1. The wildtype gene sequence.
  2. The list of nucleotide mutations for each variant as determined in the consensus calling above.

Read "wildtype" gene sequence to which we made the alignments (in order to do this, initialize an `alignparse.Targets` and get the gene sequence from it):


```python
targets = alignparse.targets.Targets(seqsfile=config['amplicons_' + background],
                                     feature_parse_specs=config['feature_parse_specs_' + background])
geneseq = targets.get_target(background).get_feature('gene').seq

print(f"Read gene of {len(geneseq)} nts for {background} from {config['amplicons_' + background]}")
```

    Read gene of 603 nts for XBB15 from data/PacBio_amplicon_XBB15.gb


Now initialize the codon variant table using this wildtype sequence and our list of nucleotide mutations for each variant:


```python
variants = dms_variants.codonvarianttable.CodonVariantTable(
                barcode_variant_file=config['nt_variant_table_file_' + background],
                geneseq=geneseq,
                primary_target=background,
                allowgaps=True,
                )

variants.barcode_variant_df.head(n=12)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
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
      <th>0</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAAATGATT</td>
      <td>5</td>
      <td>GCT192ATG</td>
      <td>A192M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAACAGCGA</td>
      <td>30</td>
      <td>AAG114GCT</td>
      <td>K114A</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAATACGTG</td>
      <td>3</td>
      <td>AAC109GAT</td>
      <td>N109D</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAACACGAAT</td>
      <td>8</td>
      <td>TCT29AGA</td>
      <td>S29R</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAGACAACT</td>
      <td>3</td>
      <td>ATC142GAA</td>
      <td>I142E</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAACAAGATAA</td>
      <td>10</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAACCACAAAT</td>
      <td>1</td>
      <td>CCA197TGG</td>
      <td>P197W</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAACCGACAGC</td>
      <td>3</td>
      <td>GTT180TTG</td>
      <td>V180L</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAACGATCCCA</td>
      <td>5</td>
      <td>GGT196TGT</td>
      <td>G196C</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>9</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAACTATTCGT</td>
      <td>3</td>
      <td>GGT86TAT</td>
      <td>G86Y</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAGTTGATCA</td>
      <td>5</td>
      <td>GTT153TCT</td>
      <td>V153S</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAGTTTAAAC</td>
      <td>4</td>
      <td>TTC17GCT</td>
      <td>F17A</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



### Basic stats on variants
We now will analyze the variants.
In this call and in the plots below, we set `samples=None` as we aren't looking at variant counts in specific samples, but are simply looking at properties of the variants in the table.

Here are the number of variants for each target:


```python
display(HTML(
    variants
    .n_variants_df(samples=None)
    .pivot_table(index=['target'],
                 columns='library',
                 values='count')
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>pool1</th>
      <th>pool2</th>
      <th>all libraries</th>
    </tr>
    <tr>
      <th>target</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>XBB15</th>
      <td>92846</td>
      <td>109743</td>
      <td>202589</td>
    </tr>
  </tbody>
</table>


Plot the number of variants supported by each number of CCSs:


```python
max_support = 10  # group variants with >= this much support

p = variants.plotVariantSupportHistogram(max_support=max_support,
                                         widthscale=1.1,
                                         heightscale=0.9)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_103_0.png)
    


### Mutations per variant
Plot the number of barcoded variants with each number of amino-acid and codon mutations.
This is for the primary target only, and doesn't include the spiked-in secondary targets:


```python
max_muts = 7  # group all variants with >= this many mutations

for mut_type in ['aa', 'codon']:
    p = variants.plotNumMutsHistogram(mut_type, samples=None, max_muts=max_muts,
                                      widthscale=1.1,
                                      heightscale=0.9)
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    _ = p.draw()
    plotfile = os.path.join(config['figs_dir'], f"n_{mut_type}_muts_per_variant_"+background+".pdf")
    print(f"Saving plot to {plotfile}")
    p.save(plotfile)
```

    Saving plot to results/figures/n_aa_muts_per_variant_XBB15.pdf
    Saving plot to results/figures/n_codon_muts_per_variant_XBB15.pdf



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_105_1.png)
    



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_105_2.png)
    


Plot the frequencies of different codon mutation types among **all** variants (any number of mutations), again only for primary target:


```python
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant',
                                    heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
plotfile = os.path.join(config['figs_dir'], f"avg_muts_per_variant_"+background+".pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile)
```

    Saving plot to results/figures/avg_muts_per_variant_XBB15.pdf



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_107_1.png)
    


Variants supported by multiple PacBio CCSs should have fewer spurious mutations since sequencing errors are very unlikely to occur on two CCSs.
Below we plot the number of codon mutations per variant among variants with at least two CCSs supporting their call.
The difference in mutation rates here and in the plot above (that does not apply the `min_support=2` filter) gives some estimate of the frequency of mutations in our variants our spurious.
In fact, we see the numbers are very similar, indicating that few of the mutations are spurious:


```python
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant', 
                                    min_support=2, heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_109_0.png)
    


### Completeness of mutation sampling
We examine how completely amino-acid mutations are sampled by the variants for the primary target, looking at single-mutant variants only and all variants.
The plot below shows that virtually every mutation is found in a variant in each library, even if we just look among the single mutants.
Things look especially good if we aggregate across libraries:


```python
for variant_type in ['all', 'single']:
    p = variants.plotCumulMutCoverage(variant_type, mut_type='aa', samples=None)
    _ = p.draw()
    plotfile = os.path.join(config['figs_dir'],
                            f"variant_cumul_{variant_type}_mut_coverage_"+background+".pdf")
    print(f"Saving plot to {plotfile}")
    p.save(plotfile)
```

    Saving plot to results/figures/variant_cumul_all_mut_coverage_XBB15.pdf
    Saving plot to results/figures/variant_cumul_single_mut_coverage_XBB15.pdf



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_111_1.png)
    



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_111_2.png)
    


To get more quantitative information like that plotted above, we determine how many mutations are found 0, 1, or >1 times both among single and all mutants for the primary target:


```python
count_dfs = []
for variant_type in ['all', 'single']:
    i_counts = (variants.mutCounts(variant_type, mut_type='aa', samples=None)
                .assign(variant_type=variant_type)
                )
    count_dfs += [i_counts.assign(include_stops=True),
                  i_counts
                  .query('not mutation.str.contains("\*")', engine='python')
                  .assign(include_stops=False)
                  ]
    
display(HTML(
    pd.concat(count_dfs)
    .assign(count=lambda x: (numpy.clip(x['count'], None, 2)
                             .map({0: '0', 1: '1', 2:'>1'}))
            )
    .groupby(['variant_type', 'include_stops', 'library', 'count'])
    .aggregate(number_of_mutations=pd.NamedAgg(column='mutation', aggfunc='count'))
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th>number_of_mutations</th>
    </tr>
    <tr>
      <th>variant_type</th>
      <th>include_stops</th>
      <th>library</th>
      <th>count</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="18" valign="top">all</th>
      <th rowspan="9" valign="top">False</th>
      <th rowspan="3" valign="top">pool1</th>
      <th>0</th>
      <td>10</td>
    </tr>
    <tr>
      <th>1</th>
      <td>5</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4005</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">pool2</th>
      <th>0</th>
      <td>6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>7</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4007</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4011</td>
    </tr>
    <tr>
      <th rowspan="9" valign="top">True</th>
      <th rowspan="3" valign="top">pool1</th>
      <th>0</th>
      <td>150</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4060</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">pool2</th>
      <th>0</th>
      <td>141</td>
    </tr>
    <tr>
      <th>1</th>
      <td>15</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4065</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>135</td>
    </tr>
    <tr>
      <th>1</th>
      <td>14</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4072</td>
    </tr>
    <tr>
      <th rowspan="18" valign="top">single</th>
      <th rowspan="9" valign="top">False</th>
      <th rowspan="3" valign="top">pool1</th>
      <th>0</th>
      <td>12</td>
    </tr>
    <tr>
      <th>1</th>
      <td>6</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4002</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">pool2</th>
      <th>0</th>
      <td>9</td>
    </tr>
    <tr>
      <th>1</th>
      <td>5</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4006</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>8</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4010</td>
    </tr>
    <tr>
      <th rowspan="9" valign="top">True</th>
      <th rowspan="3" valign="top">pool1</th>
      <th>0</th>
      <td>173</td>
    </tr>
    <tr>
      <th>1</th>
      <td>8</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4040</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">pool2</th>
      <th>0</th>
      <td>162</td>
    </tr>
    <tr>
      <th>1</th>
      <td>14</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4045</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>158</td>
    </tr>
    <tr>
      <th>1</th>
      <td>13</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>4050</td>
    </tr>
  </tbody>
</table>


### Mutation frequencies along gene
We plot the frequencies of mutations along the gene among the variants for the primary target.
Ideally, this would be uniform.
We make the plot for both all variants and single-mutant / wildtype variants:


```python
for variant_type in ['all', 'single']:
    p = variants.plotMutFreqs(variant_type, mut_type='codon', samples=None)
    p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_115_0.png)
    



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_115_1.png)
    


We can also use heat maps to examine the extent to which specific amino-acid or codon mutations are over-represented.
These heat maps are large, so we make them just for all variants and the merge of all libraries:


```python
for mut_type in ['aa', 'codon']:
    p = variants.plotMutHeatmap('all', mut_type, samples=None, #libraries='all_only',
                                widthscale=2)
    p.draw()
```


    
![png](process_ccs_XBB15_files/process_ccs_XBB15_117_0.png)
    



    
![png](process_ccs_XBB15_files/process_ccs_XBB15_117_1.png)
    


### Write codon-variant table
We write the codon variant table to a CSV file.
This table looks like this:


```python
display(HTML(
    variants.barcode_variant_df
    .head()
    .to_html(index=False)
    ))
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
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAAATGATT</td>
      <td>5</td>
      <td>GCT192ATG</td>
      <td>A192M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAACAGCGA</td>
      <td>30</td>
      <td>AAG114GCT</td>
      <td>K114A</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAATACGTG</td>
      <td>3</td>
      <td>AAC109GAT</td>
      <td>N109D</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAACACGAAT</td>
      <td>8</td>
      <td>TCT29AGA</td>
      <td>S29R</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>XBB15</td>
      <td>pool1</td>
      <td>AAAAAAAAAGACAACT</td>
      <td>3</td>
      <td>ATC142GAA</td>
      <td>I142E</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Note how this table differs from the nucleotide variant table we generated above and used to initialize the [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) in that it gives **codon** substitutions and associated amino-acid substitutions.

Write it to CSV file:


```python
print(f"Writing codon-variant table to {config['codon_variant_table_file_'+background]}")

variants.barcode_variant_df.to_csv(config['codon_variant_table_file_'+background], index=False)
```

    Writing codon-variant table to results/variants/codon_variant_table_XBB15.csv



```python

```
