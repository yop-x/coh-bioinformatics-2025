# iMARGI RNA–DNA Co-Expression Network Visualization

This repository contains Python functions and utilities to analyze **iMARGI datasets** and visualize **RNA–DNA co-expression networks**.  
The tool enables researchers to explore how RNAs interact with different DNA loci, with customizable thresholds and normalization methods.

---

## Features

- Load and process raw iMARGI datasets
- Normalize interaction frequencies using **Edges Per Million (EPM)**
- Support for both **coordinate-based** and **gene name-based** queries
- Generate network graphs with clear, interpretable labels
- Export results to CSV/Excel for further analysis

---

## Installation

Clone this repository and install the required dependencies:

```bash
git clone https://github.com/yop-x/coh-bioinformatics-2025
cd coh-bioinformatics-2025
pip install -r requirements.txt
```
## Usage 
The plotting function accepts either a gene name (e.g., LINC00607) or an RNA coordinate (e.g., chr2:215,611,563-215,843,535).
Edge weights are EPM (Edges Per Million) and nodes are highlighted when EPM > freq.

### Inputs

**query**: 
- Could be either a gene name or a RNA/DNA coordinate. 

**iMARGI file**: 
- .bedpe/.tsv (optionally gzipped).

- Works with basic 10-column files (no gene names) or richer files with 15 columns that already contain gene info.

**EPM/freq**: 
- Edge weights are EPM (Edges Per Million) and nodes are highlighted when EPM > freq.


```
from coexp_update import get_coexp_network_updated

# 1) Query by gene name
get_coexp_network_updated(
    query="LINC00607",
    iMARGI_files="data/GSM4006842_HUVEC_control_iMARGI.bedpe.gz",
    freq=0.5  # EPM threshold for highlighting
)

# 2) Query by RNA coordinate (commas/en-dash ok)
get_coexp_network_updated(
    query="chr2:215,611,563–215,843,535",
    iMARGI_files="data/GSM4006842_HUVEC_control_iMARGI.bedpe.gz",
    freq=1.0
)
```

*Note: The function expects the gene annotation file at **./genes_df.csv** (i.e., in the same directory where you run the script/notebook).
If your file is elsewhere, either move/copy it into the working directory:*






@article{calandrelli2020stress,
  title={Stress-induced RNA--chromatin interactions promote endothelial dysfunction},
  author={Calandrelli, R. and Xu, L. and Luo, Y. and others},
  journal={Nature Communications},
  volume={11},
  number={1},
  pages={5211},
  year={2020},
  publisher={Nature Publishing Group},
  doi={10.1038/s41467-020-18957-w}
}