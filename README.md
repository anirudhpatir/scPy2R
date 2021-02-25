# scPy2R
Single-cell python packages for R

## Installation
```
devtools::install_github("anirudhpatir/scrublet@main")
```

## Current packages

1. scrublet: adapted from [scfurl](https://github.com/scfurl/m3addon)

## Scrublet  
### Input & output values  
**Input**: count matrix, rows = genes, columns = cells.  
**Output**: list, predicted doublets with their scores from real and simulated data, as well as the threshold applied.  

### Usage
```
doublets = scrublet_R(count_matrix)
```

## Citation
[Wolock, Samuel L., Romain Lopez, and Allon M. Klein. "Scrublet: computational identification of cell doublets in single-cell transcriptomic data." Cell systems 8.4 (2019): 281-291.](https://scholar.google.com/scholar_url?url=https://www.sciencedirect.com/science/article/pii/S2405471218304745&hl=en&sa=T&oi=gsb-ggp&ct=res&cd=0&d=6988818207536851142&ei=IcAhYJfkPJLGmgHPtIiYAw&scisig=AAGBfm3DpZ_-V11xqpwh8SZ2g5hyXDkTPQ)
