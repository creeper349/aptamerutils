# aptamerutils

## Introduction

`aptamerutils` is a simple Python package for short nucleic acid sequences (for example, aptamers) arrangement, sorting, clustering. 

This package can be installed from PyPI:
```
pip install aptamerutils
```
For more about classes and functions of this package, see notebook files in `.docs/`.

## Example

This is a figure of aptamer clustering result. `.pdf`, `.png`, `.svg` files are all supported as output. The colored bases are feature fragments of each cluster (family), as the bolded bases are unpaired predicted by Viennarna when a nucleic acid molecule forms its secondary structure.

![figure](example/results/save.png)

## versions

### 0.1.0

- Initial release
- Clustering of short nucleic acid sequences
- Search desired nucleic acid sequences (logical operations supported)

### 0.2.0

- Webapp for clustering is developed
- Fuzzy matching is supported when trimming
- Kmer library of a `SeqList` object can be directly extracted

### 0.3.0

- Output figure is shown can be shown in a different way, now it is more similar to a clustering result published in a paper of aptamers. Style in 0.1.0 -- 0.2.0 is also preserved
- Fuzzy matching is added to webapp
- Suspicious sequences ( contain characters other than ATCGU ) in a `SeqList` object will be deleted when trimming