# stickleback
Library QC by nanopore using levenshtein distance to map insertion sites on target sequence

Stickleback maps nanopore reads containing insertions (such as molecular handles) and then identifies the insertion site on the template molecule. 

```
python stickleback.0.0.py path/to/samfile.sam queryString path/to/template.fasta [minimumReadLength] [maximumReadLength]
```

