# stickleback
![image](https://user-images.githubusercontent.com/10180619/177040107-ba9fc8ca-4571-42a1-a81c-9f52d4a48cb7.png)
Library QC by nanopore using levenshtein distance to map insertion sites on target sequence

Stickleback maps nanopore reads containing insertions (such as molecular handles) and then identifies the insertion site on the template molecule. 

```
python stickleback.0.0.py path/to/samfile.sam queryString path/to/template.fasta [minimumReadLength] [maximumReadLength]
```

