```
----------------=============------------------
--==--==--==--==   ><```º>   ==--==--==--==--==
==--==--==--==-- stickleback --==--==--==--==--
----------------=============------------------
```
![image](https://user-images.githubusercontent.com/10180619/177040107-ba9fc8ca-4571-42a1-a81c-9f52d4a48cb7.png)

## stickleback v0.1
### [SPINE](https://github.com/QVEU/SPINE_Q) Library QC by nanopore using Levenshtein distance to map insertion sites on target sequence.

`stickleback` maps nanopore reads containing insertions (such as molecular handles) and then identifies the insertion site on the template molecule. 

```
python stickleback.0.1.py path/to/samfile.sam queryString path/to/template.fasta [minimumReadLength] [maximumReadLength]
```

