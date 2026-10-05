```
----------------=============----------------
--==--==--==--==   ><```º>   ==--==--==--==--
==--==--==--==-- stickleback --==--==--==--==
----------------=============----------------
```
![image](https://user-images.githubusercontent.com/10180619/177040107-ba9fc8ca-4571-42a1-a81c-9f52d4a48cb7.png)

## stickleback v0.1
### [SPINE](https://github.com/QVEU/SPINE_Q) Library QC by nanopore using Levenshtein distance to map insertion sites on target sequence.

`stickleback` maps nanopore reads containing insertions (such as molecular handles) and then identifies the insertion site on the template molecule. Requires a merged sam file generated from mapping (e.g. minimap2). 

```
python stickleback.0.3.py path/to/samfile.sam queryString path/to/template.fasta [minimumReadLength] [maximumReadLength]
```

## Illumina (paired-end) data: `stickleback_illumina.py`
For high-accuracy, high-depth short reads. Reads FASTQs directly (no minimap2/SAM needed), finds the query by exact match on either strand, and places the junction by looking up the flanking 16 bp in a unique k-mer index of the template. Each read pair is counted once. Roughly 4 µs per pair (~7 min per 100M pairs on one core, plus decompression).

```
python stickleback_illumina.py -1 R1.fastq.gz -2 R2.fastq.gz -q queryString -t template.fasta -o out/sample [--circular] [--max-dist 2] [--flank 16]
```
Leave out `-2` for single-end data. On Slurm (create `OUT/` first for the log):
```
sbatch stickleback_illumina.sbatch R1.fastq.gz R2.fastq.gz queryString template.fasta out/sample --circular
```

Outputs:
- `out/sample_sites.csv`: `insPos_v,orientation,minD,count`. `insPos_v` uses the v0.3 convention (1-based position of the first template base after the insert). `orientation` is `+` if the insert is in the template's orientation, `-` if reverse-complemented. `count` is the number of read pairs.
- `out/sample_summary.txt`: pairs placed / unplaced (flank not unique or not found) / discordant (the two flanks or the two mates disagree) / no query, plus per-mate counts.

Options:
- `--circular`: plasmid templates, so junctions across the origin are placed.
- `--max-dist N`: allow up to N edits in the query (needs `pip install edlib`). This recovers reads with a sequencing error in the query, at ~2.5× the run time. Reads matched this way have `minD > 0`; a mismatch at the query's edge can shift the call by 1 bp.
- `--flank K`: flank length used for placement (default 16).

In R, weight the histograms by `count`, e.g. `geom_histogram(aes(insPos_v, weight = count))`.

Test: `python tests/test_stickleback_illumina.py`

