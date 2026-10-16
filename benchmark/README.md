# Benchmark: synthetic reads with known insertion sites

`simulate_reads.py` builds a paired-end Illumina dataset in which every read pair's
insertion site is known, and `compare_sites.py` scores a caller's output against that
truth table. The point is twofold: a reproducible dataset for timing runs, and a way to
answer "is this site real?" when two versions of the pipeline disagree.

## Recipe used for the October 2026 comparison

```bash
# 50k pairs, all read categories, plus a SAM so the old nanopore script can be run on it
python benchmark/simulate_reads.py --out bench/mixed --n-pairs 50000 --seed 11 --sam

Q=$(grep -m1 -P '^query\t' bench/mixed/recipe.txt | cut -f2)

# current script
python stickleback_illumina.py -1 bench/mixed/R1.fastq.gz -2 bench/mixed/R2.fastq.gz \
    -q "$Q" -t bench/mixed/template.fasta -o bench/mixed/new

# previous version (SAM-based, nanopore)
python stickleback.0.3.py bench/mixed/reads.sam "$Q" bench/mixed/template.fasta

python benchmark/compare_sites.py --truth bench/mixed/truth.csv \
    --sites bench/mixed/new_sites.csv --old bench/mixed/reads_stickleback.csv
```

Defaults: 5,000 bp template (4,700 bp of random sequence with a 300 bp tandem duplication
at 2,000, so some flanks are genuinely ambiguous), a random 34 bp query that does not occur
in the template, 150 bp reads, 250-450 bp fragments, 0.2% substitution error, realistic
quality strings, gzip level 6. `--seed` makes a run reproducible; `recipe.txt` in the output
directory records every parameter and the query that was drawn.

For a pure timing dataset, 500k pairs at these defaults is about 45 MB of gzipped FASTQ.
Use `--weights` to drop the edge cases (e.g. `--weights clean=0.7,no_insert=0.3`).

## Read categories

Each pair is drawn from one category, recorded in `truth.csv`, so a caller can be scored
on each kind of read separately (`--only <category>` makes a dataset of just one).

| category | what it is | expected call |
|---|---|---|
| `clean` | insert well inside a read, ample unique flank | placed at the true site |
| `no_insert` | plain template | noQuery |
| `query_error` | 1-2 substitutions inside the query | noQuery in exact mode, placed with `--max-dist` |
| `low_flank` | fewer than `--flank` template bases on one side of the query | placed from the other side, else unplaced |
| `repeat_flank` | insert inside the tandem duplication | unplaced (flank k-mers are not unique) |
| `terminal` | insert at the first or last template position | placed at 1 or L |
| `n_bases` | a run of 8 Ns in the flank | placed, or noQuery if the Ns hit the query |
| `chimeric` | the two mates carry different sites | discordant |
| `empty_read` | zero-length R1 | not placed, and must not crash |

`truth.csv` columns: `read_id,category,insPos_v,orientation,insPos_v2,orientation2`.
`insPos_v2` is only set for chimeric pairs. `insPos_v` is the 1-based template position of
the first template base after the insert, the same convention both scripts use.

## Scoring

`compare_sites.py` reports, against the truth table:

- for `stickleback_illumina.py`, a site-level comparison (its `_sites.csv` is aggregated, so
  there are no read names to join on): how many placed pairs sit on a true
  (position, orientation), how many sites were called that do not exist, how many true sites
  were never called, and how far the spurious calls sit from the nearest true position;
- for `stickleback.0.3.py`, a per-read comparison (its output carries read names): the
  histogram of `insPos_v` offsets from truth, broken down by category.

`--categories clean,terminal` restricts the truth to those categories.

## Notes

- The SAM written by `--sam` is minimal: QNAME, RNAME and SEQ in template orientation, which
  is all `stickleback.0.3.py` reads. Empty reads are omitted, as an aligner would omit them.
- `stickleback.0.3.py` takes about 2 minutes for 50k pairs on 4 cores, so keep the datasets
  you run it on small; `stickleback_illumina.py` does the same data in under a second.
