#!/usr/bin/env python
'''
stickleback_illumina.py
v1.0
  ><```º>
Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit

Maps insertion sites of a query sequence (e.g. a molecular handle) in a template from
Illumina reads. Built for high-accuracy, high-depth short-read data, where the nanopore
scripts (stickleback.0.x.py) are far too slow.

How it works (no alignment / SAM needed):
    1. Index every k-mer (k = --flank) of the template on both strands, keeping only k-mers
       that occur once.
    2. Stream the FASTQ(s). In each read, find the query (exact match, either strand; optional
       fuzzy fallback with edlib).
    3. Look up the k bases on either side of the query in the index. The hit gives the
       template position of the junction and the orientation of the insert.
    4. In paired-end mode both mates are examined and each pair is counted once.

Output:
    <prefix>_sites.csv    one row per (insPos_v, orientation, minD) with the number of pairs
    <prefix>_summary.txt  read/pair accounting

insPos_v uses the same convention as stickleback.0.3.py: the 1-based template position of
the first template base after the insert (the insert sits between insPos_v-1 and insPos_v).
orientation is "+" if the query is inserted in the same orientation as the template, "-" if
reverse-complemented.

USAGE:
    paired-end:  python stickleback_illumina.py -1 R1.fastq.gz -2 R2.fastq.gz -q QUERY -t template.fasta -o out/sample
    single-end:  python stickleback_illumina.py -1 R1.fastq.gz -q QUERY -t template.fasta -o out/sample
'''

##### Imports #####
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
from collections import Counter

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(s):
    return s.translate(COMP)[::-1]

##### Input #####
def readFasta(path):
    with open(path) as fasta:
        return "".join(l.strip() for l in fasta if not l.startswith(">")).upper()

def openFastq(path, threads=2):
    '''Open plain or gzipped FASTQ as a text stream. Uses pigz for decompression if available
    (much faster than Python's gzip, which otherwise becomes the bottleneck).'''
    if not path.endswith(".gz"):
        return open(path)
    pigz = shutil.which("pigz")
    if pigz:
        proc = subprocess.Popen([pigz, "-dc", "-p", str(threads), path], stdout=subprocess.PIPE,
                                text=True, bufsize=1 << 20)
        return proc.stdout
    return gzip.open(path, "rt")

def fastqRecords(handle):
    '''Yields (name, seq) from a FASTQ stream.'''
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline().rstrip()
        handle.readline()
        handle.readline()
        yield header, seq

def readID(header):
    '''Read name without '@', comment, or /1 /2 suffix, for pairing checks.'''
    name = header[1:].split(None, 1)[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name

##### Mapping #####
def buildIndex(templateSeq, k, circular):
    '''
    Maps each k-mer that occurs exactly once across both strands of the template to
    (position, strand). Position is the 0-based start of the k-mer on that strand.
    '''
    L = len(templateSeq)
    counts = Counter()
    entries = {}
    for strand, seq in (("+", templateSeq), ("-", revcomp(templateSeq))):
        if circular:
            seq = seq + seq[:k - 1]
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            counts[kmer] += 1
            entries[kmer] = (i, strand)
    index = {kmer: entries[kmer] for kmer, n in counts.items() if n == 1}
    return index, len(entries) - len(index)

class Mapper:
    def __init__(self, query, templateSeq, k, circular, maxDist):
        self.q = query
        self.qrc = revcomp(query)
        self.qL = len(query)
        self.L = len(templateSeq)
        self.k = k
        self.circular = circular
        self.maxDist = maxDist
        self.index, self.nAmbiguous = buildIndex(templateSeq, k, circular)
        if maxDist > 0:
            import edlib  # only needed for fuzzy matching
            self.edlib = edlib

    def findQuery(self, read):
        '''Returns (read oriented so the query is forward, start, end, editDistance) or None.'''
        p = read.find(self.q)
        if p >= 0:
            return read, p, p + self.qL, 0
        p = read.find(self.qrc)
        if p >= 0:
            read = revcomp(read)
            p = len(read) - p - self.qL
            return read, p, p + self.qL, 0
        if self.maxDist == 0:
            return None
        best = None
        for seq in (read, revcomp(read)):
            hit = self.edlib.align(self.q, seq, mode="HW", task="locations", k=self.maxDist)
            if hit["editDistance"] >= 0 and (best is None or hit["editDistance"] < best[3]):
                s, e = hit["locations"][0]
                best = (seq, s, e + 1, hit["editDistance"])
        return best

    def junction(self, kmer, side):
        '''Template junction (0-based index of the first template base after the insert) and
        insert orientation from a flank k-mer. side is "left" or "right" of the query.'''
        hit = self.index.get(kmer)
        if hit is None:
            return None
        pos, strand = hit
        j = pos + self.k if side == "left" else pos  # junction on the strand the read matched
        if strand == "-":
            j = self.L - j  # rc-strand coordinate -> template coordinate
        if self.circular:
            j %= self.L
        return j, strand

    def mapRead(self, read):
        '''
        Returns (status, site, minD). status is one of "placed", "noQuery", "unplaced",
        "discordant". site is (insPos_v, orientation) when placed.
        '''
        found = self.findQuery(read)
        if found is None:
            return "noQuery", None, None
        seq, s, e, d = found
        k = self.k
        left = self.junction(seq[s - k:s], "left") if s >= k else None
        right = self.junction(seq[e:e + k], "right") if len(seq) - e >= k else None
        if left and right and left != right:
            return "discordant", None, d
        hit = left or right
        if hit is None:
            return "unplaced", None, d
        return "placed", (hit[0] + 1, hit[1]), d

def combinePair(m1, m2):
    '''Combines mate calls into one call per pair.'''
    s1, site1, d1 = m1
    s2, site2, d2 = m2
    if "discordant" in (s1, s2) or (site1 and site2 and site1 != site2):
        return "discordant", None, None
    if site1 and site2:
        return "placed", site1, min(d1, d2)
    if site1:
        return m1
    if site2:
        return m2
    if "unplaced" in (s1, s2):
        return "unplaced", None, None
    return "noQuery", None, None

##### Main #####
def parseArgs(argv):
    p = argparse.ArgumentParser(description="Map query insertion sites in a template from Illumina reads.")
    p.add_argument("-1", "--r1", required=True, help="R1 FASTQ (.gz ok)")
    p.add_argument("-2", "--r2", help="R2 FASTQ (.gz ok); enables paired-end mode")
    p.add_argument("-q", "--query", required=True, help="inserted sequence to locate")
    p.add_argument("-t", "--template", required=True, help="template FASTA")
    p.add_argument("-o", "--out", required=True, help="output prefix")
    p.add_argument("--flank", type=int, default=16,
                   help="flank length used to place the junction (default 16)")
    p.add_argument("--max-dist", type=int, default=0,
                   help="max edit distance for query matching; >0 enables edlib fallback (default 0 = exact)")
    p.add_argument("--circular", action="store_true", help="template is circular (e.g. plasmid)")
    p.add_argument("--threads", type=int, default=2, help="pigz threads per input file (default 2)")
    p.add_argument("--progress", type=int, default=5_000_000, help="report every N reads/pairs")
    return p.parse_args(argv)

def main(argv=None):
    args = parseArgs(argv)
    t0 = time.time()
    print("\n----------------=============-----------------")
    print("--==--==--==--==   ><```º>   ==--==--==--==--=")
    print("==--==-- stickleback (illumina) --==--==--==-")
    print("----------------=============-----------------\n")

    # Fail fast, before the reads are streamed: a missing input or an unwritable output
    # directory should not surface only at the final write, hours into a run. Inputs are
    # checked first so a run that cannot start leaves no empty output directory behind.
    for path in (args.r1, args.r2, args.template):
        if path and not os.path.isfile(path):
            sys.exit("ERROR: input file not found: {}".format(path))
    outDir = os.path.dirname(os.path.abspath(args.out))
    try:
        os.makedirs(outDir, exist_ok=True)
    except OSError as e:
        sys.exit("ERROR: cannot create output directory {}: {}".format(outDir, e))
    if not os.access(outDir, os.W_OK):
        sys.exit("ERROR: output directory is not writable: {}".format(outDir))

    query = args.query.upper()
    templateSeq = readFasta(args.template)
    mapper = Mapper(query, templateSeq, args.flank, args.circular, args.max_dist)
    paired = args.r2 is not None
    print("Query length: {}".format(len(query)))
    print("Template length: {} ({})".format(len(templateSeq), "circular" if args.circular else "linear"))
    print("Flank k-mers: {} unique, {} non-unique (reads flanked only by these are unplaced)".format(
        len(mapper.index), mapper.nAmbiguous))
    if query in templateSeq or revcomp(query) in templateSeq:
        print("WARNING: query occurs in the template itself; reads from uninserted template will be counted.")
    print("Mode: {}".format("paired-end" if paired else "single-end"))

    stats = Counter()
    sites = Counter()
    r1 = fastqRecords(openFastq(args.r1, args.threads))
    if paired:
        r2 = fastqRecords(openFastq(args.r2, args.threads))
    n = 0
    for h1, s1 in r1:
        call = mapper.mapRead(s1)
        if paired:
            try:
                h2, s2 = next(r2)
            except StopIteration:
                sys.exit("ERROR: R2 has fewer reads than R1.")
            if readID(h1) != readID(h2):
                sys.exit("ERROR: R1/R2 out of sync at read {}: {} vs {}".format(n + 1, readID(h1), readID(h2)))
            m2 = mapper.mapRead(s2)
            stats["R1_" + call[0]] += 1
            stats["R2_" + m2[0]] += 1
            call = combinePair(call, m2)
        status, site, d = call
        stats[status] += 1
        if site:
            sites[site + (d,)] += 1
        n += 1
        if n % args.progress == 0:
            print("  {:,} {} processed, {:,} placed ({:.1f} min)".format(
                n, "pairs" if paired else "reads", stats["placed"], (time.time() - t0) / 60), flush=True)
    if paired and next(r2, None) is not None:
        sys.exit("ERROR: R2 has more reads than R1.")

    sitesFile = args.out + "_sites.csv"
    with open(sitesFile, "w") as out:
        out.write("insPos_v,orientation,minD,count\n")
        for (pos, strand, d), c in sorted(sites.items()):
            out.write("{},{},{},{}\n".format(pos, strand, d, c))

    unit = "pairs" if paired else "reads"
    lines = ["mode\t{}".format("paired-end" if paired else "single-end"),
             "total_{}\t{}".format(unit, n)]
    for status in ("placed", "unplaced", "discordant", "noQuery"):
        lines.append("{}_{}\t{}\t{:.2%}".format(status, unit, stats[status], stats[status] / max(n, 1)))
    if paired:
        for mate in ("R1", "R2"):
            for status in ("placed", "unplaced", "discordant", "noQuery"):
                lines.append("{}_{}\t{}".format(mate, status, stats[mate + "_" + status]))
    lines.append("distinct_sites\t{}".format(len({s[:2] for s in sites})))
    lines.append("minutes\t{:.2f}".format((time.time() - t0) / 60))
    with open(args.out + "_summary.txt", "w") as out:
        out.write("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print("\nWrote {} and {}_summary.txt".format(sitesFile, args.out))

if __name__ == "__main__":
    main()
