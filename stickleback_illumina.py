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

Reads are parsed in blocks and mapped in worker processes (--procs), so a run scales with
the cores it is given. With --max-dist > 0 a seed prefilter (pigeonhole: a query of length
Q matching with at most d edits must contain one of its d+1 non-overlapping Q/(d+1)-mers
exactly) skips edlib for the reads that cannot match, at no cost in sensitivity.

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
import collections
import gzip
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from multiprocessing import Pool

# Reads are handed to the workers in blocks of this many records.
BLOCK_RECORDS = 20_000

STATUSES = ("placed", "unplaced", "discordant", "noQuery")

class InputError(Exception):
    '''A malformed or inconsistent input: reported to the user, not as a traceback.'''

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

def fastqBlocks(handle, path, nRecords=None):
    '''Yields (headers, seqs) lists of up to nRecords records each. The blocks are the unit of
    work handed to the mapping processes; reading them with readline() keeps the decompressor
    busy, since it can go on filling the pipe while the block is mapped.'''
    nRecords = nRecords or BLOCK_RECORDS
    readline = handle.readline
    headers, seqs = [], []
    while True:
        header = readline()
        if not header:
            break
        seq = readline()
        readline()
        if not readline():
            if not header.strip() and not seq.strip():
                break                           # stray blank line(s) at the end of the file
            raise InputError("ERROR: truncated FASTQ record at the end of {}".format(path))
        headers.append(header)
        seqs.append(seq.rstrip())
        if len(seqs) == nRecords:
            yield headers, seqs
            headers, seqs = [], []
    if seqs:
        yield headers, seqs

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
        self.seeds = ()
        if maxDist > 0:
            import edlib  # only needed for fuzzy matching
            self.edlib = edlib
            # Pigeonhole prefilter: maxDist edits can touch at most maxDist of the
            # maxDist+1 non-overlapping seeds below, so a read that matches the query
            # within maxDist edits contains at least one of them exactly. A read with
            # none of them cannot match, and is skipped without calling edlib.
            seedLen = self.qL // (maxDist + 1)
            if seedLen > 1:
                seeds = [query[i * seedLen:(i + 1) * seedLen] for i in range(maxDist + 1)]
                self.seeds = tuple(seeds + [revcomp(s) for s in seeds])

    def findQuery(self, read):
        '''Returns (read oriented so the query is forward, start, end, editDistance) or None.'''
        qL = self.qL
        p = read.find(self.q)
        if p >= 0:
            return read, p, p + qL, 0
        p = read.find(self.qrc)
        if p >= 0:
            read = revcomp(read)
            p = len(read) - p - qL
            return read, p, p + qL, 0
        maxDist = self.maxDist
        if maxDist == 0:
            return None
        if self.seeds and not any(seed in read for seed in self.seeds):
            return None
        best = None
        for seq in (read, revcomp(read)):
            hit = self.edlib.align(self.q, seq, mode="HW", task="locations", k=maxDist)
            d = hit["editDistance"]
            if d < 0 or d > maxDist:   # no hit; an empty read reports d = len(query) > maxDist
                continue
            if best is not None and d >= best[3]:
                continue
            locations = hit["locations"]
            if not locations or locations[0][0] is None or locations[0][1] is None:
                continue
            s, e = locations[0]
            best = (seq, s, e + 1, d)
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
        junction = self.junction
        left = junction(seq[s - k:s], "left") if s >= k else None
        right = junction(seq[e:e + k], "right") if len(seq) - e >= k else None
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

##### Blocks of reads #####
# One Mapper per worker process, built once by the pool initializer: the k-mer index is
# the same for every block and too big to ship with each one.
WORKER = {}
R1KEY = {status: "R1_" + status for status in STATUSES}
R2KEY = {status: "R2_" + status for status in STATUSES}

def initWorker(query, templateSeq, k, circular, maxDist):
    WORKER["mapper"] = Mapper(query, templateSeq, k, circular, maxDist)

def checkPairing(first, headers1, headers2):
    '''Raises unless the two mate files hold the same reads in the same order. The common case
    is checked in one pass through zip(); the offending read is only located if one is found.'''
    if all(readID(h1) == readID(h2) for h1, h2 in zip(headers1, headers2)):
        return
    for i, (h1, h2) in enumerate(zip(headers1, headers2)):
        if readID(h1) != readID(h2):
            raise InputError("ERROR: R1/R2 out of sync at read {}: {} vs {}".format(
                first + i + 1, readID(h1), readID(h2)))

def mapBlock(block):
    '''Maps one block of reads, or of read pairs, to status and site counts. Runs in a worker
    process, so it returns Counters to merge rather than touching shared state.'''
    first, headers1, seqs1, headers2, seqs2 = block
    mapRead = WORKER["mapper"].mapRead
    stats = Counter()
    sites = Counter()
    if seqs2 is None:
        for seq in seqs1:
            status, site, d = mapRead(seq)
            stats[status] += 1
            if site:
                sites[site + (d,)] += 1
        return len(seqs1), stats, sites
    checkPairing(first, headers1, headers2)
    for s1, s2 in zip(seqs1, seqs2):
        m1 = mapRead(s1)
        m2 = mapRead(s2)
        stats[R1KEY[m1[0]]] += 1
        stats[R2KEY[m2[0]]] += 1
        status, site, d = combinePair(m1, m2)
        stats[status] += 1
        if site:
            sites[site + (d,)] += 1
    return len(seqs1), stats, sites

def readBlocks(r1Path, r2Path, threads):
    '''Yields blocks of reads (single-end) or of mate-paired reads, in file order.'''
    blocks1 = fastqBlocks(openFastq(r1Path, threads), r1Path)
    blocks2 = fastqBlocks(openFastq(r2Path, threads), r2Path) if r2Path else None
    first = 0
    for headers1, seqs1 in blocks1:
        if blocks2 is None:
            yield first, None, seqs1, None, None     # headers are only used for pairing
        else:
            headers2, seqs2 = next(blocks2, (None, None))
            if seqs2 is None or len(seqs2) != len(seqs1):
                # Mates that have gone out of step make for a more useful error than a count,
                # and are the usual reason for one file running out first.
                checkPairing(first, headers1, headers2 or [])
                raise InputError("ERROR: R2 has {} reads than R1.".format(
                    "fewer" if seqs2 is None or len(seqs2) < len(seqs1) else "more"))
            yield first, headers1, seqs1, headers2, seqs2
        first += len(seqs1)
    if blocks2 is not None and next(blocks2, None) is not None:
        raise InputError("ERROR: R2 has more reads than R1.")

def mapInPool(pool, blocks, inFlight):
    '''Maps blocks in the pool, keeping at most inFlight of them outstanding. Pool.imap feeds
    its whole input to the workers as fast as the queue takes it, which for a large FASTQ
    means reading the file into memory faster than it can be mapped.'''
    pending = collections.deque()
    for block in blocks:
        pending.append(pool.apply_async(mapBlock, (block,)))
        while len(pending) >= inFlight:
            yield pending.popleft().get()
    while pending:
        yield pending.popleft().get()

def workerCount(requested):
    '''Worker processes to map with: one per core, less the one the parent needs to read and
    decompress, unless the user asked for a specific number. More than a handful does not
    help, since the parent has to parse the FASTQ for all of them.'''
    if requested > 0:
        return requested
    cores = os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1
    return max(1, min(8, int(cores) - 1))

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
    p.add_argument("--procs", type=int, default=0,
                   help="worker processes used to map reads (default: one per core less one, "
                        "from SLURM_CPUS_PER_TASK when set, at most 8; 1 maps in this process)")
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

    procs = workerCount(args.procs)
    print("Mapping processes: {}".format(procs))

    stats = Counter()
    sites = Counter()
    n = 0
    reported = 0
    blocks = readBlocks(args.r1, args.r2 if paired else None, args.threads)
    pool = None
    try:
        if procs > 1:
            pool = Pool(procs, initializer=initWorker,
                        initargs=(query, templateSeq, args.flank, args.circular, args.max_dist))
            results = mapInPool(pool, blocks, procs + 2)
        else:
            WORKER["mapper"] = mapper
            results = map(mapBlock, blocks)
        for count, blockStats, blockSites in results:
            n += count
            stats += blockStats
            sites += blockSites
            if args.progress and n - reported >= args.progress:
                reported = n
                print("  {:,} {} processed, {:,} placed ({:.1f} min)".format(
                    n, "pairs" if paired else "reads", stats["placed"],
                    (time.time() - t0) / 60), flush=True)
    except InputError as e:
        sys.exit(str(e))
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()

    sitesFile = args.out + "_sites.csv"
    with open(sitesFile, "w") as out:
        out.write("insPos_v,orientation,minD,count\n")
        for (pos, strand, d), c in sorted(sites.items()):
            out.write("{},{},{},{}\n".format(pos, strand, d, c))

    unit = "pairs" if paired else "reads"
    lines = ["mode\t{}".format("paired-end" if paired else "single-end"),
             "total_{}\t{}".format(unit, n)]
    for status in STATUSES:
        lines.append("{}_{}\t{}\t{:.2%}".format(status, unit, stats[status], stats[status] / max(n, 1)))
    if paired:
        for mate in ("R1", "R2"):
            for status in STATUSES:
                lines.append("{}_{}\t{}".format(mate, status, stats[mate + "_" + status]))
    lines.append("distinct_sites\t{}".format(len({s[:2] for s in sites})))
    lines.append("minutes\t{:.2f}".format((time.time() - t0) / 60))
    with open(args.out + "_summary.txt", "w") as out:
        out.write("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print("\nWrote {} and {}_summary.txt".format(sitesFile, args.out))

if __name__ == "__main__":
    main()
