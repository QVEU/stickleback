#!/usr/bin/env python
'''
simulate_reads.py
  ><```º>

Builds a synthetic paired-end Illumina dataset with known insertion sites, for benchmarking
and for checking that stickleback_illumina.py (and the older SAM-based stickleback.0.3.py)
call the sites we think they call.

Every read pair is drawn from one of a set of named categories, so the truth table says not
just where the insert is but what kind of read it is (clean, error in the query, too little
flank, inside a template repeat, empty, ...). compare_sites.py scores a caller against it.

Outputs, under --out:
    template.fasta          the template (random sequence + one tandem duplication)
    R1.fastq.gz, R2.fastq.gz
    truth.csv               read_id,category,insPos_v,orientation,insPos_v2,orientation2
                            (insPos_v2 is only set for chimeric pairs, whose two mates carry
                            different sites and which a caller should report as discordant)
    reads.sam               (with --sam) the same reads as a minimal SAM, for stickleback.0.3.py
    recipe.txt              the exact parameters this dataset was built with

insPos_v follows the convention of both scripts: the 1-based template position of the first
template base after the insert, and orientation is "+" when the insert is in the template's
orientation.

USAGE:
    python benchmark/simulate_reads.py --out bench/data --n-pairs 500000
    python benchmark/simulate_reads.py --out bench/clean --n-pairs 20000 --only clean --sam
'''

import argparse
import gzip
import os
import random
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
B = "ACGT"

# Fraction of pairs drawn from each category. Everything but "clean" and "no_insert" is an
# edge case we want represented but not dominant.
DEFAULT_WEIGHTS = {
    "clean":        0.55,   # insert well inside a read, ample unique flank on both sides
    "no_insert":    0.20,   # plain template, no query at all
    "query_error":  0.08,   # 1-2 substitutions inside the query (exact matching must miss these)
    "low_flank":    0.07,   # query present but < flank template bases beside it in the read
    "repeat_flank": 0.04,   # insert inside the duplicated segment, so flank k-mers are not unique
    "terminal":     0.03,   # insert at the very start or end of the template
    "n_bases":      0.015,  # a run of Ns in the flank
    "chimeric":     0.01,   # the two mates carry different sites
    "empty_read":   0.005,  # zero-length read (adapter dimer trimmed to nothing)
}


def revcomp(s):
    return s.translate(COMP)[::-1]


def makeTemplate(rng, length, dupStart, dupLen):
    '''Random template carrying one tandem duplication, so some k-mers are genuinely ambiguous.'''
    seq = [rng.choice(B) for _ in range(length)]
    seq = "".join(seq)
    dup = seq[dupStart:dupStart + dupLen]
    return seq[:dupStart + dupLen] + dup + seq[dupStart + dupLen:]


def mutate(rng, s, n):
    '''n substitutions at distinct positions.'''
    if not s or n <= 0:
        return s
    out = list(s)
    for i in rng.sample(range(len(s)), min(n, len(s))):
        out[i] = rng.choice([b for b in B if b != out[i]])
    return "".join(out)


def addErrors(rng, s, err):
    if err <= 0:
        return s
    return "".join(rng.choice(B) if rng.random() < err else c for c in s)


def qualString(rng, n):
    '''Plausible Illumina quality line: mostly high, decaying towards the 3' end.'''
    out = []
    for i in range(n):
        q = 38 - int(10 * (i / max(n - 1, 1)) ** 3) - rng.randint(0, 3)
        out.append(chr(33 + max(2, min(41, q))))
    return "".join(out)


def pickSite(rng, args, category, template):
    '''Returns (junction, orientation). junction is the 0-based template index of the first
    base after the insert, so insPos_v == junction + 1.'''
    L = len(template)
    orient = rng.choice("+-")
    if category == "terminal":
        return (1 if rng.random() < 0.5 else L - 1), orient
    if category == "repeat_flank":
        # inside the second copy of the duplication, where flank k-mers are not unique
        lo = args.dup_start + args.dup_len + args.flank
        hi = args.dup_start + 2 * args.dup_len - args.flank
        return rng.randint(lo, hi), orient
    return rng.randint(args.flank + 5, L - args.flank - 5), orient


def buildFragment(rng, args, template, category, query):
    '''Returns (fragment, truth) where truth is (insPos_v, orientation) or None.'''
    L = len(template)
    fragLen = rng.randint(args.frag_min, args.frag_max)
    if category == "no_insert":
        start = rng.randint(0, max(0, L - fragLen))
        return template[start:start + fragLen], None

    j, orient = pickSite(rng, args, category, template)
    insert = query if orient == "+" else revcomp(query)
    if category == "query_error":
        insert = mutate(rng, insert, rng.randint(1, 2))

    # How much template sits upstream of the insert inside the read that will cover it.
    if category == "low_flank":
        up = rng.randint(1, args.flank - 1)
        if rng.random() < 0.5:                       # ... or too little on the downstream side
            up = args.read_len - len(insert) - rng.randint(1, args.flank - 1)
    else:
        up = rng.randint(args.flank + 2, args.read_len - len(insert) - args.flank - 2)
    up = max(0, min(up, args.read_len - len(insert)))

    if rng.random() < 0.5:                           # put the junction in the R2 end instead
        up = fragLen - len(insert) - up
        up = max(0, min(up, fragLen - len(insert)))

    left = template[max(0, j - up):j]
    right = template[j:j + (fragLen - len(insert) - len(left))]
    frag = left + insert + right
    if category == "n_bases" and len(frag) > 20:
        p = rng.randint(0, len(frag) - 10)
        frag = frag[:p] + "N" * 8 + frag[p + 8:]
    return frag, (j + 1, orient)


def simulate(args):
    rng = random.Random(args.seed)
    template = makeTemplate(rng, args.template_len, args.dup_start, args.dup_len)
    query = args.query or "".join(rng.choice(B) for _ in range(args.query_len))
    while query in template or revcomp(query) in template:   # a query that occurs in the template
        query = "".join(rng.choice(B) for _ in range(args.query_len))

    weights = dict(DEFAULT_WEIGHTS)
    if args.only:
        weights = {args.only: 1.0}
    if args.weights:
        weights = {}
        for item in args.weights.split(","):
            name, _, frac = item.partition("=")
            weights[name.strip()] = float(frac)
    for name in weights:
        if name not in DEFAULT_WEIGHTS:
            sys.exit("ERROR: unknown category: {} (known: {})".format(name, ", ".join(DEFAULT_WEIGHTS)))
    names = list(weights)
    cum = list(weights.values())

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "template.fasta"), "w") as fh:
        fh.write(">benchmark_template\n")
        for i in range(0, len(template), 60):
            fh.write(template[i:i + 60] + "\n")

    f1 = gzip.open(os.path.join(args.out, "R1.fastq.gz"), "wt", compresslevel=args.gzip_level)
    f2 = gzip.open(os.path.join(args.out, "R2.fastq.gz"), "wt", compresslevel=args.gzip_level)
    truth = open(os.path.join(args.out, "truth.csv"), "w")
    truth.write("read_id,category,insPos_v,orientation,insPos_v2,orientation2\n")
    sam = None
    if args.sam:
        sam = open(os.path.join(args.out, "reads.sam"), "w")
        sam.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:benchmark_template\tLN:{}\n".format(len(template)))

    for i in range(args.n_pairs):
        category = rng.choices(names, weights=cum)[0]
        frag, site = buildFragment(rng, args, template, category, query)
        site2 = None
        if category == "chimeric":                   # second mate from an unrelated site
            frag2, site2 = buildFragment(rng, args, template, "clean", query)
            r1 = addErrors(rng, frag[:args.read_len], args.error_rate)
            r2 = addErrors(rng, revcomp(frag2)[:args.read_len], args.error_rate)
        else:
            r1 = addErrors(rng, frag[:args.read_len], args.error_rate)
            r2 = addErrors(rng, revcomp(frag)[:args.read_len], args.error_rate)
        if category == "empty_read":
            r1 = ""

        name = "sim{:08d}".format(i)
        f1.write("@{} 1:N:0:1\n{}\n+\n{}\n".format(name, r1, qualString(rng, len(r1))))
        f2.write("@{} 2:N:0:1\n{}\n+\n{}\n".format(name, r2, qualString(rng, len(r2))))
        truth.write("{},{},{},{},{},{}\n".format(
            name, category,
            site[0] if site else "", site[1] if site else "",
            site2[0] if site2 else "", site2[1] if site2 else ""))
        if sam:
            # Minimal SAM for stickleback.0.3.py: it only reads QNAME, RNAME and SEQ. Reads are
            # written in template orientation, which is what an aligner would emit.
            for suffix, flag, seq in (("/1", 65, r1), ("/2", 129, revcomp(r2))):
                if not seq:
                    continue
                sam.write("{}{}\t{}\tbenchmark_template\t1\t60\t{}M\t=\t1\t0\t{}\t{}\n".format(
                    name, suffix, flag, len(seq), seq, "I" * len(seq)))

    f1.close(); f2.close(); truth.close()
    if sam:
        sam.close()

    with open(os.path.join(args.out, "recipe.txt"), "w") as fh:
        fh.write("query\t{}\n".format(query))
        fh.write("template_len\t{}\n".format(len(template)))
        fh.write("duplication\t{}..{} ({} bp, tandem)\n".format(
            args.dup_start, args.dup_start + args.dup_len, args.dup_len))
        for k, v in sorted(vars(args).items()):
            fh.write("arg_{}\t{}\n".format(k, v))
        fh.write("categories\t{}\n".format(weights))
    print("query: {}".format(query))
    print("template: {} bp (tandem duplication of {} bp at {})".format(
        len(template), args.dup_len, args.dup_start))
    print("wrote {} pairs to {}".format(args.n_pairs, args.out))


def parseArgs(argv):
    p = argparse.ArgumentParser(description="Simulate paired-end reads with known insertion sites.")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--n-pairs", type=int, default=500_000)
    p.add_argument("--seed", type=int, default=20261016)
    p.add_argument("--template-len", type=int, default=4700, help="before the duplication is added")
    p.add_argument("--dup-start", type=int, default=2000)
    p.add_argument("--dup-len", type=int, default=300)
    p.add_argument("--query", help="query sequence (default: random --query-len mer)")
    p.add_argument("--query-len", type=int, default=34)
    p.add_argument("--read-len", type=int, default=150)
    p.add_argument("--frag-min", type=int, default=250)
    p.add_argument("--frag-max", type=int, default=450)
    p.add_argument("--error-rate", type=float, default=0.002, help="substitution rate per base")
    p.add_argument("--flank", type=int, default=16, help="flank the caller will use, for sizing reads")
    p.add_argument("--gzip-level", type=int, default=6)
    p.add_argument("--only", help="draw every pair from this one category")
    p.add_argument("--weights", help="e.g. clean=0.8,no_insert=0.2")
    p.add_argument("--sam", action="store_true", help="also write reads.sam for stickleback.0.3.py")
    return p.parse_args(argv)


if __name__ == "__main__":
    simulate(parseArgs(sys.argv[1:]))
