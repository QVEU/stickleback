'''
Simulation test for stickleback_illumina.py: builds paired-end reads from fragments carrying the
query at known sites (both orientations, linear and circular templates, 0.3% substitution errors)
and checks every placed pair lands on its true site.
Run: python tests/test_stickleback_illumina.py   (or pytest)
'''
import gzip, os, random, sys, tempfile
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import stickleback_illumina as sb

QUERY = "AGCGGGAGACCGGGGTCTCTGAGCG"
B = "ACGT"

def simulate(rng, template, n, circular, readLen=150, err=0.003):
    L = len(template)
    ext = template * 2 if circular else template
    def noisy(s):
        return "".join(rng.choice(B) if rng.random() < err else c for c in s)
    pairs = []
    for i in range(n):
        fragLen = rng.randint(200, 450)
        j = rng.randint(1, L - 1)                     # junction: insert precedes template[j]
        orient = rng.choice("+-")
        if rng.random() < 0.1:                        # fragment without an insert
            start = rng.randint(0, (L if circular else L - fragLen) - 1)
            frag, truth = ext[start:start + fragLen], None
        else:
            # keep the insert inside the sequenced part of the fragment (first or last readLen bases)
            up = rng.randint(5, readLen - len(QUERY) - 5)
            if rng.random() < 0.5:
                up = fragLen - len(QUERY) - up
            down = fragLen - len(QUERY) - up
            if circular:
                left, right = ext[(j - up) % L:(j - up) % L + up], ext[j:j + down]
            else:
                left, right = template[max(0, j - up):j], template[j:j + down]
            frag = left + (QUERY if orient == "+" else sb.revcomp(QUERY)) + right
            truth = (j + 1, orient)
        if rng.random() < 0.5:
            frag = sb.revcomp(frag)
        r1 = noisy(frag[:readLen])
        r2 = noisy(sb.revcomp(frag)[:readLen])
        pairs.append(("read{}".format(i), r1, r2, truth))
    return pairs

def checkPairs(circular, seed):
    rng = random.Random(seed)
    template = "".join(rng.choice(B) for _ in range(7471))
    pairs = simulate(rng, template, 20000, circular)
    m = sb.Mapper(QUERY, template, 16, circular, 0)
    placed = wrong = withInsert = 0
    for _, r1, r2, truth in pairs:
        status, site, _ = sb.combinePair(m.mapRead(r1), m.mapRead(r2))
        withInsert += truth is not None
        if status == "placed":
            placed += 1
            wrong += site != truth
    assert wrong == 0, "{} pairs placed at the wrong site".format(wrong)
    assert placed / withInsert > 0.9, "only {}/{} inserts placed".format(placed, withInsert)
    return placed, withInsert

def test_linear():
    checkPairs(False, 1)

def test_circular():
    checkPairs(True, 2)

def test_cli_paired():
    rng = random.Random(3)
    template = "".join(rng.choice(B) for _ in range(3000))
    pairs = simulate(rng, template, 2000, False)
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "t.fa"), "w") as f:
            f.write(">t\n" + template + "\n")
        for mate, idx in (("R1", 1), ("R2", 2)):
            with gzip.open(os.path.join(d, mate + ".fastq.gz"), "wt") as f:
                for p in pairs:
                    f.write("@{} {}:N:0\n{}\n+\n{}\n".format(p[0], idx, p[idx], "I" * len(p[idx])))
        out = os.path.join(d, "s")
        sb.main(["-1", os.path.join(d, "R1.fastq.gz"), "-2", os.path.join(d, "R2.fastq.gz"),
                 "-q", QUERY.lower(), "-t", os.path.join(d, "t.fa"), "-o", out])
        rows = open(out + "_sites.csv").read().split("\n")[1:-1]
        counted = {}
        for r in rows:
            pos, o, _, c = r.split(",")
            counted[(int(pos), o)] = counted.get((int(pos), o), 0) + int(c)
        truthSites = {p[3] for p in pairs if p[3]}
        assert set(counted) <= truthSites
        assert sum(counted.values()) > 0.9 * sum(1 for p in pairs if p[3])

def writeFastq(path, records):
    with gzip.open(path, "wt") as f:
        for name, seq in records:
            f.write("@{} 1:N:0\n{}\n+\n{}\n".format(name, seq, "I" * len(seq)))

def test_empty_read_with_max_dist():
    '''A zero-length read (adapter-trimmed to nothing) must not crash the fuzzy path.'''
    try:
        import edlib  # noqa: F401
    except ImportError:
        print("test_empty_read_with_max_dist skipped (no edlib)")
        return
    rng = random.Random(4)
    template = "".join(rng.choice(B) for _ in range(2000))
    m = sb.Mapper(QUERY, template, 16, False, 2)
    assert m.mapRead("")[0] == "noQuery"
    assert m.mapRead("ACGT")[0] == "noQuery"
    # the query with one substitution is still found, and placed at the right site
    j = 700
    read = template[j - 60:j] + QUERY[:12] + "T" + QUERY[13:] + template[j:j + 60]
    status, site, d = m.mapRead(read)
    assert (status, site) == ("placed", (j + 1, "+")), (status, site, d)
    assert d == 1, d

def test_seed_prefilter_keeps_every_hit():
    '''The pigeonhole prefilter may not cost a single read: same calls with it and without.'''
    try:
        import edlib  # noqa: F401
    except ImportError:
        print("test_seed_prefilter_keeps_every_hit skipped (no edlib)")
        return
    rng = random.Random(5)
    template = "".join(rng.choice(B) for _ in range(3000))
    pairs = simulate(rng, template, 2000, False, err=0.02)   # enough errors to need edlib
    filtered = sb.Mapper(QUERY, template, 16, False, 2)
    assert filtered.seeds, "prefilter is not active"
    unfiltered = sb.Mapper(QUERY, template, 16, False, 2)
    unfiltered.seeds = ()
    fuzzy = 0
    for _, r1, r2, _ in pairs:
        for read in (r1, r2):
            a, b = filtered.mapRead(read), unfiltered.mapRead(read)
            assert a == b, (read, a, b)
            fuzzy += a[2] not in (None, 0)
    assert fuzzy > 50, "{} reads needed the fuzzy path; the test is not exercising it".format(fuzzy)

def test_procs_do_not_change_the_counts():
    '''Mapping in worker processes must give exactly the same output as mapping in one.'''
    rng = random.Random(6)
    template = "".join(rng.choice(B) for _ in range(3000))
    pairs = simulate(rng, template, 3000, False)
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "t.fa"), "w") as f:
            f.write(">t\n" + template + "\n")
        for mate, idx in (("R1", 1), ("R2", 2)):
            writeFastq(os.path.join(d, mate + ".fastq.gz"), [(p[0], p[idx]) for p in pairs])
        outputs = []
        for procs in ("1", "3"):
            out = os.path.join(d, "s" + procs)
            sb.main(["-1", os.path.join(d, "R1.fastq.gz"), "-2", os.path.join(d, "R2.fastq.gz"),
                     "-q", QUERY, "-t", os.path.join(d, "t.fa"), "-o", out,
                     "--procs", procs, "--progress", "500"])
            summary = [l for l in open(out + "_summary.txt") if not l.startswith("minutes")]
            outputs.append((open(out + "_sites.csv").read(), summary))
        assert outputs[0] == outputs[1], "--procs changed the result"

def test_out_of_sync_mates_are_caught():
    rng = random.Random(7)
    template = "".join(rng.choice(B) for _ in range(2000))
    pairs = simulate(rng, template, 60, False)
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "t.fa"), "w") as f:
            f.write(">t\n" + template + "\n")
        writeFastq(os.path.join(d, "R1.fastq.gz"), [(p[0], p[1]) for p in pairs])
        writeFastq(os.path.join(d, "R2.fastq.gz"), [(p[0], p[2]) for p in pairs[:30] + pairs[31:]])
        args = ["-1", os.path.join(d, "R1.fastq.gz"), "-2", os.path.join(d, "R2.fastq.gz"),
                "-q", QUERY, "-t", os.path.join(d, "t.fa"), "-o", os.path.join(d, "s"),
                "--procs", "1"]
        try:
            sb.main(args)
        except SystemExit as e:
            assert "out of sync at read 31" in str(e), e
        else:
            assert False, "a dropped R2 record went unnoticed"

if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print(name, "ok")
