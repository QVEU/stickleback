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

if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print(name, "ok")
