#!/usr/bin/env python
'''
compare_sites.py
  ><```º>

Scores a caller's output against the truth table written by simulate_reads.py, and
compares two callers with each other.

Handles both output formats:
    --sites  <prefix>_sites.csv from stickleback_illumina.py (insPos_v,orientation,minD,count)
    --old    <input>_stickleback.csv from stickleback.0.3.py (per-read, has read + insPos_v)

With --old the comparison is per read, because that file carries read names; with --sites it is
per site, because the aggregate file does not. Both are reported against the same truth.

USAGE:
    python benchmark/compare_sites.py --truth bench/data/truth.csv \
        --sites bench/out/new_sites.csv --old bench/data/reads_stickleback.csv
'''

import argparse
import csv
import sys
from collections import Counter, defaultdict


def readTruth(path):
    '''read_id -> list of true (insPos_v, orientation) for the pair, plus read_id -> category.
    A chimeric pair has two, and a caller that reports either one (or calls it discordant) is
    not wrong; a pair with no insert has none.'''
    sites, categories = {}, {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            categories[row["read_id"]] = row["category"]
            s = []
            if row["insPos_v"]:
                s.append((int(row["insPos_v"]), row["orientation"]))
            if row.get("insPos_v2"):
                s.append((int(row["insPos_v2"]), row["orientation2"]))
            sites[row["read_id"]] = s
    return sites, categories


def readSites(path):
    '''(insPos_v, orientation) -> pairs, from stickleback_illumina.py.'''
    hist = Counter()
    with open(path) as fh:
        for row in csv.DictReader(fh):
            hist[(int(row["insPos_v"]), row["orientation"])] += int(row["count"])
    return hist


def readOld(path):
    '''read_id -> insPos_v, from stickleback.0.3.py. Keeps minD so we can report its filter.'''
    calls = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            name = row["read"].rsplit("/", 1)[0]
            pos = int(round(float(row["insPos_v"])))
            # a pair is one observation; keep the first mate seen, note disagreements
            calls.setdefault(name, []).append(pos)
    return calls


def truthHistogram(truthSites, categories, keep=None):
    hist = Counter()
    for read, sites in truthSites.items():
        if not sites:
            continue
        if keep and categories[read] not in keep:
            continue
        for site in sites:
            hist[site] += 1
    return hist


def offsetTable(called, truth):
    '''How far off each call is, as a histogram of signed offsets.'''
    return Counter(c - t for c, t in zip(called, truth))


def reportSites(hist, truth, label):
    out = ["### {} (site level)".format(label)]
    tot, ttot = sum(hist.values()), sum(truth.values())
    exact = sum(min(hist[s], truth[s]) for s in hist if s in truth)
    out.append("pairs called placed: {:,}   pairs with a true insert: {:,}".format(tot, ttot))
    out.append("pairs on a true (position, orientation): {:,} ({:.2%} of calls)".format(
        exact, exact / max(tot, 1)))
    extraSites = {s: hist[s] for s in hist if s not in truth}
    missedSites = {s: truth[s] for s in truth if s not in hist}
    out.append("sites called: {:,}   true sites: {:,}   called but not true: {:,} ({:,} pairs)"
               "   true but never called: {:,} ({:,} pairs)".format(
                   len(hist), len(truth), len(extraSites), sum(extraSites.values()),
                   len(missedSites), sum(missedSites.values())))
    # how far the spurious calls sit from the nearest true site of either orientation
    truthPos = {p for p, _ in truth}
    near = Counter()
    for (p, o), c in extraSites.items():
        d = min((abs(p - t) for t in truthPos), default=None)
        near[d if d is not None and d <= 3 else ">3"] += c
    if extraSites:
        out.append("spurious pairs by distance to nearest true position: {}".format(
            ", ".join("{}:{:,}".format(k, v) for k, v in sorted(near.items(), key=str))))
    if extraSites:
        worst = sorted(extraSites.items(), key=lambda kv: -kv[1])[:10]
        out.append("largest spurious sites: {}".format(
            ", ".join("{}{} x{}".format(p, o, c) for (p, o), c in worst)))
    return out


def reportOld(calls, truthSites, categories, label):
    out = ["### {} (per read pair)".format(label)]
    offs, byCat = Counter(), defaultdict(Counter)
    noTruth = 0
    for name, positions in calls.items():
        trueSites = truthSites.get(name) or []
        cat = categories.get(name, "?")
        for pos in positions:
            site = min(trueSites, key=lambda s: abs(s[0] - pos)) if trueSites else None
            if site is None:
                noTruth += 1
                byCat[cat]["called_on_insert_free_read"] += 1
                continue
            d = pos - site[0]
            offs[d if abs(d) <= 3 else ("far" if d > 0 else "far")] += 1
            byCat[cat]["exact" if d == 0 else ("off_by_{}".format(d) if abs(d) <= 3 else "far")] += 1
            byCat[cat]["orient_" + site[1]] += 1
    total = sum(offs.values()) + noTruth
    out.append("reads in output: {:,}   of which the read has no insert at all: {:,} ({:.1%})".format(
        total, noTruth, noTruth / max(total, 1)))
    out.append("offset of called insPos_v from truth: {}".format(
        ", ".join("{}:{:,}".format(k, v) for k, v in sorted(offs.items(), key=str))))
    for cat in sorted(byCat):
        out.append("  {:14s} {}".format(cat, dict(byCat[cat])))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Score caller output against simulate_reads.py truth.")
    p.add_argument("--truth", required=True)
    p.add_argument("--sites", help="stickleback_illumina.py <prefix>_sites.csv")
    p.add_argument("--old", help="stickleback.0.3.py per-read csv")
    p.add_argument("--categories", help="only score these truth categories, comma separated")
    args = p.parse_args(argv)

    truthSites, categories = readTruth(args.truth)
    keep = set(args.categories.split(",")) if args.categories else None
    truth = truthHistogram(truthSites, categories, keep)

    lines = ["## Caller vs truth",
             "truth: {:,} pairs, {:,} carry an insert, {:,} distinct (position, orientation)".format(
                 len(truthSites), sum(1 for s in truthSites.values() if s), len(truth)),
             "category mix: {}".format(dict(Counter(categories.values()))), ""]
    if args.sites:
        lines += reportSites(readSites(args.sites), truth, "stickleback_illumina.py") + [""]
    if args.old:
        lines += reportOld(readOld(args.old), truthSites, categories, "stickleback.0.3.py") + [""]
    print("\n".join(lines))


if __name__ == "__main__":
    main()
