#!/usr/bin/env python3
"""Build a taxonomy metadata TSV from an NCBI taxdump, offline.

    make_taxonomy_metadata.py <input> <taxdump_dir> <metadata.tsv>

<input> may be any of:
  * a sample sheet written by `defensome.py inspect` (uses its species column)
  * a directory of proteomes (species taken from the file names)
  * a text file with one species or file name per line

<taxdump_dir> needs names.dmp and nodes.dmp from
ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz. No API key, no rate limit.

Writes sample_id, species, matched, common_name and superkingdom..family.
Unmatched species are kept with blank ranks so they stay visible rather than
vanishing.
"""
import csv, glob, os, re, sys

RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family"]
FASTA = ("*.fa", "*.faa", "*.fasta", "*.fa.gz", "*.faa.gz", "*.fasta.gz", "*.pep.all.fa.gz")


def die(msg):
    sys.exit("ERROR: " + msg)


def load(taxdir):
    """One pass over names.dmp; it is ~250 MB and reading it twice is wasteful."""
    for f in ("names.dmp", "nodes.dmp"):
        p = os.path.join(taxdir, f)
        if not os.path.exists(p):
            die(f"{p} not found.\n"
                f"  <taxdump_dir> must be the directory holding names.dmp and nodes.dmp.\n"
                f"  If you unpacked into ./taxdump but are running from a subdirectory,\n"
                f"  give the path relative to where you are, e.g. ../taxdump")
    name2id, id2name, common = {}, {}, {}
    # "genbank common name" is the curated one NCBI shows on a taxon page and
    # there is at most one per node; "common name" is a free-for-all with
    # several per node, so it is only a fallback and the first wins.
    # "blast name" is deliberately ignored: it labels whole clades ("insects"),
    # which beside a binomial reads as an error.
    PREF = {"genbank common name": 2, "common name": 1}
    with open(os.path.join(taxdir, "names.dmp"), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) < 4:
                continue
            if f[3] == "scientific name":
                tid = int(f[0]); name2id[f[1]] = tid; id2name[tid] = f[1]
            elif f[3] in PREF:
                tid, rank = int(f[0]), PREF[f[3]]
                if rank > common.get(tid, (0, ""))[0]:
                    common[tid] = (rank, f[1])
    common = {k: v for k, (_, v) in common.items()}
    nodes = {}
    with open(os.path.join(taxdir, "nodes.dmp"), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            nodes[int(f[0])] = (int(f[1]), f[2])
    return name2id, nodes, id2name, common


def lineage(taxid, nodes, id2name):
    out, seen = {}, set()
    while taxid and taxid not in seen:
        seen.add(taxid)
        parent, rank = nodes.get(taxid, (None, None))
        if rank in RANKS:
            out[rank] = id2name.get(taxid, "")
        if taxid == 1 or parent is None:
            break
        taxid = parent
    return out


def clean(name):
    """File name or sample id to a binomial."""
    n = re.sub(r"\.(fa|faa|fasta)(\.gz)?$", "", str(name).strip())
    n = re.sub(r"(_pep|_protein|_prot|_aa|_cds|_longest|_primary|\.pep\.all)$",
               "", n, flags=re.I)
    return n.replace("_", " ").strip()


def read_input(src):
    """-> list of (sample_id, species)"""
    if os.path.isdir(src):
        files = sorted({f for p in FASTA for f in glob.glob(os.path.join(src, p))})
        if not files:
            die(f"no FASTA files in the directory {src}")
        return [(re.sub(r"\.(fa|faa|fasta)(\.gz)?$", "", os.path.basename(f)),
                 clean(os.path.basename(f))) for f in files]
    if not os.path.exists(src):
        die(f"{src} not found.\n"
            "  Give a sample sheet, a directory of proteomes, or a list of species.\n"
            "  A sample sheet is written by:  python3 defensome.py inspect "
            "--proteomes <dir> --out <results>")
    with open(src) as fh:
        head = fh.readline()
    if "\t" in head and "sample_id" in head:
        rows = list(csv.DictReader(open(src), delimiter="\t"))
        if not rows:
            die(f"{src} has a header but no rows")
        return [(r["sample_id"], clean(r.get("species") or r["sample_id"])) for r in rows]
    return [(clean(l), clean(l)) for l in open(src) if l.strip()]


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    src, taxdir, out = sys.argv[1:4]
    pairs = read_input(src)
    name2id, nodes, id2name, common = load(taxdir)
    n_ok, n_common, unmatched = 0, 0, []
    with open(out, "w") as f:
        f.write("\t".join(["sample_id", "species", "matched", "common_name"] + RANKS) + "\n")
        for sid, sp in pairs:
            tid = name2id.get(sp)
            if tid is None and sp.count(" ") > 1:
                tid = name2id.get(" ".join(sp.split()[:2]))      # drop subspecies
            genus_only = False
            if tid is None and sp.count(" ") == 1:
                # Fall back to the genus. Everything at family rank and above is
                # still correct, which is all this file is used for.
                tid = name2id.get(sp.split()[0])
                genus_only = bool(tid)
            lin = lineage(tid, nodes, id2name) if tid else {}
            n_ok += bool(tid)
            if not tid:
                unmatched.append(sp)
            mark = "genus" if genus_only else ("yes" if tid else "no")
            # A genus-level match gives a genus-level vernacular, which beside a
            # binomial would be wrong, so only exact matches carry one.
            cn = "" if (tid is None or genus_only) else common.get(tid, "")
            n_common += bool(cn)
            f.write("\t".join([sid, sp.replace(" ", "_"), mark, cn] +
                              [lin.get(k, "") for k in RANKS]) + "\n")
    print(f"wrote {out}: {n_ok}/{len(pairs)} matched to NCBI taxonomy, "
          f"{n_common} with a common name")
    if unmatched:
        print(f"  {len(unmatched)} unmatched, kept with blank ranks: "
              + ", ".join(unmatched[:8]) + (" ..." if len(unmatched) > 8 else ""))
    print(f"\nThen: python3 defensome.py dashboard --out <results> "
          f"--metadata {out} --group-by class")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:      # piping into head
        try:
            sys.stdout.close()
        except Exception:
            pass
