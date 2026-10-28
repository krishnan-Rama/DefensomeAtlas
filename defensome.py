#!/usr/bin/env python3
"""
defensome - comparative chemical defensome annotation from proteomes.

  defensome.py run --proteomes DIR --out DIR                 everything, to the dashboard
  defensome.py run --proteomes DIR --out DIR --species "Danio rerio,Bombyx mori"
  defensome.py run ... --taxdump taxdump/ | --metadata traits.tsv --group-by habitat
  defensome.py doctor                                         what can run here

Every step is also its own subcommand (scan, annotate, qc, domains, trees,
dashboard, structure, ...); `defensome.py <step> --help`.

Only hard dependencies: HMMER 3 (for `scan`), Python 3.9+, pandas, numpy.
matplotlib is optional; figures are skipped with a warning if it is missing.
"""
import argparse, glob, gzip, math, os, re, shutil, subprocess, sys

# numpy/pandas are needed by annotate and report, but NOT by scan. Importing
# them lazily means `scan` runs in a bare HMMER module environment, which
# matters on clusters where the HMMER and SciPy toolchains cannot coexist.
try:
    import numpy as np
    import pandas as pd
except ImportError:  # pragma: no cover
    np = pd = None


def need_pandas():
    if pd is None:
        sys.exit("ERROR: this step needs pandas and numpy.\n"
                 "  module load SciPy-bundle    (or: pip install --user pandas numpy)")

__version__ = "5.2.0"

DOMTBL_COLS = ["target","t_acc","tlen","query","q_acc","qlen","e_full","score_full",
               "bias_full","dom_n","dom_of","c_evalue","i_evalue","score_dom","bias_dom",
               "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","acc"]

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MAP = os.path.join(_HERE, "defensome_map.tsv")

# The dashboard template, its JavaScript and the default map are embedded below
# so that defensome.py works on its own from any directory. A file of the same
# name sitting next to the script wins, so you can still edit any of them:
#   python3 defensome.py assets --write .
_ASSETS = {
    "defensome_map.tsv": (
        "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"
        "OsiTaHkxWT6t22qf465N/3dSzTT/B6wiPv8="
    ),
    "dashboard.html": (
        "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"
        "H/NN7QY="
    ),
    "dashboard.js": (
        "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"
        "PRSZmylA3lv2YENw0MOCXChtgB4EkrL39JSFt3wjHAFM6QMkU57M5Zaktp2b6SRzKItY0wWv3of84eiVFe2AfsCg+szd00JlBIlB"
        "NW7nYAHe9qsQO2eFkC2OwqkfrEqGQQkKw2oZ/ErB4M7Bwii4/v8s/nmnDW0rCvA53H/bMngfHArSkMzpSxH0iKt1sBCepLT6AJyk"
        "sD3Os0VfMUgHa+00ZcvlbBbnpCjFZkOogrQuCwj4XduXtMu1jJo3cGfswplyKsevMndOdUmjXXggbGl9XthSZeWKObmsL5fimL4O"
        "ifLh+KCnCgD57LcBgHwaA0A+CwCQJg/RIoCP37XjxDhMvGW5xIldqvAREZ2OkCRHpRbDRZxNwt0ZT5Rb5zsQ7wqJuQovV2jh34qY"
        "mxiYk+mCpJqPCO1sNVROFOb01XFZF8vJjai9bpzTvuJizojIqgiW0gqmC8QwPYPZ1MkJzsu4SovpeOBCrFN1mNO8ECekYXXhFuxP"
        "weAi5uYc8uccJHSpUVWjs2MdJeTKM5cg5zcubiFqBncq2pDY91BDVTecd85JyKlbPSn2lJ5aF/TL8YLrzfMiuubozVk4Fq5gBFf8"
        "K+RU5grkc44FwNh5dWaDdqTsQOmM5lKhMBeJesuW585zjYNKUKEfo7NfU/sd1GzRi/qENzB3dZXjaBmrtPHx1xMLmHlSfBTQEfRA"
        "to91TZ6wOU2ic+99TL5G4U8+DYlw/tqnjmnI9TKVvOigy87pbX7PI//vy4uKmPZZ+TFOL+va42PRKTWm+fUJlOW67P8yRW4CHDAm"
        "lVrtL1OTvCCYyOlwKCXIa9GvoLCiqefsg8XcbCE7oZTQg8O5S/1lWkvDxVCsXEHRU+x/uuyQ9m1IldPRVDCzNLzCtHs+u3fDAY1/"
        "QshOPYLxYmdKnwHR5USkD3K5feDtBl9ufL3P92OS+xiEr2Bw5gpf7PbG/dPRmKhxZuxoQ0nIG9K8eYqbF1WD9dxm55uo9vM6TYYt"
        "hs2q839RUT7mVN/2xBiNy+cFIssAoR60fB3gJa8kGZgdheAK0V6T+xeP7lDRS1TX4uwUtZ+bux922NPVAv9djjrZpQZ8Jl+nRJZ+"
        "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"
    ),
}


def asset(name):
    """Sibling file if present, otherwise the copy embedded in this script."""
    p = os.path.join(_HERE, name)
    if os.path.exists(p):
        return open(p).read()
    import base64, zlib
    if name not in _ASSETS:
        sys.exit(f"ERROR: no embedded asset '{name}' and no file at {p}")
    return zlib.decompress(base64.b64decode(_ASSETS[name])).decode()


def cmd_assets(a):
    """Write the embedded assets out so they can be edited."""
    d = a.write or "."
    os.makedirs(d, exist_ok=True)
    for name in ("defensome_map.tsv", "dashboard.html", "dashboard.js"):
        dest = os.path.join(d, name)
        if os.path.exists(dest) and not a.force:
            log(f"exists, skipping (use --force): {dest}"); continue
        import base64, zlib
        with open(dest, "w") as fh:
            fh.write(zlib.decompress(base64.b64decode(_ASSETS[name])).decode())
        log(f"wrote {dest}")


# --------------------------------------------------------------------------- #
def log(msg):
    print(f"[defensome] {msg}", file=sys.stderr, flush=True)


def species_list(proteome_dir):
    fa = sorted(glob.glob(os.path.join(proteome_dir, "*.faa")) +
                glob.glob(os.path.join(proteome_dir, "*.fasta")) +
                glob.glob(os.path.join(proteome_dir, "*.fa")))
    if not fa:
        sys.exit(f"ERROR: no .faa/.fa/.fasta files in {proteome_dir}")
    return [(os.path.splitext(os.path.basename(f))[0], f) for f in fa]


def read_map(path):
    need_pandas()
    import io
    src = path if os.path.exists(path) else None
    if src is None and os.path.abspath(path) == os.path.abspath(DEFAULT_MAP):
        df = pd.read_csv(io.StringIO(asset("defensome_map.tsv")), sep="\t", comment="#")
    elif src is None:
        sys.exit(f"ERROR: map file not found: {path}")
    else:
        df = pd.read_csv(src, sep="\t", comment="#")
    need = {"category","family","pfam_ids","rule","min_cov","min_len","tier"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path} is missing columns: {sorted(missing)}")
    df["pfam_ids"] = df.pfam_ids.str.replace(" ", "").str.split(",")
    # optional columns (older maps do not have them)
    if "rescue" not in df.columns:
        df["rescue"] = "no"
    df["rescue"] = df.rescue.fillna("no").astype(str).str.strip().str.lower()
    if "validate" not in df.columns:
        df["validate"] = ""
    df["validate"] = df.validate.fillna("").astype(str).str.strip()
    bad = df[~df.rule.isin(["ANY", "ALL"])]
    if len(bad):
        sys.exit(f"ERROR: rule must be ANY or ALL; bad rows: {list(bad.family)}")
    return df


# =========================================================================== #
#  Pfam HMM file utilities (pure Python, no HMMER needed)
# =========================================================================== #
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
# The defensome HMMs, prebuilt from one Pfam release (Pfam is CC0), so a user
# does not need the 300 MB Pfam-A. Published as a GitHub release asset; see
# docs/RELEASING.md. DEFENSOME_DB_URL overrides it (a mirror, or a file path).
DB_REPO = "krishnan-Rama/DEFENSOME"
DB_TAG = "pfam-db"
DB_ASSET = "defensome_db.tar.gz"
DB_URL = os.environ.get("DEFENSOME_DB_URL",
                        f"https://github.com/{DB_REPO}/releases/download/{DB_TAG}/{DB_ASSET}")
PFAM_VERSION_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz"


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def iter_hmm_records(path):
    """Yield (header_dict, record_text) for each HMM in a HMMER3 text file."""
    buf, hdr = [], {}
    with _open_text(path) as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("//"):
                yield hdr, "".join(buf)
                buf, hdr = [], {}
                continue
            if len(line) > 6 and line[:6].rstrip() in ("NAME", "ACC", "LENG", "GA", "DESC"):
                key = line[:6].rstrip()
                if key not in hdr:
                    hdr[key] = line[6:].strip()


def hmm_index(path):
    """accession (no version) -> dict(name, leng, ga_seq, ga_dom, desc)"""
    idx = {}
    for h, _ in iter_hmm_records(path):
        acc = h.get("ACC", "").split(".")[0]
        if not acc:
            continue
        ga = h.get("GA", "").rstrip(";").split()
        idx[acc] = {"name": h.get("NAME", ""), "leng": int(h.get("LENG", "0") or 0),
                    "ga_seq": float(ga[0]) if len(ga) > 0 else None,
                    "ga_dom": float(ga[1]) if len(ga) > 1 else None,
                    "desc": h.get("DESC", "")}
    return idx


def write_hmm_subset(src, accessions, dest):
    """Copy only the HMMs whose accession is in `accessions`. Atomic: parallel
    array tasks may race to build the same file; identical content, os.replace."""
    want = {a.split(".")[0] for a in accessions}
    found = set()
    d = os.path.dirname(os.path.abspath(dest))
    os.makedirs(d, exist_ok=True)
    tmp = f"{dest}.{os.getpid()}.tmp"
    with open(tmp, "w") as out:
        for h, rec in iter_hmm_records(src):
            acc = h.get("ACC", "").split(".")[0]
            if acc in want:
                out.write(rec)
                found.add(acc)
    os.replace(tmp, dest)
    return found


def _map_rows_plain(mapfile):
    """Minimal map reader for steps that must run without pandas (scan)."""
    import csv, io
    if os.path.exists(mapfile):
        text = open(mapfile).read()
    else:
        text = asset("defensome_map.tsv")
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    rows = list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))
    for r in rows:
        r["pfam_ids"] = [x.strip() for x in (r.get("pfam_ids") or "").split(",") if x.strip()]
        r["rescue"] = (r.get("rescue") or "no").strip().lower()
    return rows


def _download(url, dest):
    import urllib.request
    log(f"downloading {url}")
    tmp = dest + ".part"
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as fh:
        total = int(r.headers.get("Content-Length") or 0)
        got, last = 0, 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            if total and got - last > total // 20:
                last = got
                print(f"   {got/1e6:8.1f} / {total/1e6:.1f} MB", file=sys.stderr)
    os.replace(tmp, dest)


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_db_bundle(db_dir):
    """defensome.hmm + MANIFEST.tsv + PFAM_VERSION as one tarball with a
    checksum: the file to attach to the GitHub release."""
    import tarfile
    out = os.path.join(db_dir, DB_ASSET)
    with tarfile.open(out, "w:gz") as tf:
        for f in ("defensome.hmm", "MANIFEST.tsv", "PFAM_VERSION"):
            p_ = os.path.join(db_dir, f)
            if os.path.exists(p_):
                tf.add(p_, arcname=f)
    open(out + ".sha256", "w").write(f"{_sha256(out)}  {DB_ASSET}\n")
    return out


def fetch_db(db_dir, url=None, force=False):
    """Get the prebuilt defensome HMMs into db_dir: a bundle beside the tool
    first (offline), then the release asset. Returns the path of the HMM file."""
    import tarfile
    hmm = os.path.join(db_dir, "defensome.hmm")
    if os.path.exists(hmm) and not force:
        return hmm
    os.makedirs(db_dir, exist_ok=True)
    local = [p_ for p_ in (os.path.join(_HERE, "db", DB_ASSET), os.path.join(_HERE, DB_ASSET))
             if os.path.exists(p_)]
    url = url or DB_URL
    tgz = os.path.join(db_dir, DB_ASSET)
    if local:
        shutil.copy(local[0], tgz)
        log(f"using the bundled HMM database {local[0]}")
    elif os.path.exists(url):
        shutil.copy(url, tgz)
    else:
        try:
            _download(url, tgz)
        except Exception as e:
            sys.exit(f"ERROR: could not download the defensome HMMs from\n  {url}\n  ({e})\n"
                     "  Offline or behind a firewall? Either download that file yourself and\n"
                     "  pass --db-url /path/to/defensome_db.tar.gz, or build the database from\n"
                     "  a Pfam-A you have:  --pfam /path/to/Pfam-A.hmm.gz")
        try:
            _download(url + ".sha256", tgz + ".sha256")
            want = open(tgz + ".sha256").read().split()[0]
            if want != _sha256(tgz):
                sys.exit(f"ERROR: checksum mismatch for {tgz}; delete it and retry")
            log("checksum OK")
        except SystemExit:
            raise
        except Exception:
            log("no checksum published beside the database; skipping verification")
    with tarfile.open(tgz) as tf:
        for m in tf.getmembers():
            if m.isfile() and os.path.basename(m.name) == m.name:
                tf.extract(m, db_dir)
    if not os.path.exists(hmm):
        sys.exit(f"ERROR: {tgz} does not contain defensome.hmm")
    v = os.path.join(db_dir, "PFAM_VERSION")
    log(f"HMM database ready: {hmm}" + (f" ({open(v).read().split(chr(10))[0].strip()})"
                                         if os.path.exists(v) else ""))
    return hmm


def check_db_covers_map(hmm, mapfile):
    """A map edited to use an accession the prebuilt database lacks would make
    that family silently empty."""
    rows = _map_rows_plain(mapfile)
    want = {x for r in rows for x in r["pfam_ids"]}
    have = set(hmm_index(hmm))
    lost = sorted(want - have)
    if lost:
        log(f"WARNING: the HMM database has no model for {', '.join(lost)}, used by the map; "
            "those families will be empty. Rebuild from a full Pfam-A with --pfam.")
    return lost


def cmd_setup(a):
    """Prepare the HMM database: extract the defensome families from Pfam-A.

    Scanning 30,000 proteins against all ~25,000 Pfam families takes hours per
    proteome. The map needs about 45 of them. Family calls are identical either
    way, because hmmsearch scores each HMM independently and --cut_ga uses
    per-family bit-score thresholds, so the subset gives the same counts in a
    fraction of the time. The full Pfam-A is only needed for the optional
    architecture pass, which reports every other domain on defensome proteins.

    Writes <db-dir>/defensome.hmm and <db-dir>/MANIFEST.tsv. The manifest says
    which accessions your Pfam release actually contains: an accession missing
    from it is the first thing to rule out when a family comes back empty."""
    rows = _map_rows_plain(a.map)
    accs = sorted({x for r in rows for x in r["pfam_ids"]})
    os.makedirs(a.db_dir, exist_ok=True)
    if getattr(a, "fetch", False) or (not a.pfam and not a.download):
        hmm = fetch_db(a.db_dir, getattr(a, "db_url", None), force=a.force)
        check_db_covers_map(hmm, a.map)
        print(f"\nscan with:  python3 defensome.py scan --pfam {hmm} ...")
        return
    src = a.pfam
    if a.download:
        src = os.path.join(a.db_dir, "Pfam-A.hmm.gz")
        if not os.path.exists(src) or a.force:
            _download(a.url or PFAM_URL, src)
        try:
            vgz = os.path.join(a.db_dir, "Pfam.version.gz")
            _download(PFAM_VERSION_URL, vgz)
            with gzip.open(vgz, "rt") as fh:
                open(os.path.join(a.db_dir, "PFAM_VERSION"), "w").write(fh.read())
        except Exception as e:
            log(f"could not fetch Pfam.version ({e}); release not recorded")
    if not src or not os.path.exists(src):
        sys.exit("ERROR: give --pfam /path/to/Pfam-A.hmm[.gz], or --download")
    # which Pfam release this is: stated, or the Pfam.version file EBI ships beside it
    ver_f = os.path.join(a.db_dir, "PFAM_VERSION")
    if getattr(a, "pfam_version", None):
        open(ver_f, "w").write(a.pfam_version.strip() + "\n")
    elif not a.download:
        for c in ("Pfam.version", "Pfam.version.gz"):
            sib = os.path.join(os.path.dirname(os.path.abspath(src)), c)
            if os.path.exists(sib):
                op_ = gzip.open if sib.endswith(".gz") else open
                with op_(sib, "rt") as fh:
                    open(ver_f, "w").write(fh.read())
                break

    log(f"indexing {src} (reads the whole file once; a few minutes for full Pfam-A)")
    dest = os.path.join(a.db_dir, "defensome.hmm")
    found = write_hmm_subset(src, accs, dest)
    idx = hmm_index(dest)
    fam_of = {}
    for r in rows:
        for x in r["pfam_ids"]:
            fam_of.setdefault(x, []).append(r["family"])
    man = os.path.join(a.db_dir, "MANIFEST.tsv")
    with open(man, "w") as fh:
        fh.write("pfam_id\tfamilies\tfound\tname\tlength\tga_seq\tga_dom\tshort_domain\tdescription\n")
        for x in accs:
            i = idx.get(x, {})
            short = "yes" if i and i["leng"] and i["leng"] < 60 else ("no" if i else "")
            fh.write(f"{x}\t{','.join(fam_of.get(x, []))}\t{'yes' if x in found else 'NO'}\t"
                     f"{i.get('name','')}\t{i.get('leng','')}\t{i.get('ga_seq','')}\t"
                     f"{i.get('ga_dom','')}\t{short}\t{i.get('desc','')}\n")
    missing = [x for x in accs if x not in found]
    print(f"\n=== {len(found)} of {len(accs)} map accessions found in {os.path.basename(src)} ===")
    for x in accs:
        i = idx.get(x)
        if i:
            flag = "  SHORT DOMAIN" if i["leng"] < 60 else ""
            print(f"  {x}  {i['name']:18s} length {i['leng']:4d}  GA {i['ga_seq']}{flag}")
    if missing:
        print("\nMISSING from this Pfam release (these families cannot be called):")
        for x in missing:
            print(f"  {x}  used by {', '.join(fam_of.get(x, []))}")
        print("Check the accession on https://www.ebi.ac.uk/interpro/ ; Pfam retires")
        print("and merges families between releases.")
    if a.download and a.keep_full:
        full = os.path.join(a.db_dir, "Pfam-A.hmm")
        with gzip.open(src, "rb") as fi, open(full, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        log(f"kept full Pfam-A for the optional architecture pass: {full}")
    ver = os.path.join(a.db_dir, "PFAM_VERSION")
    if os.path.exists(ver):
        print("\n" + open(ver).read().strip())
    print(f"\nwrote {dest}  ({os.path.getsize(dest)/1e6:.1f} MB, {len(found)} HMMs)")
    print(f"      {man}")
    if getattr(a, "bundle", False):
        if not os.path.exists(ver):
            log("WARNING: the bundle will not say which Pfam release it is. Put Pfam.version "
                "beside the Pfam file, or pass --pfam-version 'Pfam release 37.4'.")
        b = write_db_bundle(a.db_dir)
        print(f"\nbundle for the GitHub release: {b}\n                               {b}.sha256"
              f"\n  gh release create {DB_TAG} {b} {b}.sha256 --repo {DB_REPO} "
              f"--title 'Defensome HMMs' --notes-file {os.path.join(a.db_dir, 'PFAM_VERSION')}")
    print(f"\nscan with:  python3 defensome.py scan --pfam {dest} ...")


# --------------------------------------------------------------------------- #
def _gz_move(tmp, dest):
    with open(tmp, "rb") as fi, gzip.open(dest, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(tmp)


# HMMER's comparison pipeline aborts (SIGABRT, "Target sequence length > 100K")
# on any target over this length, killing the whole search rather than skipping
# the sequence. Genome annotations occasionally emit one: a BRAKER model that
# merged a scaffold produced a single 100K+ residue "protein" in one DToL
# sawfly proteome and took the entire 1,575-species scan down with it.
HMMER_MAX_TARGET = 100_000


def filter_oversized(faa, workdir, sp):
    """Return a FASTA path safe for hmmsearch, plus the records it dropped.

    If nothing is oversized the original path is returned untouched, so the
    common case costs one pass over the file and no extra disk."""
    big = [(h.split()[0], len(q)) for h, q in read_fasta(faa) if len(q) > HMMER_MAX_TARGET]
    if not big:
        return faa, []
    os.makedirs(workdir, exist_ok=True)
    clean = os.path.join(workdir, f"{sp}.filtered.faa")
    drop = {p for p, _ in big}
    with open(clean, "w") as fh:
        for h, q in read_fasta(faa):
            if h.split()[0] not in drop:
                fh.write(f">{h}\n{q}\n")
    return clean, big


def cmd_scan(a):
    """hmmsearch per proteome, in up to three passes.

    1. Main pass, --cut_ga: Pfam's curated per-family bit-score thresholds.
       An E-value cutoff would depend on proteome size and so would not be
       comparable across the proteomes being compared. -> hmmsearch/
    2. Rescue pass, for families marked rescue=yes in the map: the same HMMs
       with an E-value cutoff instead of --cut_ga. Very short domains have few
       positions to accumulate score, so a genuine but divergent copy can sit
       below a threshold tuned on other taxa. These hits are kept apart and
       must pass the map's validate rules before they count. -> rescue/
    3. Architecture pass, only with --full-pfam: every Pfam domain on the
       proteins that hit a defensome family, for the architecture strings.
       Restricted to those proteins, so it takes minutes. -> hmmsearch_full/
    """
    hs = tool("hmmsearch")
    if not hs:
        sys.exit("ERROR: hmmsearch not found.\n"
                 "  module load HMMER    (or set DEFENSOME_HMMSEARCH=/path/to/hmmsearch)\n"
                 "  python3 defensome.py doctor   shows everything that is missing")
    outdir = os.path.join(a.out, "hmmsearch"); os.makedirs(outdir, exist_ok=True)
    # record which database was searched, so qc can find its MANIFEST later
    with open(os.path.join(a.out, "scan_db.txt"), "w") as fh:
        fh.write(os.path.abspath(a.pfam) + "\n")
    todo = species_list(a.proteomes)
    if a.species:
        todo = [(s_, f) for s_, f in todo if s_ == a.species]
        if not todo:
            sys.exit(f"ERROR: species '{a.species}' not found in {a.proteomes}")

    # rescue HMMs: built once from the scan database, atomically
    rescue_hmm, rescue_accs = None, set()
    if not a.no_rescue:
        rows = _map_rows_plain(a.map)
        rescue_accs = {x for r in rows if r["rescue"] == "yes" for x in r["pfam_ids"]}
        if rescue_accs:
            rdir = os.path.join(a.out, "rescue"); os.makedirs(rdir, exist_ok=True)
            rescue_hmm = os.path.join(rdir, "rescue_families.hmm")
            if not os.path.exists(rescue_hmm) or a.force:
                got = write_hmm_subset(a.pfam, rescue_accs, rescue_hmm)
                lost = sorted(rescue_accs - got)
                if lost:
                    log(f"WARNING: rescue accessions not in {a.pfam}: {lost}")
                if not got:
                    rescue_hmm = None

    oversized_log = []
    for sp, faa in todo:
        # Must happen before any pass: all three call hmmsearch on this file.
        faa, big = filter_oversized(faa, os.path.join(a.out, "oversized"), sp)
        if big:
            log(f"{sp}: {len(big)} sequence(s) over {HMMER_MAX_TARGET:,} aa excluded; "
                "hmmsearch aborts on these rather than skipping them")
            for pid, L in big:
                log(f"    {pid}  {L:,} aa")
                oversized_log.append((sp, pid, L))
        dest = os.path.join(outdir, f"{sp}.domtblout.gz")
        if os.path.exists(dest) and not a.force:
            log(f"{sp}: main pass exists, skipping (use --force to redo)")
        else:
            log(f"{sp}: hmmsearch --cut_ga")
            tmp = dest[:-3]
            subprocess.run([hs, "--cut_ga", "--cpu", str(a.threads),
                            "--domtblout", tmp, "-o", os.devnull, a.pfam, faa], check=True)
            _gz_move(tmp, dest)

        if rescue_hmm:
            rdest = os.path.join(a.out, "rescue", f"{sp}.domtblout.gz")
            if not os.path.exists(rdest) or a.force:
                log(f"{sp}: rescue pass, E <= {a.rescue_evalue}, no --cut_ga")
                tmp = rdest[:-3]
                subprocess.run([hs, "-E", str(a.rescue_evalue), "--domE", str(a.rescue_evalue),
                                "--cpu", str(a.threads), "--domtblout", tmp,
                                "-o", os.devnull, rescue_hmm, faa], check=True)
                _gz_move(tmp, rdest)

        if a.full_pfam:
            fdir = os.path.join(a.out, "hmmsearch_full"); os.makedirs(fdir, exist_ok=True)
            fdest = os.path.join(fdir, f"{sp}.domtblout.gz")
            if not os.path.exists(fdest) or a.force:
                cand = set()
                with gzip.open(dest, "rt") as fh:
                    for line in fh:
                        if not line.startswith("#") and line.strip():
                            cand.add(line.split(None, 1)[0])
                sub = os.path.join(fdir, f"{sp}.candidates.faa")
                n = 0
                with open(sub, "w") as fo:
                    for h, seq in read_fasta(faa):
                        if h.split()[0] in cand:
                            fo.write(f">{h}\n{seq}\n"); n += 1
                log(f"{sp}: architecture pass on {n} candidate proteins against full Pfam")
                tmp = fdest[:-3]
                subprocess.run([hs, "--cut_ga", "--cpu", str(a.threads), "--domtblout", tmp,
                                "-o", os.devnull, a.full_pfam, sub], check=True)
                _gz_move(tmp, fdest)
                os.remove(sub)
    if oversized_log:
        p = os.path.join(a.out, "oversized_proteins.tsv")
        new_file = not os.path.exists(p)
        with open(p, "a") as fh:
            if new_file:
                fh.write("species\tprotein\tlength_aa\n")
            for sp_, pid, L in oversized_log:
                fh.write(f"{sp_}\t{pid}\t{L}\n")
        log(f"{len(oversized_log)} oversized sequence(s) recorded in {p}")
    log("scan done")


# --------------------------------------------------------------------------- #
def read_domtbl(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 22)
            if len(f) >= 22:
                rows.append(f[:22])
    if not rows:
        return pd.DataFrame(columns=DOMTBL_COLS + ["pfam_id"])
    df = pd.DataFrame(rows, columns=DOMTBL_COLS)
    for c in ("tlen", "qlen", "hmm_from", "hmm_to"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["pfam_id"] = df.q_acc.str.split(".").str[0]
    return df.dropna(subset=["tlen", "qlen"])


def domain_coverage(df):
    """Per (protein, pfam), sum non-overlapping HMM coverage.
    Repeat domains and split hits are merged rather than double counted."""
    out = []
    for (prot, pf), g in df.groupby(["target", "pfam_id"], sort=False):
        iv = sorted(zip(g.hmm_from.astype(int), g.hmm_to.astype(int)))
        merged, cs, ce = [], None, None
        for s, e in iv:
            if cs is None:
                cs, ce = s, e
            elif s <= ce + 1:
                ce = max(ce, e)
            else:
                merged.append((cs, ce)); cs, ce = s, e
        if cs is not None:
            merged.append((cs, ce))
        qlen = float(g.qlen.iloc[0])
        out.append((prot, pf, int(g.tlen.iloc[0]),
                    sum(e - s + 1 for s, e in merged) / qlen if qlen else 0.0))
    return pd.DataFrame(out, columns=["protein", "pfam_id", "prot_len", "hmm_cov"])


def count_proteins(faa):
    n = 0
    op = gzip.open if faa.endswith(".gz") else open
    with op(faa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


_VAL_RE = re.compile(r"^(len|cys|copies)\s*(<=|>=|<|>|==)\s*([0-9.]+)$")


def _parse_validate(spec, family=""):
    rules = []
    for part in str(spec or "").split(";"):
        part = part.strip()
        if not part:
            continue
        m = _VAL_RE.match(part)
        if not m:
            sys.exit(f"ERROR: bad validate rule '{part}' for {family}. "
                     "Use len<=N, len>=N, cys>=F or copies>=N, separated by ';'")
        rules.append((m.group(1), m.group(2), float(m.group(3))))
    return rules


def _check_rules(rules, vals):
    """vals: dict key -> number or None. Returns (ok, reason)."""
    ops = {"<=": lambda x, y: x <= y, ">=": lambda x, y: x >= y, "<": lambda x, y: x < y,
           ">": lambda x, y: x > y, "==": lambda x, y: x == y}
    for k, op, v in rules:
        x = vals.get(k)
        if x is None:
            return False, f"{k} unavailable (proteome not given)"
        if not ops[op](x, v):
            return False, f"{k}={x:.3g} fails {k}{op}{v:g}"
    return True, "passes"


MT_RULE = {"min_len": 25, "max_len": 120, "min_cys": 8, "min_cys_frac": 0.18,
           "min_motifs": 4, "max_aromatic": 1}


def mt_features(seq):
    u = seq.upper().rstrip("*")
    L = len(u)
    c = u.count("C")
    return {"length": L, "n_cys": c, "cys_fraction": round(c / L, 3) if L else 0.0,
            "n_motifs": len(re.findall(r"(?=C.C|C..C)", u)),
            "n_aromatic": sum(u.count(x) for x in "FWY")}


def _is_mt_like(seq):
    """Composition screen for metallothionein-like proteins.

    No Pfam model covers insect metallothioneins (InterPro IPR000966, family 5,
    has no Pfam member), so they are found by what defines them: short, rich
    in cysteine arranged in C-x-C and C-x-x-C motifs, and devoid of aromatic
    residues. The aromatic rule does the discriminating work. Cysteine alone
    admits defensins, knottins and chitin-binding peptides, which are also
    short and cysteine-rich but carry F, W or Y; an earlier version of this
    screen used cysteine only and accepted such a decoy while rejecting a real
    metallothionein at 23% cysteine. Heuristic: reported apart from HMM calls."""
    f = mt_features(seq)
    r = MT_RULE
    return (r["min_len"] <= f["length"] <= r["max_len"] and f["n_cys"] >= r["min_cys"]
            and f["cys_fraction"] >= r["min_cys_frac"] and f["n_motifs"] >= r["min_motifs"]
            and f["n_aromatic"] <= r["max_aromatic"])


def _copies(df):
    """distinct non-overlapping target-sequence segments hit by one HMM"""
    iv = sorted(zip(pd.to_numeric(df.ali_from).astype(int), pd.to_numeric(df.ali_to).astype(int)))
    n, ce = 0, -1
    for s_, e_ in iv:
        if s_ > ce:
            n += 1
            ce = e_
        else:
            ce = max(ce, e_)
    return n


def cmd_annotate(a):
    need_pandas()
    dm = read_map(a.map)
    scan_dir = os.path.join(a.out, "hmmsearch")
    files = sorted(glob.glob(os.path.join(scan_dir, "*.domtblout*")))
    if not files:
        sys.exit(f"ERROR: no domtblout files in {scan_dir}. Run `scan` first.")
    calls_dir = os.path.join(a.out, "gene_calls"); os.makedirs(calls_dir, exist_ok=True)
    rescue_fams = dm[dm.rescue == "yes"]
    has_mt = "Metallothionein" in set(dm.family)
    faa_of = dict(species_list(a.proteomes)) if a.proteomes and os.path.isdir(a.proteomes) else {}

    prot_totals, stat_rows, comp_rows, resc_rows, all_calls = {}, [], [], [], []
    for path in files:
        sp = re.sub(r"\.domtblout(\.gz)?$", "", os.path.basename(path))
        cov = domain_coverage(read_domtbl(path))
        rows = []
        if not cov.empty:
            has = cov.groupby("protein").pfam_id.apply(set).to_dict()
            for _, fam in dm.iterrows():
                sel = cov[(cov.pfam_id.isin(fam.pfam_ids)) &
                          (cov.hmm_cov >= fam.min_cov) & (cov.prot_len >= fam.min_len)]
                for prot, g in sel.groupby("protein"):
                    if fam.rule == "ALL" and not set(fam.pfam_ids).issubset(has.get(prot, set())):
                        continue
                    rows.append((sp, fam.category, fam.family, fam.tier, prot,
                                 int(g.prot_len.iloc[0]), round(float(g.hmm_cov.max()), 3), "GA"))
        else:
            log(f"WARNING: {sp} produced no Pfam hits")
        c = pd.DataFrame(rows, columns=["species", "category", "family", "tier",
                                        "protein", "prot_len", "hmm_cov", "evidence"])
        called = set(zip(c.family, c.protein))

        # rescue candidates from the below-threshold pass
        cands = []
        rpath = os.path.join(a.out, "rescue", f"{sp}.domtblout.gz")
        if len(rescue_fams) and os.path.exists(rpath):
            rd = read_domtbl(rpath)
            if not rd.empty:
                rd["i_evalue"] = pd.to_numeric(rd.i_evalue, errors="coerce")
                rcov = domain_coverage(rd)
                for _, fam in rescue_fams.iterrows():
                    sel = rcov[rcov.pfam_id.isin(fam.pfam_ids)]
                    for prot, g in sel.groupby("protein"):
                        if (fam.family, prot) in called:
                            continue
                        sub = rd[(rd.target == prot) & (rd.pfam_id.isin(fam.pfam_ids))]
                        cands.append({"fam": fam, "protein": prot,
                                      "prot_len": int(g.prot_len.iloc[0]),
                                      "hmm_cov": float(g.hmm_cov.max()),
                                      "best_evalue": float(sub.i_evalue.min()),
                                      "copies": _copies(sub)})

        # one pass over the proteome: size, short-protein capacity,
        # composition screen, and sequences for rescue candidates
        need = {x["protein"] for x in cands}
        seqs = {}
        faa = faa_of.get(sp)
        if faa:
            lens = []
            for h, seq in read_fasta(faa):
                pid = h.split()[0]
                L = len(seq.rstrip("*"))
                lens.append(L)
                if pid in need:
                    seqs[pid] = seq
                if has_mt and L <= 120 and _is_mt_like(seq):
                    f_ = mt_features(seq)
                    comp_rows.append((sp, pid, f_["length"], f_["n_cys"], f_["cys_fraction"],
                                      f_["n_motifs"], f_["n_aromatic"]))
            prot_totals[sp] = len(lens)
            lens.sort()
            stat_rows.append((sp, len(lens), sum(1 for x in lens if x <= 60),
                              sum(1 for x in lens if x <= 100),
                              lens[len(lens) // 2] if lens else 0,
                              lens[0] if lens else 0))

        for x in cands:
            fam = x["fam"]
            s_ = seqs.get(x["protein"])
            vals = {"len": x["prot_len"], "copies": x["copies"],
                    "cys": (s_.upper().count("C") / max(len(s_), 1)) if s_ else None}
            rules = _parse_validate(fam.validate, fam.family)
            ok, why = (x["hmm_cov"] >= fam.min_cov and x["prot_len"] >= fam.min_len,
                       "below map min_cov/min_len")
            if ok:
                ok, why = _check_rules(rules, vals)
            resc_rows.append((sp, fam.family, x["protein"], x["prot_len"],
                              round(x["hmm_cov"], 3), x["best_evalue"], x["copies"],
                              None if vals["cys"] is None else round(vals["cys"], 3),
                              "PASS" if ok else "FAIL", why))

        c.to_csv(os.path.join(calls_dir, f"{sp}.tsv"), sep="\t", index=False)
        all_calls.append(c)
        log(f"{sp}: {c.protein.nunique()} defensome proteins in {c.family.nunique()} families"
            + (f"; {sum(1 for r in resc_rows if r[0]==sp and r[8]=='PASS')} rescued"
               if len(rescue_fams) else ""))

    calls = pd.concat(all_calls, ignore_index=True)
    calls.to_csv(os.path.join(a.out, "gene_calls_all.tsv"), sep="\t", index=False)
    samples = [re.sub(r"\.domtblout(\.gz)?$", "", os.path.basename(p)) for p in files]
    counts = (calls.groupby(["species", "family"]).protein.nunique()
              .unstack(fill_value=0).reindex(index=samples, columns=dm.family, fill_value=0))
    counts.index.name = "species"
    counts.to_csv(os.path.join(a.out, "counts.tsv"), sep="\t")

    rc = pd.DataFrame(resc_rows, columns=["species", "family", "protein", "prot_len",
                                          "hmm_cov", "best_evalue", "copies",
                                          "cys_fraction", "validation", "reason"])
    rc.to_csv(os.path.join(a.out, "rescue_calls.tsv"), sep="\t", index=False)
    ok = rc[rc.validation == "PASS"]
    add = (ok.groupby(["species", "family"]).protein.nunique().unstack(fill_value=0)
             .reindex(index=counts.index, columns=counts.columns, fill_value=0))
    (counts + add).to_csv(os.path.join(a.out, "counts_rescued.tsv"), sep="\t")
    if len(rescue_fams):
        for f in rescue_fams.family:
            n_ga, n_r = int(counts[f].sum()), int(add[f].sum()) if f in add else 0
            log(f"{f}: {n_ga} at gathering threshold, {n_r} more rescued below it "
                f"({int((rc.family == f).sum())} candidates examined)")

    pd.DataFrame(comp_rows, columns=["species", "protein", "length", "n_cys", "cys_fraction",
                                     "n_cys_motifs", "n_aromatic"]) \
      .to_csv(os.path.join(a.out, "composition_candidates.tsv"), sep="\t", index=False)
    if stat_rows:
        pd.DataFrame(stat_rows, columns=["species", "n_proteins", "n_le60aa", "n_le100aa",
                                         "median_length", "min_length"]) \
          .to_csv(os.path.join(a.out, "proteome_stats.tsv"), sep="\t", index=False)

    if prot_totals:
        tot = pd.Series(prot_totals).reindex(counts.index)
        norm = counts.div(tot, axis=0) * 10000
        norm.round(3).to_csv(os.path.join(a.out, "counts_per10k.tsv"), sep="\t")
        pd.DataFrame({"species": tot.index, "n_proteins": tot.values}).to_csv(
            os.path.join(a.out, "proteome_sizes.tsv"), sep="\t", index=False)
        log("wrote counts.tsv, counts_rescued.tsv and counts_per10k.tsv")
    else:
        log("wrote counts.tsv and counts_rescued.tsv (no --proteomes, so no normalisation "
            "and rescue rules needing sequence cannot be checked)")
    log(f"annotate done: {counts.shape[0]} samples x {counts.shape[1]} families")


# --------------------------------------------------------------------------- #
def _dim(n, per, lo=6.0, hi=38.0):
    """Figure size in inches for n rows or columns, clamped.

    Matplotlib allocates the whole canvas up front, so an unclamped "inches
    per species" rule becomes a 630-inch figure at 1,575 species, roughly
    1.5 billion pixels, and the process is killed. Clamping keeps large runs
    readable as a block even when individual labels no longer fit."""
    return float(min(max(lo, per * max(n, 1)), hi))


def _tick_step(n, max_labels=90):
    """Label every k-th tick so a large axis stays legible."""
    return max(1, int(np.ceil(n / max_labels)))


def cmd_report(a):
    need_pandas()
    counts = pd.read_csv(os.path.join(a.out, "counts.tsv"), sep="\t", index_col=0)
    npath = os.path.join(a.out, "counts_per10k.tsv")
    norm = pd.read_csv(npath, sep="\t", index_col=0) if os.path.exists(npath) else None
    dm = read_map(a.map)
    core = list(dm.loc[dm.tier == "CORE", "family"])

    summary = pd.DataFrame({
        "family": counts.columns,
        "tier": dm.set_index("family").tier.reindex(counts.columns).values,
        "min": counts.min(), "median": counts.median(),
        "max": counts.max(), "mean": counts.mean().round(2),
        "cv": (counts.std() / counts.mean().replace(0, np.nan)).round(3),
        "n_species_absent": (counts == 0).sum(),
    }).reset_index(drop=True).sort_values("cv", ascending=False)
    summary.to_csv(os.path.join(a.out, "family_summary.tsv"), sep="\t", index=False)
    print("\n=== Most variable families (CORE tier) ===")
    print(summary[summary.tier == "CORE"].head(10).to_string(index=False))

    if a.metadata and a.group_by:
        md = pd.read_csv(a.metadata, sep="\t", comment="#")
        if a.group_by not in md.columns:
            sys.exit(f"ERROR: column '{a.group_by}' not in {a.metadata}.\n"
                     f"Available columns: {list(md.columns)}")
        # Find whichever column actually matches the proteome names. Metadata
        # often stores "Genus species" while the file is "Genus_species.faa".
        want, key, best = set(counts.index), None, 0
        for c in md.columns:
            v = md[c].astype(str).str.strip().str.replace(" ", "_", regex=False)
            n = len(want & set(v))
            if n > best:
                key, best = c, n
        if key is None or best == 0:
            sys.exit(f"ERROR: no column in {a.metadata} matches the proteome names.\n"
                     f"Expected values like: {sorted(want)[:3]}\n"
                     f"Columns checked: {list(md.columns)}")
        log(f"metadata key column: '{key}' ({best}/{len(want)} species matched)")
        md = md.assign(**{key: md[key].astype(str).str.strip()
                          .str.replace(" ", "_", regex=False)}).set_index(key)
        use = norm if norm is not None else counts
        g = md[a.group_by].reindex(use.index)
        # A key column can match every row while the grouping column is blank.
        # Saying "35/35 matched" and then "no metadata for 35 species" in the
        # same breath is not a diagnosis; say which of the two actually failed.
        g = g.replace(r"^\s*$", np.nan, regex=True)
        miss = g[g.isna()].index.tolist()
        if len(miss) == len(g):
            sys.exit(f"ERROR: every row matched on '{key}', but column "
                     f"'{a.group_by}' is empty for all {len(g)} of them.\n"
                     f"  {a.metadata} has the right species; it just has no values in\n"
                     f"  '{a.group_by}'. Fill that column in, then rerun. Nothing here can\n"
                     f"  guess how you want the species grouped.")
        if miss:
            log(f"WARNING: '{a.group_by}' is blank for {len(miss)} of {len(g)} species: "
                f"{miss[:5]}{' ...' if len(miss) > 5 else ''}")
            log("         They are dropped from the group means and tests.")
        gm = use.groupby(g).mean().round(2).T
        gm["n_groups"] = gm.notna().sum(axis=1)
        gm.to_csv(os.path.join(a.out, f"group_means_{a.group_by}.tsv"), sep="\t")
        print(f"\n=== Mean copies per 10k {size_unit(a.out)} by {a.group_by} (CORE) ===")
        print(gm.loc[[f for f in core if f in gm.index]].to_string())
        try:
            import warnings
            from scipy import stats
            rows = []
            for fam in use.columns:
                grp = [use.loc[g == lv, fam].values for lv in g.dropna().unique()]
                grp = [x for x in grp if len(x) > 1]
                if len(grp) >= 2 and any(len(set(x)) > 1 for x in grp):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        h, p = stats.kruskal(*grp)
                    if np.isfinite(p):
                        rows.append((fam, h, p))
            kw = pd.DataFrame(rows, columns=["family", "H", "p"]).sort_values("p")
            kw["p_bh"] = kw.p * len(kw) / kw.p.rank()   # Benjamini-Hochberg
            kw["p_bh"] = kw.p_bh[::-1].cummin()[::-1].clip(upper=1)
            kw.to_csv(os.path.join(a.out, f"kruskal_{a.group_by}.tsv"), sep="\t", index=False)
            if kw.empty:
                print(f"\nNo group test run: '{a.group_by}' needs at least two "
                      "groups with more than one species each.")
            else:
                print(f"\nKruskal-Wallis by {a.group_by}, top 5 (see kruskal_*.tsv):")
                print(kw.head(5).to_string(index=False))
            print("\nNOTE: species are not independent. These tests ignore phylogeny\n"
                  "      and are a screen, not a result.")
        except ImportError:
            log("scipy not available, skipping group tests")

    meta_df = None
    if a.metadata and a.group_by:
        meta_df = md
    make_figures(a.out, counts, norm, dm, meta_df, a.group_by)
    log("report done")


# =========================================================================== #
VERDICTS = {
    "ACCESSION_NOT_IN_DATABASE":
        "the HMM is not in the database that was searched; check the accession "
        "and the Pfam release in MANIFEST.tsv",
    "FILTERED_BY_MAP":
        "hits passed Pfam's threshold but every one failed this family's min_cov, "
        "min_len or ALL rule; inspect them before loosening the map",
    "FOUND_BELOW_GA":
        "no hits at Pfam's gathering threshold, but validated hits below it; "
        "reported in counts_rescued.tsv and rescue_calls.tsv, not in counts.tsv",
    "BELOW_GA_FAILED_VALIDATION":
        "below-threshold hits exist but none passed the map's validate rules; "
        "read rescue_calls.tsv to see why",
    "RECOVERED_FROM_TRANSCRIPTS":
        "absent from the proteomes, but MT-like ORFs were recovered from the transcripts by "
        "`short-orfs`; see short_orf_mt.tsv. The proteomes lost them to the ORF caller's length cut-off",
    "INPUT_LACKS_SHORT_PROTEINS":
        "a short family, and most proteomes contain almost no proteins of that "
        "length; the gene predictor's minimum length has probably removed them "
        "(TransDecoder's default is 100 aa), so absence cannot be concluded",
    "COMPOSITION_CANDIDATES_ONLY":
        "no HMM evidence, but proteins with metallothionein-like composition; "
        "see composition_candidates.tsv; heuristic, verify by alignment",
    "NOT_DETECTED":
        "searched with the right HMM, below threshold too, nothing found; "
        "consistent with genuine absence or extreme divergence",
    "NOT_DETECTED_AT_GA_ONLY":
        "nothing at Pfam's threshold, and the below-threshold rescue was not run "
        "for this family (rescue=no in the map, or scan --no-rescue); absence "
        "has not been tested as thoroughly as it could be",
}


def _manifest_for(out):
    p = os.path.join(out, "scan_db.txt")
    if not os.path.exists(p):
        return None
    db = open(p).read().strip()
    m = os.path.join(os.path.dirname(db), "MANIFEST.tsv")
    return pd.read_csv(m, sep="\t") if os.path.exists(m) else None


def diagnose_zero_families(a, c):
    """Turn every family that came back empty into a labelled verdict.

    A zero can mean the HMM was never searched, that a map rule filtered
    everything, that real copies score just below Pfam's threshold, that the
    input proteomes cannot contain proteins that short, or genuine absence.
    Those need different responses, and a count table cannot tell them apart."""
    dmap = read_map(a.map).set_index("family")
    dead = [f for f in c.columns if f in dmap.index and c[f].sum() == 0]
    if not dead:
        return
    man = _manifest_for(a.out)
    in_db = set(man[man.found == "yes"].pfam_id) if man is not None else None
    leng = dict(zip(man.pfam_id, man.length)) if man is not None else {}

    def raw_hits(sub, accs):
        n = {x: 0 for x in accs}
        for p in glob.glob(os.path.join(a.out, sub, "*.domtblout*")):
            with _open_text(p) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split(None, 6)
                    if len(f) > 4 and f[4].split(".")[0] in n:
                        n[f[4].split(".")[0]] += 1
        return n
    accs = {x for f in dead for x in dmap.loc[f, "pfam_ids"]}
    ga = raw_hits("hmmsearch", accs)
    rs = raw_hits("rescue", accs)
    rescue_ran = bool(glob.glob(os.path.join(a.out, "rescue", "*.domtblout*")))
    rc = _read_tsv(os.path.join(a.out, "rescue_calls.tsv"))
    cc = _read_tsv(os.path.join(a.out, "composition_candidates.tsv"))
    ps = _read_tsv(os.path.join(a.out, "proteome_stats.tsv"))
    so = _read_tsv(os.path.join(a.out, "short_orf_summary.tsv"))
    n_so_total = int(so.mt_like_genes.sum()) if so is not None and len(so) else 0

    rows = []
    for f in dead:
        ids = dmap.loc[f, "pfam_ids"]
        n_ga = sum(ga[x] for x in ids)
        n_rs = sum(rs[x] for x in ids)
        n_pass = int(((rc.family == f) & (rc.validation == "PASS")).sum()) if rc is not None and len(rc) else 0
        n_comp = len(cc) if (cc is not None and f == "Metallothionein") else 0
        L = min([int(leng[x]) for x in ids if x in leng and str(leng[x]).isdigit()] or [999])
        short = L < 60 or int(dmap.loc[f, "min_len"]) < 60
        lacking = None
        if ps is not None and len(ps):
            lacking = int((ps.n_le60aa <= 5).sum())
        missing = [x for x in ids if in_db is not None and x not in in_db]
        # Evidence outranks a missing sibling accession: a family with two
        # accessions, one absent from the release and one producing rescued
        # hits, has been detected, and saying otherwise sends the user to fix
        # the database when the interesting result is the hits.
        if n_ga > 0:
            v = "FILTERED_BY_MAP"
        elif n_pass > 0:
            v = "FOUND_BELOW_GA"
        elif f == "Metallothionein" and n_so_total > 0:
            v = "RECOVERED_FROM_TRANSCRIPTS"
        elif n_rs > 0:
            v = "BELOW_GA_FAILED_VALIDATION"
        elif missing and len(missing) == len(ids):
            v = "ACCESSION_NOT_IN_DATABASE"
        elif short and lacking is not None and lacking >= max(1, len(ps) // 2):
            v = "INPUT_LACKS_SHORT_PROTEINS"
        elif n_comp > 0:
            v = "COMPOSITION_CANDIDATES_ONLY"
        elif not (rescue_ran and dmap.loc[f, "rescue"] == "yes"):
            v = "NOT_DETECTED_AT_GA_ONLY"
        else:
            v = "NOT_DETECTED"
        note = VERDICTS[v]
        if missing and v != "ACCESSION_NOT_IN_DATABASE":
            note += (f"; note {','.join(missing)} is not in the database, so this family "
                     "was only partly searched")
        rows.append({"family": f, "pfam_ids": ",".join(ids), "verdict": v,
                     "accessions_missing": ",".join(missing),
                     "in_database": "unknown" if in_db is None else ("no" if missing else "yes"),
                     "hmm_length": "" if L == 999 else L,
                     "raw_hits_at_GA": n_ga, "raw_hits_below_GA": n_rs,
                     "rescued_validated": n_pass, "composition_candidates": n_comp,
                     "short_orf_candidates": n_so_total if f == "Metallothionein" else 0,
                     "proteomes_with_<=5_proteins_under_60aa":
                         "" if lacking is None else f"{lacking}/{len(ps)}",
                     "what_it_means": note})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "qc_zero_families.tsv"), sep="\t", index=False)
    print("\n=== Families with no calls at Pfam's gathering threshold ===")
    for r in rows:
        print(f"  {r['family']:18s} {r['pfam_ids']:10s} {r['verdict']}")
        print(f"  {'':18s} {r['what_it_means']}")
    if man is None:
        print("  (no MANIFEST.tsv found next to the scanned database; run `setup` so "
              "missing accessions can be ruled out)")




def cmd_qc(a):
    """Annotation-quality report from low-copy control families.

    Two separate questions, reported separately:

    (1) Is a FAMILY behaving? Compares the median count across species against
        the literature expectation. If catalase has a median of 3 in every
        species, the map rule is catching something extra; that is a problem
        with the map, not with any one genome.

    (2) Is a SPECIES behaving? Compares each species against the observed
        median for that family, not against my prior. A genome with seven times
        the median catalase count has duplicated gene models. This is
        self-calibrating, so it stays valid on taxa where my priors do not.
    """
    need_pandas()
    check_fresh(a.out, read_map(a.map), a.map)
    c = pd.read_csv(os.path.join(a.out, "counts.tsv"), sep="\t", index_col=0)
    diagnose_zero_families(a, c)
    sizes_f = os.path.join(a.out, "proteome_sizes.tsv")
    sizes = (pd.read_csv(sizes_f, sep="\t").set_index("species").n_proteins
             if os.path.exists(sizes_f) else None)
    # literature expectation for insect genomes, used only for question (1)
    # Literature expectation for insect genomes. Where the observed median
    # stays high after tightening the thresholds, the domain is simply shared
    # by a wider set of proteins than the family name implies, and the family
    # needs phylogeny to subdivide rather than a stricter cutoff.
    expect = {"GCL": 1, "Catalase": 1, "Hsp90": 2, "SOD_Fe": 1, "SOD_CuZn": 3,
              "NR_xeno": 21, "GST_MAPEG": 4, "BTB_BACK_Kelch": 7, "Keap1": 1}
    expect_src = "insect genomes (built-in)"
    if getattr(a, "expect", None):
        if not os.path.exists(a.expect):
            sys.exit(f"ERROR: --expect file not found: {a.expect}")
        user = {}
        with open(a.expect) as fh:
            for line in fh:
                line = line.split("#")[0].strip()
                if not line:
                    continue
                f_, _, v = line.partition("\t")
                if not v.strip():
                    f_, _, v = line.partition(None)
                try:
                    user[f_.strip()] = float(v)
                except ValueError:
                    sys.exit(f"ERROR: cannot read '{line}' in {a.expect}; "
                             "use 'family<TAB>expected_count' per line")
        if not user:
            sys.exit(f"ERROR: {a.expect} has no usable rows")
        expect = user
        expect_src = a.expect
    ctrl = [k for k in expect if k in c.columns]
    if not ctrl:
        sys.exit("ERROR: no low-copy control families found in counts.tsv")

    fam = pd.DataFrame({
        "expected": pd.Series({k: expect[k] for k in ctrl}),
        "observed_median": c[ctrl].median(),
        "observed_min": c[ctrl].min(), "observed_max": c[ctrl].max(),
    })
    fam["median_vs_expected"] = (fam.observed_median / fam.expected).round(1)
    # A median of zero has several possible causes (never searched, input too
    # short, genuinely absent) and the zero-family diagnosis above names the
    # right one. Calling every zero MAP_RULE_TOO_STRICT sent users to loosen a
    # rule for families whose HMM had never been in the database at all.
    zdiag = _read_tsv(os.path.join(a.out, "qc_zero_families.tsv"))
    zv = dict(zip(zdiag.family, zdiag.verdict)) if zdiag is not None and len(zdiag) else {}
    def _verdict(f, r):
        if r.observed_median == 0:
            return ("NOT_SEARCHED" if zv.get(f) == "ACCESSION_NOT_IN_DATABASE"
                    else "ZERO_SEE_DIAGNOSIS")
        if r.median_vs_expected >= 2:
            return "MAP_RULE_TOO_LOOSE"
        if r.median_vs_expected <= 0.5:
            return "MAP_RULE_TOO_STRICT"
        return "ok"
    fam["verdict"] = [_verdict(f, r) for f, r in fam.iterrows()]
    fam.sort_values("median_vs_expected", ascending=False) \
       .to_csv(os.path.join(a.out, "qc_families.tsv"), sep="\t")
    print(f"\n=== (1) Are the control FAMILIES behaving? ===")
    print(f"expectations from: {expect_src}")
    if expect_src.endswith("(built-in)"):
        print("These counts are priors for INSECT genomes. Vertebrates legitimately carry")
        print("many more nuclear receptors and BTB-Kelch proteins, so on a non-insect or")
        print("mixed-lineage set a 'MAP_RULE_TOO_LOOSE' verdict here may just mean the")
        print("prior is for the wrong clade. Supply your own with --expect family<TAB>count.")
        print("Question (2) below is self-calibrating and unaffected.")
    print(fam.sort_values("median_vs_expected", ascending=False).to_string())
    rule = fam[fam.verdict.isin(["MAP_RULE_TOO_LOOSE", "MAP_RULE_TOO_STRICT"])]
    if len(rule):
        print(f"\n-> {list(rule.index)}: counts are present but off the expectation, so the")
        print("   map rule is the likely cause. Adjust min_cov/min_len or require a partner domain.")
    unsearched = fam[fam.verdict == "NOT_SEARCHED"]
    if len(unsearched):
        print(f"\n-> {list(unsearched.index)}: never searched; the HMM is not in the database.")
        print("   Fix the database (run `setup` against a full Pfam-A), not the map.")
    zero = fam[fam.verdict == "ZERO_SEE_DIAGNOSIS"]
    if len(zero):
        print(f"\n-> {list(zero.index)}: zero everywhere; the zero-family table above says why.")

    # question (2): species outliers against the observed median
    base = c[ctrl].median().replace(0, np.nan)
    usable = [k for k in ctrl if pd.notna(base[k])]
    if not usable:
        print("\n=== (2) Species outliers: not assessable ===")
        print("Every control family has a median of zero across samples, so no")
        print("species can be compared against it. Either the controls are genuinely")
        print("absent from these proteomes or the map rules for them are too strict;")
        print("the table under (1) shows which.")
        log("qc done")
        return
    ratio = c[usable].div(base[usable], axis=1)
    # pandas >= 3 raises on idxmax over an all-NA row; compute it only where
    # at least one control family is informative for that sample
    has = ratio.notna().any(axis=1)
    worst = pd.Series("", index=ratio.index, dtype=object)
    if has.any():
        worst[has] = ratio[has].idxmax(axis=1)
    qc = pd.DataFrame({
        "n_proteins": sizes.reindex(c.index) if sizes is not None else np.nan,
        "max_ratio": ratio.max(axis=1).round(1),
        "worst_family": worst,
        "n_families_inflated": (ratio >= 2).sum(axis=1),
        "n_families_missing": (c[ctrl] == 0).sum(axis=1),
    })
    for k in ctrl:
        qc[f"n_{k}"] = c[k]
    qc["flag"] = np.where((qc.max_ratio >= 3) | (qc.n_families_inflated >= 2), "INFLATED",
                  np.where(qc.n_families_missing >= 2, "SPARSE", "ok"))
    qc = qc.sort_values("max_ratio", ascending=False)
    qc.to_csv(os.path.join(a.out, "qc_annotation.tsv"), sep="\t")
    print("\n=== (2) Are any SPECIES outliers? (ratio to observed median) ===")
    print(qc.head(10).to_string())
    n_bad = (qc.flag != "ok").sum()
    print(f"\n{n_bad} of {len(qc)} species flagged. Use max_ratio as an "
          "annotation-quality covariate,\nor exclude the flagged species. "
          "This replaces a BUSCO duplication run for most purposes.")
    log("qc done")


# =========================================================================== #
# =========================================================================== #
#  Figure style, shared by every figure the tool writes
# =========================================================================== #
# Same colours as the dashboard, so a method means one colour everywhere.
METHOD_COLOURS = {"DToL": "#1d3557", "genome_guided": "#457b9d", "denovo": "#e63946"}
_METHOD_RANK = {"DToL": 0, "genome": 0, "genome_guided": 1, "denovo": 2}


def method_order(methods):
    """Genome annotation first, then genome-guided, then de novo: the order of
    decreasing reliance on a genome, so lines between them read as a gradient
    rather than zig-zagging alphabetically."""
    return sorted(methods, key=lambda m: (_METHOD_RANK.get(m, 9), str(m)))


def method_colour(m, i=0):
    fallback = ["#6a994e", "#bc6c25", "#7b2cbf", "#8d99ae"]
    return METHOD_COLOURS.get(m, fallback[i % len(fallback)])


def _pub_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": 9,
        "axes.titlesize": 11, "axes.titleweight": "bold", "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e9ecef", "grid.linewidth": .6,
        "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        # editable text in Illustrator/Inkscape rather than outlined glyphs
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.bbox": "tight", "savefig.dpi": 300,
    })


class FigureWriter:
    """Every figure as PDF (vector, for the thesis) and PNG (300 dpi, for
    slides), plus one combined PDF and a captions file so each figure says
    what it shows without the reader having to open this code."""
    def __init__(self, fig_dir, summary_name):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        self.plt, self.dir = plt, fig_dir
        os.makedirs(fig_dir, exist_ok=True)
        self.pages = PdfPages(os.path.join(fig_dir, summary_name))
        self.captions = []

    def save(self, fig, name, caption=""):
        base = os.path.join(self.dir, name)
        fig.savefig(base + ".pdf")
        fig.savefig(base + ".png")
        self.pages.savefig(fig)
        self.plt.close(fig)
        if caption:
            self.captions.append((name, caption))

    def close(self, caption_file="FIGURES.md"):
        self.pages.close()
        with open(os.path.join(self.dir, caption_file), "a") as fh:
            for n, c in self.captions:
                fh.write(f"**{n}**  {c}\n\n")



def _split_families(dm, mat):
    """CORE families that can be z-scored, and those that cannot.

    An all-zero or constant family has no z-score. Dropping it silently is
    wrong: 'GCL is single copy in every species' and 'MATE was never detected'
    are both findings and both belong in the output.
    """
    core = [f for f in dm.loc[dm.tier == "CORE", "family"] if f in mat.columns]
    var, absent, const = [], [], []
    for f in core:
        if mat[f].sum() == 0:
            absent.append(f)
        elif mat[f].std() == 0:
            const.append(f)
        else:
            var.append(f)
    return core, var, absent, const


def make_figures(out, counts, norm, dm, meta=None, group=None):
    np.seterr(divide="ignore", invalid="ignore")
    try:
        _pub_style()
    except Exception:
        pass
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        log("matplotlib not available, skipping figures (pip install --user matplotlib)")
        return
    fig_dir = os.path.join(out, "figures"); os.makedirs(fig_dir, exist_ok=True)
    use = norm if norm is not None else counts
    unit = size_unit(out)
    lab = f"copies per 10k {unit}" if norm is not None else "copies"
    cats = dm.set_index("family").category

    # every count figure is produced twice: normalised and raw
    variants = [("", use, lab)]
    if norm is not None:
        variants.append(("_raw", counts, "raw copies"))

    core, var, _, _ = _split_families(dm, use)
    # "absent" and "constant" are statements about the genomes, so judge them on
    # RAW counts. GCL is exactly 1 in every species but varies after dividing by
    # proteome size, and that single-copy result is the calibration finding.
    _, _, absent, const = _split_families(dm, counts)
    if absent or const:
        log(f"invariant families held out of z-scored panels: "
            f"absent={absent} constant={const} (see 00_invariant_families.pdf)")

    lv, cmap = [], {}
    if meta is not None and group is not None and group in meta.columns:
        lv = sorted(meta[group].dropna().unique())
        cmap = dict(zip(lv, plt.cm.viridis(np.linspace(0, .85, len(lv)))))

    if os.path.exists(os.path.join(fig_dir, "FIGURES.md")):
        os.remove(os.path.join(fig_dir, "FIGURES.md"))
    W = FigureWriter(fig_dir, "report_all_figures.pdf")
    _caps = {'00_invariant': 'Families that cannot be z-scored because they are absent or constant in every sample, with their raw counts; judged on raw counts.', '01_heatmap': 'Samples by CORE family, colour is the z-score of each family across samples. Colour is relative to the other samples in this run, not an absolute amount.', '02_total': 'Total CORE defensome per sample.', '03_pca': 'PCA of z-scored CORE families: sample scores (left) and family loadings (right). Read the loadings to see which families drive each axis.', '04_category': 'Share of the CORE defensome in each functional tier, per sample.', '05_family': 'Pearson correlation between families across samples: blocks of red are families that expand together.', '06_normalisation': 'Evidence that dividing by proteome size removed the size confound: correlation with proteome size before (raw) and after (per 10k).', '07_family': 'Every CORE family by group, boxes with points, n per group in the axis labels.', '08_domain': 'Share of called proteins whose required domains are complete, truncated (PARTIAL) or unusually short (FRAGMENT), per family.'}

    def save(fig, name):
        fig.tight_layout()
        stem = name[:-4] if name.endswith(".pdf") else name
        cap = next((v for k, v in _caps.items() if stem.startswith(k)), "")
        if "_raw" in stem:
            cap = (cap + " Raw counts.").strip()
        W.save(fig, stem, cap)

    # ---- 00. the families a heatmap cannot show ----------------------------
    if absent or const:
        show = absent + const
        fig, axes = plt.subplots(1, len(show), figsize=(3.6*len(show), 4.6),
                                 squeeze=False)
        for ax, f in zip(axes[0], show):
            vals = counts[f]
            ax.bar(range(len(vals)), vals.values,
                   color="#adb5bd" if f in absent else "#2a6f97")
            ax.set_title(f, fontsize=11)
            ax.set_xticks([]); ax.set_xlabel(f"{len(vals)} species", fontsize=8)
            ax.set_ylabel("raw copies", fontsize=9)
            if f in absent:
                ax.set_ylim(0, 1)
                ax.text(.5, .55, "never detected\nin any species", ha="center",
                        va="center", transform=ax.transAxes, fontsize=9, color="#c1121f")
                ax.text(.5, .3, "run `qc` to see whether the\nHMM produced raw hits at all",
                        ha="center", va="center", transform=ax.transAxes, fontsize=7.5,
                        color="#666")
            else:
                ax.set_ylim(0, max(2, vals.max()*1.6))
                ax.text(.5, .8, f"constant at {int(vals.iloc[0])}\nin all species",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=9, color="#2a6f97")
                ax.text(.5, .6, "single-copy control:\nthe pipeline is calibrated",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=7.5, color="#666")
        fig.suptitle("Families excluded from the z-scored panels, and why", y=1.02)
        save(fig, "00_invariant_families.pdf")

    # ---- 01. heatmaps, normalised and raw ----------------------------------
    for suffix, mat_, lab_ in variants:
        _, v_, _, _ = _split_families(dm, mat_)
        mm = mat_[v_]
        zz = ((mm - mm.mean()) / mm.std().replace(0, np.nan)).fillna(0)
        os_ = list(zz.index[np.argsort(zz.values @ zz.values.sum(0))])
        of_ = list(zz.columns[np.argsort(zz.values.sum(0))])
        fig, ax = plt.subplots(figsize=(_dim(len(v_), .40, 9), _dim(len(mm), .30, 7)))
        im = ax.imshow(zz.loc[os_, of_].values, aspect="auto", cmap="RdBu_r",
                       vmin=-3, vmax=3)
        ax.set_xticks(range(len(of_))); ax.set_xticklabels(of_, rotation=90, fontsize=8)
        ax.set_yticks(range(len(os_)))
        ax.set_yticklabels([s.replace("_", " ") for s in os_], fontsize=7, style="italic")
        if lv:
            for t_, s_ in zip(ax.get_yticklabels(), os_):
                t_.set_color(cmap.get(meta[group].get(s_), "#222"))
        fig.colorbar(im, ax=ax, label="z-score", shrink=.6)
        note = ""
        if absent or const:
            note = "\nnot shown: " + ", ".join(absent + const) + " (see 00_invariant_families.pdf)"
        ax.set_title(f"Defensome, CORE families ({lab_}, z-scored per family){note}",
                     fontsize=11)
        save(fig, f"01_heatmap_core{suffix}.pdf")

    # ---- 02. totals, normalised and raw ------------------------------------
    for suffix, mat_, lab_ in variants:
        t_ = mat_[var].sum(axis=1).sort_values()
        cc = [cmap.get(meta[group].get(s_), "#999999") for s_ in t_.index] if lv else "#2a6f97"
        fig, ax = plt.subplots(figsize=(7.5, _dim(len(t_), .26, 6)))
        ax.barh(range(len(t_)), t_.values, color=cc)
        ax.set_yticks(range(len(t_)))
        ax.set_yticklabels([s.replace("_", " ") for s in t_.index], fontsize=7.5,
                           style="italic")
        ax.set_xlabel(f"total CORE defensome, {lab_}")
        if lv:
            ax.legend(handles=[Patch(color=cmap[k], label=str(k)) for k in lv],
                      fontsize=8, title=group)
        save(fig, f"02_total_defensome{suffix}.pdf")

    # ---- 03. PCA -----------------------------------------------------------
    # needs at least two varying families and three samples, or there is no
    # second axis to plot; small or narrow runs skip it rather than crash
    if len(var) >= 2 and len(use) >= 3:
        m = use[var]
        z = ((m - m.mean()) / m.std().replace(0, np.nan)).fillna(0)
        x = z.values - z.values.mean(0)
        u, s, vt = np.linalg.svd(x, full_matrices=False)
        pc = u[:, :2] * s[:2]; ev = (s**2 / (s**2).sum() * 100)[:2]
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6.5))
        if lv:
            for k in lv:
                i = [j for j, sp in enumerate(z.index) if meta[group].get(sp) == k]
                ax.scatter(pc[i, 0], pc[i, 1], s=42, color=cmap[k], label=str(k),
                           edgecolor="white", linewidth=.6)
            ax.legend(fontsize=9, title=group)
        else:
            ax.scatter(pc[:, 0], pc[:, 1], s=34, color="#2a6f97")
        for i, sp in enumerate(z.index):
            ax.annotate(sp.replace("_", " "), (pc[i, 0], pc[i, 1]), fontsize=5.5,
                        xytext=(4, 2), textcoords="offset points", style="italic")
        ax.set_xlabel(f"PC1 ({ev[0]:.1f}% of variance)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
        ax.axhline(0, lw=.4, c="grey"); ax.axvline(0, lw=.4, c="grey")
        ax.set_title(f"Species scores ({lab}, {len(var)} CORE families z-scored)")
        ld = vt[:2].T * s[:2] / np.sqrt(len(z))
        ax2.axhline(0, lw=.5, c="grey"); ax2.axvline(0, lw=.5, c="grey")
        for i, f in enumerate(z.columns):
            ax2.annotate("", xy=(ld[i, 0], ld[i, 1]), xytext=(0, 0),
                         arrowprops=dict(arrowstyle="->", color="#c1121f", lw=.7))
            ax2.annotate(f, (ld[i, 0], ld[i, 1]), fontsize=7)
        ax2.set_xlabel("PC1 loading"); ax2.set_ylabel("PC2 loading")
        ax2.set_title("Family loadings: which families drive the separation")
        save(fig, "03_pca.pdf")

    else:
        log(f"PCA skipped: needs >=2 varying CORE families and >=3 samples "
            f"(have {len(var)} and {len(use)})")

    # ---- 04. category composition, normalised and raw ----------------------
    for suffix, mat_, lab_ in variants:
        cs = mat_[var].T.groupby(cats.reindex(var).values).sum().T
        cs = cs.loc[cs.sum(axis=1).sort_values().index]
        frac = cs.div(cs.sum(axis=1), axis=0)
        fig, ax = plt.subplots(figsize=(9.5, _dim(len(frac), .26, 6)))
        left = np.zeros(len(frac))
        for c_ in frac.columns:
            ax.barh(range(len(frac)), frac[c_].values, left=left, label=c_)
            left += frac[c_].values
        ax.set_yticks(range(len(frac)))
        ax.set_yticklabels([s.replace("_", " ") for s in frac.index], fontsize=7.5,
                           style="italic")
        ax.set_xlabel(f"fraction of CORE defensome ({lab_})")
        ax.legend(fontsize=8, ncol=3, loc="lower right")
        save(fig, f"04_category_composition{suffix}.pdf")

    # ---- 05. family co-variation ------------------------------------------
    if len(var) >= 2 and len(use) >= 3:
        with np.errstate(divide="ignore", invalid="ignore"):
            cm = np.nan_to_num(np.corrcoef(z.values.T))
        o = np.argsort(cm.sum(0))
        fig, ax = plt.subplots(figsize=(_dim(len(var), .36, 8), _dim(len(var), .36, 7)))
        im = ax.imshow(cm[np.ix_(o, o)], cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(var))); ax.set_xticklabels([var[i] for i in o], rotation=90, fontsize=8)
        ax.set_yticks(range(len(var))); ax.set_yticklabels([var[i] for i in o], fontsize=8)
        fig.colorbar(im, ax=ax, label="Pearson r", shrink=.7)
        ax.set_title(f"Do defensome families expand together? ({lab})")
        save(fig, "05_family_correlation.pdf")


    # ---- 06. normalisation QC ----------------------------------------------
    sizes_f = os.path.join(out, "proteome_sizes.tsv")
    if os.path.exists(sizes_f) and norm is not None:
        sz = pd.read_csv(sizes_f, sep="\t").set_index("species").n_proteins.reindex(counts.index)
        def _r(xx, yy):
            return 0.0 if np.std(xx) == 0 or np.std(yy) == 0 else float(np.corrcoef(xx, yy)[0, 1])
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
        axes[0].scatter(sz, counts[var].sum(axis=1), s=30, c="#c1121f")
        axes[0].set_xlabel(f"{unit} in proteome"); axes[0].set_ylabel("raw CORE defensome")
        axes[0].set_title(f"raw counts track proteome size (r={_r(counts[var].sum(axis=1), sz):.2f})")
        axes[1].scatter(sz, norm[var].sum(axis=1), s=30, c="#2a6f97")
        axes[1].set_xlabel(f"{unit} in proteome"); axes[1].set_ylabel(f"CORE per 10k {unit}")
        axes[1].set_title(f"normalised counts should not (r={_r(norm[var].sum(axis=1), sz):.2f})")
        rr = sorted([(f, _r(counts[f], sz), _r(norm[f], sz)) for f in var], key=lambda t: -t[1])
        axes[2].barh(range(len(rr)), [t[1] for t in rr], .4, label="raw", color="#c1121f")
        axes[2].barh([i+.4 for i in range(len(rr))], [t[2] for t in rr], .4,
                     label="per 10k", color="#2a6f97")
        axes[2].set_yticks([i+.2 for i in range(len(rr))])
        axes[2].set_yticklabels([t[0] for t in rr], fontsize=7)
        axes[2].axvline(0, c="k", lw=.6); axes[2].set_xlabel(f"r with number of {unit}")
        axes[2].legend(fontsize=8)
        save(fig, "06_normalisation_qc.pdf")

    # ---- 07. every family by group, paginated, full labels -----------------
    if lv:
        g = meta[group].reindex(use.index)
        counts_per_group = {k: int((g == k).sum()) for k in lv}
        for suffix, mat_, lab_ in variants:
            sel = core                      # ALL core families, including invariant
            per_page, nc = 20, 4
            pages = [sel[i:i+per_page] for i in range(0, len(sel), per_page)]
            for pi, page in enumerate(pages, 1):
                nr = int(np.ceil(len(page)/nc))
                fig, axes = plt.subplots(nr, nc, figsize=(4.3*nc, 3.3*nr), squeeze=False)
                for ax_, f in zip(np.ravel(axes), page):
                    data = [mat_.loc[g == k, f].values for k in lv]
                    bp = ax_.boxplot(data, positions=range(1, len(lv)+1),
                                     widths=.55, showfliers=False, patch_artist=True)
                    for patch, k in zip(bp["boxes"], lv):
                        patch.set_facecolor(cmap[k]); patch.set_alpha(.35)
                    for i, dv in enumerate(data, 1):
                        ax_.scatter(np.random.normal(i, .07, len(dv)), dv, s=11,
                                    color=cmap[lv[i-1]], alpha=.85, zorder=3,
                                    edgecolor="white", linewidth=.3)
                    ax_.set_xticks(range(1, len(lv)+1))
                    ax_.set_xticklabels([f"{k}\n(n={counts_per_group[k]})" for k in lv],
                                        fontsize=7.5)
                    ax_.set_title(f, fontsize=10)
                    ax_.tick_params(axis="y", labelsize=7.5)
                    if mat_[f].sum() == 0:
                        ax_.text(.5, .5, "never detected", ha="center", va="center",
                                 transform=ax_.transAxes, fontsize=9, color="#c1121f")
                    elif mat_[f].std() == 0:
                        ax_.text(.5, .85, "constant", ha="center", va="center",
                                 transform=ax_.transAxes, fontsize=8, color="#2a6f97")
                for ax_ in np.ravel(axes)[len(page):]:
                    ax_.axis("off")
                fig.suptitle(f"CORE families by {group} ({lab_})  "
                             f"page {pi} of {len(pages)}", y=1.005, fontsize=13)
                save(fig, f"07_family_by_group{suffix}_p{pi}.pdf")

    # ---- 08. domain completeness -------------------------------------------
    comp_f = os.path.join(out, "domains", "completeness.tsv")
    if os.path.exists(comp_f):
        comp = pd.read_csv(comp_f, sep="\t")
        t = comp.groupby(["family", "status"]).size().unstack(fill_value=0)
        for s_ in ("COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"):
            if s_ not in t:
                t[s_] = 0
        n_calls = t.sum(axis=1)
        t = t[["COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"]].div(n_calls, axis=0)
        t = t.loc[t.COMPLETE.sort_values().index]
        fig, ax = plt.subplots(figsize=(10, _dim(len(t), .28, 6)))
        left = np.zeros(len(t))
        for c_, col in zip(t.columns, ["#2a6f97", "#f4a261", "#e76f51", "#adb5bd"]):
            ax.barh(range(len(t)), t[c_].values, left=left, label=c_, color=col)
            left += t[c_].values
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels([f"{f}  (n={int(n_calls[f])})" for f in t.index], fontsize=8)
        ax.set_xlabel("fraction of called proteins"); ax.set_xlim(0, 1)
        ax.legend(fontsize=8, ncol=4, loc="lower left")
        ax.set_title("Domain architecture completeness per family")
        save(fig, "08_domain_completeness.pdf")


    # ---- R1. what the below-threshold rescue added -----------------------
    rc_f = os.path.join(out, "rescue_calls.tsv")
    cr_f = os.path.join(out, "counts_rescued.tsv")
    if os.path.exists(rc_f) and os.path.exists(cr_f):
        rc = pd.read_csv(rc_f, sep="\t")
        if len(rc):
            cr = pd.read_csv(cr_f, sep="\t", index_col=0)
            fams_r = sorted(rc.family.unique())
            fig, axes = plt.subplots(1, len(fams_r), figsize=(4.4 * len(fams_r), _dim(len(counts), .22, 3.2)),
                                     squeeze=False)
            for ax, f in zip(axes[0], fams_r):
                ga = counts[f] if f in counts else pd.Series(0, index=counts.index)
                tot_r = cr[f].reindex(counts.index).fillna(0) if f in cr else ga
                order = tot_r.sort_values().index
                y = np.arange(len(order))
                ax.barh(y, ga.reindex(order).values, color="#1d3557", label="at gathering threshold")
                ax.barh(y, (tot_r - ga).reindex(order).values, left=ga.reindex(order).values,
                        color="#f4a261", label="rescued below it")
                ax.set_yticks(y); ax.set_yticklabels([s_.replace("_", " ") for s_ in order], fontsize=7)
                sub = rc[rc.family == f]
                ax.set_title(f"{f}\n{int((sub.validation=='PASS').sum())} rescued, "
                             f"{int((sub.validation=='FAIL').sum())} rejected", loc="left", fontsize=10)
                ax.set_xlabel("genes"); ax.grid(axis="y", visible=False)
            from matplotlib.patches import Patch
            fig.legend(handles=[Patch(color="#1d3557", label="at Pfam's gathering threshold"),
                                Patch(color="#f4a261", label="rescued below it, validated")],
                       loc="upper center", ncol=2, bbox_to_anchor=(.5, 1.06), fontsize=8.5)
            save(fig, "R1_rescued_families.pdf")
            W.captions[-1] = ("R1_rescued_families",
                "Genes called at Pfam's gathering threshold (dark) and the additional genes found "
                "by the below-threshold rescue that passed the map's validation rules (orange). "
                "Rescued genes are reported in counts_rescued.tsv only, never in counts.tsv.")

    # ---- D1. why each empty family is empty --------------------------------
    z_f = os.path.join(out, "qc_zero_families.tsv")
    if os.path.exists(z_f):
        z0 = pd.read_csv(z_f, sep="\t")
        if len(z0):
            vcol = {"ACCESSION_NOT_IN_DATABASE": "#6a040f", "FILTERED_BY_MAP": "#e9a13b",
                    "FOUND_BELOW_GA": "#2a9d8f", "BELOW_GA_FAILED_VALIDATION": "#f4a261",
                    "INPUT_LACKS_SHORT_PROTEINS": "#7b2cbf", "COMPOSITION_CANDIDATES_ONLY": "#8d99ae",
                    "NOT_DETECTED": "#adb5bd", "NOT_DETECTED_AT_GA_ONLY": "#ced4da"}
            z0 = z0.sort_values(["verdict", "family"])
            fig, ax = plt.subplots(figsize=(9.5, _dim(len(z0), .38, 4)))
            for yi, r in enumerate(z0.itertuples(index=False)):
                ax.add_patch(plt.Rectangle((0, yi - .4), 1, .8, color=vcol.get(r.verdict, "#adb5bd")))
                ax.text(-.03, yi, f"{r.family}  ({r.pfam_ids})", ha="right", va="center", fontsize=8.5)
                label = (r.verdict.replace("_", " ").lower()
                         .replace(" ga", " GA").replace("map", "map rules"))
                light = r.verdict in ("NOT_DETECTED", "NOT_DETECTED_AT_GA_ONLY",
                                      "COMPOSITION_CANDIDATES_ONLY")
                # text is always dark enough to read; the swatch carries the colour
                ax.text(1.05, yi - .1, label, va="center", fontsize=9, fontweight="bold",
                        color="#495057" if light else vcol.get(r.verdict, "#495057"))
                detail = (f"GA hits {r.raw_hits_at_GA} · below-GA hits {r.raw_hits_below_GA} · "
                          f"rescued {r.rescued_validated}")
                ax.text(1.05, yi + .24, detail, va="center", fontsize=7.2, color="#6c757d")
            ax.set_xlim(-2.2, 5.2); ax.set_ylim(len(z0) - .5, -.7)
            ax.axis("off")
            ax.set_title("Why each empty family is empty", loc="left")
            save(fig, "D1_zero_family_verdicts.pdf")
            W.captions[-1] = ("D1_zero_family_verdicts",
                "Every family with no calls at Pfam's gathering threshold, with the pipeline's verdict "
                "on why: missing from the database, filtered by the map, found only below the "
                "threshold, absent because the input lacks proteins that short, or not detected. "
                "Full detail and advice in qc_zero_families.tsv.")

    W.close()
    log(f"figures written to {fig_dir} (PDF, PNG, report_all_figures.pdf, FIGURES.md)")


# =========================================================================== #
#  Sequence and domain utilities
# =========================================================================== #
def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def tool(name):
    """Absolute path to an external tool, or None.
    DEFENSOME_<NAME> (e.g. DEFENSOME_HMMSEARCH=/opt/hmmer/bin/hmmsearch) wins
    over PATH, for clusters where HMMER and Python cannot share a module set."""
    env = os.environ.get("DEFENSOME_" + name.upper().replace("-", "_"))
    if env and os.path.exists(env) and os.access(env, os.X_OK):
        return env
    return shutil.which(name)


def run(cmd, **kw):
    cmd = [str(c) for c in cmd]
    resolved = tool(cmd[0])
    if resolved:
        cmd[0] = resolved
    log("$ " + " ".join(cmd))
    return subprocess.run(cmd, check=True, **kw)


def all_domains(path):
    """Every Pfam domain on every protein, with coordinates. Unlike
    domain_coverage() this keeps each hit separately so architectures and
    per-domain completeness can be reconstructed."""
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 22)
            if len(f) < 22:
                continue
            rows.append((f[0], int(f[2]), f[3], f[4].split(".")[0], int(f[5]),
                         int(f[15]), int(f[16]), int(f[17]), int(f[18]), float(f[13])))
    return pd.DataFrame(rows, columns=["protein","prot_len","pfam_name","pfam_id",
                                       "hmm_len","hmm_from","hmm_to","ali_from",
                                       "ali_to","score"])


def check_fresh(out, dm, mapfile):
    """Stop if gene_calls_all.tsv was built with a different version of the map.

    Editing the map and rerunning a downstream command without redoing
    `annotate` silently mixes two family vocabularies. `annotate` is cheap
    because it reads the cached domtblout, so the fix is always the same.
    """
    calls_f = os.path.join(out, "gene_calls_all.tsv")
    if not os.path.exists(calls_f):
        sys.exit(f"ERROR: {calls_f} not found. Run `annotate` first.")
    fams = set(pd.read_csv(calls_f, sep="\t", usecols=["family"]).family.unique())
    orphan = sorted(fams - set(dm.family))
    if orphan:
        sys.exit(
            f"\nERROR: gene_calls_all.tsv is stale.\n"
            f"  Families present in the calls but not in {os.path.basename(mapfile)}: "
            f"{orphan}\n"
            f"  The map was edited after `annotate` last ran.\n\n"
            f"  Fix (seconds, reuses the cached hmmsearch output):\n"
            f"    python3 defensome.py annotate --proteomes <DIR> --out {out}\n")
    # mapfile may not exist on disk at all when the embedded default is in use
    if os.path.exists(mapfile) and os.path.getmtime(mapfile) > os.path.getmtime(calls_f):
        log(f"WARNING: {os.path.basename(mapfile)} is newer than gene_calls_all.tsv. "
            "Family names still match, but thresholds may have changed. "
            "Rerun `annotate` to be sure.")


def scan_files(out):
    fs = sorted(glob.glob(os.path.join(out, "hmmsearch", "*.domtblout*")))
    if not fs:
        sys.exit(f"ERROR: no domtblout files in {out}/hmmsearch. Run `scan` first.")
    return [(re.sub(r"\.domtblout(\.gz)?$", "", os.path.basename(f)), f) for f in fs]


# =========================================================================== #
def cmd_extract(a):
    """Per-family protein FASTA. Needed by `trees` and `cyp`."""
    need_pandas()
    check_fresh(a.out, read_map(a.map), a.map)
    calls_f = os.path.join(a.out, "gene_calls_all.tsv")
    calls = pd.read_csv(calls_f, sep="\t")
    if a.families:
        want = set(a.families.split(","))
        calls = calls[calls.family.isin(want)]
        unknown = want - set(calls.family)
        if unknown:
            log(f"WARNING: no calls for {sorted(unknown)}")
    fdir = os.path.join(a.out, "fasta"); os.makedirs(fdir, exist_ok=True)

    wanted = calls.groupby("species").protein.apply(set).to_dict()
    seqs = {}
    for sp, faa in species_list(a.proteomes):
        if sp not in wanted:
            continue
        keep = wanted[sp]
        for hdr, seq in read_fasta(faa):
            pid = hdr.split()[0]
            if pid in keep:
                seqs[(sp, pid)] = seq
    missing = sum(1 for sp, g in calls.groupby("species")
                  for p in g.protein.unique() if (sp, p) not in seqs)
    if missing:
        log(f"WARNING: {missing} called proteins not found in the proteome FASTAs "
            "(header/ID mismatch?)")

    n = 0
    for fam, g in calls.groupby("family"):
        path = os.path.join(fdir, f"{fam}.faa")
        with open(path, "w") as fh:
            for sp, prot in g[["species","protein"]].drop_duplicates().itertuples(index=False):
                s = seqs.get((sp, prot))
                if s:
                    fh.write(f">{sp}|{prot}\n{s}\n")
        n += 1
    log(f"extract done: {n} family FASTA files in {fdir}")


# =========================================================================== #
def cmd_domains(a):
    """Per-protein domain architecture and completeness.

    For every defensome protein this records each Pfam domain present, where it
    sits, and how much of the HMM it covers. A protein is COMPLETE when every
    domain its family requires is present above the family coverage threshold.
    Anything else is PARTIAL (domain present but truncated) or MISSING_DOMAIN.
    FRAGMENT flags proteins far shorter than the family norm, which is usually a
    split gene model rather than a real short paralogue."""
    need_pandas()
    dm_raw = read_map(a.map)
    check_fresh(a.out, dm_raw, a.map)
    dm = dm_raw.set_index("family")
    calls = pd.read_csv(os.path.join(a.out, "gene_calls_all.tsv"), sep="\t")
    ddir = os.path.join(a.out, "domains"); os.makedirs(ddir, exist_ok=True)

    want = calls.groupby("species").protein.apply(set).to_dict()
    dom_rows = []
    for sp, path in scan_files(a.out):
        if sp not in want:
            continue
        d = all_domains(path)
        d = d[d.protein.isin(want[sp])].copy()
        d.insert(0, "species", sp)
        dom_rows.append(d)
    dom = pd.concat(dom_rows, ignore_index=True)
    dom["dom_cov"] = ((dom.hmm_to - dom.hmm_from + 1) / dom.hmm_len).round(3)
    dom.to_csv(os.path.join(ddir, "domain_table.tsv"), sep="\t", index=False)
    log(f"domain_table.tsv: {len(dom)} domain hits on {dom.protein.nunique()} proteins")

    # Ordered architecture string per protein, e.g. GST_N-GST_C
    arch = (dom.sort_values(["species","protein","ali_from"])
              .groupby(["species","protein"])
              .agg(architecture=("pfam_name", lambda s: "-".join(s)),
                   n_domains=("pfam_name", "size"),
                   prot_len=("prot_len", "first")).reset_index())
    arch.to_csv(os.path.join(ddir, "architectures.tsv"), sep="\t", index=False)

    # Per-protein completeness against the family's required domains
    # Merge non-overlapping HMM segments per (protein, pfam) exactly as
    # annotate does. Many Pfam models match in two pieces (FMO's paired
    # Rossmann folds, ABC_membrane, thioredoxin); scoring on the best single
    # segment understates coverage and mislabels intact proteins as PARTIAL.
    best = {}
    for (sp_, pr_, pf_), g_ in dom.groupby(["species","protein","pfam_id"], sort=False):
        iv = sorted(zip(g_.hmm_from.astype(int), g_.hmm_to.astype(int)))
        cs = ce = None; tot = 0
        for s_, e_ in iv:
            if cs is None: cs, ce = s_, e_
            elif s_ <= ce + 1: ce = max(ce, e_)
            else: tot += ce - cs + 1; cs, ce = s_, e_
        if cs is not None: tot += ce - cs + 1
        hl = float(g_.hmm_len.iloc[0])
        best[(sp_, pr_, pf_)] = tot / hl if hl else 0.0
    med_len = calls.groupby("family").prot_len.median().to_dict()

    rows = []
    skipped = set()
    for r in calls.itertuples(index=False):
        if r.family not in dm.index:
            skipped.add(r.family); continue
        fam = dm.loc[r.family]
        req = list(fam.pfam_ids) if fam.rule == "ALL" else [
            p for p in fam.pfam_ids if (r.species, r.protein, p) in best]
        covs = {p: best.get((r.species, r.protein, p), 0.0) for p in req}
        n_ok = sum(1 for v in covs.values() if v >= fam.min_cov)
        if not req:
            status = "MISSING_DOMAIN"
        elif n_ok == len(req):
            status = "COMPLETE"
        elif any(v > 0 for v in covs.values()):
            status = "PARTIAL"
        else:
            status = "MISSING_DOMAIN"
        lr = r.prot_len / med_len.get(r.family, r.prot_len)
        if status == "COMPLETE" and lr < 0.5:
            status = "FRAGMENT"
        rows.append((r.species, r.protein, r.family, r.category, status,
                     len(req), n_ok, round(min(covs.values()) if covs else 0.0, 3),
                     round(sum(covs.values())/len(covs), 3) if covs else 0.0,
                     r.prot_len, round(lr, 2),
                     ";".join(f"{k}:{v:.2f}" for k, v in covs.items())))
    comp = pd.DataFrame(rows, columns=["species","protein","family","category",
        "status","n_required","n_complete","min_domain_cov","mean_domain_cov",
        "prot_len","length_ratio","domain_coverages"])
    if skipped:
        log(f"WARNING: skipped families absent from the map: {sorted(skipped)}")
    comp.to_csv(os.path.join(ddir, "completeness.tsv"), sep="\t", index=False)

    summ = (comp.groupby(["species","family"]).status
              .value_counts().unstack(fill_value=0).reset_index())
    for s in ("COMPLETE","PARTIAL","FRAGMENT","MISSING_DOMAIN"):
        if s not in summ:
            summ[s] = 0
    summ["n"] = summ[["COMPLETE","PARTIAL","FRAGMENT","MISSING_DOMAIN"]].sum(axis=1)
    summ["pct_complete"] = (summ.COMPLETE / summ.n * 100).round(1)
    summ.to_csv(os.path.join(ddir, "completeness_by_species_family.tsv"),
                sep="\t", index=False)

    byfam = comp.groupby("family").status.value_counts().unstack(fill_value=0)
    byfam["pct_complete"] = (byfam.get("COMPLETE", 0) /
                             byfam.sum(axis=1) * 100).round(1)
    print("\n=== Domain completeness by family ===")
    print(byfam.sort_values("pct_complete").to_string())
    log("domains done")


# =========================================================================== #
CLAN_RE = re.compile(r"clan=(\S+)")


def cmd_cyp(a):
    """Assign CYPs to clans using a labelled reference set.

    Builds one HMM per clan from the reference sequences, then scores every
    CYP against all clan HMMs and takes the best. `margin` is the bitscore gap
    to the runner-up: small margins are assignments you should not trust, and
    the deep clan splits mean a Drosophila-only reference transfers imperfectly
    to Lepidoptera. Check the tree before believing a marginal call."""
    need_pandas()
    for exe in ("mafft", "hmmbuild", "hmmsearch"):
        if not tool(exe):
            sys.exit(f"ERROR: {exe} not found. module load HMMER MAFFT\n"
                     f"  (or set DEFENSOME_{exe.upper()}=/path/to/{exe})")
    if not os.path.exists(a.refs):
        sys.exit(f"ERROR: reference FASTA not found: {a.refs}")
    cdir = os.path.join(a.out, "cyp"); os.makedirs(os.path.join(cdir, "clan_hmms"), exist_ok=True)

    fam_faa = os.path.join(a.out, "fasta", f"{a.family}.faa")
    if not os.path.exists(fam_faa):
        sys.exit(f"ERROR: {fam_faa} not found. Run `extract` first.")
    if getattr(a, "samples", None):
        keep = {x.strip() for x in a.samples.split(",") if x.strip()}
        sub = os.path.join(cdir, f"{a.family}.subset.faa")
        n_in = n_out = 0
        with open(sub, "w") as fh:
            for hdr, seq in read_fasta(fam_faa):
                n_in += 1
                if hdr.split("|", 1)[0] in keep:
                    fh.write(f">{hdr}\n{seq}\n"); n_out += 1
        if not n_out:
            sys.exit(f"ERROR: --samples matched none of the {n_in} sequences in {fam_faa}.\n"
                     f"  Sample ids are the part before '|' in each header.")
        log(f"--samples: {n_out} of {n_in} {a.family} sequences, from {len(keep)} samples")
        fam_faa = sub

    by_clan = {}
    for hdr, seq in read_fasta(a.refs):
        m = CLAN_RE.search(hdr)
        if m:
            by_clan.setdefault(m.group(1), []).append((hdr.split()[0], seq))
    if not by_clan:
        sys.exit(f"ERROR: no 'clan=' tags in {a.refs} headers")
    log("reference clans: " + ", ".join(f"{k}({len(v)})" for k, v in sorted(by_clan.items())))
    # Best-hit assignment has no "none of these" outcome, so a reference built
    # for one lineage will confidently mislabel every P450 from another.
    print("\nNOTE: each query goes to its best-scoring reference clan; there is no\n"
          "      'unassigned' outcome. Clan labels mean something only for the lineage\n"
          "      the references came from: an insect CYP2/3/4/MITO reference will place\n"
          "      a fish or mollusc P450 into an insect clan with a large and meaningless\n"
          "      margin. Use --samples to restrict the run, or references for that lineage.")

    hmms = []
    for clan, seqs in sorted(by_clan.items()):
        if len(seqs) < 2:
            log(f"skipping clan {clan}: only {len(seqs)} reference sequence"); continue
        base = os.path.join(cdir, "clan_hmms", clan)
        with open(base + ".faa", "w") as fh:
            for n, s in seqs:
                fh.write(f">{n}\n{s}\n")
        with open(base + ".aln", "w") as fh:
            run(["mafft", "--auto", "--quiet", "--thread", a.threads, base + ".faa"], stdout=fh)
        run(["hmmbuild", "--amino", "-n", clan, base + ".hmm", base + ".aln"],
            stdout=subprocess.DEVNULL)
        hmms.append((clan, base + ".hmm"))

    scores = {}
    for clan, hmm in hmms:
        tbl = os.path.join(cdir, f"hits_{clan}.tbl")
        run(["hmmsearch", "--noali", "-E", "1e-5", "--cpu", a.threads,
             "--tblout", tbl, "-o", os.devnull, hmm, fam_faa])
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                scores.setdefault(f[0], {})[clan] = float(f[5])

    rows = []
    for tgt, sc in scores.items():
        ranked = sorted(sc.items(), key=lambda x: -x[1])
        top, tops = ranked[0]
        margin = tops - ranked[1][1] if len(ranked) > 1 else tops
        sp, prot = tgt.split("|", 1) if "|" in tgt else ("NA", tgt)
        rows.append((sp, prot, top, round(tops, 1), round(margin, 1),
                     "LOW" if margin < 20 else "OK"))
    calls = pd.DataFrame(rows, columns=["species","protein","clan","score",
                                        "margin","confidence"])
    unassigned = set()
    for hdr, _ in read_fasta(fam_faa):
        if hdr.split()[0] not in scores:
            unassigned.add(hdr.split()[0])
    calls.to_csv(os.path.join(cdir, f"{a.family}_clan_calls.tsv"), sep="\t", index=False)

    counts = calls.groupby(["species","clan"]).size().unstack(fill_value=0)
    counts["UNASSIGNED"] = pd.Series(
        {s.split("|")[0]: 1 for s in unassigned}).groupby(level=0).sum() \
        if unassigned else 0
    counts = counts.fillna(0).astype(int)
    counts.to_csv(os.path.join(cdir, f"{a.family}_clan_counts.tsv"), sep="\t")
    print(f"\n=== {a.family} clan counts per species ===")
    print(counts.to_string())
    print(f"\nlow-confidence assignments (margin < 20 bits): "
          f"{(calls.confidence == 'LOW').sum()} of {len(calls)}")
    if unassigned:
        log(f"{len(unassigned)} sequences matched no clan HMM above E=1e-5")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fdir = os.path.join(a.out, "figures"); os.makedirs(fdir, exist_ok=True)
        cl = counts.drop(columns=[c for c in ["UNASSIGNED"] if c in counts])
        frac = cl.div(cl.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        frac = frac.loc[frac.sort_values(list(frac.columns)[0]).index]
        fig, ax = plt.subplots(figsize=(9, _dim(len(frac), .24, 6)))
        left = np.zeros(len(frac))
        for c_ in frac.columns:
            ax.barh(range(len(frac)), frac[c_].values, left=left, label=c_)
            left += frac[c_].values
        ax.set_yticks(range(len(frac))); ax.set_yticklabels(frac.index, fontsize=6)
        ax.set_xlabel(f"fraction of {a.family}s"); ax.legend(fontsize=7, ncol=len(frac.columns))
        ax.set_title(f"{a.family} clan composition")
        fig.tight_layout(); fig.savefig(os.path.join(fdir, f"{a.family}_clans.pdf")); plt.close(fig)
        log(f"wrote {fdir}/{a.family}_clans.pdf")
    except ImportError:
        pass
    log("cyp done")


# =========================================================================== #
def metadata_lookup(metadata, ids):
    """Metadata rows keyed by sample id. The key column is whichever column
    matches the most ids (spaces read as underscores), which is how the
    dashboard matches it too, so a filter means the same thing in both."""
    need_pandas()
    if not metadata or not os.path.exists(metadata):
        sys.exit(f"ERROR: metadata file not found: {metadata}")
    md = pd.read_csv(metadata, sep="\t", comment="#", dtype=str)
    want, key, best = set(ids), None, 0
    for c in md.columns:
        v = set(md[c].astype(str).str.strip().str.replace(" ", "_", regex=False))
        if len(want & v) > best:
            key, best = c, len(want & v)
    if not key:
        sys.exit(f"ERROR: no column of {metadata} matches the sample IDs "
                 f"(e.g. {sorted(want)[:3]})")
    md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
    return md.drop_duplicates(subset=[key]).set_index(key), key, best


def parse_where(items):
    """--where COLUMN=VALUE[,VALUE...]; case-insensitive, repeatable (AND)."""
    out = []
    for w in items or []:
        if "=" not in w:
            sys.exit(f"ERROR: --where needs COLUMN=VALUE, got {w!r}")
        c, v = w.split("=", 1)
        vals = [x.strip().lower() for x in v.split(",") if x.strip()]
        if not c.strip() or not vals:
            sys.exit(f"ERROR: --where needs COLUMN=VALUE, got {w!r}")
        out.append((c.strip(), vals))
    return out


def select_samples(metadata, where, ids):
    """Sample IDs whose metadata satisfies every --where condition."""
    md, key, n = metadata_lookup(metadata, ids)
    cols = {c.lower(): c for c in md.columns}
    keep = set(md.index) & set(ids)
    for c, vals in where:
        if c.lower() not in cols:
            sys.exit(f"ERROR: --where column '{c}' is not in {metadata}. "
                     f"Columns: {', '.join(md.columns)}")
        col = md[cols[c.lower()]].astype(str).str.strip().str.lower()
        keep &= set(col[col.isin(vals)].index)
    return keep


def tree_label(fam, where, complete_only, clan, name, thin=None):
    """File stem for a tree: the family, then each restriction, so a subset
    tree can never overwrite the whole-family one and its name says what it is."""
    if name:
        return f"{fam}__{re.sub(r'[^A-Za-z0-9._+-]+', '_', name)}"
    parts = [fam]
    for c, vals in where:
        parts.append(f"{c}-{'+'.join(vals)}")
    if clan:
        parts.append("clan-" + "+".join(clan))
    if complete_only:
        parts.append("complete")
    if thin:
        c_, _, k_ = thin.partition(":")
        parts.append(f"thin-{c_}{k_ or 1}")
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", "__".join(parts))


def cmd_trees(a):
    """MAFFT alignment and FastTree per family, optionally on a subset:
    a taxon (--where order=Lepidoptera), complete gene models only, or one
    CYP clan. Each tree gets a JSON sidecar recording exactly what went in."""
    import json
    for exe in ("mafft", "FastTree"):
        if not tool(exe) and not tool(exe.lower()):
            sys.exit(f"ERROR: {exe} not found. module load MAFFT FastTree")
    ft = "FastTree" if tool("FastTree") else "fasttree"
    fdir = os.path.join(a.out, "fasta")
    tdir = os.path.join(a.out, "trees"); os.makedirs(tdir, exist_ok=True)
    fams = a.families.split(",") if a.families else \
        [os.path.basename(f)[:-4] for f in sorted(glob.glob(os.path.join(fdir, "*.faa")))
         if "__" not in os.path.basename(f)]
    where = parse_where(getattr(a, "where", None))
    if (where or getattr(a, "thin", None)) and not a.metadata:
        sys.exit("ERROR: --where and --thin need --metadata")
    if getattr(a, "thin", None) and not re.fullmatch(r"[^:]+(:\d+)?", a.thin):
        sys.exit(f"ERROR: --thin needs COLUMN or COLUMN:N, got {a.thin!r}")
    clan = [c.strip().upper() for c in (getattr(a, "clan", None) or "").split(",") if c.strip()]
    complete_only = bool(getattr(a, "complete_only", False))
    comp = None
    if complete_only:
        cf = os.path.join(a.out, "domains", "completeness.tsv")
        if not os.path.exists(cf):
            sys.exit(f"ERROR: --complete-only needs {cf}. Run `domains` first.")
        need_pandas()
        c_ = pd.read_csv(cf, sep="\t", usecols=["species", "protein", "family", "status"])
        comp = set(zip(c_.family, c_.species, c_.protein, c_.status))

    for fam in fams:
        faa = os.path.join(fdir, f"{fam}.faa")
        if not os.path.exists(faa):
            log(f"skipping {fam}: no {faa}"); continue
        label = tree_label(fam, where, complete_only, clan, getattr(a, "name", None),
                           getattr(a, "thin", None))
        recs = [(h, s_) for h, s_ in read_fasta(faa)]
        n_all = len(recs)
        why = []
        if where:
            ids = {h.split("|", 1)[0] for h, _ in recs}
            keep = select_samples(a.metadata, where, ids)
            recs = [r for r in recs if r[0].split("|", 1)[0] in keep]
            why.append(f"{len(keep)} of {len(ids)} species match "
                       + " and ".join(f"{c} in {{{', '.join(v)}}}" for c, v in where))
        if comp is not None:
            recs = [r for r in recs
                    if (fam, *r[0].split()[0].split("|", 1), "COMPLETE") in comp]
            why.append("complete gene models only")
        if clan:
            cf = os.path.join(a.out, "cyp", f"{fam}_clan_calls.tsv")
            if not os.path.exists(cf):
                log(f"skipping {fam}: --clan needs {cf} (run `cyp`)"); continue
            c_ = pd.read_csv(cf, sep="\t", usecols=["species", "protein", "clan"])
            ok = {(s_, p_) for s_, p_, c2 in zip(c_.species.astype(str), c_.protein.astype(str),
                                                 c_.clan.astype(str).str.upper()) if c2 in clan}
            recs = [r for r in recs if tuple(r[0].split()[0].split("|", 1)) in ok]
            why.append(f"clan {'+'.join(clan)}")
        # Over-long models are read-through or fused gene predictions. They do
        # not belong in a family alignment, they open gaps across every other
        # sequence, and MAFFT's memory grows with the square of the longest one.
        n_long, longest = 0, 0
        ratio = getattr(a, "max_len_ratio", 2.0)
        if ratio and recs:
            med = sorted(len(s_) for _, s_ in recs)[len(recs) // 2]
            keep_ = [r for r in recs if len(r[1]) <= ratio * med]
            n_long, longest = len(recs) - len(keep_), max(len(s_) for _, s_ in recs)
            if n_long:
                why.append(f"{n_long} over {ratio:g}x the median length of {med} aa "
                           f"dropped (longest {longest:,} aa)")
            recs = keep_
        # Thin species, not genes: N representative species per value of a
        # metadata column, each keeping every copy, so the tree still shows the
        # family's full diversity across the span of the taxon.
        thin = getattr(a, "thin", None)
        if thin:
            col_, _, k_ = thin.partition(":")
            k_ = int(k_ or 1)
            ids = sorted({h.split("|", 1)[0] for h, _ in recs})
            md, _, _ = metadata_lookup(a.metadata, ids)
            cols = {c.lower(): c for c in md.columns}
            if col_.lower() not in cols:
                sys.exit(f"ERROR: --thin column '{col_}' is not in {a.metadata}. "
                         f"Columns: {', '.join(md.columns)}")
            grp = md[cols[col_.lower()]].reindex(ids).fillna("(none)").astype(str)
            per = {}
            for h, _ in recs:
                sp_ = h.split("|", 1)[0]; per[sp_] = per.get(sp_, 0) + 1
            chosen = set()
            for g_, sps in grp.groupby(grp).groups.items():
                sps = list(sps)
                mid = sorted(per[x] for x in sps)[len(sps) // 2]
                # the most typical copy number in the group: not the inflated
                # assembly, not the sparse one
                sps.sort(key=lambda x: (abs(per[x] - mid), x))
                chosen.update(sps[:k_])
            recs = [r for r in recs if r[0].split("|", 1)[0] in chosen]
            why.append(f"{len(chosen)} representative species, {k_} per {col_} "
                       f"(of {len(ids)})")
        n = len(recs)
        if why:
            log(f"{label}: {n:,} of {n_all:,} sequences ({'; '.join(why)})")
        if n < 4:
            log(f"skipping {label}: only {n} sequences"); continue
        # Thinning genes within species is the last resort: it keeps the
        # taxonomic span but drops most of each species' copies, so the tree no
        # longer shows the family's diversity. The copies kept are those
        # closest to the median length, never the longest, which are the fused
        # and read-through models.
        sub_n = 0
        if a.per_species and n > a.max_seqs:
            med = sorted(len(s_) for _, s_ in recs)[n // 2]
            best = {}
            for hdr, seq in recs:
                best.setdefault(hdr.split("|", 1)[0], []).append((abs(len(seq) - med), hdr, seq))
            recs = []
            for rows in best.values():
                rows.sort()
                recs.extend((h, s_) for _, h, s_ in rows[:a.per_species])
            log(f"{label}: {n} sequences over --max-seqs {a.max_seqs}; keeping {a.per_species} "
                f"per species closest to {med} aa -> {len(recs)} from {len(best)} species. "
                "Most copies are now missing; --thin COLUMN:N keeps them all.")
            n, sub_n = len(recs), a.per_species
        if n > a.max_seqs:
            log(f"skipping {label}: {n} sequences exceeds --max-seqs {a.max_seqs}"
                + ("" if a.per_species else "; --per-species N would subsample instead"))
            continue
        aln, tre = os.path.join(tdir, f"{label}.aln"), os.path.join(tdir, f"{label}.tre")
        if os.path.exists(tre) and not a.force:
            log(f"{label}: tree exists, skipping (--force rebuilds)"); continue
        src = faa
        if why or sub_n:
            src = os.path.join(tdir, f"{label}.input.faa")
            with open(src, "w") as fh:
                for h, s_ in recs:
                    fh.write(f">{h}\n{s_}\n")
        # --retree 1 for large families: --auto picks an O(N^2) method that is
        # impractical past a few thousand sequences.
        mode = ["--auto"] if n <= 1000 else ["--retree", "1", "--maxiterate", "0"]
        with open(aln, "w") as fh:
            cmd = [tool("mafft") or "mafft", *mode, "--quiet", "--anysymbol",
                   "--thread", str(a.threads), src]
            log("$ " + " ".join(cmd))
            pr = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
        if pr.returncode != 0:
            os.remove(aln)
            mx = max(len(s_) for _, s_ in recs)
            killed = "Killed" in (pr.stderr or "") or pr.returncode in (-9, 137)
            log(f"{label}: MAFFT failed on {n:,} sequences (longest {mx:,} aa)"
                + (": killed, almost certainly out of memory." if killed else
                   f" (exit {pr.returncode}): {(pr.stderr or '').strip()[-300:]}"))
            if killed:
                log("  Run it as a job with memory to spare (sbatch/srun --mem=32G), "
                    "or make the input smaller: --thin, --clan, a narrower --where.")
            continue
        with open(tre, "w") as fh:
            run([ft, "-lg", "-quiet", aln], stdout=fh)
        meta = {"label": label, "family": fam, "n_tips": n, "thin": thin or "",
                "max_len_ratio": ratio, "dropped_long": n_long,
                "n_species": len({h.split("|", 1)[0] for h, _ in recs}),
                "where": {c: v for c, v in where}, "complete_only": complete_only,
                "clan": clan, "per_species": sub_n, "of_family": n_all,
                "built": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                "tool_version": __version__}
        with open(os.path.join(tdir, f"{label}.json"), "w") as fh:
            json.dump(meta, fh, indent=1)
        log(f"{label}: {n} sequences -> {tre}")
    log("trees done")


# =========================================================================== #
#  Newick parsing and tree layout. Self-contained: no ete3, no ggtree.
# =========================================================================== #
class Node:
    __slots__ = ("name","length","children","parent","x","y","angle","r","_sz")

    def __init__(self, name="", length=0.0):
        self.name, self.length, self.children, self.parent = name, length, [], None
        self.x = self.y = self.angle = self.r = 0.0
        self._sz = 0

    def is_leaf(self):
        return not self.children


def parse_newick(text):
    text = re.sub(r"\[[^\]]*\]", "", text.strip())     # drop comments
    text = text.rstrip(";").strip()
    pos = [0]

    def parse_node():
        node = Node()
        if text[pos[0]] == "(":
            pos[0] += 1
            while True:
                node.children.append(parse_node())
                node.children[-1].parent = node
                if text[pos[0]] == ",":
                    pos[0] += 1
                else:
                    break
            assert text[pos[0]] == ")", f"malformed Newick at {pos[0]}"
            pos[0] += 1
        m = re.match(r"([^,():;]*)", text[pos[0]:])
        label = m.group(1).strip()
        pos[0] += m.end()
        if text[pos[0]:pos[0]+1] == ":":
            pos[0] += 1
            m2 = re.match(r"[-\d.eE+]+", text[pos[0]:])
            node.length = float(m2.group(0)); pos[0] += m2.end()
        # A label on an internal node is support, not a name.
        if node.is_leaf():
            node.name = label.strip("'\"")
        return node
    return parse_node()


def tree_nodes(root):
    stack, out = [root], []
    while stack:
        n = stack.pop(); out.append(n); stack.extend(n.children)
    return out


def tree_leaves(root):
    return [n for n in tree_nodes(root) if n.is_leaf()]


def prune_tree(root, keep):
    """Drop tips not in `keep`, then collapse the resulting single-child nodes."""
    def rec(n):
        if n.is_leaf():
            return n if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length += n.length
            return kids[0]
        n.children = kids
        for k in kids:
            k.parent = n
        return n
    r = rec(root)
    if r is None:
        sys.exit("ERROR: no tips in the tree match the count matrix")
    r.length = 0.0
    return r


def ladderize(root):
    def rec(n):
        n._sz = 1 if n.is_leaf() else sum(rec(c) for c in n.children)
        n.children.sort(key=lambda c: c._sz)
        return n._sz
    rec(root)
    return root


def set_depths(root, cladogram=True):
    """Radial/horizontal coordinate for every node."""
    if cladogram:
        def depth(n, d=0):
            n.r = d
            for c in n.children:
                depth(c, d + 1)
        depth(root)
        maxd = max(n.r for n in tree_nodes(root))
        # push all tips to the rim so labels line up
        for n in tree_nodes(root):
            n.r = n.r / maxd if maxd else 0.0
        for n in tree_leaves(root):
            n.r = 1.0
    else:
        def dist(n, d=0.0):
            n.r = d
            for c in n.children:
                dist(c, d + max(c.length, 0.0))
        dist(root)
        maxd = max(n.r for n in tree_nodes(root)) or 1.0
        for n in tree_nodes(root):
            n.r /= maxd
    return root


def set_angles(root, span=350.0, start=0.0):
    lv = tree_leaves(root)
    step = span / max(len(lv), 1)
    for i, n in enumerate(lv):
        n.angle = math.radians(start + i * step)
    def rec(n):
        if n.is_leaf():
            return n.angle
        a = [rec(c) for c in n.children]
        n.angle = (min(a) + max(a)) / 2
        return n.angle
    rec(root)
    return root, step


# =========================================================================== #
# Distinct palettes so concentric rings never share a colour. Sharing one made
# CYP clan "CYP3" and diet "Polyphagous" the same orange, which is unreadable.
_PALETTES = [
    ["#1b4965","#2a9d8f","#8ecae6","#023047","#5fa8d3","#0b6e4f","#89c2d9"],   # blues/teals
    ["#bc4749","#f4a261","#e9c46a","#6a040f","#dda15e","#9c6644","#ffb703"],   # warm
    ["#7b2cbf","#c77dff","#5a189a","#e0aaff","#3c096c","#9d4edd","#240046"],   # purples
]

def _palette(levels, ring=0):
    base = _PALETTES[ring % len(_PALETTES)]
    return {lv: base[i % len(base)] for i, lv in enumerate(levels)}


def cmd_tree(a):
    """Radial and rectangular species cladograms with trait and defensome rings.

    Reads any Newick (OrthoFinder's SpeciesTree_rooted.txt works directly),
    prunes it to the species in counts.tsv, then draws:
      - tip labels coloured by the metadata column given to --group-by
      - concentric rings, one per defensome family, z-scored across species
      - an outer bar ring for total CORE defensome
    """
    need_pandas()
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge, Rectangle
    except ImportError:
        sys.exit("ERROR: matplotlib required. pip install --user matplotlib")

    counts = pd.read_csv(os.path.join(a.out, "counts.tsv"), sep="\t", index_col=0)
    npath = os.path.join(a.out, "counts_per10k.tsv")
    norm = pd.read_csv(npath, sep="\t", index_col=0) if os.path.exists(npath) else None
    use = counts if a.raw or norm is None else norm
    lab = "copies" if use is counts else f"copies per 10k {size_unit(a.out)}"
    dm = read_map(a.map)

    root = ladderize(prune_tree(parse_newick(open(a.tree).read()), set(counts.index)))
    tips = [n.name for n in tree_leaves(root)]
    log(f"tree: {len(tips)} tips after pruning to counts.tsv")
    missing = sorted(set(counts.index) - set(tips))
    if missing:
        log(f"WARNING: not in the tree: {missing}")

    grp, cmap = None, {}
    if a.metadata and a.group_by:
        md = pd.read_csv(a.metadata, sep="\t", comment="#")
        want, key, bestn = set(counts.index), None, 0
        for c_ in md.columns:
            v = md[c_].astype(str).str.strip().str.replace(" ", "_", regex=False)
            if len(want & set(v)) > bestn:
                key, bestn = c_, len(want & set(v))
        if key and a.group_by in md.columns:
            md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
            grp = md.set_index(key)[a.group_by]
            cmap = _palette(sorted(grp.dropna().unique()))
            log(f"trait '{a.group_by}' from column '{key}' ({bestn} matched)")

    fams = (a.families.split(",") if a.families
            else [f for f in dm.loc[dm.tier == "CORE", "family"]
                  if f in use.columns and use[f].sum() > 0])
    fams = [f for f in fams if f in use.columns]
    m = use.loc[tips, fams]
    z = ((m - m.mean()) / m.std().replace(0, np.nan)).fillna(0).clip(-2.5, 2.5)
    total = use.loc[tips, fams].sum(axis=1)

    fig_dir = os.path.join(a.out, "figures"); os.makedirs(fig_dir, exist_ok=True)

    # ------------------------------------------------------------------ radial
    set_depths(root, cladogram=not a.phylogram)
    # Leave a wedge open at the top. Ring labels sit in it: at 90 degrees the
    # radial direction is vertical, so successive rings stack downward-to-upward
    # and horizontal label text never overlaps.
    GAP_DEG = 26.0
    root, step = set_angles(root, span=360.0 - GAP_DEG, start=90.0 + GAP_DEG/2)
    R_TREE, RING_W, GAP = 1.0, 0.085, 0.30
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.set_aspect("equal"); ax.axis("off")

    def pol(r, th):
        return r*math.cos(th), r*math.sin(th)

    for n in tree_nodes(root):                       # branches
        for c in n.children:
            r0 = n.r * R_TREE
            x0, y0 = pol(r0, c.angle); x1, y1 = pol(c.r*R_TREE, c.angle)
            ax.plot([x0, x1], [y0, y1], c="#444", lw=.9, zorder=1)
        if n.children:                                # connector arc
            a0 = min(c.angle for c in n.children); a1 = max(c.angle for c in n.children)
            arc = np.linspace(a0, a1, 40)
            ax.plot(n.r*R_TREE*np.cos(arc), n.r*R_TREE*np.sin(arc),
                    c="#444", lw=.9, zorder=1)

    half = math.radians(step) / 2 * 0.92
    cm_ = plt.cm.RdBu_r
    for k, fam in enumerate(fams):                    # family rings
        r0 = R_TREE + GAP + k*RING_W
        for n in tree_leaves(root):
            v = (z.loc[n.name, fam] + 2.5) / 5.0
            ax.add_patch(Wedge((0, 0), r0+RING_W*.9, math.degrees(n.angle-half),
                               math.degrees(n.angle+half), width=RING_W*.9,
                               facecolor=cm_(v), edgecolor="none", zorder=2))
        ax.text(0.0, r0 + RING_W*.45, fam, fontsize=7, ha="center",
                va="center", color="#333",
                bbox=dict(fc="white", ec="none", pad=.6, alpha=.85), zorder=4)

    r_bar = R_TREE + GAP + len(fams)*RING_W + 0.10    # outer total bar ring
    bmax = total.max() or 1
    bmin = total.min()
    ax.add_patch(Wedge((0, 0), r_bar, 0, 360, width=.004,
                       facecolor="#ccc", edgecolor="none", zorder=1))
    ax.text(0.0, r_bar + 0.16, f"total CORE ({lab})", fontsize=7.5,
            ha="center", va="center", color="#333",
            bbox=dict(fc="white", ec="none", pad=.6, alpha=.85), zorder=4)
    for n in tree_leaves(root):
        h = 0.30 * (total[n.name] - bmin*.85) / (bmax - bmin*.85)
        col = cmap.get(grp.get(n.name), "#888") if grp is not None else "#2a6f97"
        ax.add_patch(Wedge((0, 0), r_bar+h, math.degrees(n.angle-half),
                           math.degrees(n.angle+half), width=h,
                           facecolor=col, edgecolor="none", alpha=.85, zorder=2))
    for n in tree_leaves(root):                       # tip labels
        col = cmap.get(grp.get(n.name), "#222") if grp is not None else "#222"
        deg = math.degrees(n.angle) % 360
        flip = 90 < deg < 270
        ax.text(*pol(r_bar + 0.36, n.angle), n.name.replace("_", " "),
                rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                ha="right" if flip else "left", va="center",
                fontsize=7, color=col, style="italic", zorder=3)

    lim = r_bar + 1.05
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    if grp is not None:
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=v, label=str(k_)) for k_, v in cmap.items()],
                  loc="upper left", frameon=False, fontsize=10, title=a.group_by)
    sm = plt.cm.ScalarMappable(cmap=cm_, norm=plt.Normalize(-2.5, 2.5))
    cb = fig.colorbar(sm, ax=ax, shrink=.25, pad=.02, location="right")
    cb.set_label(f"z-score, {lab}", fontsize=9)
    ax.set_title(f"Defensome across the species tree  ({lab})", fontsize=14, pad=18)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "10_tree_radial.pdf"))
    plt.close(fig)

    # ------------------------------------------------------------- rectangular
    lv = tree_leaves(root)
    for i, n in enumerate(lv):
        n.y = i
    def sety(n):
        if n.is_leaf():
            return n.y
        ys = [sety(c) for c in n.children]
        n.y = (min(ys) + max(ys)) / 2
        return n.y
    sety(root)
    fig, (axt, axh) = plt.subplots(1, 2, figsize=(17, _dim(len(lv), .30, 8)),
                                   gridspec_kw={"width_ratios": [1, 2.1]})
    for n in tree_nodes(root):
        for c in n.children:
            axt.plot([n.r, c.r], [c.y, c.y], c="#444", lw=.9)
        if n.children:
            axt.plot([n.r, n.r], [min(c.y for c in n.children),
                                  max(c.y for c in n.children)], c="#444", lw=.9)
    for n in lv:
        col = cmap.get(grp.get(n.name), "#222") if grp is not None else "#222"
        axt.text(n.r + .02, n.y, n.name.replace("_", " "), va="center",
                 fontsize=7, color=col, style="italic")
    axt.set_xlim(-.02, 1.55); axt.set_ylim(-1, len(lv)); axt.axis("off")
    axt.set_title("species tree", fontsize=11)
    im = axh.imshow(z.values, aspect="auto", cmap=cm_, vmin=-2.5, vmax=2.5,
                    extent=[0, len(fams), len(lv)-.5, -.5])
    axh.set_xticks(np.arange(len(fams))+.5)
    axh.set_xticklabels(fams, rotation=90, fontsize=7.5)
    axh.set_yticks([]); axh.set_title(f"defensome, z-scored ({lab})", fontsize=11)
    fig.colorbar(im, ax=axh, shrink=.4, label="z-score")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "11_tree_heatmap.pdf"))
    plt.close(fig)
    log(f"wrote {fig_dir}/10_tree_radial.pdf and 11_tree_heatmap.pdf")


# =========================================================================== #
def cmd_genetree(a):
    """Radial gene tree for one family, annotated with clan and species trait.

    A 4000-tip tree with unlabelled dots tells you nothing. This draws two
    concentric rings (clan, then species trait), labels clan blocks that form
    contiguous arcs, and writes species labels when the tree is small enough
    for them to be legible. It also writes a clan-by-trait contingency table,
    which is usually the number you actually wanted from the picture.
    """
    need_pandas()
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge, Patch
    except ImportError:
        sys.exit("ERROR: matplotlib required")
    tre = a.tree or os.path.join(a.out, "trees", f"{a.family}.tre")
    if not os.path.exists(tre):
        sys.exit(f"ERROR: {tre} not found. Run `trees` first.")
    root = ladderize(parse_newick(open(tre).read()))
    lv = tree_leaves(root)
    n = len(lv)
    log(f"{a.family} gene tree: {n} tips")

    sp_of = {x.name: x.name.split("|")[0] for x in lv}

    # clan, if it has been computed
    clan_of, clan_cmap = {}, {}
    clan_f = os.path.join(a.out, "cyp", f"{a.family}_clan_calls.tsv")
    if os.path.exists(clan_f):
        cl = pd.read_csv(clan_f, sep="\t")
        clan_of = {f"{r.species}|{r.protein}": r.clan for r in cl.itertuples()}
        clan_cmap = _palette(sorted(set(clan_of.values())), ring=0)

    # species trait; with a sample sheet and no other trait, colour by method,
    # which is the variable a method-comparison run is actually about
    trait_of, trait_cmap, trait_name = {}, {}, None
    if getattr(a, "samplesheet", None) and not (a.metadata and a.group_by):
        sh_ = pd.read_csv(a.samplesheet, sep="\t")
        mth = dict(zip(sh_.sample_id.astype(str), sh_.method.astype(str)))
        trait_of = {t: mth.get(sp_of[t]) for t in sp_of}
        if any(trait_of.values()):
            trait_cmap = _palette(sorted({v for v in trait_of.values() if v}), ring=1)
            trait_name = "method"
            log(f"colouring by method from {a.samplesheet}")
    elif a.metadata and a.group_by:
        md = pd.read_csv(a.metadata, sep="\t", comment="#")
        sp_all = set(sp_of.values())
        key, bestn = None, 0
        for c_ in md.columns:
            v = md[c_].astype(str).str.strip().str.replace(" ", "_", regex=False)
            if len(sp_all & set(v)) > bestn:
                key, bestn = c_, len(sp_all & set(v))
        if key and a.group_by in md.columns:
            md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
            g = md.set_index(key)[a.group_by]
            trait_of = {t: g.get(sp_of[t]) for t in sp_of}
            trait_cmap = _palette(sorted({v for v in trait_of.values() if pd.notna(v)}), ring=1)
            trait_name = a.group_by
            log(f"trait '{a.group_by}' from column '{key}' ({bestn} species matched)")

    rings = []
    if clan_of:
        rings.append(("clan", clan_of, clan_cmap))
    if trait_of:
        rings.append((trait_name, trait_of, trait_cmap))
    if not rings:                      # fall back to colouring by species
        cm_ = _palette(sorted(set(sp_of.values())), ring=2)
        rings.append(("species", sp_of, cm_))

    set_depths(root, cladogram=False)
    GAP_DEG = 16.0
    root, step = set_angles(root, span=360.0 - GAP_DEG, start=90.0 + GAP_DEG/2)

    R, RW = 1.0, 0.055
    fig, ax = plt.subplots(figsize=(16, 16)); ax.set_aspect("equal"); ax.axis("off")
    for nd in tree_nodes(root):
        for c in nd.children:
            ax.plot([nd.r*R*math.cos(c.angle), c.r*R*math.cos(c.angle)],
                    [nd.r*R*math.sin(c.angle), c.r*R*math.sin(c.angle)],
                    c="#666", lw=.3, zorder=1)
        if nd.children:
            arc = np.linspace(min(c.angle for c in nd.children),
                              max(c.angle for c in nd.children), 24)
            ax.plot(nd.r*R*np.cos(arc), nd.r*R*np.sin(arc), c="#666", lw=.3, zorder=1)

    half = math.radians(step)/2 * 1.0
    for k, (name, mapping, cm_) in enumerate(rings):
        r0 = R + 0.06 + k*RW
        for nd in lv:
            col = cm_.get(mapping.get(nd.name), "#e9ecef")
            ax.add_patch(Wedge((0, 0), r0+RW*.92, math.degrees(nd.angle-half),
                               math.degrees(nd.angle+half), width=RW*.92,
                               facecolor=col, edgecolor="none", zorder=2))
        ax.text(0.0, r0+RW*.45, name, fontsize=9, ha="center", va="center",
                bbox=dict(fc="white", ec="#ddd", pad=.9, alpha=.95,
                          boxstyle="round,pad=0.3"), zorder=5)

    r_out = R + 0.06 + len(rings)*RW

    # label clan blocks that form contiguous arcs
    if clan_of:
        seq = [clan_of.get(nd.name) for nd in lv]
        i = 0
        min_run = max(12, n // 60)
        while i < len(seq):
            j = i
            while j+1 < len(seq) and seq[j+1] == seq[i]:
                j += 1
            if seq[i] is not None and (j - i + 1) >= min_run:
                a0, a1 = lv[i].angle, lv[j].angle
                mid = (a0 + a1) / 2
                arc = np.linspace(a0-half, a1+half, 40)
                rr = r_out + 0.05
                ax.plot(rr*np.cos(arc), rr*np.sin(arc),
                        c=clan_cmap.get(seq[i], "#333"), lw=3.0, solid_capstyle="butt",
                        zorder=3)
                deg = math.degrees(mid) % 360
                flip = 90 < deg < 270
                ax.text((rr+0.05)*math.cos(mid), (rr+0.05)*math.sin(mid),
                        f"{seq[i]} ({j-i+1})",
                        rotation=deg+180 if flip else deg, rotation_mode="anchor",
                        ha="right" if flip else "left", va="center", fontsize=8,
                        color=clan_cmap.get(seq[i], "#333"), zorder=4)
            i = j + 1
        r_out += 0.30

    # species labels only when they can actually be read
    label_limit = a.max_labels
    if n <= label_limit:
        for nd in lv:
            deg = math.degrees(nd.angle) % 360
            flip = 90 < deg < 270
            col = trait_cmap.get(trait_of.get(nd.name), "#333") if trait_of else "#333"
            ax.text((r_out+0.02)*math.cos(nd.angle), (r_out+0.02)*math.sin(nd.angle),
                    sp_of[nd.name].replace("_", " "),
                    rotation=deg+180 if flip else deg, rotation_mode="anchor",
                    ha="right" if flip else "left", va="center",
                    fontsize=max(2.0, min(6.0, 260.0/n)), color=col, style="italic",
                    zorder=4)
        r_out += 0.42
    else:
        ax.text(0, -r_out-0.16,
                f"{n} tips: species labels suppressed above --max-labels {label_limit}. "
                f"Ring colours and {a.family}_tip_table.tsv carry the same information.",
                ha="center", fontsize=9, color="#666")
        r_out += 0.22

    handles = []
    for name, mapping, cm_ in rings:
        handles.append(Patch(color="white", label=f"— {name} —"))
        handles += [Patch(color=v, label=str(k_)) for k_, v in cm_.items()][:12]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.5,
              bbox_to_anchor=(-0.02, 1.0))
    ax.set_xlim(-r_out-.12, r_out+.12); ax.set_ylim(-r_out-.12, r_out+.12)
    ttl = f"{a.family} gene tree: {n} sequences from {len(set(sp_of.values()))} species"
    ax.set_title(ttl, fontsize=15, pad=16)
    fig_dir = os.path.join(a.out, "figures"); os.makedirs(fig_dir, exist_ok=True)
    tag = "_".join(r[0] for r in rings).replace(" ", "")
    out = os.path.join(fig_dir, f"12_{a.family}_genetree_{tag}.pdf")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    log(f"wrote {out}")

    # the numbers behind the picture
    tips = pd.DataFrame({
        "tip": [nd.name for nd in lv],
        "species": [sp_of[nd.name] for nd in lv],
        "clan": [clan_of.get(nd.name) for nd in lv],
        "trait": [trait_of.get(nd.name) for nd in lv],
        "tip_order": range(len(lv)),
    })
    tips.to_csv(os.path.join(a.out, "trees", f"{a.family}_tip_table.tsv"),
                sep="\t", index=False)
    if clan_of:
        ct = pd.crosstab(tips.species, tips.clan)
        ct["total"] = ct.sum(axis=1)
        ct.to_csv(os.path.join(a.out, "trees", f"{a.family}_clan_by_species.tsv"), sep="\t")
        print(f"\n=== {a.family} clan counts per species ===")
        print(ct.to_string())
        if trait_of:
            ct2 = pd.crosstab(tips.trait, tips.clan)
            print(f"\n=== {a.family} clan by {trait_name} (totals across species) ===")
            print(ct2.to_string())
            frac = ct2.div(ct2.sum(axis=1), axis=0).round(3)
            print("\nas fractions within each group:"); print(frac.to_string())
            ct2.to_csv(os.path.join(a.out, "trees", f"{a.family}_clan_by_{trait_name}.tsv"),
                       sep="\t")

            # per-species clan composition, which is the comparable version
            per = pd.crosstab(tips.species, tips.clan)
            sz = os.path.join(a.out, "proteome_sizes.tsv")
            if os.path.exists(sz):
                nprot = pd.read_csv(sz, sep="\t").set_index("species").n_proteins
                per10k = per.div(nprot.reindex(per.index), axis=0) * 10000
                per10k.round(2).to_csv(
                    os.path.join(a.out, "trees", f"{a.family}_clan_per10k.tsv"), sep="\t")
                tmap = tips.drop_duplicates("species").set_index("species").trait
                fig, ax = plt.subplots(figsize=(10, _dim(len(per10k), .26, 6)))
                order = per10k.sum(axis=1).sort_values().index
                left = np.zeros(len(order))
                for c_ in per10k.columns:
                    ax.barh(range(len(order)), per10k.loc[order, c_].values, left=left,
                            label=c_, color=clan_cmap.get(c_, "#999"))
                    left += per10k.loc[order, c_].values
                ax.set_yticks(range(len(order)))
                ax.set_yticklabels([s.replace("_", " ") for s in order], fontsize=7.5,
                                   style="italic")
                if trait_cmap:
                    for t_, s_ in zip(ax.get_yticklabels(), order):
                        t_.set_color(trait_cmap.get(tmap.get(s_), "#222"))
                ax.set_xlabel(f"{a.family} copies per 10k {size_unit(a.out)}")
                ax.legend(fontsize=8, ncol=len(per10k.columns), loc="lower right")
                ax.set_title(f"{a.family} clan composition per species "
                             f"(tip labels coloured by {trait_name})")
                fig.tight_layout()
                fig.savefig(os.path.join(fig_dir, f"13_{a.family}_clans_per_species.pdf"))
                plt.close(fig)
                log(f"wrote {fig_dir}/13_{a.family}_clans_per_species.pdf")


# =========================================================================== #
def size_unit(out):
    """What the per-10k denominator actually counts.

    After `collapse` the proteomes hold one sequence per gene, so the
    denominator is genes. Without that step it is proteins, and with several
    isoforms per gene those are not the same number. Labelling both "proteins"
    overstated what the ratio means."""
    return "genes" if os.path.exists(os.path.join(out, "collapse_stats.tsv")) else "proteins"


def _read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", **kw) if os.path.exists(path) else None


def build_payload(out, mapfile, metadata=None, group_by=None, full=True,
                  samplesheet=None, ga_only=False):
    """Collect every output file that exists into one JSON-serialisable dict."""
    dm = read_map(mapfile)
    counts = _read_tsv(os.path.join(out, "counts.tsv"), index_col=0)
    if counts is None:
        sys.exit(f"ERROR: {out}/counts.tsv not found. Run `annotate` first.")
    per10k = _read_tsv(os.path.join(out, "counts_per10k.tsv"), index_col=0)
    sizes = _read_tsv(os.path.join(out, "proteome_sizes.tsv"))
    qc = _read_tsv(os.path.join(out, "qc_annotation.tsv"), index_col=0)
    qcf = _read_tsv(os.path.join(out, "qc_families.tsv"), index_col=0)
    qcz = _read_tsv(os.path.join(out, "qc_zero_families.tsv"))
    fsum = _read_tsv(os.path.join(out, "family_summary.tsv"))
    comp = _read_tsv(os.path.join(out, "domains", "completeness.tsv"))
    arch = _read_tsv(os.path.join(out, "domains", "architectures.tsv"))
    calls = _read_tsv(os.path.join(out, "gene_calls_all.tsv"))

    species = list(counts.index)
    fams = [f for f in dm.family if f in counts.columns]

    trait, trait_name = {}, None
    if metadata and group_by and os.path.exists(metadata):
        md = pd.read_csv(metadata, sep="\t", comment="#")
        if group_by in md.columns:
            key, best = None, 0
            want = set(species)
            for c in md.columns:
                v = md[c].astype(str).str.strip().str.replace(" ", "_", regex=False)
                if len(want & set(v)) > best:
                    key, best = c, len(want & set(v))
            if key:
                md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
                g = md.set_index(key)[group_by]
                trait = {s: (None if pd.isna(g.get(s)) else str(g.get(s))) for s in species}
                trait_name = group_by
                log(f"dashboard trait '{group_by}' from column '{key}' ({best} matched)")

    # Every categorical column of the metadata, not just --group-by, so the
    # dashboard can filter and regroup by any of them without regenerating.
    # Taxonomic ranks come first and in rank order; anything else follows.
    tax, tax_cols, tax_all = {}, [], []
    if metadata and os.path.exists(metadata):
        md2 = pd.read_csv(metadata, sep="\t", comment="#")
        key2, best2 = None, 0
        want2 = set(species)
        for c in md2.columns:
            v = md2[c].astype(str).str.strip().str.replace(" ", "_", regex=False)
            if len(want2 & set(v)) > best2:
                key2, best2 = c, len(want2 & set(v))
        if key2:
            md2[key2] = md2[key2].astype(str).str.strip().str.replace(" ", "_", regex=False)
            md2 = md2.drop_duplicates(subset=[key2]).set_index(key2)
            RANK_ORDER = ["superkingdom", "kingdom", "phylum", "subphylum", "class",
                          "subclass", "order", "suborder", "superfamily", "family",
                          "subfamily", "genus"]
            # Per-species labels (a vernacular name, a taxid) are carried
            # through for display but are useless as filters or grouping
            # variables: one value per row makes every group a singleton.
            LABEL_ONLY = {"common_name", "common", "vernacular_name", "vernacular",
                          "english_name", "taxid", "tax_id", "ncbi_taxid",
                          "accession", "assembly", "notes", "comment"}
            label_cols = [c for c in md2.columns
                          if c != key2 and c.lower() in LABEL_ONLY]
            ranked = [c for c in RANK_ORDER if c in md2.columns and c not in label_cols]
            # Do NOT test `dtype == object` here. Pandas 3 reads text columns as
            # the "str" extension dtype, not object, so that test silently
            # rejected every non-rank column (habitat included) on newer
            # pandas while working on older. Cardinality is the real criterion:
            # a column with a handful of repeated values is a grouping
            # variable whatever it is stored as, and a continuous measurement
            # fails it by having one value per species.
            def _groupable(c):
                s_ = md2[c]
                if pd.api.types.is_float_dtype(s_):
                    return False
                n_ = s_.nunique(dropna=True)
                return 1 < n_ <= max(60, len(species) // 4)
            others = [c for c in md2.columns
                      if c not in ranked and c not in label_cols and c != key2
                      and _groupable(c)]
            tax_cols = ranked + others
            tax_all = tax_cols + label_cols
            for c in tax_all:
                col_ = md2[c]
                tax[c] = {s: (None if s not in col_.index or pd.isna(col_.get(s))
                              or str(col_.get(s)).strip() in ("", "nan")
                              else str(col_.get(s)).strip()) for s in species}
            if tax_cols:
                log(f"dashboard taxonomy columns: {', '.join(tax_cols)}")

    samples = {}
    if samplesheet and os.path.exists(samplesheet):
        sh_ = pd.read_csv(samplesheet, sep="\t")
        for r in sh_.itertuples(index=False):
            d = r._asdict()
            samples[str(d["sample_id"])] = {
                "species": str(d.get("species", "") or ""),
                "method": str(d.get("method", "") or ""),
                "pipeline": str(d.get("pipeline", "") or "")}
        if not trait and group_by and group_by in sh_.columns:
            m_ = dict(zip(sh_.sample_id.astype(str), sh_[group_by].astype(str)))
            trait = {s_: (m_.get(s_) or None) for s_ in species}
            trait_name = group_by
            log(f"dashboard trait '{group_by}' taken from the sample sheet")
        elif not trait and any(v["method"] for v in samples.values()):
            trait = {s_: samples.get(s_, {}).get("method") or None for s_ in species}
            trait_name = "method"
            log("dashboard trait defaults to 'method' from the sample sheet")

    nprot = {}
    if sizes is not None:
        nprot = dict(zip(sizes.species, sizes.n_proteins.astype(int)))

    sp_rows = []
    for s in species:
        row = {"name": s, "n_proteins": int(nprot.get(s, 0)),
               "trait": trait.get(s),
               "tax": {c: tax[c].get(s) for c in tax_all} if tax_all else {},
               "species_name": samples.get(s, {}).get("species", ""),
               "method": samples.get(s, {}).get("method", ""),
               "pipeline": samples.get(s, {}).get("pipeline", "")}
        if qc is not None and s in qc.index:
            row["max_ratio"] = float(qc.loc[s].get("max_ratio", float("nan")))
            row["flag"] = str(qc.loc[s].get("flag", ""))
            row["worst_family"] = str(qc.loc[s].get("worst_family", ""))
        sp_rows.append(row)

    P = {
        "size_unit": size_unit(out),
        "generated": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "out_dir": os.path.abspath(out),
        "trait_name": trait_name,
        "species": sp_rows,
        "families": [{"family": r.family, "category": r.category, "tier": r.tier,
                      "pfam": ",".join(r.pfam_ids), "rule": r.rule,
                      "min_cov": float(r.min_cov), "min_len": int(r.min_len),
                      "notes": ("" if pd.isna(r.notes) else str(r.notes))}
                     for r in dm.itertuples() if r.family in counts.columns],
        "counts": {f: [int(v) for v in counts[f]] for f in fams},
        "per10k": ({f: [round(float(v), 4) for v in per10k[f]] for f in fams}
                   if per10k is not None else None),
    }

    if fsum is not None:
        P["family_summary"] = fsum.where(pd.notna(fsum), None).to_dict("records")
    if qcf is not None:
        P["qc_families"] = [dict(family=i, **{k: (None if pd.isna(v) else v)
                                              for k, v in r.items()})
                            for i, r in qcf.iterrows()]
    if qcz is not None:
        P["qc_zero"] = qcz.where(pd.notna(qcz), None).to_dict("records")
    cr = _read_tsv(os.path.join(out, "counts_rescued.tsv"), index_col=0)
    rcalls = _read_tsv(os.path.join(out, "rescue_calls.tsv"))
    P["rescued"] = {}
    if cr is not None and rcalls is not None and len(rcalls):
        for f in sorted(rcalls[rcalls.validation == "PASS"].family.unique()):
            if f in cr.columns and f in counts.columns:
                extra = (cr[f] - counts[f]).reindex(counts.index).fillna(0)
                P["rescued"][f] = [int(v) for v in extra]
    # For families the map marks `rescue`, the copies found below Pfam's
    # gathering threshold that passed the family's validation are real genes
    # the threshold missed (short or divergent domains). The dashboard counts
    # them, marked as such; counts.tsv stays GA-only for the statistics.
    P["counts_basis"] = "GA"
    if P["rescued"] and not ga_only:
        nprot_ = {r["name"]: r["n_proteins"] for r in sp_rows}
        for f, extra in P["rescued"].items():
            if f in P["counts"]:
                P["counts"][f] = [a_ + b_ for a_, b_ in zip(P["counts"][f], extra)]
                if P.get("per10k"):
                    P["per10k"][f] = [round(c_ / nprot_[s_] * 1e4, 4) if nprot_.get(s_) else 0.0
                                      for c_, s_ in zip(P["counts"][f], species)]
        P["counts_basis"] = "GA + validated below-GA for " + ", ".join(sorted(P["rescued"]))
        log(f"dashboard counts include validated below-threshold genes for "
            f"{', '.join(sorted(P['rescued']))} (--ga-only to leave them out)")

    if comp is not None:
        t = comp.groupby(["family", "status"]).size().unstack(fill_value=0)
        for s_ in ("COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"):
            if s_ not in t:
                t[s_] = 0
        P["completeness"] = {f: {k: int(t.loc[f, k]) for k in
                                 ["COMPLETE","PARTIAL","FRAGMENT","MISSING_DOMAIN"]}
                             for f in t.index}
        sf = comp.groupby(["species", "family"]).status.value_counts().unstack(fill_value=0)
        for s_ in ("COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"):
            if s_ not in sf:
                sf[s_] = 0
        P["comp_sf"] = [[i[0], i[1], int(r.COMPLETE), int(r.PARTIAL),
                         int(r.FRAGMENT), int(r.MISSING_DOMAIN)]
                        for i, r in sf.iterrows()]
        # the non-complete proteins are the ones worth eyeballing
        bad = comp[comp.status != "COMPLETE"]
        P["flagged_proteins"] = bad[["species","protein","family","status","prot_len",
                                     "length_ratio","min_domain_cov","domain_coverages"]] \
            .head(4000).values.tolist()

    if arch is not None and calls is not None:
        a2 = arch.merge(calls[["species","protein","family"]].drop_duplicates(),
                        on=["species","protein"], how="left")
        ag = (a2.groupby(["family","architecture"]).size()
                .reset_index(name="n").sort_values(["family","n"], ascending=[True, False]))
        P["architectures"] = ag.values.tolist()

    P["stats"] = {}
    for f in glob.glob(os.path.join(out, "kruskal_*.tsv")):
        t = pd.read_csv(f, sep="\t")
        P["stats"][os.path.basename(f)[8:-4]] = t.where(pd.notna(t), None).to_dict("records")
    P["group_means"] = {}
    for f in glob.glob(os.path.join(out, "group_means_*.tsv")):
        t = pd.read_csv(f, sep="\t", index_col=0)
        P["group_means"][os.path.basename(f)[12:-4]] = {
            "groups": [c for c in t.columns if c != "n_groups"],
            "data": {i: [None if pd.isna(v) else float(v)
                         for k, v in r.items() if k != "n_groups"]
                     for i, r in t.iterrows()}}

    # provenance: what was found, what was not, so a reader can see the gaps
    prov = []
    for label, rel in [("hmmsearch domtblout", "hmmsearch"),
                       ("gene calls", "gene_calls_all.tsv"),
                       ("count matrix", "counts.tsv"),
                       ("normalised counts", "counts_per10k.tsv"),
                       ("proteome sizes", "proteome_sizes.tsv"),
                       ("annotation QC", "qc_annotation.tsv"),
                       ("domain table", "domains/domain_table.tsv"),
                       ("completeness", "domains/completeness.tsv"),
                       ("architectures", "domains/architectures.tsv"),
                       ("per-family FASTA", "fasta"),
                       ("gene trees", "trees"),
                       ("CYP clan calls", "cyp"),
                       ("figures", "figures")]:
        p_ = os.path.join(out, rel)
        found = os.path.exists(p_)
        n = len(glob.glob(os.path.join(p_, "*"))) if found and os.path.isdir(p_) else None
        prov.append({"item": label, "path": rel, "found": bool(found),
                     "n_files": n,
                     "modified": (__import__("datetime").datetime.fromtimestamp(
                         os.path.getmtime(p_)).isoformat(timespec="seconds")
                         if found else None)})
    P["provenance"] = prov

    clan_dir = os.path.join(out, "trees")
    P["clans"] = {}
    for f in glob.glob(os.path.join(clan_dir, "*_clan_by_species.tsv")):
        fam = os.path.basename(f).replace("_clan_by_species.tsv", "")
        t = pd.read_csv(f, sep="\t", index_col=0)
        t = t.drop(columns=[c for c in ["total"] if c in t.columns])
        P["clans"][fam] = {"clans": list(t.columns),
                           "data": {i: [int(v) for v in r] for i, r in t.iterrows()}}

    P["tax_cols"] = tax_cols
    P["samples"] = samples
    cmpd = os.path.join(out, "comparison")
    def _recs(p, **kw):
        t = _read_tsv(p, **kw)
        return None if t is None else t.where(pd.notna(t), None).to_dict("records")
    P["cmp"] = {
        "long": _recs(os.path.join(cmpd, "long_counts.tsv")),
        "recovery": _recs(os.path.join(cmpd, "recovery.tsv")),
        "tests": _recs(os.path.join(cmpd, "paired_tests.tsv")),
        "collapse": _recs(os.path.join(out, "collapse_stats.tsv")),
    }
    so = _read_tsv(os.path.join(out, "short_orf_summary.tsv"))
    P["short_orfs"] = None if so is None else so.to_dict("records")
    cc = _read_tsv(os.path.join(out, "composition_candidates.tsv"))
    P["mt_composition"] = ({} if cc is None or not len(cc)
                           else cc.groupby("species").size().astype(int).to_dict())
    # Per-tip annotations for the tree viewer. They are attached to each tree
    # (embedded or in its own data file) by write_site, never globally: at
    # catalogue scale a global map of every non-complete protein is megabytes
    # the landing page does not need.
    P["_tip_clan"] = {}
    for f in glob.glob(os.path.join(out, "cyp", "*_clan_calls.tsv")):
        t = pd.read_csv(f, sep="\t")
        for r in t.itertuples(index=False):
            P["_tip_clan"][f"{r.species}|{r.protein}"] = str(r.clan)
    P["_tip_status"] = {}
    if comp is not None:
        for r in comp[comp.status != "COMPLETE"].itertuples(index=False):
            P["_tip_status"][f"{r.species}|{r.protein}"] = str(r.status)

    P["_tree_files"] = {}
    sp_tree = os.path.join(out, "species_tree.nwk")
    if os.path.exists(sp_tree):
        P["_tree_files"]["species"] = sp_tree
    for f in sorted(glob.glob(os.path.join(out, "trees", "*.tre"))):
        P["_tree_files"][os.path.basename(f)[:-4]] = f
    return P


SITE_FORMAT = 1
TREE_EMBED_MAX = int(os.environ.get("DEFENSOME_TREE_EMBED_MAX", 300_000))  # bytes per tree in index.html
TREE_EMBED_TOTAL = 3_000_000   # and in total; anything larger loads on demand


def site_key(name, used):
    """A file-safe key for a species, unique within the site."""
    k = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("._") or "x"
    base, n = k, 2
    while k in used:
        k, n = f"{base}_{n}", n + 1
    used.add(k)
    return k


def _js_put(path, kind, key, obj):
    """A data file the page loads with a <script> tag. Script tags work from
    file://, where fetch() does not, so the site opens straight off disk."""
    import json
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(f"DEF.put({json.dumps(kind)},{json.dumps(key)},")
        json.dump(obj, fh, separators=(",", ":"))
        fh.write(");\n")
    os.replace(tmp, path)


def _tips_of(nwk):
    return set(m.group(1).strip("'\"") for m in re.finditer(r"[(,]\s*([^(),:;\[\]]+)", nwk))


def _site_trees(P, data):
    """Small trees ride in the page; large ones get a data file and load when
    opened. Tip annotations travel with their tree."""
    import json
    trees, meta, embedded = {}, {}, 0
    tip_clan, tip_status = {}, {}
    for key, path in P.pop("_tree_files", {}).items():
        nwk = open(path).read().strip()
        tips = _tips_of(nwk)
        ann_c = {t: P["_tip_clan"][t] for t in tips if t in P["_tip_clan"]}
        ann_s = {t: P["_tip_status"][t] for t in tips if t in P["_tip_status"]}
        side = path[:-4] + ".json"
        m = json.load(open(side)) if os.path.exists(side) else {}
        m.update({"n_tips": len(tips), "bytes": len(nwk)})
        m.setdefault("family", "" if key == "species" else key.split("__")[0])
        meta[key] = m
        if len(nwk) <= TREE_EMBED_MAX and embedded + len(nwk) <= TREE_EMBED_TOTAL:
            trees[key] = nwk; embedded += len(nwk)
            tip_clan.update(ann_c); tip_status.update(ann_s)
        else:
            trees[key] = None
            _js_put(os.path.join(data, "trees", re.sub(r"[^A-Za-z0-9._-]+", "_", key) + ".js"),
                    "tree", key, {"nwk": nwk, "clan": ann_c, "status": ann_s})
    P["trees"], P["tree_meta"] = trees, meta
    P["tip_clan"], P["tip_status"] = tip_clan, tip_status
    P.pop("_tip_clan", None); P.pop("_tip_status", None)


def write_site(P, out, site, proteomes=None, sequences=True, force=False):
    """Everything below the landing page: one data file per species (genes,
    domain coordinates, sequences), one per family (members), one per large
    tree, and gzipped family FASTA. Rewritten only when an input changed."""
    import json, hashlib
    data = os.path.join(site, "data")
    for d in ("species", "families", "trees", "fasta"):
        os.makedirs(os.path.join(data, d), exist_ok=True)
    species = [r["name"] for r in P["species"]]
    used = set()
    keys = [site_key(s, used) for s in species]
    for r, k in zip(P["species"], keys):
        r["key"] = k

    calls = _read_tsv(os.path.join(out, "gene_calls_all.tsv"))
    if calls is None:
        log("no gene_calls_all.tsv: species and gene pages will be summary-only")
        P["site"] = {"genes": False, "sequences": False, "family_fasta": {}}
        _site_trees(P, data)
        return
    calls = calls[calls.species.isin(set(species))]
    calls["evidence"] = "GA" if "evidence" not in calls.columns else calls["evidence"].fillna("GA")
    rc = _read_tsv(os.path.join(out, "rescue_calls.tsv"))
    if rc is not None and len(rc) and P.get("counts_basis", "GA") != "GA":
        rc = rc[(rc.validation == "PASS") & rc.species.isin(set(species))]
        if len(rc):
            calls = pd.concat([calls, pd.DataFrame({"species": rc.species, "family": rc.family,
                "protein": rc.protein, "prot_len": rc.prot_len, "hmm_cov": rc.hmm_cov,
                "evidence": "below_GA"})], ignore_index=True)
    comp = _read_tsv(os.path.join(out, "domains", "completeness.tsv"))
    status = ({} if comp is None else
              dict(zip(zip(comp.species, comp.protein, comp.family), comp.status)))
    clan = P.get("_tip_clan", {})

    # Inputs that decide the data files. If none changed and every file is
    # there, the expensive pass over the proteomes is skipped.
    dom_f = os.path.join(out, "domains", "domain_table.tsv")
    inputs = [os.path.join(out, "gene_calls_all.tsv"), dom_f,
              os.path.join(out, "domains", "completeness.tsv")]
    inputs += sorted(glob.glob(os.path.join(out, "cyp", "*_clan_calls.tsv")))
    src = {}
    if sequences and proteomes:
        src = dict(species_list(proteomes))
        inputs += [src[s] for s in species if s in src]
    elif sequences:
        inputs += sorted(glob.glob(os.path.join(out, "fasta", "*.faa")))
    sig = hashlib.sha1()
    sig.update(f"{SITE_FORMAT}|{sequences}|{bool(proteomes)}|{P.get('counts_basis')}|{','.join(keys)}".encode())
    inputs.append(os.path.join(out, "rescue_calls.tsv"))
    for f in inputs:
        if os.path.exists(f):
            st_ = os.stat(f); sig.update(f"{f}|{st_.st_size}|{int(st_.st_mtime)}".encode())
    sig = sig.hexdigest()
    man_f = os.path.join(data, "manifest.json")
    man = json.load(open(man_f)) if os.path.exists(man_f) else {}
    fams = [f["family"] for f in P["families"]]
    have = all(os.path.exists(os.path.join(data, "species", k + ".js")) for k in keys)
    fresh = not force and man.get("signature") == sig and have

    by_sp = {s: g for s, g in calls.groupby("species")}
    if not fresh:
        if os.path.exists(dom_f):
            dom = pd.read_csv(dom_f, sep="\t", usecols=["species", "protein", "pfam_name",
                              "pfam_id", "ali_from", "ali_to", "hmm_from", "hmm_to",
                              "hmm_len", "score"])
            dom_by = {s: g for s, g in dom.groupby("species")}
        else:
            dom_by = None
            log("no domains/domain_table.tsv; reading domain coordinates from the "
                "hmmsearch tables (run `domains` to make this faster)")
            scans = dict(scan_files(out)) if os.path.isdir(os.path.join(out, "hmmsearch")) else {}
        pooled = {}
        if sequences and not proteomes:
            n_ = 0
            for f in glob.glob(os.path.join(out, "fasta", "*.faa")):
                if "__" in os.path.basename(f):
                    continue
                for h, sq in read_fasta(f):
                    sp_, _, pid = h.split()[0].partition("|")
                    pooled[(sp_, pid)] = sq; n_ += 1
            if pooled:
                log(f"sequences from {out}/fasta ({len(pooled):,} proteins); "
                    "--proteomes reads them per species with less memory")
            else:
                log("no --proteomes and no fasta/*.faa: gene pages will have no sequence. "
                    "Pass --proteomes, or run `extract` first.")
        for old in glob.glob(os.path.join(data, "fasta", "*.faa.gz")):
            os.remove(old)          # never leave FASTA from an earlier export behind
        fa_out = {f: gzip.open(os.path.join(data, "fasta", f + ".faa.gz"), "wt",
                               compresslevel=6) for f in fams} if sequences else {}
        n_seq = n_miss = n_miss_below = 0
        for k_i, (s, key) in enumerate(zip(species, keys)):
            g = by_sp.get(s)
            genes, doms, seqs = [], {}, {}
            if g is not None:
                want = set(g.protein.astype(str))
                for r in g.itertuples(index=False):
                    pid = str(r.protein)
                    genes.append([pid, r.family, int(r.prot_len), round(float(r.hmm_cov), 3),
                                  status.get((s, r.protein, r.family), "COMPLETE" if comp is None else ""),
                                  clan.get(f"{s}|{pid}", ""), "" if r.evidence == "GA" else r.evidence])
                if dom_by is not None:
                    d_ = dom_by.get(s)
                else:
                    d_ = all_domains(scans[s]) if s in scans else None
                # domains of below-threshold genes are only in the rescue table
                below = set(g.protein[g.evidence != "GA"].astype(str))
                rp = os.path.join(out, "rescue", f"{s}.domtblout.gz")
                if below and os.path.exists(rp):
                    rd_ = all_domains(rp)
                    rd_ = rd_[rd_.protein.astype(str).isin(below)].copy()
                    rd_.insert(0, "species", s)
                    d_ = rd_ if d_ is None else pd.concat([d_, rd_[d_.columns]], ignore_index=True)
                if d_ is not None:
                    d_ = d_[d_.protein.astype(str).isin(want)].sort_values(["protein", "ali_from"])
                    for r in d_.itertuples(index=False):
                        doms.setdefault(str(r.protein), []).append(
                            [r.pfam_name, r.pfam_id, int(r.ali_from), int(r.ali_to),
                             int(r.hmm_from), int(r.hmm_to), int(r.hmm_len), round(float(r.score), 1)])
                if sequences:
                    if proteomes and s in src:
                        for h, sq in read_fasta(src[s]):
                            pid = h.split()[0]
                            if pid in want:
                                seqs[pid] = sq.rstrip("*")
                    elif pooled:
                        for pid in want:
                            if (s, pid) in pooled:
                                seqs[pid] = pooled[(s, pid)].rstrip("*")
                    n_seq += len(seqs)
                    miss = want - set(seqs)
                    n_miss += len(miss & set(g.protein[g.evidence == "GA"].astype(str)))
                    n_miss_below += len(miss) - len(miss & set(g.protein[g.evidence == "GA"].astype(str)))
                    for pid, fam in zip(g.protein.astype(str), g.family):
                        if pid in seqs and fam in fa_out:
                            fa_out[fam].write(f">{s}|{pid}\n{seqs[pid]}\n")
            _js_put(os.path.join(data, "species", key + ".js"), "species", key,
                    {"name": s, "genes": genes, "domains": doms, "seq": seqs})
            if (k_i + 1) % 200 == 0:
                log(f"  species data {k_i + 1:,}/{len(species):,}")
        for fh in fa_out.values():
            fh.close()
        if sequences:
            log(f"sequences: {n_seq:,} written, {n_miss:,} not found in the proteomes"
                + (" (header/ID mismatch?)" if n_miss else ""))
            if n_miss_below:
                log(f"  {n_miss_below:,} below-threshold genes have no sequence: `extract` only "
                    "writes GA calls, so pass --proteomes to include them")
        json.dump({"signature": sig, "format": SITE_FORMAT,
                   "generated": P["generated"]}, open(man_f, "w"))
    else:
        log("species data unchanged since the last export; reusing it (--force rewrites)")

    # Family member lists are cheap and index into this run's species order,
    # so they are always rewritten.
    sidx = {s: i for i, s in enumerate(species)}
    fsize = {}
    for fam, g in calls.groupby("family"):
        rows = [[sidx[s], str(pid), int(ln),
                 status.get((s, pid, fam), "COMPLETE" if comp is None else ""),
                 clan.get(f"{s}|{pid}", ""), "" if ev == "GA" else ev]
                for s, pid, ln, ev in zip(g.species, g.protein, g.prot_len, g.evidence)]
        fp = os.path.join(data, "families", re.sub(r"[^A-Za-z0-9._-]+", "_", fam) + ".js")
        _js_put(fp, "family", fam, {"members": rows})
        fsize[fam] = os.path.getsize(fp)

    _site_trees(P, data)
    fa_size = {f: os.path.getsize(os.path.join(data, "fasta", f + ".faa.gz"))
               for f in fams if sequences and os.path.exists(os.path.join(data, "fasta", f + ".faa.gz"))}
    P["site"] = {"genes": True, "sequences": bool(sequences), "family_bytes": fsize,
                 "family_fasta": fa_size}


# =========================================================================== #
#  run: everything, for a folder of proteomes or a chosen subset of it
# =========================================================================== #
def _norm_name(x):
    x = os.path.splitext(os.path.basename(str(x).strip()))[0]
    return re.sub(r"[\s_]+", "_", x).lower()


FASTA_EXT = (".faa", ".fa", ".fasta")


def list_fasta(d):
    """(stem, path) for every protein FASTA in d, gzipped or not."""
    out = []
    for f in sorted(os.listdir(d)):
        g = f[:-3] if f.endswith(".gz") else f
        ext = os.path.splitext(g)[1].lower()
        if ext in FASTA_EXT:
            out.append((os.path.splitext(g)[0], os.path.join(d, f)))
    return out


def select_proteomes(proteome_dir, wanted):
    """Proteome files matching the requested species. A name matches a file by
    its stem, or as a prefix of it, with spaces and underscores equivalent and
    case ignored: 'Danio rerio' finds Danio_rerio_pep.faa."""
    files = list_fasta(proteome_dir)
    if not wanted:
        return files, []
    keep, missing = [], []
    for w in wanted:
        n = _norm_name(w)
        # exact stem first; then a binomial or a genus followed by the file's suffix
        hit = [f for f in files if _norm_name(f[0]) == n] or \
              [f for f in files if _norm_name(f[0]).startswith(n + "_")]
        if hit:
            keep.extend(h for h in hit if h not in keep)
        else:
            missing.append(w)
    return keep, missing


def _read_species_arg(x):
    if not x:
        return []
    if os.path.exists(x):
        out = []
        for line in open(x):
            line = line.split("#")[0].strip()
            if line:
                out.append(line.split("\t")[0])
        return out
    return [y.strip() for y in x.split(",") if y.strip()]


def cmd_run(a):
    """inspect -> collapse -> HMMs -> scan -> annotate -> qc -> extract ->
    domains -> [cyp] -> [trees] -> [compare] -> report -> dashboard, for every
    proteome in a folder or the species you name. Each step is the ordinary
    subcommand, so a failure names the command to rerun by hand."""
    import time
    me = [sys.executable, os.path.abspath(__file__)]
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    failed = []

    def step(name, args, must=True):
        log(f"==== {name} ====")
        t0 = time.time()
        rc = subprocess.run(me + args).returncode
        if rc == 0:
            log(f"---- {name} OK ({time.time() - t0:.0f} s)")
            return True
        log(f"---- {name} FAILED (exit {rc})")
        if must:
            sys.exit(f"\nSTOPPED at {name}: every later step needs it.\n"
                     f"  rerun alone:  python3 {os.path.basename(__file__)} {' '.join(args)}\n"
                     f"  what is missing:  python3 {os.path.basename(__file__)} doctor")
        failed.append(name)
        return False

    # 1. which proteomes
    if not a.proteomes or not os.path.isdir(a.proteomes):
        sys.exit(f"ERROR: --proteomes must be a directory of .faa/.fa/.fasta files")
    wanted = _read_species_arg(a.species)
    files, missing = select_proteomes(a.proteomes, wanted)
    if missing:
        sys.exit(f"ERROR: no proteome matches {', '.join(missing)} in {a.proteomes}.\n"
                 f"  Available, e.g.: {', '.join(s_ for s_, _ in list_fasta(a.proteomes)[:5])}")
    if not files:
        sys.exit(f"ERROR: no .faa/.fa/.fasta files (gzipped or not) in {a.proteomes}")
    src = os.path.abspath(a.proteomes)
    if wanted or any(f.endswith(".gz") for _, f in files):
        # a subset, or gzipped files: stage plain-text proteomes in <out>/input
        src = os.path.join(out, "input")
        if os.path.isdir(src):
            shutil.rmtree(src)
        os.makedirs(src)
        for sp, f in files:
            if f.endswith(".gz"):
                with gzip.open(f, "rb") as fi, open(os.path.join(src, sp + ".faa"), "wb") as fo:
                    shutil.copyfileobj(fi, fo)
                continue
            dst = os.path.join(src, os.path.basename(f))
            try:
                os.symlink(os.path.abspath(f), dst)
            except OSError:
                shutil.copy(f, dst)
        log(f"{len(files)} of {len(list_fasta(a.proteomes))} proteomes -> {src}")
    else:
        log(f"{len(files)} proteomes in {src}")

    # 2. tools, before any work
    step("doctor", ["doctor", "--require", "core,scan"])

    # taxonomy from an NCBI taxdump, when no metadata was given
    if not a.metadata and getattr(a, "taxdump", None):
        tx = os.path.join(out, "taxonomy.tsv")
        script = os.path.join(_HERE, "scripts", "make_taxonomy_metadata.py")
        log("==== taxonomy ====")
        if subprocess.run([sys.executable, script, src, a.taxdump, tx]).returncode == 0:
            a.metadata = tx
            a.group_by = a.group_by or "class"
            log(f"---- taxonomy OK: {tx}, grouped by {a.group_by}")
        else:
            failed.append("taxonomy")

    # 3. one protein per gene
    sheet = os.path.join(out, "samplesheet.tsv")
    prot = src
    if not a.no_collapse:
        step("inspect", ["inspect", "--proteomes", src, "--out", out, "--samplesheet", sheet])
        if a.review:
            print(f"\nCheck {sheet} (species, method and the gene-ID rule), then rerun "
                  "the same command without --review.")
            return
        step("collapse", ["collapse", "--samplesheet", sheet, "--out", out]
             + (["--drop-utrorf"] if a.drop_utrorf else []))
        prot = os.path.join(out, "proteomes")

    # 4. the HMMs
    db = os.path.join(out, "db")
    if a.hmm:
        hmm = os.path.abspath(a.hmm)
        if hmm.endswith(".gz"):
            os.makedirs(db, exist_ok=True)
            plain = os.path.join(db, "defensome.hmm")
            with gzip.open(hmm, "rb") as fi, open(plain, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            hmm = plain
    elif a.pfam:
        hmm = os.path.join(db, "defensome.hmm")
        if not os.path.exists(hmm) or os.path.getmtime(a.pfam) > os.path.getmtime(hmm):
            step("setup", ["setup", "--pfam", os.path.abspath(a.pfam), "--db-dir", db, "--map", a.map])
    else:
        hmm = fetch_db(db, a.db_url)
    check_db_covers_map(hmm, a.map)

    # 5. scan and call
    T = str(a.threads)
    step("scan", ["scan", "--proteomes", prot, "--pfam", hmm, "--out", out, "--threads", T,
                  "--map", a.map])
    step("annotate", ["annotate", "--proteomes", prot, "--out", out, "--map", a.map])
    step("qc", ["qc", "--out", out, "--map", a.map], must=False)
    step("extract", ["extract", "--out", out, "--proteomes", prot, "--map", a.map], must=False)
    step("domains", ["domains", "--out", out, "--map", a.map], must=False)

    # 6. optional analyses
    if a.cyp_refs:
        step("cyp", ["cyp", "--out", out, "--refs", os.path.abspath(a.cyp_refs), "--family", "CYP",
                     "--threads", T], must=False)
    if a.trees:
        targs = ["trees", "--out", out, "--families", a.trees, "--threads", T, "--complete-only"]
        if a.metadata and a.thin:
            targs += ["--metadata", os.path.abspath(a.metadata), "--thin", a.thin]
        if a.per_species:
            targs += ["--per-species", str(a.per_species)]
        step("trees", targs, must=False)
    try:
        meth = {(r.get("method") or "").strip() for r in csv_rows(sheet)} - {""} if os.path.exists(sheet) else set()
    except Exception:
        meth = set()
    if len(meth) >= 2:
        step("compare", ["compare", "--samplesheet", sheet, "--out", out, "--complete-only"], must=False)

    # 7. figures and the dashboard
    md = ["--metadata", os.path.abspath(a.metadata)] if a.metadata else []
    gb = ["--group-by", a.group_by] if a.group_by else []
    if not a.no_figures:
        step("report", ["report", "--out", out, "--map", a.map] + (md + gb if md and gb else []), must=False)
    dargs = ["dashboard", "--out", out, "--proteomes", prot, "--map", a.map] + md + gb
    if os.path.exists(sheet):
        dargs += ["--samplesheet", sheet]
    if a.species_tree:
        dargs += ["--species-tree", os.path.abspath(a.species_tree)]
    step("dashboard", dargs, must=False)

    print("\n" + ("all steps OK" if not failed else "FAILED: " + ", ".join(failed)))
    print(f"results:    {out}")
    print(f"dashboard:  {os.path.join(out, 'dashboard', 'index.html')}")
    if failed:
        sys.exit(1)


def csv_rows(path):
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# =========================================================================== #
#  Structure and binding for a chosen gene (Boltz-2)
# =========================================================================== #
# Cofactors a family's fold needs before a substrate pose means anything: a P450
# without its haem is an empty pocket. PDB Chemical Component codes.
FAMILY_COFACTORS = {"CYP": ["HEM"], "Cyt_b5": ["HEM"], "Catalase": ["HEM"],
                    "FMO": ["FAD", "NAP"], "TrxR_GR": ["FAD"], "GST": ["GSH"],
                    "GST_MAPEG": ["GSH"], "UGT": ["UPG"], "ALDH": ["NAD"], "AKR": ["NAP"],
                    "ABC_full": ["ATP"], "ABC_NBD": ["ATP"], "Hsp70": ["ATP"], "Hsp90": ["ATP"]}


def _safe(x):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(x)).strip("._") or "x"


def structure_dir(out, sample, pid, ligand_name=""):
    stem = f"{_safe(sample)}__{_safe(pid)}" + (f"__{_safe(ligand_name)}" if ligand_name else "")
    return os.path.join(out, "structures", stem), stem


def boltz_yaml(seq, cofactors, smiles=None):
    """Boltz input, version 1: protein A, cofactors B.., the ligand last, and
    an affinity request on the ligand when there is one."""
    ids = iter("BCDEFGHIJKLMNOPQRSTUVWXYZ")
    lines = ["version: 1", "sequences:", "  - protein:", "      id: A",
             f"      sequence: {seq}"]
    for c in cofactors:
        lines += ["  - ligand:", f"      id: {next(ids)}", f"      ccd: {c}"]
    lig = None
    if smiles:
        lig = next(ids)
        lines += ["  - ligand:", f"      id: {lig}", f"      smiles: '{smiles}'"]
        lines += ["properties:", "  - affinity:", f"      binder: {lig}"]
    return "\n".join(lines) + "\n"


def cmd_structure(a):
    """Write a Boltz-2 input for a gene (with its family's cofactors and, if
    given, a ligand with an affinity request) and run it when asked."""
    import json
    need_pandas()
    calls = _read_tsv(os.path.join(a.out, "gene_calls_all.tsv"))
    if calls is None:
        sys.exit(f"ERROR: {a.out}/gene_calls_all.tsv not found. Run `annotate` first.")
    src = dict(species_list(a.proteomes)) if a.proteomes else {}
    if a.ligand and not a.ligand_name:
        sys.exit("ERROR: --ligand needs --ligand-name, which names the output folder")
    made = []
    for g in a.gene:
        sample, sep, pid = g.replace("/", "|").partition("|")
        if not sep:
            sys.exit(f"ERROR: --gene takes SAMPLE|PROTEIN, got {g!r}")
        rows = calls[(calls.species.astype(str) == sample) & (calls.protein.astype(str) == pid)]
        if not len(rows):
            sys.exit(f"ERROR: {pid} is not a defensome gene call in {sample}")
        fams = list(dict.fromkeys(rows.family))
        seq = None
        if sample in src:
            seq = next((s_ for h, s_ in read_fasta(src[sample]) if h.split()[0] == pid), None)
        if seq is None:
            for f in fams:
                faa = os.path.join(a.out, "fasta", f + ".faa")
                if os.path.exists(faa):
                    seq = next((s_ for h, s_ in read_fasta(faa) if h.split()[0] == f"{sample}|{pid}"), None)
                    if seq:
                        break
        if not seq:
            sys.exit(f"ERROR: no sequence for {sample}|{pid}. Pass --proteomes, or run `extract`.")
        seq = seq.rstrip("*").upper()
        if "U" in seq:
            log(f"{pid}: selenocysteine read as cysteine for the prediction")
            seq = seq.replace("U", "C")
        seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYX]", "X", seq)
        if a.cofactor in (None, "auto"):
            cof = list(dict.fromkeys(c for f in fams for c in FAMILY_COFACTORS.get(f, [])))
        elif a.cofactor == "none":
            cof = []
        else:
            cof = [c.strip().upper() for c in a.cofactor.split(",") if c.strip()]
        d, stem = structure_dir(a.out, sample, pid, a.ligand_name or "")
        os.makedirs(d, exist_ok=True)
        yml = os.path.join(d, stem + ".yaml")
        open(yml, "w").write(boltz_yaml(seq, cof, a.ligand))
        meta = {"sample": sample, "protein": pid, "family": fams, "length": len(seq),
                "cofactors": cof, "ligand": a.ligand or "", "ligand_name": a.ligand_name or "",
                "created": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                "tool_version": __version__}
        json.dump(meta, open(os.path.join(d, "meta.json"), "w"), indent=1)
        cmd = ["boltz", "predict", yml, "--use_msa_server", "--out_dir", d,
               "--diffusion_samples", str(a.samples), "--output_format", "mmcif"]
        if a.accelerator:
            cmd += ["--accelerator", a.accelerator]
        # the fast kernels need cuequivariance_torch, which `pip install boltz`
        # does not bring; without it the run dies at the first pairformer layer
        if not getattr(a, "kernels", False):
            cmd += ["--no_kernels"]
        with open(os.path.join(d, "run_boltz.sh"), "w") as fh:
            import shlex
            fh.write("#!/usr/bin/env bash\nset -euo pipefail\n" + " ".join(shlex.quote(c) for c in cmd) + "\n")
        os.chmod(os.path.join(d, "run_boltz.sh"), 0o755)
        log(f"{sample}|{pid}: {', '.join(fams)}, {len(seq)} aa, cofactors "
            f"{'+'.join(cof) or 'none'}" + (f", ligand {a.ligand_name}" if a.ligand else "")
            + f" -> {yml}")
        made.append((d, cmd))
    if not a.run:
        print("\nInputs written. To predict, on a GPU node:")
        for d, _ in made:
            print(f"  bash {os.path.join(d, 'run_boltz.sh')}")
        print("Then rerun `dashboard` to show the structures on the gene pages.")
        return
    if not tool("boltz"):
        sys.exit("ERROR: boltz not found. pip install boltz  (Boltz-2; a GPU is strongly "
                 "advised), or run the run_boltz.sh scripts where it is installed.")
    for d, cmd in made:
        run(cmd)
    log("structure done; rerun `dashboard` to show the predictions")


def _site_structures(P, out, data):
    """Predicted structures for the gene pages: one data file per prediction,
    and an index in the page."""
    import json
    idx = {}
    for m_f in sorted(glob.glob(os.path.join(out, "structures", "*", "meta.json"))):
        d = os.path.dirname(m_f)
        try:
            m = json.load(open(m_f))
        except ValueError:
            continue
        sid = os.path.basename(d)
        cif = sorted(glob.glob(os.path.join(d, "boltz_results_*", "predictions", "*", "*_model_0.cif")))
        conf = sorted(glob.glob(os.path.join(d, "boltz_results_*", "predictions", "*", "confidence_*_model_0.json")))
        aff = sorted(glob.glob(os.path.join(d, "boltz_results_*", "predictions", "*", "affinity_*.json")))
        e = {"id": sid, "ligand_name": m.get("ligand_name", ""), "ligand": m.get("ligand", ""),
             "cofactors": m.get("cofactors", []), "predicted": bool(cif)}
        obj = {"meta": m}
        if conf:
            c = json.load(open(conf[0]))
            obj["confidence"] = c
            e.update({k: round(float(c[k]), 3) for k in ("confidence_score", "ptm", "iptm", "complex_plddt")
                      if c.get(k) is not None})
        if aff:
            f_ = json.load(open(aff[0]))
            obj["affinity"] = f_
            e["affinity_probability"] = round(float(f_.get("affinity_probability_binary", 0)), 3)
            e["affinity_log10_ic50_uM"] = round(float(f_.get("affinity_pred_value", 0)), 3)
        if cif:
            obj["cif"] = open(cif[0]).read()
            os.makedirs(os.path.join(data, "structures"), exist_ok=True)
            _js_put(os.path.join(data, "structures", sid + ".js"), "structure", sid, obj)
        idx.setdefault(f"{m.get('sample')}|{m.get('protein')}", []).append(e)
    P["structures"] = idx
    if idx:
        n = sum(x["predicted"] for v in idx.values() for x in v)
        log(f"structures: {n} predicted, {sum(len(v) for v in idx.values()) - n} prepared but not run")



def cmd_dashboard(a):
    """A static site: index.html carries the comparative views; species, gene
    and family pages load their data from data/ on demand. No server, no CDN:
    it opens from disk, and the folder can be copied or hosted as it is."""
    need_pandas()
    import json
    P = build_payload(a.out, a.map, a.metadata, a.group_by, full=False,
                      samplesheet=a.samplesheet, ga_only=bool(getattr(a, "ga_only", False)))
    if a.species_tree and os.path.exists(a.species_tree):
        P["_tree_files"]["species"] = a.species_tree
    site = a.dashboard_out or os.path.join(a.out, "dashboard")
    if site.endswith(".html"):
        site = os.path.dirname(site) or "."
    os.makedirs(site, exist_ok=True)
    write_site(P, a.out, site, proteomes=a.proteomes,
               sequences=not getattr(a, "no_sequences", False), force=bool(a.force))
    _site_structures(P, a.out, os.path.join(site, "data"))
    P["cofactors"] = FAMILY_COFACTORS
    if not P.get("tip_clan") and any(m.get("family") == "CYP" for m in P.get("tree_meta", {}).values()):
        log("NOTE: no CYP clan calls found, so the tree clan ring will be empty. "
            "Run `cyp --refs <clan-labelled FASTA>` first.")
    P["tool_version"] = __version__
    html = (asset("dashboard.html")
            .replace("/*__JS__*/", asset("dashboard.js"))
            .replace("/*__PAYLOAD__*/null",
                     json.dumps(P, separators=(",", ":"), default=str)))
    dest = os.path.join(site, "index.html")
    with open(dest, "w") as fh:
        fh.write(html)
    old = os.path.join(a.out, "dashboard.html")
    if os.path.exists(old) and not a.dashboard_out:
        os.remove(old)
        log(f"removed the old single-file {old}; the dashboard is now {site}/")
    tot = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(site) for f in fs)
    log(f"wrote {dest} ({os.path.getsize(dest) / 1e6:.1f} MB; site {tot / 1e6:,.0f} MB)")
    print(f"\nOpen in a browser (no server needed):\n  {os.path.abspath(dest)}\n"
          f"Keep the data/ folder beside it: species, gene and family pages load from it.")


# =========================================================================== #
#  CYP clan reference building and external-dataset triage
# =========================================================================== #
def norm_clan(v):
    if v is None:
        return None
    s = str(v).strip().upper().replace("CLAN", "").replace("_", "").replace("-", "").strip()
    if s.startswith("MITO"):
        return "MITO"
    m = re.search(r"([234])", s)
    return {"2": "CYP2", "3": "CYP3", "4": "CYP4"}.get(m.group(1)) if m else None


def clan_from_name(name):
    """Fallback: infer the clade from the CYP family number using Feyereisen's
    assignment of insect families to the four clades."""
    m = re.search(r"CYP(\d+)", str(name).upper())
    if not m:
        return None
    f = int(m.group(1))
    if f in (12, 301, 302, 306, 314, 315):
        return "MITO"
    if f in (4, 311, 312, 313, 316, 318, 325, 340, 341, 349, 351, 353):
        return "CYP4"
    if f in (6, 9, 28, 308, 309, 317, 321, 324, 327, 332, 337, 345, 346, 347):
        return "CYP3"
    if f in (15, 18, 303, 304, 305, 307) or 300 <= f <= 307:
        return "CYP2"
    return None


def build_icpd_refs(fasta, table, orders, species, min_len, max_per_clan, out,
                    id_col=None, include_fragments=False):
    """Clan-labelled reference FASTA from the ICPD protein set.

    Uses a clan column from the info table when one exists, and otherwise
    infers the clade from the CYP family number. Caps sequences per clan and
    samples evenly across species, because ten near-identical CYP6 paralogues
    from one over-studied pest would otherwise pull every query toward it."""
    meta = {}
    fasta_ids = {h.split()[0] for h, _ in read_fasta(fasta)} if fasta else set()
    if table:
        need_pandas()
        t = (pd.read_excel(table) if str(table).endswith((".xlsx", ".xls"))
             else pd.read_csv(table, sep=None, engine="python"))
        print(f"info table: {t.shape[0]} rows; columns: {list(t.columns)}")
        # Pick the ID column by testing which one actually matches the FASTA
        # headers. ICPD's FASTA is keyed on "Database ID" (Abter001) while the
        # first ID-looking column is "Protein ID"
        # (Abscondita_terminalisg020.t1); guessing from the column name picks
        # the wrong one and every join silently fails.
        idc, best = None, 0
        if fasta_ids:
            for c in t.columns:
                n = len(fasta_ids & set(t[c].astype(str).str.strip()))
                if n > best:
                    idc, best = c, n
            if idc:
                print(f"id column '{idc}' matches {best}/{len(fasta_ids)} FASTA headers")
        if not idc:
            idc = next((c for c in t.columns if re.search(
                r"^(id|gene.?id|protein.?id|accession|seq.?id)", str(c), re.I)), t.columns[0])
            if fasta_ids:
                print(f"WARNING: no column matched the FASTA headers; falling back to '{idc}'.\n"
                      "         Check --id-col, or drop --fasta and build from the table's "
                      "protein sequence column.")
        spc = next((c for c in t.columns if re.search(r"species|organism", str(c), re.I)), None)
        ordc = next((c for c in t.columns if re.search(r"order", str(c), re.I)), None)
        clc = next((c for c in t.columns if re.search(r"clan|clade", str(c), re.I)), None)
        nmc = next((c for c in t.columns if re.search(r"genesymbol|description", str(c), re.I)),
                   None) or next((c for c in t.columns if re.search(r"symbol|cyp", str(c), re.I)), None)
        cmp_ = next((c for c in t.columns if re.search(r"complete", str(c), re.I)), None)
        seqc = next((c for c in t.columns if re.search(r"protein.*seq", str(c), re.I)), None)
        if id_col:
            idc = id_col
        print(f"using id='{idc}' species='{spc}' order='{ordc}' clan='{clc}' "
              f"name='{nmc}' completeness='{cmp_}' sequence='{seqc}'")
        for r in t.itertuples(index=False):
            d = dict(zip(t.columns, r))
            meta[str(d[idc]).strip()] = {
                "species": str(d.get(spc, "")).strip() if spc else "",
                "order": str(d.get(ordc, "")).strip() if ordc else "",
                "clan": norm_clan(d.get(clc)) if clc else None,
                "name": str(d.get(nmc, "")).strip() if nmc else "",
                "complete": str(d.get(cmp_, "")).strip().lower() if cmp_ else "",
                "seq": str(d.get(seqc, "")).strip() if seqc else ""}
        if cmp_:
            vc = t[cmp_].astype(str).str.strip().str.lower().value_counts()
            print(f"\ncompleteness in the table: {vc.to_dict()}")

    keep_ord = {o.strip().lower() for o in (orders or "").split(",") if o.strip()}
    keep_sp = {sp.strip().lower() for sp in (species or "").split(",") if sp.strip()}
    by_clan, no_clan, filtered, frag, bad_seq = {}, 0, 0, 0, 0
    if fasta:
        src = ((h.split()[0], s_) for h, s_ in read_fasta(fasta))
    else:
        src = ((k, v.get("seq", "")) for k, v in meta.items())
        print("no --fasta given: building from the table's protein sequence column")
    for sid, seq in src:
        m = meta.get(sid) or {"species": sid, "order": sid, "clan": None, "name": sid}
        # ICPD marks a large fraction of its predictions "fragmented". A
        # truncated P450 in a clan reference degrades the profile it goes into,
        # and these are exactly the models that lack the diagnostic motifs.
        # ICPD marks entries "full-length", "fragmented" or "artificially
        # fused". Fragments lack the diagnostic motifs; fused models are two
        # adjacent genes merged into one chimera, which is worse. Neither
        # belongs in a profile HMM.
        cmpl = m.get("complete", "")
        if not include_fragments and (cmpl.startswith("frag") or "fused" in cmpl):
            frag += 1; continue
        # ---- BUG 3: some ICPD "full-length" sequences carry internal stop
        # codons. cd-hit silently discards them; hmmbuild would keep them and
        # corrupt the alignment column statistics.
        core = seq.rstrip("*")
        if "*" in core or "X" * 10 in core.upper():
            bad_seq += 1; continue
        if keep_ord and not any(o in str(m.get("order", "")).lower()
                                or o in str(m.get("species", "")).lower() for o in keep_ord):
            filtered += 1; continue
        if keep_sp and not any(x in str(m.get("species", "")).lower() for x in keep_sp):
            filtered += 1; continue
        if len(seq) < min_len:
            filtered += 1; continue
        hm = CLAN_RE.search(sid)
        clan = (m.get("clan") or (norm_clan(hm.group(0)) if hm else None)
                or clan_from_name(m.get("name") or m.get("description") or sid))
        if not clan:
            no_clan += 1; continue
        by_clan.setdefault(clan, []).append((sid, m.get("species", ""), m.get("name", ""), seq))

    if not by_clan:
        sys.exit("ERROR: nothing could be assigned to a clan. Check the column "
                 "detection printed above, or supply a table with a clan column.")
    d = os.path.dirname(os.path.abspath(out))
    if d:
        os.makedirs(d, exist_ok=True)
    n = 0
    with open(out, "w") as fh:
        for clan, items in sorted(by_clan.items()):
            items.sort(key=lambda x: (x[1], x[0]))
            step = max(1, len(items) // max(max_per_clan, 1))
            for sid, sp_, nm, seq in items[::step][:max_per_clan]:
                fh.write(f">{sid} clan={clan} species={sp_.replace(' ', '_')} name={nm}\n"
                         f"{seq.rstrip('*')}\n")
                n += 1
    print("\n=== reference composition ===")
    for clan, items in sorted(by_clan.items()):
        print(f"  {clan:6s} available {len(items):6d}  written {min(len(items), max_per_clan):5d}")
    print(f"\nskipped: {filtered} filtered out, {frag} fragmented or fused, "
          f"{bad_seq} with internal stops or long X runs, "
          f"{no_clan} with no assignable clan")
    log(f"wrote {n} sequences to {out}")
    print("\nNext, reduce redundancy so one over-sampled pest cannot dominate:")
    print(f"  cd-hit -i {out} -o {out}.c80.faa -c 0.8 -n 5")
    print(f"  python3 defensome.py cyp --out <results> --refs {out}.c80.faa")


def cmd_cyp_refs(a):
    """Build a clan-labelled CYP reference FASTA for `cyp --refs`.

    --source flybase  uses the Drosophila set you already have. Fast, but
                      ~350 Myr from Lepidoptera, so expect LOW-confidence calls.
    --source icpd     uses the Insect Cytochrome P450 Database (66,477 P450s
                      from 682 insect species). Restrict with --orders
                      Lepidoptera and the clan calls get much better.
    """
    if a.source == "flybase":
        if not a.fasta or not os.path.exists(a.fasta):
            sys.exit("ERROR: --fasta is required for --source flybase\n"
                     "  e.g. db/flybase/flybase_dmel_clan_ref.canonical.faa")
        n = 0
        d = os.path.dirname(os.path.abspath(a.out))
        if d:
            os.makedirs(d, exist_ok=True)
        seen = {}
        with open(a.out, "w") as fh:
            for hdr, seq in read_fasta(a.fasta):
                m = re.search(r"clan=(\S+)", hdr)
                clan = norm_clan(m.group(1)) if m else clan_from_name(hdr)
                if not clan or len(seq) < a.min_len:
                    continue
                seen[clan] = seen.get(clan, 0) + 1
                fh.write(f">{hdr.split()[0]} clan={clan} species=Drosophila_melanogaster\n{seq}\n")
                n += 1
        if not n:
            sys.exit(f"ERROR: no clan= tags or recognisable CYP names in {a.fasta}")
        print("=== reference composition ===")
        for k, v in sorted(seen.items()):
            print(f"  {k:6s} {v:5d}")
        log(f"wrote {n} sequences to {a.out}")
        print("\nThis is a Drosophila-only reference. Marginal clan calls will be "
              "flagged LOW.\nFor Lepidoptera, download ICPD and rerun with "
              "--source icpd (see `dataset-help`).")
        return

    if a.fasta and not os.path.exists(a.fasta):
        sys.exit(f"ERROR: --fasta not found: {a.fasta}\n"
                 "  python3 defensome.py dataset-help")
    if not a.fasta and not a.table:
        sys.exit("ERROR: give --fasta, --table, or both.\n"
                 "  python3 defensome.py dataset-help")
    build_icpd_refs(a.fasta, a.table, a.orders, a.species, a.min_len,
                    a.max_per_clan, a.out, a.id_col, a.include_fragments)


INSECT_HINTS = ["insect", "arthropod", "drosophila", "bombyx", "helicoverpa",
                "spodoptera", "plutella", "anopheles", "aedes", "culex",
                "tribolium", "apis", "myzus", "nilaparvata", "papilio",
                "musca", "leptinotarsa", "bemisia", "tetranychus", "acyrthosiphon",
                "animal", "lepidoptera", "diptera", "coleoptera", "hymenoptera"]


def triage_p450rdb(indir, out):
    os.makedirs(out, exist_ok=True)
    def find(pat):
        for f in sorted(os.listdir(indir)):
            if re.search(pat, f, re.I):
                return os.path.join(indir, f)
        return None
    fr = find(r"Reactions")
    if not fr:
        sys.exit(f"ERROR: no Reactions CSV in {indir}")
    rx = pd.read_csv(fr, low_memory=False)
    print(f"Reactions: {len(rx):,} rows x {rx.shape[1]} columns")
    if "Species" in rx.columns:
        k = rx.Species.fillna("unknown").value_counts()
        print("\n=== Rows by kingdom ==="); print(k.to_string())
        k.to_csv(os.path.join(out, "rows_by_kingdom.tsv"), sep="\t")
    sp = rx["Species name"].fillna("unknown").value_counts()
    print(f"\n=== Distinct source organisms: {len(sp)} ==="); print(sp.head(15).to_string())
    sp.to_csv(os.path.join(out, "rows_by_species.tsv"), sep="\t")
    blob = (rx["Species name"].fillna("") + " " +
            rx.get("Species", pd.Series([""] * len(rx))).fillna("")).str.lower()
    ins = rx[blob.apply(lambda s: any(h in s for h in INSECT_HINTS))]
    print(f"\n=== Rows matching an insect/arthropod hint: {len(ins):,} of {len(rx):,} "
          f"({len(ins)/max(len(rx),1)*100:.1f}%) ===")
    if len(ins):
        print(ins["Species name"].value_counts().head(12).to_string())
        ins.to_csv(os.path.join(out, "insect_rows.tsv"), sep="\t", index=False)
    else:
        print("NONE. P450RDB cannot supply insect enzyme-substrate pairs.\n"
              "Use it for host-plant chemistry; get insect P450 data from ICPD.")
    m = melt_reactions(rx)
    m = m[m.sub_SMILES.notna() & ~m.sub_SMILES.astype(str).str.strip().isin(["", "/"])]
    if "sequence" in m.columns:
        m = m[m.sequence.notna() & (m.sequence.astype(str).str.len() > 100)]
    m.to_csv(os.path.join(out, "enzyme_substrate_pairs.tsv"), sep="\t", index=False)
    print(f"\n=== Pairs with both a sequence and a SMILES: {len(m):,} ===")
    print(f"    distinct enzymes: {m['Uniprot ID'].nunique() if 'Uniprot ID' in m else '?'}"
          f"  distinct substrates: {m.substrate.nunique()}")
    print("    POSITIVES ONLY. There are no curated non-substrates here, and")
    print("    absence of a record is not evidence of a non-substrate.")
    if "Transformations" in rx.columns:
        t = rx.Transformations.fillna("unknown").value_counts()
        print("\n=== Transformation types ==="); print(t.head(12).to_string())
        t.to_csv(os.path.join(out, "transformation_types.tsv"), sep="\t")
    if "Species" in rx.columns:
        pl = melt_reactions(rx[rx.Species.fillna("").str.lower().eq("plant")])
        cols = [c for c in ["Species name", "Symbol", "substrate", "sub_CID",
                            "sub_SMILES", "Transformations", "PMID"] if c in pl.columns]
        pl[cols].drop_duplicates().to_csv(
            os.path.join(out, "plant_metabolite_inventory.tsv"), sep="\t", index=False)
        print(f"\n=== Plant metabolite inventory: {pl.substrate.nunique():,} compounds "
              f"from {pl['Species name'].nunique()} plants ===")
        print("    Join this to your host-plant table: it is the chemistry your")
        print("    herbivores actually face, with structures attached.")
    log(f"wrote {out}/")


def melt_reactions(rx):
    """One row per (enzyme, substrate). P450RDB stores up to four substrates
    per reaction in parallel column blocks."""
    keep = ["Symbol", "Name", "Uniprot ID", "EC number", "Species name",
            "Species", "Txid", "Transformations", "bond", "PMID", "DOI", "sequence"]
    keep = [c for c in keep if c in rx.columns]
    out = []
    for i in (1, 2, 3, 4):
        sc, cc = f"Substrate{i}", f"sub_CID{i}"
        sm = f"sub_Smiles{i}"
        if sc not in rx.columns:
            continue
        blk = rx[keep + [c for c in (sc, cc, sm) if c in rx.columns]].copy()
        blk = blk.rename(columns={sc: "substrate", cc: "sub_CID", sm: "sub_SMILES"})
        blk["substrate_slot"] = i
        out.append(blk)
    m = pd.concat(out, ignore_index=True)
    m = m[m.substrate.notna() & (m.substrate.astype(str).str.strip() != "")]
    # cofactors and co-substrates are not the chemistry you are modelling
    cof = {"o2", "h2o", "h+", "nadph", "nadp+", "nadh", "nad+", "fmn", "fmnh2",
           "fad", "fadh2", "co2", "glutathione", "electron"}
    m = m[~m.substrate.astype(str).str.strip().str.lower().isin(cof)]
    return m


def cmd_chem_triage(a):
    """Triage P450RDB before committing to it."""
    need_pandas()
    triage_p450rdb(a.dir, a.out)


# --------------------------------------------------------------------------- #
#  ICPD reference-table evidence triage
# --------------------------------------------------------------------------- #
# Evidence tiers, strongest first. A gene inherits the strongest tier of any
# paper that mentions it. The distinction is not cosmetic: overexpression in a
# resistant strain says nothing about which compound an enzyme turns over, and
# a substrate label built from correlative records will not survive review.
EVIDENCE_TIERS = [
    ("FUNCTIONAL_METABOLISM", [
        r"baculovirus", r"heterologously\s+express", r"recombinant\s+(protein|enzyme|P450|CYP)",
        r"metaboli[sz]ed?\s+(by|the|six|two|imperatorin|xanthotoxin)", r"metabolic\s+assay",
        r"pmol/min", r"\bin\s+vitro\s+metabolis", r"\bSf9\b", r"expressed\s+in\s+(E\.\s?coli|yeast|insect cells)",
        r"catalyz\w+\s+the", r"enzyme\s+converts", r"substrate\s+specificit",
        r"functional\s+characteri[sz]", r"enzyme\s+activit", r"kinetic\s+parameter",
        r"\bK\s?m\b", r"\bkcat\b", r"purified\s+(protein|enzyme)",
        r"microsome[s]?\s+(from|of)", r"turnover", r"depletion\s+assay",
        r"metaboli[sz]ing\s+activity", r"oxidat\w+\s+of\s+\w+"]),
    ("GENETIC_PERTURBATION", [
        r"RNA\s?i\b", r"RNA\s+interference", r"knock-?down", r"silenc", r"CRISPR",
        r"transgenic", r"piggyBac", r"mutant", r"loss[- ]of[- ]function",
        r"gene\s+editing", r"dsRNA", r"antisense"]),
    ("CORRELATIVE_EXPRESSION", [
        r"over-?express", r"up-?regulat", r"induc", r"expression\s+(profile|level|analysis)",
        r"qRT-?PCR", r"resistant\s+strain", r"microarray", r"transcriptom"]),
]
# Compounds worth pulling out of abstracts: these become the candidate ligand
# list for anything downstream.
# The functional insect-P450 literature is dominated by insecticide resistance
# in crop pests. In ICPD the insecticide:allelochemical ratio of gene-records is
# roughly 12:1. If the question is host-plant chemistry in non-pest herbivores,
# that is a second domain shift on top of any human-to-insect one, and it has to
# be stated rather than absorbed silently.
COMPOUND_CLASS = {
    "insecticide": ["imidacloprid", "thiamethoxam", "permethrin", "deltamethrin",
                    "cypermethrin", "cyhalothrin", "chlorantraniliprole", "phoxim",
                    "abamectin", "indoxacarb", "metaflumizone", "diazinon", "DDT",
                    "aldrin", "pyrethroid", "neonicotinoid", "organophosphate"],
    "synergist": ["piperonyl butoxide"],
    "allelochemical": ["imperatorin", "xanthotoxin", "furanocoumarin", "quercetin",
                       "flavone", "chlorogenic acid", "indole-3-carbinol", "rutin",
                       "tannin", "gossypol", "caffeine", "nicotine", "capsaicin",
                       "myristicin", "aflatoxin"],
    "endogenous": ["juvenile hormone", "ecdysone", "20-hydroxyecdysone",
                   "methyl farnesoate", "theobromine"],
}
COMPOUND_HINTS = [
    "imperatorin", "xanthotoxin", "furanocoumarin", "quercetin", "flavone",
    "chlorogenic acid", "indole-3-carbinol", "rutin", "tannin", "gossypol",
    "caffeine", "theobromine", "nicotine", "capsaicin", "myristicin",
    "piperonyl butoxide", "diazinon", "cypermethrin", "aldrin", "permethrin",
    "deltamethrin", "cyhalothrin", "phoxim", "imidacloprid", "thiamethoxam",
    "chlorantraniliprole", "abamectin", "indoxacarb", "metaflumizone",
    "DDT", "pyrethroid", "neonicotinoid", "organophosphate", "aflatoxin",
    "juvenile hormone", "ecdysone", "20-hydroxyecdysone", "methyl farnesoate",
]
GENE_SPLIT = re.compile(r"[;,]|\s{2,}")
CYP_NAME = re.compile(r"\bCYP\s?[0-9]+[A-Z]*[0-9]*(?:[-_v][0-9A-Za-z]+)?\b", re.I)


def cmd_icpd_evidence(a):
    """Sort ICPD reference records into evidence tiers.

    The References tables list a paper per row and often dozens of genes per
    paper. This explodes them to one row per (gene, paper), classifies each
    paper by the strongest evidence its abstract describes, and reports how
    many genes survive at each tier. Only FUNCTIONAL_METABOLISM records carry
    a defensible substrate label."""
    need_pandas()
    # Detect columns PER FILE. References1 names its gene column "Gene name or
    # symbol" and References2 names it "Gene symbol"; concatenating first and
    # detecting once picks whichever file came first and silently discards
    # every row of the other.
    tables = []
    for f in a.tables:
        if not os.path.exists(f):
            sys.exit(f"ERROR: not found: {f}")
        t = (pd.read_excel(f) if f.endswith((".xlsx", ".xls"))
             else pd.read_csv(f, sep=None, engine="python"))
        t["__source"] = os.path.basename(f)
        cols = {
            "gene": next((c for c in t.columns if re.search(r"gene", str(c), re.I)), None),
            "abstract": next((c for c in t.columns if re.search(r"abstract", str(c), re.I)), None),
            "title": next((c for c in t.columns if re.search(r"title", str(c), re.I)), None),
            "doi": next((c for c in t.columns if re.search(r"doi|reference", str(c), re.I)), None),
            "phenotype": next((c for c in t.columns if re.search(r"phenotype", str(c), re.I)), None),
            "effect": next((c for c in t.columns if re.search(r"variant.*effect|effect", str(c), re.I)), None),
            "species": next((c for c in t.columns if re.search(r"species", str(c), re.I)), None),
            "desc": next((c for c in t.columns if re.search(r"description", str(c), re.I)), None),
        }
        print(f"{os.path.basename(f)}: {t.shape[0]} rows")
        print("   " + "  ".join(f"{k}='{v}'" for k, v in cols.items() if v))
        if not cols["gene"]:
            print(f"   WARNING: no gene column; skipping this file")
            continue
        tables.append((t, cols))
    if not tables:
        sys.exit("ERROR: no usable table")

    rows = []
    for t, cols in tables:
      gcol, acol, tcol = cols["gene"], cols["abstract"], cols["title"]
      dcol, pcol = cols["doi"], cols["phenotype"]
      for r in t.itertuples(index=False):
        d = dict(zip(t.columns, r))
        text = " ".join(str(d.get(c, "")) for c in
                        (acol, tcol, cols["phenotype"], cols["effect"], cols["desc"]) if c)
        tier, hits = "UNCLASSIFIED", []
        for name, pats in EVIDENCE_TIERS:
            h = [p for p in pats if re.search(p, text, re.I)]
            if h:
                tier, hits = name, h
                break
        cpds = sorted({c for c in COMPOUND_HINTS if re.search(re.escape(c), text, re.I)})
        klass = sorted({k for k, v in COMPOUND_CLASS.items() if any(c in v for c in cpds)})
        genes = set()
        for chunk in GENE_SPLIT.split(str(d.get(gcol, ""))):
            genes.update(m.group(0).replace(" ", "").upper() for m in CYP_NAME.finditer(chunk))
            for w in chunk.split():
                w = w.strip("()[]").strip()
                if w and not CYP_NAME.match(w) and re.match(r"^[A-Za-z][A-Za-z0-9]{3,}$", w) \
                        and w.lower() not in ("none", "gene", "name", "symbol"):
                    genes.add(w)
        for g in sorted(genes):
            rows.append({"gene": g, "evidence_tier": tier,
                         "n_genes_in_paper": len(genes),
                         # A compound named in a paper covering 45 genes is
                         # associated with the PAPER, not with any one gene.
                         # Only single- or two-gene papers give a usable
                         # gene-compound pair without reading the paper.
                         "compound_confidence": ("HIGH" if len(genes) <= 2 else
                                                 "MEDIUM" if len(genes) <= 5 else "LOW"),
                         "species": str(d.get(cols["species"], "")) if cols["species"] else "",
                         "compounds_mentioned": "; ".join(cpds),
                         "compound_classes": "; ".join(klass),
                         "matched_patterns": "; ".join(hits[:3]),
                         "title": str(d.get(tcol, ""))[:200] if tcol else "",
                         "doi": str(d.get(dcol, "")) if dcol else "",
                         "phenotype": str(d.get(pcol, "")) if pcol else "",
                         "source_file": d.get("__source", "")})
      # end per-table loop
    ev = pd.DataFrame(rows)
    if ev.empty:
        sys.exit("ERROR: no gene names could be parsed")
    os.makedirs(a.out, exist_ok=True)
    ev.to_csv(os.path.join(a.out, "gene_paper_evidence.tsv"), sep="\t", index=False)

    order = {n: i for i, (n, _) in enumerate(EVIDENCE_TIERS)}
    order["UNCLASSIFIED"] = 99
    best = (ev.assign(_o=ev.evidence_tier.map(order))
              .sort_values("_o").groupby("gene").first().reset_index()
              .drop(columns=["_o"]))
    best.to_csv(os.path.join(a.out, "gene_best_evidence.tsv"), sep="\t", index=False)

    print("\n=== Papers by strongest evidence described ===")
    print(ev.drop_duplicates(["doi", "title"]).evidence_tier.value_counts().to_string())
    print("\n=== Distinct genes by strongest evidence anywhere ===")
    print(best.evidence_tier.value_counts().to_string())

    fn = best[best.evidence_tier == "FUNCTIONAL_METABOLISM"]
    print(f"\n=== {len(fn)} genes with a functional metabolism record ===")
    gold = fn[(fn.compound_confidence == "HIGH") & (fn.compounds_mentioned != "")]
    print(f"    of which {len(gold)} come from a paper covering at most two genes,")
    print(f"    so the gene-compound pairing can be trusted without reading it.")
    if len(gold):
        print("\n=== The set worth using as an external test set ===")
        print(gold[["gene", "species", "compounds_mentioned"]].head(30).to_string(index=False))
        gold.to_csv(os.path.join(a.out, "gene_compound_highconfidence.tsv"),
                    sep="\t", index=False)
    print("\nCompounds are matched at PAPER level. A paper on 45 genes naming five")
    print("compounds yields 225 gene-compound pairs and almost all are spurious:")
    print("that is why compound_confidence exists. Anything below HIGH needs the")
    print("paper read before it can carry a substrate label.")
    print("\nEverything else is correlative or genetic. Overexpression in a")
    print("resistant strain is not evidence that an enzyme metabolises a compound,")
    print("and RNAi knockdown is genetic evidence of involvement, not of substrate")
    print("specificity. Only the functional set can carry a substrate label.")

    big = ev.drop_duplicates(["doi", "title"]).nlargest(5, "n_genes_in_paper")
    if len(big):
        print("\n=== Papers listing the most genes (these are almost always expression screens) ===")
        print(big[["n_genes_in_paper", "evidence_tier", "title"]].to_string(index=False))

    kc = {}
    for ks in ev.compound_classes.dropna():
        for k in [x.strip() for x in ks.split(";") if x.strip()]:
            kc[k] = kc.get(k, 0) + 1
    if kc:
        print("\n=== What kind of chemistry has the literature actually tested? ===")
        for k, v in sorted(kc.items(), key=lambda x: -x[1]):
            print(f"  {k:16s} {v:6d} gene-records")
        if kc.get("insecticide", 0) > 3 * kc.get("allelochemical", 1):
            print("\n  The functional record is dominated by synthetic insecticides in")
            print("  crop pests. If your question is host-plant allelochemistry in")
            print("  non-pest herbivores, that is a domain shift you must report, not")
            print("  a set you can train on and quietly generalise from.")

    cnt = {}
    for cs in ev.compounds_mentioned.dropna():
        for c in [x.strip() for x in cs.split(";") if x.strip()]:
            cnt[c] = cnt.get(c, 0) + 1
    if cnt:
        cp = pd.DataFrame(sorted(cnt.items(), key=lambda x: -x[1]),
                          columns=["compound", "n_gene_records"])
        cp.to_csv(os.path.join(a.out, "compounds_mentioned.tsv"), sep="\t", index=False)
        print(f"\n=== Candidate compounds named across the abstracts: {len(cp)} ===")
        print(cp.head(15).to_string(index=False))
    log(f"wrote {a.out}/")


def cmd_cyp_benchmark(a):
    """Measure clan-assignment accuracy against ICPD's own curated clade calls.

    A large bitscore margin says an assignment is unambiguous given the four
    HMMs. It says nothing about whether the assignment is right. ICPD carries a
    curated Clade for every one of its 66,477 P450s, and the reference used only
    a few hundred per clan, so the remainder is a large held-out labelled test
    set that costs nothing to use.

    Sequences that went into the reference are excluded, so this is a genuine
    hold-out and not an accuracy estimate on the training data."""
    need_pandas()
    for exe in ("hmmsearch",):
        if not shutil.which(exe):
            sys.exit(f"ERROR: {exe} not found. module load HMMER")
    hmms = sorted(glob.glob(os.path.join(a.clan_hmms, "*.hmm")))
    if not hmms:
        sys.exit(f"ERROR: no .hmm files in {a.clan_hmms}\n"
                 "  build them first: defensome.py cyp --out <results> --refs <ref.faa>")
    os.makedirs(a.out, exist_ok=True)

    t = (pd.read_excel(a.table) if str(a.table).endswith((".xlsx", ".xls"))
         else pd.read_csv(a.table, sep=None, engine="python"))
    fasta_ids = {h.split()[0] for h, _ in read_fasta(a.fasta)}
    idc, best = None, 0
    for c in t.columns:
        n = len(fasta_ids & set(t[c].astype(str).str.strip()))
        if n > best:
            idc, best = c, n
    clc = next((c for c in t.columns if re.search(r"clan|clade", str(c), re.I)), None)
    ordc = next((c for c in t.columns if re.search(r"order", str(c), re.I)), None)
    cmp_ = next((c for c in t.columns if re.search(r"complete", str(c), re.I)), None)
    spc = next((c for c in t.columns if re.search(r"species", str(c), re.I)), None)
    if not (idc and clc):
        sys.exit(f"ERROR: need an ID and a clade column. Columns: {list(t.columns)}")
    log(f"truth from '{clc}', keyed on '{idc}' ({best}/{len(fasta_ids)} matched)")

    truth, meta = {}, {}
    keep_ord = {o.strip().lower() for o in (a.orders or "").split(",") if o.strip()}
    for r in t.itertuples(index=False):
        d = dict(zip(t.columns, r))
        sid = str(d[idc]).strip()
        cl = norm_clan(d.get(clc))
        if not cl:
            continue
        if keep_ord and ordc and not any(o in str(d.get(ordc, "")).lower() for o in keep_ord):
            continue
        if cmp_ and not str(d.get(cmp_, "")).strip().lower().startswith("full"):
            continue
        truth[sid] = cl
        meta[sid] = {"species": str(d.get(spc, "")) if spc else "",
                     "order": str(d.get(ordc, "")) if ordc else ""}

    used = set()
    if a.refs and os.path.exists(a.refs):
        used = {h.split()[0] for h, _ in read_fasta(a.refs)}
        log(f"excluding {len(used)} sequences that went into the reference")

    test = os.path.join(a.out, "heldout.faa")
    n_test, tested_ids, n_short, n_stop = 0, set(), 0, 0
    with open(test, "w") as fh:
        for hdr, seq in read_fasta(a.fasta):
            sid = hdr.split()[0]
            if sid not in truth or sid in used:
                continue
            if len(seq) < a.min_len:
                n_short += 1; continue
            if "*" in seq.rstrip("*"):
                n_stop += 1; continue
            fh.write(f">{sid}\n{seq.rstrip('*')}\n")
            tested_ids.add(sid); n_test += 1
    if not n_test:
        sys.exit("ERROR: held-out set is empty")
    log(f"held-out labelled sequences: {n_test} tested "
        f"({n_short} below --min-len {a.min_len}, {n_stop} with internal stops, excluded)")

    scores = {}
    for hmm in hmms:
        clan = os.path.basename(hmm)[:-4]
        tbl = os.path.join(a.out, f"hits_{clan}.tbl")
        run(["hmmsearch", "--noali", "-E", "1e-5", "--cpu", str(a.threads),
             "--tblout", tbl, "-o", os.devnull, hmm, test])
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                scores.setdefault(f[0], {})[clan] = float(f[5])

    rows = []
    for sid in truth:
        if sid in used or sid not in scores:
            continue
        sc = sorted(scores[sid].items(), key=lambda x: -x[1])
        pred, top = sc[0]
        margin = top - sc[1][1] if len(sc) > 1 else top
        rows.append((sid, meta[sid]["species"], truth[sid], pred,
                     round(top, 1), round(margin, 1), truth[sid] == pred))
    res = pd.DataFrame(rows, columns=["id", "species", "icpd_clade", "predicted",
                                      "score", "margin", "correct"])
    res.to_csv(os.path.join(a.out, "benchmark_calls.tsv"), sep="\t", index=False)
    unscored = len(tested_ids - set(scores))

    acc = res.correct.mean()
    print(f"\n=== Hold-out accuracy against ICPD's curated clade ===")
    print(f"  tested {len(res)} sequences; {unscored} of the tested set matched "
          f"no clan HMM at E<1e-5")
    print(f"  accuracy {acc*100:.2f}%  ({int(res.correct.sum())} correct, "
          f"{int((~res.correct).sum())} wrong)")
    cm = pd.crosstab(res.icpd_clade, res.predicted)
    print("\n=== Confusion matrix (rows = ICPD, columns = predicted) ===")
    print(cm.to_string())
    cm.to_csv(os.path.join(a.out, "confusion_matrix.tsv"), sep="\t")
    print("\n=== Per-clan recall and precision ===")
    for c in sorted(set(res.icpd_clade) | set(res.predicted)):
        tp = int(((res.icpd_clade == c) & (res.predicted == c)).sum())
        fn_ = int(((res.icpd_clade == c) & (res.predicted != c)).sum())
        fp = int(((res.icpd_clade != c) & (res.predicted == c)).sum())
        rec = tp / (tp + fn_) if tp + fn_ else float("nan")
        pre = tp / (tp + fp) if tp + fp else float("nan")
        print(f"  {c:6s} recall {rec*100:6.2f}%  precision {pre*100:6.2f}%  (n={tp+fn_})")
    print("\n=== Margin distribution, correct vs incorrect ===")
    for lbl, sub in (("correct", res[res.correct]), ("incorrect", res[~res.correct])):
        if len(sub):
            print(f"  {lbl:9s} n={len(sub):6d}  min {sub.margin.min():8.1f}  "
                  f"median {sub.margin.median():8.1f}  max {sub.margin.max():8.1f}")
    print("\n=== Accuracy against retention, by margin threshold ===")
    print("  threshold   kept        accuracy   errors kept")
    for th in (0, 10, 20, 50, 100, 117, 200, 300):
        k = res[res.margin >= th]
        if not len(k):
            continue
        print(f"  {th:9d}  {len(k):6d} ({len(k)/len(res)*100:5.1f}%)  "
              f"{k.correct.mean()*100:7.3f}%  {int((~k.correct).sum()):5d}")
    print("\n  Pick a threshold here rather than assuming one. If the minimum")
    print("  margin in your own dataset sits above the point where errors")
    print("  disappear, that is quantitative grounds for trusting every call.")

    print("\n=== Caveat that belongs in the results ===")
    print("  This compares two automated assignments. ICPD's clades come from")
    print("  its own pipeline, not from manual curation of every gene, so a")
    print("  disagreement identifies a gene where two methods differ, not a")
    print("  gene where yours is wrong. Report it as agreement, not accuracy.")

    worst = None
    try:
        per = [(c, ((res.icpd_clade == c) & (res.predicted == c)).sum() /
                max(int((res.icpd_clade == c).sum()), 1))
               for c in sorted(set(res.icpd_clade))]
        worst = min(per, key=lambda x: x[1])
    except Exception:
        pass
    if worst and worst[1] < 0.999:
        print(f"\n  Weakest clan: {worst[0]} at {worst[1]*100:.2f}% recall. If that is")
        print(f"  also your smallest reference, rerun cyp-refs with a larger")
        print(f"  --max-per-clan and rebuild; the two are usually connected.")

    wrong = res[~res.correct]
    if len(wrong):
        wrong.to_csv(os.path.join(a.out, "misassigned.tsv"), sep="\t", index=False)
        print(f"\n=== The {len(wrong)} disagreements ===")
        print(wrong.nlargest(min(15, len(wrong)), "margin")[
            ["id", "species", "icpd_clade", "predicted", "margin"]].to_string(index=False))
        print("\nHigh-margin disagreements are the interesting ones: either ICPD's")
        print("curated clade is wrong for that gene, or the clan boundary is not")
        print("where one of you thinks it is. Neither is a pipeline bug, and both")
        print("are worth a sentence in the results.")
    log(f"wrote {a.out}/")


# =========================================================================== #
#  Sample sheets, isoform collapsing, and method comparison
# =========================================================================== #
# A transcriptome-derived proteome is not comparable to a genome annotation
# until isoforms are collapsed. Trinity/EvidentialGene IDs carry an explicit
# isoform field (c0_g1_i21 is the 21st isoform of one gene), so counting
# proteins counts isoforms. Genome annotations vary: some headers carry a
# gene: tag, some use .t1 transcript suffixes, some are already gene-level.
# Getting this wrong inflates one arm of a comparison and invalidates it.
GENE_PATTERNS = [
    ("trinity_isoform", r"_i\d+\w*$",
     "Trinity/EvidentialGene isoform suffix (c0_g1_i1)"),
    ("transcript_t", r"\.t\d+$", "AUGUSTUS/BRAKER transcript suffix (.t1)"),
    ("transcript_p", r"\.p\d+$", "TransDecoder ORF suffix (.p1)"),
    ("transcript_dash_R", r"-R[A-Z]+\d*$", "FlyBase-style transcript suffix (-RA)"),
    ("uniprot_isoform", r"^((?:sp|tr)\|[A-Z0-9]+)-\d+(\|.*)?$",
     "UniProt accession with an isoform suffix (sp|P12345-2|)"),
    ("version_dot", r"\.\d+$", "trailing version number, usually NOT an isoform"),
]
# Ensembl writes "gene:ENSG000..."; UniProt writes "GN=xyz" in the same slot.
# Without the second pattern a UniProt proteome collapses not at all: zebrafish
# UP000000437 is 68,855 proteins for roughly 26,000 genes.
GENE_TAG = re.compile(r"\bgene[:=](\S+)|\bGN=(\S+)")


def gene_of(header, regex=None):
    """Gene ID for a protein header. A gene: tag in the description always
    wins; otherwise the supplied regex is stripped from the first token."""
    m = GENE_TAG.search(header)
    if m:
        return m.group(1) or m.group(2)
    tok = header.split()[0]
    return re.sub(regex, "", tok) if regex else tok


# "DToL" is not one annotation method. In practice a DToL release is BRAKER
# for some species, Ensembl for others, and occasionally raw AUGUSTUS. They
# differ systematically: on one ten-species set, BRAKER gave a median of about
# 18,900 genes and Ensembl about 14,000 after collapsing. Lumping them makes the
# reference arm of any comparison heterogeneous without anyone noticing.
_PIPELINE_HINTS = [("BRAKER", r"^BRAKER"), ("UniProt", r"^(?:sp|tr)\|"),
                   ("Ensembl", r"^ENS[A-Z]*P\d"),
                   ("AUGUSTUS", r"^AUGUSTUS"), ("Trinity", r"_c\d+_g\d+_i\d+"),
                   ("FlyBase", r"^FBpp"), ("RefSeq", r"^[XN]P_\d"),
                   ("TransDecoder", r"\.p\d+$")]


def _guess_pipeline(ids):
    for name, rx in _PIPELINE_HINTS:
        if sum(1 for i in ids if re.search(rx, i)) > len(ids) * 0.5:
            return name
    return "unknown"


_METHOD_HINTS = [("DToL", r"dtol|_pep_"),
                 ("genome_guided", r"_GG_|_GG$|genome.?guided"),
                 ("denovo", r"denovo|de.?novo|dnovo")]


def _guess_species_method(sample, pipeline=""):
    """Method from the file name, or "" when nothing supports a guess.

    This used to fall back to "denovo", which quietly labelled 35 one-per-species
    genome annotations as de novo transcriptomes. A wrong guess that looks like
    an answer is worse than a blank: the blank is visible in the sheet and the
    comparison step refuses to run on it. "denovo" is now inferred only with
    positive evidence, namely Trinity-style IDs and no genome-guided or
    reference-annotation hint in the name."""
    meth = ""
    for name, rx in _METHOD_HINTS:
        if re.search(rx, sample, re.I):
            meth = name
            break
    if not meth and pipeline == "Trinity":
        meth = "denovo"
    sp = re.sub(r"_GG(_\d{6,8})?", "", sample)
    sp = re.sub(r"(_\d{6,8})?(_okay.*|_pep.*|_DToL.*)$", "", sp, flags=re.I)
    sp = re.sub(r"_(GG|dtol|denovo|dnovo)$", "", sp, flags=re.I)
    return sp, meth


def cmd_inspect(a):
    """Examine every proteome, propose a gene-ID rule, write a draft sheet.

    Every failure in this pipeline has come from an ID assumption that was
    wrong for one input. This makes the assumption visible before anything is
    computed: it prints real headers, shows how many proteins collapse to how
    many genes under each candidate rule, and writes a sample sheet you edit
    rather than a guess you never see."""
    need_pandas()
    pats = ("*.fa", "*.faa", "*.fasta", "*.fa.gz", "*.faa.gz", "*.fasta.gz")
    files = sorted({f for p in pats for f in glob.glob(os.path.join(a.proteomes, p))})
    if not files:
        sys.exit(f"ERROR: no FASTA files in {a.proteomes}")
    rows = []
    for f in files:
        # Rule DETECTION uses a sample; the collapse ESTIMATE uses every header.
        # Sampling the first N for the estimate is biased: Trinity/EvidentialGene
        # output lists isoform-rich clusters first, which overstated transcriptome
        # collapse by up to 16 percentage points on real data.
        hdrs, allh, n, lens = [], [], 0, []
        for h, sq in read_fasta(f):
            n += 1
            allh.append(h)
            lens.append(len(sq.rstrip("*")))
            if len(hdrs) < 500:
                hdrs.append(h)
        base = os.path.basename(f)
        print(f"\n=== {base}  ({n:,} sequences) ===")
        for h in hdrs[:3]:
            print(f"    {h[:150]}")
        if sum(1 for h in hdrs if GENE_TAG.search(h)) > len(hdrs) * 0.8:
            rule, rx, desc = "gene_tag", "", "gene: or GN= tag present in the description"
            ng = len({gene_of(h) for h in hdrs})
        else:
            cands = []
            for name, crx, cdesc in GENE_PATTERNS:
                hit = sum(1 for h in hdrs if re.search(crx, h.split()[0]))
                if hit > len(hdrs) * 0.5:
                    ngx = len({gene_of(h, crx) for h in hdrs})
                    cands.append((name, crx, cdesc, hit, ngx))
                    print(f"    candidate {name:18s} matches {hit}/{len(hdrs)} -> "
                          f"{ngx} genes from {len(hdrs)} proteins")
            merging = [c for c in cands if c[4] < len(hdrs)]
            pick = (merging or cands)
            if pick:
                rule, rx, desc, _, ng = pick[0]
            else:
                rule, rx, desc = "none", "", "IDs already look gene-level"
                ng = len({h.split()[0] for h in hdrs})
        ng_all = len({gene_of(h, rx or None) for h in allh}) if rule != "none" \
            else len({h.split()[0] for h in allh})
        ratio = ng_all / n if n else 1.0
        print(f"    -> rule '{rule}': {desc}")
        print(f"    -> all {n:,} proteins collapse to {ng_all:,} genes "
              f"({(1-ratio)*100:.1f}%)")
        if ratio > 0.99 and re.search(r"_i\d", " ".join(hdrs[:50])):
            print("    WARNING: an isoform-looking suffix is present but nothing "
                  "collapsed. Check by hand.")
        utr = sum(1 for h in hdrs if "utrorf" in h.lower())
        if utr:
            print(f"    NOTE: {utr}/{len(hdrs)} headers carry 'utrorf'. EvidentialGene")
            print("          marks ORFs found in another transcript's UTR; many are")
            print("          spurious. Consider --drop-utrorf when collapsing.")
        sample = re.sub(r"\.gz$", "", base)
        sample = re.sub(r"\.(aa\.)?(fa|faa|fasta)$", "", sample)
        pipe = _guess_pipeline([h.split()[0] for h in hdrs])
        sp, meth = _guess_species_method(sample, pipe)
        est = ng_all
        # an insect genome annotation should land roughly 10k-30k genes
        flag = ""
        if meth == "DToL" and est > 35000:
            flag = "IMPLAUSIBLE_GENE_COUNT"
            print(f"    WARNING: {est:,} genes after collapsing is not plausible for an")
            print(f"             insect genome annotation (expect roughly 10k-30k).")
            if rule == "transcript_t":
                print("             The .t suffix only groups isoforms when the base is a")
                print("             GENE id (g123.t1, g123.t2). Here the base looks like a")
                print("             unique protein serial, so nothing groups. This is")
                print("             probably unfiltered ab initio output. See UPGRADING.md.")
        if pipe == "UniProt":
            tr = sum(1 for h in hdrs if h.startswith("tr|"))
            notag = sum(1 for h in hdrs if not GENE_TAG.search(h))
            print(f"    NOTE: UniProt proteome; {tr}/{len(hdrs)} sampled entries unreviewed (tr|)"
                  + (f", {notag} with no GN= tag" if notag else ""))
            print("          Counts mix curation levels with the Ensembl and BRAKER")
            print("          annotations in the same run; treat the source as a covariate.")
        if pipe == "Trinity":
            print("    NOTE: Trinity 'genes' (c_g components) are not biological genes.")
            print("          A fragmented assembly splits one real gene across several")
            print("          components, so collapsed counts still run above a genome")
            print("          annotation. Use `compare --complete-only`.")
        shortest = min(lens) if lens else 0
        n_lt100 = sum(1 for x in lens if x < 100)
        if shortest >= 95 and n > 1000:
            print(f"    NOTE: the shortest protein is {shortest} aa and none are under 95.")
            print("          An ORF caller dropped short ORFs before this file was made")
            print("          (TransDecoder's default is -m 100). Metallothioneins (40-64 aa)")
            print("          cannot be in it. Add a `nucleotide` column pointing at the")
            print("          transcripts and run `short-orfs` to recover them.")
            flag = (flag + ";" if flag else "") + "NO_SHORT_PROTEINS"
        rows.append({"sample_id": sample, "species": sp, "method": meth,
                     "pipeline": pipe, "fasta": os.path.abspath(f),
                     "gene_regex": rx, "rule": rule, "n_proteins": n,
                     "est_genes": est, "collapse_pct": round((1 - ratio) * 100, 1),
                     "shortest_aa": shortest, "n_under_100aa": n_lt100,
                     "nucleotide": "", "flag": flag})
    sh = pd.DataFrame(rows)
    out = a.samplesheet or os.path.join(a.out, "samplesheet.tsv")
    d = os.path.dirname(os.path.abspath(out))
    if d:
        os.makedirs(d, exist_ok=True)
    sh.to_csv(out, sep="\t", index=False)
    print(f"\n=== draft sample sheet: {out} ===")
    print(sh[["sample_id", "species", "method", "pipeline", "rule", "n_proteins",
              "est_genes", "collapse_pct", "shortest_aa", "flag"]].to_string(index=False))
    if "pipeline" in sh and sh[sh.method == "DToL"].pipeline.nunique() > 1:
        print("\nNOTE: the DToL arm mixes annotation pipelines:")
        print(sh[sh.method == "DToL"].pipeline.value_counts().to_string())
        print("Recovery against 'DToL' will be recovery against different pipelines")
        print("for different species. compare reports it split by pipeline.")
    n_meth = sh.method.replace("", np.nan).nunique()
    dup = sh.species.duplicated().any()
    print("\n=== design ===")
    if n_meth > 1 and dup:
        print(pd.crosstab(sh.species, sh.method).to_string())
        print("\nSeveral methods per species: `compare` can pair them within species.")
    else:
        print(f"{len(sh)} proteomes, {sh.species.nunique()} species, "
              f"{'no method column filled in' if not n_meth else str(n_meth) + ' method'}.")
        print("One proteome per species, so there is NO method comparison to make and")
        print("`compare` will be skipped. Fill the method column only if some species")
        print("really were annotated more than one way.")
        if "pipeline" in sh and sh.pipeline.nunique() > 1:
            print("\nThe annotation pipelines differ between samples, though:")
            print(sh.pipeline.value_counts().to_string())
            print("That is a covariate, not a paired comparison. Different pipelines emit")
            print("different numbers of genes, so check it is not confounded with whatever")
            print("you are actually comparing, and pass --group-by pipeline to `dashboard`")
            print("to see it. Group your species with a metadata TSV for the real contrast.")
    print("\nEDIT THIS FILE. species and method are guessed from the file name and")
    print("gene_regex from the headers. Everything downstream trusts it.")


def cmd_collapse(a):
    """Longest protein per gene, per the sample sheet, into one directory."""
    need_pandas()
    sh = pd.read_csv(a.samplesheet, sep="\t")
    for c in ("sample_id", "fasta"):
        if c not in sh.columns:
            sys.exit(f"ERROR: sample sheet needs a '{c}' column")
    missing = [r.fasta for r in sh.itertuples(index=False) if not os.path.exists(r.fasta)]
    if missing:
        sys.exit("ERROR: missing FASTA files:\n  " + "\n  ".join(missing[:5]))
    dest = os.path.join(a.out, "proteomes")
    os.makedirs(dest, exist_ok=True)
    rows = []
    for r in sh.itertuples(index=False):
        rx = getattr(r, "gene_regex", "")
        rx = "" if rx is None or (isinstance(rx, float) and np.isnan(rx)) else str(rx).strip()
        best, n_in, n_utr = {}, 0, 0
        for h, seq in read_fasta(r.fasta):
            n_in += 1
            if a.drop_utrorf and "utrorf" in h.lower():
                n_utr += 1
                continue
            g = gene_of(h, rx or None)
            if g not in best or len(seq) > len(best[g][1]):
                best[g] = (h.split()[0], seq)
        with open(os.path.join(dest, f"{r.sample_id}.faa"), "w") as fh:
            for g, (pid, seq) in best.items():
                fh.write(f">{pid}\n{seq.rstrip('*')}\n")
        rows.append({"sample_id": r.sample_id,
                     "species": getattr(r, "species", ""),
                     "method": getattr(r, "method", ""),
                     "n_proteins_in": n_in, "n_genes_out": len(best),
                     "dropped_utrorf": n_utr,
                     "collapse_pct": round((1 - len(best) / max(n_in, 1)) * 100, 1)})
        log(f"{r.sample_id}: {n_in} proteins -> {len(best)} genes"
            + (f", {n_utr} utrorf dropped" if n_utr else ""))
        if getattr(r, "method", "") == "DToL" and len(best) > 35000:
            print(f"    WARNING: {r.sample_id} has {len(best):,} genes after collapsing;")
            print(f"             an insect genome should be roughly 10k-30k. Probably")
            print(f"             unfiltered ab initio output. See UPGRADING.md, or")
            print(f"             leave it out: compare --exclude {getattr(r, 'species', r.sample_id)}")
    st = pd.DataFrame(rows)
    st.to_csv(os.path.join(a.out, "collapse_stats.tsv"), sep="\t", index=False)
    print("\n" + st.to_string(index=False))
    if st.method.nunique() > 1:
        print("\n=== median by method ===")
        print(st.groupby("method")[["n_proteins_in", "n_genes_out", "collapse_pct"]]
                .median().round(1).to_string())
        print("\nIf one method collapses far more than another, its protein count was")
        print("measuring isoforms. Compare n_genes_out, never n_proteins_in.")
    log(f"collapsed proteomes in {dest}")


def cmd_compare(a):
    """Compare defensome recovery across annotation methods, paired by species.

    The design is what makes this work: the same species assayed by several
    methods turns annotation source from a confound into a measured variable,
    and pairs every comparison within species. Cross-sectional differences
    between 30 proteomes would tell you almost nothing; ten paired triplets
    tell you how much of a genome's defensome a transcriptome recovers."""
    need_pandas()
    sh = pd.read_csv(a.samplesheet, sep="\t")
    need = {"sample_id", "species", "method"}
    if not need.issubset(sh.columns):
        sys.exit(f"ERROR: sample sheet needs {sorted(need)}")
    cf = os.path.join(a.out, "counts.tsv")
    if not os.path.exists(cf):
        sys.exit(f"ERROR: {cf} not found. `scan` and `annotate` have not completed.\n"
                 "  python3 defensome.py doctor   checks the tools they need")
    counts = pd.read_csv(cf, sep="\t", index_col=0)
    npath = os.path.join(a.out, "counts_per10k.tsv")
    norm = pd.read_csv(npath, sep="\t", index_col=0) if os.path.exists(npath) else None
    dm = read_map(a.map)
    outdir = os.path.join(a.out, "comparison"); os.makedirs(outdir, exist_ok=True)

    if a.exclude:
        drop = {x.strip() for x in a.exclude.split(",") if x.strip()}
        log(f"excluding samples: {sorted(drop)}")
        sh = sh[~sh.sample_id.isin(drop) & ~sh.species.isin(drop)]
    sh = sh[sh.sample_id.isin(counts.index)]
    if sh.empty:
        sys.exit("ERROR: no sample_id in the sheet matches counts.tsv.\n"
                 "  The scan/annotate steps must use the collapsed proteomes, whose\n"
                 "  file names are the sample_id values.")
    meta = sh.set_index("sample_id")
    methods = sorted(meta.method.unique())
    species = sorted(meta.species.unique())
    design = pd.crosstab(meta.species, meta.method)
    print("=== design ===")
    print(design.to_string())
    if len(methods) < 2:
        sys.exit("Nothing to compare: the sample sheet has fewer than two methods.\n"
                 "  `compare` pairs several annotations OF THE SAME SPECIES. With one\n"
                 "  proteome per species there is no paired contrast; use a metadata TSV\n"
                 "  and --group-by for a between-species comparison instead.")
    if not (design > 0).sum(axis=1).gt(1).any():
        sys.exit("Nothing to compare: no species has more than one method.\n"
                 "  Check the method column in the sample sheet.")
    complete = design[(design > 0).all(axis=1)].index.tolist()
    print(f"\n{len(complete)} of {len(species)} species have every method; "
          f"paired tests use those {len(complete)}.")
    if len(complete) < 3:
        print("WARNING: too few complete sets for a paired test; reporting "
              "descriptives only.")

    core = [f for f in dm.loc[dm.tier == "CORE", "family"]
            if f in counts.columns and counts[f].sum() > 0]
    # Within one species the true defensome is fixed, so for a PAIRED method
    # comparison raw counts are the right measure. Dividing by proteome size
    # re-injects the artefact: the denominator is annotation-dependent (Trinity
    # fragmentation inflates it, an unfiltered ab initio set inflates it more),
    # so per-10k would penalise exactly the arms with the noisiest gene totals.
    # Per-10k is for comparing DIFFERENT species; it is opt-in here.
    if a.complete_only:
        cf = os.path.join(a.out, "domains", "completeness.tsv")
        if not os.path.exists(cf):
            sys.exit(f"ERROR: --complete-only needs {cf}. Run `domains` first.")
        comp = pd.read_csv(cf, sep="\t")
        counts = (comp[comp.status == "COMPLETE"]
                  .groupby(["species", "family"]).protein.nunique()
                  .unstack(fill_value=0).reindex(index=counts.index,
                                                 columns=counts.columns, fill_value=0))
        log("counting COMPLETE domain architectures only (fragments excluded)")
    use = norm if (norm is not None and a.per10k and not a.complete_only) else counts
    lab = ("copies per 10k genes" if use is norm else
           "COMPLETE gene counts" if a.complete_only else "raw gene counts")

    # long table: one row per species x method x family
    rows = []
    for sid in counts.index:
        if sid not in meta.index:
            continue
        for f in core:
            rows.append({"species": meta.loc[sid, "species"],
                         "method": meta.loc[sid, "method"], "sample_id": sid,
                         "family": f, "raw": int(counts.loc[sid, f]),
                         "value": float(use.loc[sid, f])})
    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(outdir, "long_counts.tsv"), sep="\t", index=False)

    piv = long.pivot_table(index=["species", "family"], columns="method",
                           values="value").reset_index()
    piv.to_csv(os.path.join(outdir, "species_family_by_method.tsv"), sep="\t", index=False)

    print(f"\n=== Median {lab} by method, summed over CORE families ===")
    tot = long.groupby(["species", "method"]).value.sum().reset_index()
    print(tot.pivot(index="species", columns="method", values="value")
             .round(1).to_string())

    # recovery against the reference method
    ref = a.reference if a.reference in methods else (
        "DToL" if "DToL" in methods else methods[0])
    print(f"\n=== Recovery relative to '{ref}' (raw gene counts, paired within species) ===")
    rec_rows = []
    for sp_ in complete:
        base = long[(long.species == sp_) & (long.method == ref)].set_index("family").raw
        for m in methods:
            if m == ref:
                continue
            other = long[(long.species == sp_) & (long.method == m)].set_index("family").raw
            for f in core:
                b, o = base.get(f, 0), other.get(f, 0)
                rec_rows.append({"species": sp_, "method": m, "family": f,
                                 "ref_count": b, "count": o,
                                 "recovery": (o / b) if b else np.nan})
    rec = pd.DataFrame(rec_rows)
    if not rec.empty and "pipeline" in meta.columns:
        refpipe = (meta[meta.method == ref].reset_index()
                   .set_index("species").pipeline.to_dict())
        rec["ref_pipeline"] = rec.species.map(refpipe)
    if not rec.empty:
        rec.to_csv(os.path.join(outdir, "recovery.tsv"), sep="\t", index=False)
        summ = (rec.groupby("method")
                   .agg(median_recovery=("recovery", "median"),
                        families_lost=("recovery", lambda x: int((x == 0).sum())),
                        families_exceeding=("recovery", lambda x: int((x > 1).sum())),
                        n=("recovery", "size")).round(3))
        print(summ.to_string())
        if "ref_pipeline" in rec and rec.ref_pipeline.nunique() > 1:
            print(f"\n=== Same, split by the pipeline that produced the '{ref}' reference ===")
            print(rec.groupby(["ref_pipeline", "method"]).recovery.median()
                     .unstack().round(3).to_string())
            print("If these differ, the reference arm is not one thing and a single")
            print("'recovery against DToL' number averages over two denominators.")
        print("\nRecovery above 1 does not mean the transcriptome found extra genes.")
        print("It usually means residual isoform inflation, or a family where the")
        print("genome annotation missed members. Check those families by hand.")
        byfam = (rec.groupby(["method", "family"]).recovery.median()
                    .unstack(0).round(2).sort_values(rec.method.unique()[0]))
        byfam.to_csv(os.path.join(outdir, "recovery_by_family.tsv"), sep="\t")
        print(f"\n=== Worst-recovered CORE families ===")
        print(byfam.head(10).to_string())

    # paired tests
    if len(complete) >= 3:
        try:
            from scipy import stats as _st
        except ImportError:
            _st = None
        if _st:
            trows = []
            for f in core:
                for i, m1 in enumerate(methods):
                    for m2 in methods[i + 1:]:
                        x, y = [], []
                        for sp_ in complete:
                            v1 = long[(long.species == sp_) & (long.method == m1) &
                                      (long.family == f)].value
                            v2 = long[(long.species == sp_) & (long.method == m2) &
                                      (long.family == f)].value
                            if len(v1) and len(v2):
                                x.append(float(v1.iloc[0])); y.append(float(v2.iloc[0]))
                        if len(x) >= 3 and any(p != q for p, q in zip(x, y)):
                            try:
                                w = _st.wilcoxon(x, y)
                                trows.append({"family": f, "method_a": m1, "method_b": m2,
                                              "n_pairs": len(x),
                                              "median_a": float(np.median(x)),
                                              "median_b": float(np.median(y)),
                                              "W": float(w.statistic), "p": float(w.pvalue)})
                            except Exception:
                                pass
            if trows:
                tt = pd.DataFrame(trows).sort_values("p")
                tt["p_bh"] = tt.p * len(tt) / tt.p.rank()
                tt["p_bh"] = tt.p_bh[::-1].cummin()[::-1].clip(upper=1)
                tt.to_csv(os.path.join(outdir, "paired_tests.tsv"), sep="\t", index=False)
                print(f"\n=== Wilcoxon signed-rank, paired within species ({len(complete)} pairs) ===")
                print(tt.head(12).round(4).to_string(index=False))
                print(f"\n{int((tt.p_bh < 0.05).sum())} of {len(tt)} family-by-method-pair")
                print("comparisons survive a 5% false discovery rate. These are paired")
                print("within species, so unlike a cross-species diet test they are not")
                print("confounded by phylogeny: the same genome is being measured twice.")

    _compare_figures(outdir, long, rec, core, methods, complete, lab, ref)
    log(f"wrote {outdir}/")


def _compare_figures(outdir, long, rec, core, methods, complete, lab, ref):
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
    except ImportError:
        log("matplotlib unavailable, skipping comparison figures"); return
    _pub_style()
    fd = os.path.join(outdir, "figures")
    if os.path.exists(os.path.join(fd, "FIGURES.md")):
        os.remove(os.path.join(fd, "FIGURES.md"))
    W = FigureWriter(fd, "comparison_all_figures.pdf")
    methods = method_order(methods)
    cols = {m: method_colour(m, i) for i, m in enumerate(methods)}
    xs = {m: i for i, m in enumerate(methods)}

    # C1: slopegraph, one line per species across methods, y from zero
    tot = long.groupby(["species", "method"]).value.sum().unstack().reindex(columns=methods)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
    ends = []
    for sp_ in tot.index:
        v = tot.loc[sp_]
        pts = [(xs[m], v[m]) for m in methods if not np.isnan(v[m])]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#adb5bd", lw=1, zorder=1)
        for m in methods:
            if not np.isnan(v[m]):
                ax.scatter(xs[m], v[m], s=36, color=cols[m], zorder=3, edgecolor="white", lw=.6)
        last = [m for m in methods if not np.isnan(v[m])][-1]
        ends.append([v[last], xs[last], sp_])
    # push end labels apart so species with equal totals stay readable
    ends.sort()
    span = (np.nanmax(tot.values) or 1) * .045
    for k in range(1, len(ends)):
        if ends[k][0] - ends[k - 1][0] < span:
            ends[k][0] = ends[k - 1][0] + span
    for ylab, xpos, sp_ in ends:
        yreal = tot.loc[sp_, methods[int(xpos)]]
        ax.annotate(sp_.replace("_", " "), (xpos, yreal), xytext=(xpos + .08, ylab),
                    textcoords="data", va="center", fontsize=7.5, style="italic",
                    arrowprops=dict(arrowstyle="-", color="#ced4da", lw=.5)
                    if abs(ylab - yreal) > 1e-9 else None)
    ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods)
    ax.set_xlim(-.3, len(methods) - .3); ax.set_ylim(0, None)
    ax.set_ylabel(f"total CORE defensome ({lab})")
    ax.set_title("A  Same genome, three annotation routes", loc="left")
    if ref in tot.columns:
        others = [m for m in methods if m != ref]
        ratio = tot[others].div(tot[ref], axis=0)
        for i, m in enumerate(others):
            vals = ratio[m].dropna()
            ax2.scatter(np.random.default_rng(1).normal(i, .05, len(vals)), vals, s=30,
                        color=cols[m], edgecolor="white", lw=.6, zorder=3)
            if len(vals):
                ax2.hlines(vals.median(), i - .25, i + .25, color=cols[m], lw=2.5, zorder=4)
                ax2.annotate(f"median {vals.median():.2f}", (i + .28, vals.median()),
                             va="center", fontsize=8, color=cols[m])
        ax2.axhline(1, ls="--", color="#6c757d", lw=1)
        ax2.set_xticks(range(len(others))); ax2.set_xticklabels(others)
        ax2.set_xlim(-.5, len(others) - .1)
        ax2.set_ylabel(f"total relative to {ref}")
        ax2.set_title(f"B  Recovery of the {ref} total", loc="left")
    fig.tight_layout()
    W.save(fig, "C1_paired_totals",
           f"(A) Total CORE defensome per species under each annotation route, joined within "
           f"species; y-axis from zero. (B) Each route's total divided by the {ref} total for the "
           f"same species; dashed line is parity, bar is the median. Values are {lab}.")

    # C2: Cleveland dot plot, scales to any number of families
    med = long.groupby(["family", "method"]).value.median().unstack().reindex(columns=methods)
    sort_by = ref if ref in med.columns else methods[0]
    med = med.sort_values(sort_by, ascending=True)
    fig, ax = plt.subplots(figsize=(7.5, _dim(len(med), .32, 3.2)))
    dodge = {m: (i - (len(methods) - 1) / 2) * .18 for i, m in enumerate(methods)}
    for yi, (fam, r) in enumerate(med.iterrows()):
        v = r.dropna()
        if len(v):
            ax.plot([v.min(), v.max()], [yi, yi], color="#dee2e6", lw=3, zorder=1,
                    solid_capstyle="round")
        for m in methods:
            if not np.isnan(r[m]):
                # vertical dodge: equal values would otherwise hide all but one colour
                ax.scatter(r[m], yi + dodge[m], s=38, color=cols[m], zorder=3,
                           edgecolor="white", lw=.6, label=m if yi == 0 else None)
    ax.set_yticks(range(len(med))); ax.set_yticklabels(med.index)
    pos = med.values[~np.isnan(med.values)]
    if len(pos) and pos.max() > 0 and pos[pos > 0].min() > 0 and pos.max() / pos[pos > 0].min() > 25:
        ax.set_xscale("log")
    ax.set_xlabel(f"median across species ({lab})")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", ncol=len(methods))
    ax.set_title("Per-family gene counts by annotation route", loc="left")
    fig.tight_layout()
    W.save(fig, "C2_family_by_method",
           f"Median {lab} per family across species, one dot per annotation route. The grey bar "
           f"spans the routes, so a long bar is a family whose count depends heavily on how the "
           f"proteome was made. Families sorted by the {sort_by} value.")

    # C3: recovery heatmap, family x method, centred on parity
    if rec is not None and not rec.empty:
        others = [m for m in methods if m != ref]
        hm = rec.groupby(["family", "method"]).recovery.median().unstack().reindex(columns=others)
        hm = hm.loc[hm.mean(axis=1).sort_values().index]
        fig, ax = plt.subplots(figsize=(2.2 + 1.5 * len(others), _dim(len(hm), .3, 3.2)))
        data = hm.values.astype(float)
        vmax = max(2.0, np.nanmax(data) if np.isfinite(np.nanmax(data)) else 2.0)
        im = ax.imshow(np.nan_to_num(data, nan=1.0), cmap="RdBu", aspect="auto",
                       norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax))
        for yi in range(data.shape[0]):
            for xi in range(data.shape[1]):
                v = data[yi, xi]
                if np.isnan(v):
                    ax.add_patch(plt.Rectangle((xi - .5, yi - .5), 1, 1, color="#f1f3f5"))
                    ax.text(xi, yi, "n/a", ha="center", va="center", fontsize=7, color="#adb5bd")
                else:
                    ax.text(xi, yi, f"{v:.2f}", ha="center", va="center", fontsize=7.5,
                            color="white" if abs(v - 1) > .45 else "#212529")
        ax.set_xticks(range(len(others))); ax.set_xticklabels(others)
        ax.set_yticks(range(len(hm))); ax.set_yticklabels(hm.index)
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, shrink=.6, pad=.02)
        cb.set_label(f"median recovery vs {ref}")
        ax.set_title(f"Recovery relative to {ref}", loc="left")
        fig.tight_layout()
        W.save(fig, "C3_recovery_heatmap",
               f"Median over species of (gene count in route / gene count in {ref}) for each "
               f"family. 1 is parity; red below, blue above. Above 1 usually means residual "
               f"isoform or fragment inflation rather than genes {ref} missed. n/a: {ref} had none.")
    W.close()
    log(f"comparison figures in {fd} (PDF, PNG, combined PDF, FIGURES.md)")


def cmd_doctor(a):
    """Report which tools and libraries are present, and which steps can run."""
    import platform
    print(f"defensome {__version__}  ({os.path.abspath(__file__)})")
    print(f"python    {platform.python_version()}  ({sys.executable})\n")
    libs = {}
    for m in ("pandas", "numpy", "scipy", "matplotlib", "openpyxl"):
        try:
            mod = __import__(m)
            libs[m] = getattr(mod, "__version__", "yes")
        except Exception:
            libs[m] = None
    print("python libraries")
    for m, v in libs.items():
        need = "required" if m in ("pandas", "numpy") else "optional"
        print(f"  {'OK ' if v else '-- '} {m:11s} {v or 'missing':12s} {need}")
    print("\nexternal tools")
    tools = {}
    for t in ("hmmsearch", "hmmbuild", "mafft", "FastTree", "cd-hit", "node"):
        p = tool(t) or (tool("fasttree") if t == "FastTree" else None)
        ver = ""
        if p:
            try:
                o = subprocess.run([p, "-h"] if t != "node" else [p, "--version"],
                                   capture_output=True, text=True, timeout=10)
                txt = (o.stdout + o.stderr)
                m_ = re.search(r"(\d+\.\d+(?:\.\d+)?)", txt)
                ver = m_.group(1) if m_ else ""
            except Exception:
                pass
        tools[t] = p
        print(f"  {'OK ' if p else '-- '} {t:11s} {ver or ('missing' if not p else ''):12s} {p or ''}")
    steps = {
        "inspect / collapse / compare / report / dashboard": libs["pandas"] and libs["numpy"],
        "scan": tools["hmmsearch"],
        "cyp": tools["hmmbuild"] and tools["hmmsearch"] and tools["mafft"],
        "trees": tools["mafft"] and tools["FastTree"],
        "figures": libs["matplotlib"],
        "dashboard test": tools["node"],
    }
    print("\nsteps that can run here")
    for k, ok in steps.items():
        print(f"  {'OK ' if ok else '-- '} {k}")
    missing = [k for k, ok in steps.items() if not ok]
    if missing:
        print("\nOn an Lmod cluster, typically:")
        print("  module load HMMER MAFFT FastTree")
        print("If HMMER and your Python need incompatible toolchains, point at the")
        print("binary instead of loading both modules:")
        print("  module load HMMER && export DEFENSOME_HMMSEARCH=$(which hmmsearch) \\")
        print("                    && export DEFENSOME_HMMBUILD=$(which hmmbuild)")
        print("  module purge && module load Python")
    if a.require:
        req = {x.strip() for x in a.require.split(",") if x.strip()}
        bad = [r for r in req if r in ("scan",) and not tools["hmmsearch"]] + \
              [r for r in req if r in ("cyp",) and not steps["cyp"]] + \
              [r for r in req if r in ("trees",) and not steps["trees"]] + \
              [r for r in req if r in ("core",) and not steps[
                  "inspect / collapse / compare / report / dashboard"]]
        if bad:
            sys.exit(f"\nERROR: required step(s) cannot run: {', '.join(sorted(set(bad)))}")


# --------------------------------------------------------------------------- #
#  Short-ORF recovery from nucleotide sequence
# --------------------------------------------------------------------------- #
_COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")
_CODON = {}
for _i, _a in enumerate("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    _CODON["TCAG"[_i // 16] + "TCAG"[(_i // 4) % 4] + "TCAG"[_i % 4]] = _a


def _translate(nt):
    return "".join(_CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def short_orfs(seq, min_aa=25, max_aa=99):
    """ATG-to-stop ORFs of min_aa..max_aa residues on both strands.

    Stop and ATG positions are found once with regular expressions and paired
    by frame, so only the candidate ORFs are ever translated. Translating whole
    transcriptomes six ways in pure Python would be far slower."""
    import bisect
    seq = seq.upper()
    out = []
    for strand, s in (("+", seq), ("-", seq.translate(_COMP)[::-1])):
        stops = [[], [], []]
        for m in re.finditer(r"(?=(TAA|TAG|TGA))", s):
            stops[m.start() % 3].append(m.start())
        atgs = [[], [], []]
        for m in re.finditer(r"(?=ATG)", s):
            atgs[m.start() % 3].append(m.start())
        for f in range(3):
            prev = f - 3
            for st in stops[f] + [None]:
                end = st if st is not None else len(s) - ((len(s) - f) % 3)
                j = bisect.bisect_left(atgs[f], prev + 3)
                if j < len(atgs[f]) and atgs[f][j] < end and st is not None:
                    a0 = atgs[f][j]
                    n_aa = (st - a0) // 3
                    if min_aa <= n_aa <= max_aa:
                        out.append((strand, a0, n_aa, _translate(s[a0:st])))
                if st is None:
                    break
                prev = st
    return out


def cmd_short_orfs(a):
    """Recover proteins too short for the ORF caller that built the proteome.

    TransDecoder keeps ORFs of at least 100 residues by default, and most gene
    predictors skip very short genes, so an insect metallothionein (40-64
    residues) is usually absent from a transcriptome proteome before any
    search runs. No domain rule can find a protein that is not in the file.
    This goes back to the transcripts named in the sample sheet's `nucleotide`
    column, takes ATG-to-stop ORFs in the band the ORF caller discarded, and
    applies the metallothionein composition rule. Genes, not transcripts, are
    counted, using the sample sheet's gene_regex."""
    need_pandas()
    sh = pd.read_csv(a.samplesheet, sep="\t")
    col = a.nucleotide_col
    if col not in sh.columns:
        sys.exit(f"ERROR: the sample sheet has no '{col}' column.\n"
                 f"  Add one giving the transcript FASTA for each sample (Trinity.fasta etc.).\n"
                 f"  Leave it empty for samples without one; they are skipped.")
    d = os.path.join(a.out, "short_orfs"); os.makedirs(d, exist_ok=True)
    hits, summ = [], []
    for r in sh.itertuples(index=False):
        nt = getattr(r, col, None)
        if nt is None or (isinstance(nt, float) and np.isnan(nt)) or not str(nt).strip():
            continue
        nt = str(nt).strip()
        if not os.path.exists(nt):
            log(f"WARNING: {r.sample_id}: nucleotide file not found: {nt}"); continue
        if os.path.getsize(nt) > a.max_mb * 1e6:
            log(f"WARNING: {r.sample_id}: {nt} is over {a.max_mb} MB; this step is meant for "
                f"transcript sets, not whole genomes. Use --max-mb to override.")
            continue
        rx = getattr(r, "gene_regex", "")
        rx = "" if rx is None or (isinstance(rx, float) and np.isnan(rx)) else str(rx).strip()
        n_tx = n_orf = 0
        genes = {}
        with open(os.path.join(d, f"{r.sample_id}.mt_like.faa"), "w") as fo:
            for hdr, seq in read_fasta(nt):
                n_tx += 1
                for strand, pos, n_aa, prot in short_orfs(seq, a.min_aa, a.max_aa):
                    n_orf += 1
                    if not _is_mt_like(prot):
                        continue
                    tid = hdr.split()[0]
                    gene = gene_of(hdr, rx or r"_i\d+\w*$")
                    f_ = mt_features(prot)
                    fo.write(f">{tid}|{strand}{pos} gene={gene} len={n_aa}\n{prot}\n")
                    hits.append({"sample_id": r.sample_id, "species": getattr(r, "species", ""),
                                 "method": getattr(r, "method", ""), "gene": gene,
                                 "transcript": tid, "strand": strand, "start": pos,
                                 **f_, "protein": prot})
                    genes.setdefault(gene, prot)
        summ.append({"sample_id": r.sample_id, "species": getattr(r, "species", ""),
                     "method": getattr(r, "method", ""), "transcripts": n_tx,
                     "short_orfs_scanned": n_orf, "mt_like_genes": len(genes)})
        log(f"{r.sample_id}: {n_tx:,} transcripts, {n_orf:,} short ORFs, "
            f"{len(genes)} MT-like genes")
    if not summ:
        sys.exit(f"ERROR: no sample had a usable '{col}' entry")
    hd = pd.DataFrame(hits)
    hd.to_csv(os.path.join(a.out, "short_orf_mt.tsv"), sep="\t", index=False)
    sd = pd.DataFrame(summ)
    sd.to_csv(os.path.join(a.out, "short_orf_summary.tsv"), sep="\t", index=False)
    print("\n" + sd.to_string(index=False))
    print("\nThese are composition-screen candidates from ORFs the proteome never")
    print("contained. They are reported apart from every count table. Confirm a few")
    print("by alignment to a known insect metallothionein before citing numbers.")
    log(f"wrote {a.out}/short_orf_mt.tsv and short_orf_summary.tsv")


def cmd_dataset_help(a):
    print(DATASET_HELP)


DATASET_HELP = """
External datasets: what each can and cannot do
==============================================

                what it holds                    protein?  negatives?  insect?
  ICPD          66,477 P450s / 682 species       yes       no          yes
  P450RDB       enzyme -> substrate + sequence   yes       NO          check
  Ni 2025       ~2000 compounds per human CYP    none      YES         no

1. ICPD  http://www.insectp450.net/ui/#/page/download

   Once downloaded (all_prot.fa, information.xlsx, References1.xlsx,
   References2.xlsx):

       python3 defensome.py cyp-refs --source icpd \\
           --fasta db/icpd/all_prot.fa --table db/icpd/information.xlsx \\
           --orders Lepidoptera --out db/icpd/ref.faa
       cd-hit -i db/icpd/ref.faa -o db/icpd/ref.c80.faa -c 0.8 -n 5
       python3 defensome.py cyp --out results/ --refs db/icpd/ref.c80.faa

       python3 defensome.py icpd-evidence \\
           --tables db/icpd/References1.xlsx db/icpd/References2.xlsx \\
           --out db/icpd/evidence

   The FASTA is keyed on the table's "Database ID" column (Abter001), NOT
   "Protein ID"; cyp-refs finds the right column by matching against the FASTA
   headers rather than guessing from the column name. Entries marked
   "fragmented" in the Completeness column are excluded by default: a truncated
   P450 degrades any profile it goes into.

   OLD INSTRUCTIONS BELOW, superseded:
   Download 'Protein Sequence' (fasta) and 'Sequence information table' (xlsx).
   The links are behind a web UI, so fetch them in a browser and copy across:

       mkdir -p db/icpd
       # then scp/rsync the two files into db/icpd/

   Build a lepidopteran clan reference and use it:

       python3 defensome.py cyp-refs --source icpd \
           --fasta db/icpd/ICPD_protein.fasta \
           --table db/icpd/ICPD_sequence_info.xlsx \
           --orders Lepidoptera --out db/icpd/ref.faa
       cd-hit -i db/icpd/ref.faa -o db/icpd/ref.c80.faa -c 0.8 -n 5
       python3 defensome.py cyp --out results/ --refs db/icpd/ref.c80.faa

   Validate the clan assignment against ICPD's own curated clades. The
   reference uses only a few hundred sequences per clan, so the rest is a large
   held-out labelled test set:

       python3 defensome.py cyp-benchmark \\
           --fasta db/icpd/all_prot.fa --table db/icpd/information.xlsx \\
           --clan-hmms results/cyp/clan_hmms --refs db/icpd/ref.c80.faa \\
           --orders Lepidoptera --out db/icpd/benchmark

   A large bitscore margin says an assignment is unambiguous. It does not say
   it is right. This measures accuracy directly and writes a confusion matrix.

   Also download 'References2' (genes associated with phenotypes). Separating
   'overexpressed in a resistant strain' (correlative) from 'heterologously
   expressed and shown to metabolise X' (functional) is the single most
   important curation decision in this analysis. Only the second is a
   substrate label.

   Cite: Wu et al. (2025) Mol Ecol Resour 25:e14070.

2. P450RDB  https://www.cellknowledge.com.cn/p450rdb_v2/download.html
       python3 defensome.py chem-triage --dir db/p450rdb --out db/p450rdb/triage

   Expect mostly plant and microbial BIOSYNTHETIC reactions. A reaction
   database has no negatives, so it cannot train a binary classifier. Its real
   value is triage/plant_metabolite_inventory.tsv: the chemistry your
   herbivores face, with structures, joinable to your host-plant table.

3. Ni et al. (2025) Sci Data 12:1427, doi 10.1038/s41597-025-05753-8
   The only source with curated NON-substrates. Two jobs, neither of them
   predicting insect biology: reproduce their GCN as a methods benchmark
   (published MCC 0.51 to 0.72), then measure how far a human-trained model
   transfers to functionally validated insect pairs. That drop is the result.
   Licence is CC BY-NC-ND; NoDerivatives constrains merged datasets. Check the
   Figshare record separately.

Nothing you already have under db/ can substitute for ICPD. Until it is
downloaded, --source flybase keeps you running with Drosophila references and
LOW-confidence calls.
"""


# --------------------------------------------------------------------------- #


def main():
    p = argparse.ArgumentParser(prog="defensome.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version",
                   version=f"defensome {__version__}  ({os.path.abspath(__file__)})")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(q, need_proteomes=False):
        q.add_argument("--out", required=True, help="output directory")
        q.add_argument("--map", default=DEFAULT_MAP, help="defensome map TSV")
        q.add_argument("--proteomes", required=need_proteomes,
                       help="directory of .faa proteomes")

    s = sub.add_parser("scan", help="run hmmsearch per proteome")
    common(s, True)
    s.add_argument("--pfam", required=True, help="Pfam-A.hmm (hmmpress'd)")
    s.add_argument("--threads", type=int, default=8)
    s.add_argument("--species", help="run one species only (for SLURM arrays)")
    s.add_argument("--force", action="store_true")
    s.add_argument("--no-rescue", action="store_true", dest="no_rescue",
                   help="skip the below-threshold rescue pass for rescue=yes families")
    s.add_argument("--rescue-evalue", type=float, default=1e-3, dest="rescue_evalue")
    s.add_argument("--full-pfam", dest="full_pfam",
                   help="full Pfam-A.hmm for an architecture pass on candidate proteins")
    s.set_defaults(func=cmd_scan)

    n = sub.add_parser("annotate", help="apply domain rules, build count matrices")
    common(n)
    n.set_defaults(func=cmd_annotate)

    r = sub.add_parser("report", help="summaries and figures")
    common(r)
    r.add_argument("--metadata", help="TSV, first column or 'species' = proteome name")
    r.add_argument("--group-by", help="metadata column to compare groups by")
    r.set_defaults(func=cmd_report)

    e = sub.add_parser("extract", help="write per-family protein FASTA")
    common(e, True)
    e.add_argument("--families", help="comma-separated subset, default all")
    e.set_defaults(func=cmd_extract)

    d = sub.add_parser("domains", help="per-protein domain architecture and completeness")
    common(d)
    d.set_defaults(func=cmd_domains)

    q = sub.add_parser("qc", help="annotation quality from low-copy control families")
    common(q)
    q.add_argument("--expect", help="TSV of family<TAB>expected_count, overriding the "
                   "built-in insect priors used by question (1)")
    q.set_defaults(func=cmd_qc)

    y = sub.add_parser("cyp", help="assign a family to clans using labelled references")
    common(y)
    y.add_argument("--refs", required=True,
                   help="FASTA with clan=XXX in each header (e.g. FlyBase CYP refs)")
    y.add_argument("--samples",
                   help="comma-separated sample_ids to assign; use when the clan reference "
                        "covers only some of the lineages in the run")
    y.add_argument("--family", default="CYP", help="which family to classify")
    y.add_argument("--threads", type=int, default=8)
    y.set_defaults(func=cmd_cyp)

    t = sub.add_parser("trees", help="MAFFT + FastTree per family")
    common(t)
    t.add_argument("--families", help="comma-separated subset, default all extracted")
    t.add_argument("--max-seqs", type=int, default=6000, dest="max_seqs")
    t.add_argument("--per-species", type=int, default=0, dest="per_species",
                   help="when a family exceeds --max-seqs, keep at most N sequences "
                        "per species (closest to the median length) instead of "
                        "skipping it; prefer --thin, which keeps every copy")
    t.add_argument("--metadata", help="metadata TSV, needed by --where")
    t.add_argument("--where", action="append", metavar="COLUMN=VALUE",
                   help="restrict to species whose metadata matches, e.g. "
                        "order=Lepidoptera or habitat=marine,freshwater; repeat to AND")
    t.add_argument("--complete-only", action="store_true", dest="complete_only",
                   help="only gene models with every required domain intact "
                        "(needs `domains`); fragments distort alignments")
    t.add_argument("--clan", help="CYP clan(s) to keep, e.g. CYP3 (needs `cyp`)")
    t.add_argument("--name", help="label for the subset tree instead of one built "
                                  "from the restrictions")
    t.add_argument("--thin", metavar="COLUMN[:N]",
                   help="keep N representative species (default 1) per value of a "
                        "metadata column, with all their copies, e.g. family:1")
    t.add_argument("--max-len-ratio", type=float, default=2.0, dest="max_len_ratio",
                   help="drop sequences longer than this multiple of the median "
                        "length (fused or read-through models); 0 keeps all")
    t.add_argument("--threads", type=int, default=8)
    t.add_argument("--force", action="store_true")
    t.set_defaults(func=cmd_trees)

    w = sub.add_parser("tree", help="radial + rectangular species tree with defensome rings")
    common(w)
    w.add_argument("--tree", required=True, help="Newick species tree (e.g. OrthoFinder SpeciesTree_rooted.txt)")
    w.add_argument("--metadata"); w.add_argument("--group-by", dest="group_by")
    w.add_argument("--families", help="comma-separated ring families, default all CORE")
    w.add_argument("--phylogram", action="store_true", help="use branch lengths instead of a cladogram")
    w.add_argument("--raw", action="store_true", help="plot raw counts instead of per-10k")
    w.set_defaults(func=cmd_tree)

    gt = sub.add_parser("genetree", help="radial gene tree for one family")
    common(gt)
    gt.add_argument("--family", default="CYP")
    gt.add_argument("--tree", help="default results/trees/<FAMILY>.tre")
    gt.add_argument("--by", default="trait", choices=["trait","clan"])
    gt.add_argument("--samplesheet", help="colour by method from a sample sheet")
    gt.add_argument("--metadata"); gt.add_argument("--group-by", dest="group_by")
    gt.add_argument("--max-labels", type=int, default=400, dest="max_labels",
                    help="species labels only when the tree has at most this many tips")
    gt.set_defaults(func=cmd_genetree)

    cr = sub.add_parser("cyp-refs", help="build a clan-labelled CYP reference FASTA")
    cr.add_argument("--source", choices=["icpd", "flybase"], default="icpd")
    cr.add_argument("--fasta", help="protein FASTA (ICPD, or the FlyBase clan ref)")
    cr.add_argument("--table", help="ICPD sequence information table (xlsx/tsv/csv)")
    cr.add_argument("--orders", default="", help="e.g. Lepidoptera")
    cr.add_argument("--species", default="")
    cr.add_argument("--min-len", type=int, default=400, dest="min_len")
    cr.add_argument("--max-per-clan", type=int, default=300, dest="max_per_clan")
    cr.add_argument("--id-col", dest="id_col",
                    help="table column holding the FASTA header IDs; auto-detected by matching")
    cr.add_argument("--include-fragments", action="store_true", dest="include_fragments",
                    help="keep entries marked fragmented (default: exclude them)")
    cr.add_argument("--out", required=True)
    cr.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    cr.add_argument("--proteomes", help=argparse.SUPPRESS)
    cr.set_defaults(func=cmd_cyp_refs, out_is_dir=False)

    ct = sub.add_parser("chem-triage", help="triage P450RDB before committing to it")
    ct.add_argument("--dir", required=True, help="directory holding the P450RDB CSVs")
    ct.add_argument("--out", required=True)
    ct.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    ct.add_argument("--proteomes", help=argparse.SUPPRESS)
    ct.set_defaults(func=cmd_chem_triage)

    ie = sub.add_parser("icpd-evidence",
                        help="sort ICPD reference records into evidence tiers")
    ie.add_argument("--tables", nargs="+", required=True,
                    help="References1.xlsx References2.xlsx")
    ie.add_argument("--out", required=True)
    ie.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    ie.add_argument("--proteomes", help=argparse.SUPPRESS)
    ie.set_defaults(func=cmd_icpd_evidence)

    cb = sub.add_parser("cyp-benchmark",
                        help="measure clan-assignment accuracy on held-out ICPD sequences")
    cb.add_argument("--fasta", required=True, help="ICPD protein FASTA")
    cb.add_argument("--table", required=True, help="ICPD information table")
    cb.add_argument("--clan-hmms", required=True, dest="clan_hmms",
                    help="directory of clan HMMs, e.g. results/cyp/clan_hmms")
    cb.add_argument("--refs", help="the reference FASTA, so its sequences are excluded")
    cb.add_argument("--orders", default="", help="restrict the test set, e.g. Lepidoptera")
    cb.add_argument("--min-len", type=int, default=400, dest="min_len")
    cb.add_argument("--threads", type=int, default=8)
    cb.add_argument("--out", required=True)
    cb.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    cb.add_argument("--proteomes", help=argparse.SUPPRESS)
    cb.set_defaults(func=cmd_cyp_benchmark)

    su = sub.add_parser("setup", help="get the defensome HMMs: download the prebuilt set, or build it from Pfam-A")
    su.add_argument("--pfam", help="existing Pfam-A.hmm or Pfam-A.hmm.gz")
    su.add_argument("--download", action="store_true", help="fetch Pfam-A from EBI")
    su.add_argument("--url", help="override the Pfam-A download URL")
    su.add_argument("--keep-full", action="store_true", dest="keep_full",
                    help="with --download, also keep an uncompressed Pfam-A for --full-pfam")
    su.add_argument("--db-dir", default="db", dest="db_dir")
    su.add_argument("--fetch", action="store_true",
                    help="download the prebuilt defensome HMMs (the default when no --pfam)")
    su.add_argument("--db-url", dest="db_url", help=f"override the prebuilt database URL ({DB_URL})")
    su.add_argument("--pfam-version", dest="pfam_version",
                    help="record the Pfam release, e.g. 'Pfam release 37.4' (else read from Pfam.version)")
    su.add_argument("--bundle", action="store_true",
                    help="also write defensome_db.tar.gz + .sha256 to publish as a release asset")
    su.add_argument("--force", action="store_true")
    su.add_argument("--map", default=DEFAULT_MAP)
    su.add_argument("--out", default=".", help=argparse.SUPPRESS)
    su.add_argument("--proteomes", help=argparse.SUPPRESS)
    su.set_defaults(func=cmd_setup, out_is_dir=False)

    dr = sub.add_parser("doctor", help="check tools and libraries; report which steps can run")
    dr.add_argument("--require", default="",
                    help="comma list (core,scan,cyp,trees); exit non-zero if any cannot run")
    dr.add_argument("--out", default=".", help=argparse.SUPPRESS)
    dr.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    dr.add_argument("--proteomes", help=argparse.SUPPRESS)
    dr.set_defaults(func=cmd_doctor, out_is_dir=False)

    insp = sub.add_parser("inspect",
                          help="examine proteome headers, write a draft sample sheet")
    insp.add_argument("--proteomes", required=True)
    insp.add_argument("--out", default=".")
    insp.add_argument("--samplesheet", help="output path, default <out>/samplesheet.tsv")
    insp.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    insp.set_defaults(func=cmd_inspect)

    col = sub.add_parser("collapse",
                         help="longest protein per gene, per the sample sheet")
    col.add_argument("--samplesheet", required=True)
    col.add_argument("--out", required=True)
    col.add_argument("--drop-utrorf", action="store_true", dest="drop_utrorf",
                     help="discard EvidentialGene utrorf entries")
    col.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    col.add_argument("--proteomes", help=argparse.SUPPRESS)
    col.set_defaults(func=cmd_collapse)

    cmp_ = sub.add_parser("compare",
                          help="compare defensome recovery across annotation methods")
    cmp_.add_argument("--samplesheet", required=True)
    cmp_.add_argument("--out", required=True)
    cmp_.add_argument("--reference", default="DToL",
                      help="method to express recovery against (default DToL)")
    cmp_.add_argument("--per10k", action="store_true",
                      help="normalise per 10k genes (default: raw, correct for paired within-species comparison)")
    cmp_.add_argument("--complete-only", action="store_true", dest="complete_only",
                      help="count only COMPLETE domain architectures; strips Trinity fragments")
    cmp_.add_argument("--exclude", default="",
                      help="comma-separated sample_ids or species to leave out")
    cmp_.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    cmp_.add_argument("--proteomes", help=argparse.SUPPRESS)
    cmp_.set_defaults(func=cmd_compare)

    so = sub.add_parser("short-orfs",
                        help="recover proteins too short for the ORF caller (metallothioneins)")
    so.add_argument("--samplesheet", required=True)
    so.add_argument("--out", required=True)
    so.add_argument("--nucleotide-col", default="nucleotide", dest="nucleotide_col",
                    help="sample-sheet column holding transcript FASTA paths")
    so.add_argument("--min-aa", type=int, default=25, dest="min_aa")
    so.add_argument("--max-aa", type=int, default=99, dest="max_aa",
                    help="upper bound; default 99, the band TransDecoder's default -m 100 discards")
    so.add_argument("--max-mb", type=float, default=2000, dest="max_mb")
    so.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    so.add_argument("--proteomes", help=argparse.SUPPRESS)
    so.set_defaults(func=cmd_short_orfs)

    dh = sub.add_parser("dataset-help", help="how to obtain and use ICPD, P450RDB and the human CYP dataset")
    dh.add_argument("--out", default=".", help=argparse.SUPPRESS)
    dh.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    dh.add_argument("--proteomes", help=argparse.SUPPRESS)
    dh.set_defaults(func=cmd_dataset_help)

    asx = sub.add_parser("assets", help="write the embedded map/dashboard files out for editing")
    asx.add_argument("--write", default=".", help="directory to write into")
    asx.add_argument("--force", action="store_true")
    asx.add_argument("--out", default=".", help=argparse.SUPPRESS)
    asx.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    asx.add_argument("--proteomes", help=argparse.SUPPRESS)
    asx.set_defaults(func=cmd_assets)

    st = sub.add_parser("structure", help="Boltz-2 structure (and ligand binding) for chosen genes")
    st.add_argument("--out", required=True, help="output directory")
    st.add_argument("--gene", action="append", required=True, metavar="SAMPLE|PROTEIN",
                    help="gene to predict; repeat for several")
    st.add_argument("--ligand", metavar="SMILES", help="ligand to dock, with an affinity request")
    st.add_argument("--ligand-name", dest="ligand_name", help="short name for the ligand (folder name)")
    st.add_argument("--cofactor", default="auto",
                    help="auto (the family's cofactor, e.g. haem for CYP), none, or CCD codes: HEM,FAD")
    st.add_argument("--proteomes", help="directory of .faa proteomes (else results/fasta)")
    st.add_argument("--samples", type=int, default=1, help="diffusion samples (Boltz default 1)")
    st.add_argument("--accelerator", choices=["gpu", "cpu"], help="passed to boltz")
    st.add_argument("--run", action="store_true", help="run boltz now, not just write the inputs")
    st.add_argument("--kernels", action="store_true",
                    help="use Boltz's fast GPU kernels (needs pip install 'boltz[cuda]')")
    st.set_defaults(func=cmd_structure)

    db = sub.add_parser("dashboard", help="static dashboard site: comparative views, "
                                          "species, family and gene pages")
    common(db)
    db.add_argument("--no-sequences", action="store_true", dest="no_sequences",
                    help="leave protein sequences out of the gene pages and FASTA downloads")
    db.add_argument("--force", action="store_true",
                    help="rewrite the species data even if no input changed")
    db.add_argument("--ga-only", action="store_true", dest="ga_only",
                    help="count only genes at Pfam's gathering threshold, leaving out "
                         "the validated below-threshold ones for `rescue` families")
    db.add_argument("--metadata"); db.add_argument("--group-by", dest="group_by")
    db.add_argument("--species-tree", dest="species_tree",
                    help="Newick species tree to embed (OrthoFinder SpeciesTree_rooted.txt)")
    db.add_argument("--dashboard-out", dest="dashboard_out",
                    help="output folder, default <out>/dashboard")
    db.add_argument("--samplesheet",
                    help="sample sheet; adds species/method/pipeline, and makes method the default trait")
    db.set_defaults(func=cmd_dashboard)

    def run_args(q):
        q.add_argument("--proteomes", required=True, help="directory of .faa/.fa/.fasta proteomes")
        q.add_argument("--out", required=True, help="output directory")
        q.add_argument("--species", help="only these: a comma list or a file with one name per line "
                                          "(file stem or binomial; 'Danio rerio' finds Danio_rerio_pep.faa)")
        q.add_argument("--metadata", help="TSV of traits and taxonomy, one row per proteome")
        q.add_argument("--group-by", dest="group_by", help="metadata column to group by, e.g. class or habitat")
        q.add_argument("--taxdump", help="NCBI taxdump folder (names.dmp, nodes.dmp): build the taxonomy "
                                         "and common names when there is no --metadata")
        q.add_argument("--threads", type=int, default=4)
        q.add_argument("--hmm", help="defensome.hmm(.gz) to use instead of downloading it")
        q.add_argument("--pfam", help="full Pfam-A.hmm(.gz) to build the defensome HMMs from")
        q.add_argument("--db-url", dest="db_url", help="where to download the prebuilt HMMs from")
        q.add_argument("--no-collapse", action="store_true", dest="no_collapse",
                       help="proteomes already have one protein per gene; skip isoform collapse")
        q.add_argument("--drop-utrorf", action="store_true", dest="drop_utrorf")
        q.add_argument("--review", action="store_true",
                       help="stop after writing the sample sheet, to check it")
        q.add_argument("--trees", help="build gene trees for these families, e.g. CYP,GST (needs MAFFT, FastTree)")
        q.add_argument("--thin", help="with --trees and --metadata: N species per value of a column, e.g. family:1")
        q.add_argument("--per-species", type=int, default=0, dest="per_species")
        q.add_argument("--cyp-refs", dest="cyp_refs", help="clan-labelled CYP reference FASTA")
        q.add_argument("--species-tree", dest="species_tree", help="Newick species tree for the dashboard")
        q.add_argument("--no-figures", action="store_true", dest="no_figures")
        q.add_argument("--map", default=DEFAULT_MAP, help="defensome map TSV")
        q.set_defaults(func=cmd_run)

    run_args(sub.add_parser("run", help="everything, from a folder of proteomes to the dashboard"))
    run_args(sub.add_parser("all", help="the same as run"))

    a = p.parse_args()
    # Most subcommands treat --out as a results directory; cyp-refs writes a
    # single FASTA there. Creating a directory with the file's name silently
    # breaks the write, so only mkdir when --out really is a directory.
    if getattr(a, "out_is_dir", True) and getattr(a, "out", None):
        os.makedirs(a.out, exist_ok=True)
    for opt in ("metadata", "group_by", "proteomes", "species", "families",
                "threads", "force", "refs", "family", "max_seqs", "tree",
                "phylogram", "raw", "by", "max_labels", "species_tree",
                "dashboard_out", "light", "write", "force", "source",
                "table", "orders", "min_len", "max_per_clan", "dir",
                "id_col", "include_fragments", "tables", "clan_hmms",
                "samplesheet", "drop_utrorf", "reference", "per10k",
                "complete_only", "exclude", "require", "no_rescue", "full",
                "nucleotide_col", "min_aa", "max_aa", "max_mb", "samples", "expect",
                "per_species", "where", "clan", "name", "no_sequences", "thin",
                "ligand", "ligand_name", "cofactor", "accelerator", "run", "ga_only", "kernels",
                "hmm", "db_url", "pfam_version", "taxdump", "no_collapse", "review", "trees", "cyp_refs", "no_figures", "fetch", "bundle",
                "rescue_evalue", "full_pfam", "download", "url", "keep_full",
                "db_dir"):
        if not hasattr(a, opt):
            setattr(a, opt, None)
    a.func(a)


if __name__ == "__main__":
    main()
