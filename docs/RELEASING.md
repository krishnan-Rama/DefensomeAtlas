# Publishing the HMM database (maintainers)

Users never need Pfam-A: `run` and `setup` download a small prebuilt database
from this repository's GitHub release. This page is how to make and publish it.

## What lives where

| file | where | why |
|---|---|---|
| `defensome_map.tsv` | the repository root, in git | it is code: versioned with `defensome.py`, reviewed in pull requests. `package.sh` also embeds it in `defensome.py`, so a lone copy of the script still works |
| `defensome_db.tar.gz` (+ `.sha256`) | a **GitHub Release** with the tag `pfam-db`, as release assets | it is data from one Pfam release; binaries do not belong in git history, and a release gives it a fixed download URL |

The URL the tool downloads from is

```
https://github.com/krishnan-Rama/DEFENSOME/releases/download/pfam-db/defensome_db.tar.gz
```

set by `DB_REPO`, `DB_TAG` and `DB_ASSET` near the top of `defensome.py`
(`DEFENSOME_DB_URL` overrides it at run time). Keep the tag and file name
fixed and replace the assets when Pfam moves on, so the URL never breaks.

## Build it

On any machine with Python and pandas (HMMER is not needed for this step):

```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
python3 defensome.py setup --pfam Pfam-A.hmm.gz --db-dir release_db --bundle
```

or in one step, `python3 defensome.py setup --download --db-dir release_db --bundle`.

This writes `release_db/defensome_db.tar.gz`, holding `defensome.hmm` (only the
map's models), `MANIFEST.tsv` (which accessions were found, their lengths and
gathering thresholds) and `PFAM_VERSION`, plus `defensome_db.tar.gz.sha256`.
Read the list it prints: any accession reported **MISSING** has been retired
or merged in that Pfam release and its family will be empty until the map is
updated.

If you already have a `defensome.hmm` from a run (for example `results/db/`),
`setup --pfam` on the same Pfam file rebuilds it identically; there is no need
to download Pfam again if it is still on disk.

## Publish it

**On the website:** Releases → *Draft a new release* → *Choose a tag*: type
`pfam-db`, *Create new tag* → title `Defensome HMMs (Pfam 37.x)` → paste the
contents of `release_db/PFAM_VERSION` as the notes → attach
`defensome_db.tar.gz` and `defensome_db.tar.gz.sha256` → *Publish release*.
Tick *Set as a pre-release* if you do not want it shown as the latest code
release; the download URL works either way.

**With the GitHub CLI:**

```bash
gh release create pfam-db release_db/defensome_db.tar.gz release_db/defensome_db.tar.gz.sha256 \
    --repo krishnan-Rama/DEFENSOME --title "Defensome HMMs (Pfam 37.x)" \
    --notes-file release_db/PFAM_VERSION --prerelease
```

**Updating to a new Pfam release:**

```bash
gh release upload pfam-db release_db/defensome_db.tar.gz release_db/defensome_db.tar.gz.sha256 \
    --repo krishnan-Rama/DEFENSOME --clobber
gh release edit pfam-db --repo krishnan-Rama/DEFENSOME --title "Defensome HMMs (Pfam 38.0)"
```

Results record the Pfam release they used (`results/db/PFAM_VERSION`), so an
older run stays interpretable after the asset is replaced. To keep old
releases downloadable too, also publish each under its own tag (`pfam-37.4`)
and point users at it with `--db-url`.

## Check it

```bash
python3 defensome.py setup --fetch --db-dir /tmp/check_db --force
cat /tmp/check_db/PFAM_VERSION
```

## Licence

Pfam is distributed under CC0, so the models can be redistributed. Cite Pfam
or InterPro in the release notes and in papers.

## Code releases

Tag code releases separately (`v5.2.0`). GitHub attaches a source archive to
each automatically; run `bash package.sh` first so the embedded map and
dashboard in `defensome.py` match the files beside it and the tests pass.
