# db/

`defensome.py setup` writes the defensome HMM database here by default
(`defensome.hmm`, `MANIFEST.tsv`, `PFAM_VERSION`). Nothing here is tracked by
git.

For machines without internet, put a downloaded `defensome_db.tar.gz` (from
the repository's `pfam-db` release) in this folder: `run` and `setup` use it
before trying to download.
