#!/usr/bin/env python3
"""Transcripts carrying metallothioneins the ORF caller would have discarded."""
import os, random, sys
random.seed(21)
code = {"M":"ATG","V":"GTG","C":"TGC","K":"AAG","G":"GGC","T":"ACC","N":"AAC","Q":"CAG",
        "S":"AGC","A":"GCC","D":"GAC","P":"CCC","R":"CGC","L":"CTG","I":"ATC","H":"CAC",
        "Y":"TAC","E":"GAG","F":"TTC","W":"TGG"}
bt = lambda p: "".join(code[x] for x in p) + "TAA"
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
nt = lambda n: "".join(random.choice("ACGT") for _ in range(n))
MtnB = "MVCKGCGTNCQCSAQKCGDNCACNKDCQCVCKNGPKDQCCSNK"        # real Drosophila MtnB
MT2 = "MPCPCGSGCKCASQATKMSCGKCCTDLKDDLKCKCGCSGKDNK"
DEF = "MATCDLLSGTGINHSACAAHCLLRGNRGGYCNGKAVCVCRN"            # defensin decoy
LONG = "M" + "".join(random.choice("ACDEGHIKLMNPQRSTV") for _ in range(160))
root = sys.argv[1]; os.makedirs(root, exist_ok=True)
tx = {"TRINITY_DN10_c0_g1_i1": nt(120) + bt(MtnB) + nt(200),
      "TRINITY_DN10_c0_g1_i2": nt(80) + bt(MtnB) + nt(90),
      "TRINITY_DN22_c0_g1_i1": nt(150) + rc(bt(MT2)) + nt(140),
      "TRINITY_DN33_c0_g1_i1": nt(100) + bt(DEF) + nt(100),
      "TRINITY_DN44_c0_g1_i1": nt(60) + bt(LONG) + nt(60)}
for i in range(200):
    tx[f"TRINITY_DN{1000+i}_c0_g1_i1"] = nt(random.randint(300, 2500))
with open(f"{root}/Trinity.fasta", "w") as f:
    for k, v in tx.items():
        f.write(f">{k} len={len(v)}\n{v}\n")
with open(f"{root}/samplesheet.tsv", "w") as f:
    f.write("sample_id\tspecies\tmethod\tfasta\tgene_regex\tnucleotide\n")
    f.write(f"S1_denovo\tS1\tdenovo\t/x.faa\t_i\\d+\\w*$\t{os.path.abspath(root)}/Trinity.fasta\n")
    f.write("S1_DToL\tS1\tDToL\t/y.faa\t\t\n")
