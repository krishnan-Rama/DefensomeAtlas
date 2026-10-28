#!/usr/bin/env python3
"""Known-answer fixture for the rescue pass, validation rules and MT screen."""
import os, random, sys
random.seed(11)
root = sys.argv[1]; os.makedirs(f"{root}/peps", exist_ok=True)
aa = "ACDEFGHIKLMNPQRSTVWY"
bg = lambda n: "M" + "".join(random.choice(aa) for _ in range(n - 1))
P = {
 # Drosophila MtnB, 43 aa, 12 Cys, no aromatic residues (Wikipedia, MT family 5)
 "mt_real_1":  "MVCKGCGTNCQCSAQKCGDNCACNKDCQCVCKNGPKDQCCSNK",
 # synthetic MT-like: Cys-rich, no aromatics, no HMM hit -> composition only
 "mt_real_2":  "MPCPCGSGCKCASQATKMSCGKCCTDLKDDLKCKCGCSGKDNK",
 # insect defensin-like: 6 Cys / 40 aa, has Y and H -> must be rejected
 "defensin":   "ATCDLLSGTGINHSACAAHCLLRGNRGGYCNGKAVCVCRN",
 # Cys-rich, motif-dense, but aromatic: passes a Cys-only rule, must fail
 "cys_arom":   "MWCKCYCGFCKCWCGYCSCFCKWCGYCQCFCGWCKYCGNKDNK",
 # long protein with a spurious sub-GA MT hit -> must fail len<=120
 "mt_fp_long": bg(300),
 # MATE with two MatE repeats below GA -> rescued (copies>=2)
 "mate_two":   bg(560),
 # MATE-like with a single repeat below GA -> must fail copies>=2
 "mate_one":   bg(560),
 # an ordinary CYP found at GA
 "cyp_ga":     bg(510),
}
for i in range(40):
    P[f"bg{i:02d}"] = bg(random.randint(150, 600))
with open(f"{root}/peps/Sp_test.faa", "w") as f:
    for k, v in P.items():
        f.write(f">{k}\n{v}\n")
# a minimal HMMER3 text database with the accessions the map uses
accs = [("PF00067","p450",459),("PF01554","MatE",161),("PF02067","Metallothio_5",42),
        ("PF03074","Glu_cys_ligase",373),("PF00199","Catalase",382)]
with open(f"{root}/Pfam-A.hmm", "w") as f:
    for acc, name, L in accs:
        f.write(f"HMMER3/f [3.4 | Aug 2023]\nNAME  {name}\nACC   {acc}.1\nDESC  {name} fixture\n"
                f"LENG  {L}\nGA    25.00 25.00;\nHMM stub\n//\n")
print("fixture written")
