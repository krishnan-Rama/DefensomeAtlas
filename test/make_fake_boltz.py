#!/usr/bin/env python3
"""Stand-in for a finished Boltz-2 run in a `structure` folder: the files
Boltz writes (model CIF, confidence and affinity JSON) at the paths it writes
them, with known values, so the dashboard's reading of them can be checked
without a GPU. Usage: make_fake_boltz.py <structure folder>"""
import json, math, os, sys

d = sys.argv[1].rstrip("/")
stem = next(f[:-5] for f in os.listdir(d) if f.endswith(".yaml"))
seq = next(l.split(":", 1)[1].strip() for l in open(os.path.join(d, stem + ".yaml"))
           if l.strip().startswith("sequence:"))
pred = os.path.join(d, f"boltz_results_{stem}", "predictions", stem)
os.makedirs(pred, exist_ok=True)
three = dict(zip("ACDEFGHIKLMNPQRSTVWYX", "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR UNK".split()))
cols = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id",
        "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
        "auth_seq_id", "auth_asym_id"]
out = ["data_model", "loop_"] + [f"_atom_site.{c}" for c in cols]
n = 0
# an ideal alpha helix of backbone atoms (radius, phase offset in degrees, rise)
BB = [("N", "N", 1.55, -28.0, -0.86), ("C", "CA", 2.30, 0.0, 0.0),
      ("C", "C", 1.65, 24.6, 0.67), ("O", "O", 1.85, 21.0, 1.84)]
for i, aa in enumerate(seq, 1):
    b = 92 if i < len(seq) * .8 else 45     # a confident core and a floppy tail
    for el, nm, r, ph, dz in BB:
        t = math.radians(i * 100 + ph)
        n += 1
        out.append(f"ATOM {n} {el} {nm} {three.get(aa, 'UNK')} A {i} {r * math.cos(t):.3f} "
                   f"{r * math.sin(t):.3f} {1.5 * i + dz:.3f} 1.0 {b} {i} A")
for k, (el, x) in enumerate([("FE", 0.0), ("C", 1.4), ("C", 2.8)]):
    n += 1
    out.append(f"HETATM {n} {el} {el}{k} HEM B . {8 + x:.3f} 0.000 20.000 1.0 90 1 B")
open(os.path.join(pred, f"{stem}_model_0.cif"), "w").write("\n".join(out) + "\n")
json.dump({"confidence_score": 0.81, "ptm": 0.77, "iptm": 0.69, "ligand_iptm": 0.66,
           "protein_iptm": 0.0, "complex_plddt": 0.84, "complex_iplddt": 0.8,
           "complex_pde": 0.9, "complex_ipde": 1.2},
          open(os.path.join(pred, f"confidence_{stem}_model_0.json"), "w"))
json.dump({"affinity_pred_value": 0.5, "affinity_probability_binary": 0.73},
          open(os.path.join(pred, f"affinity_{stem}.json"), "w"))
print(pred)
