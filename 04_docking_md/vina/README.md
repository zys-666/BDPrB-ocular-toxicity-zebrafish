# 04_docking_md/vina — docking inputs and poses

Structures used for AutoDock Vina and for building the MD starting complexes.

| File | Content |
|---|---|
| `receptor_zf_ALDH1A2_26-518.pdb` | zebrafish ALDH1A2, AlphaFold model AF-Q90XS8-F1-model_v6, residues 26-518 (493 aa) |
| `ligand_BDPrB.pdb`, `ligand_BDPrB.sdf` | BDPrB, PubChem CID 22184028, bond orders preserved |
| `KyA33_reference_pose.pdb` | co-crystallised inhibitor A1IDQ (KyA33) extracted from PDB 9R6N, mapped onto the zebrafish model |
| `out_BDPrB_pose1.pdb`, `.sdf` | best-scoring BDPrB pose (score -7.3 kcal/mol) |
| `out_KyA33_redock_pose1.pdb` | KyA33 redocked into 9R6N, heavy-atom RMSD 0.52 A |
| `out_KyA33_crystal_vs_redock.pdb` | crystal pose superimposed on the redocked pose |
| `complex_BDPrB.pdb` | protein-ligand complex used as the MD starting structure |

Not included
------------
* Vina configuration files and `.pdbqt` conversions. Regenerate them with
  `../scripts/dock_prepare_receptor.py` using the box and scoring parameters in
  Methods 2.9 (centre -11.9, 8.9, -10.9; size 16.1 x 25.6 x 17.6 A;
  exhaustiveness 32; num_modes 20; seeds 42, 1 and 7).
* MD trajectories (`.xtc`, `.trr`, `.tpr`, `.edr`), excluded by `.gitignore`.

Third-party structures
----------------------
PDB 9R6N (human ALDH1A2 with A1IDQ) is distributed by the RCSB PDB. The AlphaFold
model is distributed by the AlphaFold Protein Structure Database under CC BY 4.0.
Both are included here as derived coordinates for reproducibility; cite the
original entries (Methods 2.9, references 31-33).
