# 04_docking_md

Configuration files and analysis scripts for Methods 2.9 (molecular docking and molecular dynamics, Figure 4 and Figures S8–S9).

Trajectories are **not** included: the raw MD data set is 9.2 GB. What is provided here is everything
needed to repeat the runs and to re-derive every reported metric from a trajectory.

## Contents

| Folder | Contents |
|---|---|
| `mdp/` | GROMACS run inputs for the three systems (`em`, `ions`, `nvt`, `npt`, `md`), one set per system: `controls_apo_*`, `controls_kya33_*`, and the BDPrB complex |
| `vina/` | AutoDock Vina docking logs and the MMPBSA input/decomposition tables |
| `scripts/` | Trajectory-analysis Python scripts (MDAnalysis-based) and shell drivers used for the production runs |

## Docking (AutoDock Vina 1.1.2)

Receptor: AlphaFold2 model AF-Q90XS8-F1-model_v6 (zebrafish ALDH1A2), residues 1–25 removed
(flexible, low-pLDDT segment clashing with the substrate channel); residues 26–518 used.
Human reference: PDB 9R6N.

Search box 16.1 × 25.6 × 17.6 Å, centre (−11.9, 8.9, −10.9); `exhaustiveness = 32`, `num_modes = 20`;
three seeds (42, 1, 7) were used to test reproducibility. KyA33 was redocked into 9R6N as a positive
control (heavy-atom RMSD < 2.0 Å criterion).

## Molecular dynamics (GROMACS 2025.1)

Three systems (apo zebrafish ALDH1A2, BDPrB complex, KyA33 complex), 100 ns each.
AMBER14SB protein force field, GAFF2/AM1-BCC ligands, TIP3P water, dodecahedral box, ≥1.0 nm from
the solute, neutralised with 0.15 mol/L NaCl. Minimisation → 500 ps NVT (V-rescale, 310 K) →
500 ps NPT (C-rescale, 1 bar) → production with a 2-fs timestep, LINCS on hydrogen bonds, PME
electrostatics (1.0 nm real-space cutoff), 1.0 nm van der Waals cutoff with long-range dispersion
correction. Frames saved every 50 ps (2,000 frames per system).

Run order per system:

```bash
gmx grompp -f mdp/<system>_ions.mdp  -c <em.gro>   -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
gmx grompp -f mdp/<system>_em.mdp    -c ions.gro   -p topol.top -o em.tpr
gmx mdrun  -v -deffnm em
gmx grompp -f mdp/<system>_nvt.mdp   -c em.gro     -r em.gro -p topol.top -o nvt.tpr
gmx mdrun  -v -deffnm nvt
gmx grompp -f mdp/<system>_npt.mdp   -c nvt.gro    -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun  -v -deffnm npt
gmx grompp -f mdp/<system>_md.mdp    -c npt.gro    -t npt.cpt -p topol.top -o md.tpr
gmx mdrun  -v -deffnm md
```

## Binding free energy (gmx_MMPBSA)

MM-PBSA over the last 50 ns of each trajectory (168 frames for the BDPrB complex, 167 for KyA33),
APBS solver, solute dielectric 2.0, solvent dielectric 78.54, mBondi radii, no entropy term.
Frames in which the ligand crossed a periodic boundary were detected from an anomalous minimum
distance and re-extracted.

## Trajectory analysis

Backbone RMSD, radius of gyration, per-residue Cα RMSF, ligand RMSD relative to the docked pose,
ligand–protein minimum distance, pocket-residue SASA and ligand–protein hydrogen bonds.
Contact occupancy is the fraction of frames in which any heavy atom of a residue lies within 4.0 Å of
any ligand heavy atom (300-ps sampling, 334 frames per system). Metrics are reported for the
50–100 ns segment.
