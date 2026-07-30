#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find residues in contact with the Center of Mass (COM) of catalytic sites 
for distance cutoffs of 5, 10, and 15 Angstroms.
"""

import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch, Selection

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────
PDBS_PATH     = "selected_pdbs"
INPUT_CSV     = "catalytic_sites_w_decoys.csv"
OUTPUT_CSV    = "contact_residues_byCOM.csv"
CONTACT_RADII = [5, 10, 15]   # Angstroms


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────
def parse_act_site_list(x):
    """Parse semicolon-separated residue IDs into a list of ints."""
    if isinstance(x, list):
        return x
    return [int(i) for i in str(x).split(';') if i.strip()]


def get_neighbor_indices(row, pdbs_path, radius):
    """
    Calculate the Center of Mass (COM) of catalytic residues and find
    neighboring residues within the specified radius.
    
    Args:
        row: Pandas Series containing 'Entry' and 'ACT_SITE_list'.
        pdbs_path: Base directory for PDB files.
        radius: Distance cutoff in Angstroms.
        
    Returns:
        Sorted 0-based indices of residues near the COM.
    """
    entry          = row['Entry']
    target_res_ids = row['ACT_SITE_list']   # 1-based PDB numbering
    pdb_path       = os.path.join(pdbs_path, f"{entry}.pdb")

    if not os.path.exists(pdb_path):
        return []

    try:
        parser    = PDBParser(QUIET=True)
        structure = parser.get_structure(entry, pdb_path)

        # 1. Collect all atoms belonging to the catalytic site
        target_atoms = [
            atom
            for model in structure
            for chain in model
            for residue in chain
            if residue.id[1] in target_res_ids
            for atom in residue.get_atoms()
        ]

        if not target_atoms:
            return []

        # 2. Calculate Center of Mass (COM)
        # We assume uniform mass for simplicity, taking the mean of coordinates
        coords = np.array([atom.coord for atom in target_atoms])
        center_of_mass = coords.mean(axis=0)

        # 3. Search for neighbors around the COM
        all_atoms = Selection.unfold_entities(structure, 'A')
        ns = NeighborSearch(all_atoms)

        # search(center, radius, level='R') returns residues
        neighbors = ns.search(center_of_mass, radius, level='R')

        contact_residues = set()
        for res in neighbors:
            # Exclude the catalytic residues themselves from the neighbor list
            if res.id[1] not in target_res_ids:
                contact_residues.add(res)

        # Return sorted 0-based indices
        return sorted(res.id[1] - 1 for res in contact_residues)

    except Exception as exc:
        print(f"  [ERROR] {entry}: {exc}")
        return []


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    print(f"Loading {INPUT_CSV} …")
    df = pd.read_csv(INPUT_CSV)
    df['ACT_SITE_list'] = df['ACT_SITE_list'].apply(parse_act_site_list)

    dataset = df[["Entry", "EC number", "ACT_SITE_list"]].copy()
    print(f"  {len(dataset):,} rows loaded.\n")

    for radius in CONTACT_RADII:
        col = f"contacts_{radius}A"
        print(f"── Cutoff {radius} Å (from COM) ────────────────")

        dataset[col] = dataset.apply(
            lambda row: get_neighbor_indices(row, PDBS_PATH, radius),
            axis=1
        )

        found = dataset[col].apply(len) > 0
        print(f"  Rows with contacts found  : {found.sum():,} / {len(dataset):,}")
        avg_n = dataset[col].apply(len).mean()
        print(f"  Avg contacts per protein  : {avg_n:.1f}\n")

        # Serialize list as semicolon-separated string for CSV
        dataset[col] = dataset[col].apply(lambda x: ';'.join(map(str, x)))

    # Also serialize ACT_SITE_list back to string
    dataset['ACT_SITE_list'] = dataset['ACT_SITE_list'].apply(
        lambda x: ';'.join(map(str, x))
    )

    print(f"Saving results to {OUTPUT_CSV} …")
    dataset.to_csv(OUTPUT_CSV, index=False)
    print("Done.")


if __name__ == "__main__":
    main()
