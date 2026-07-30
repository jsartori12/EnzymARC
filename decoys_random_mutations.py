#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate N decoy sequences per protein by randomly mutating specific residues.
Documentation in English as requested.
"""

import pandas as pd
import random
from tqdm import tqdm

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────
CONTACTS_CSV = "contact_residues_byCOM.csv"
SITES_CSV    = "catalytic_sites_w_decoys.csv"
OUTPUT_CSV   = "decoys_random.csv"
N_DECOYS     = 10
CUTOFFS      = [5, 10, 15]

STANDARD_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# ─────────────────────────────────────────────
# CORE FUNCTIONS
# ─────────────────────────────────────────────

def generate_random_decoys(sequence, contact_positions, n_decoys):
    """
    Creates decoys by replacing residues at specified positions with 
    a random amino acid from the standard pool.
    
    Args:
        sequence (str): Original protein sequence.
        contact_positions (list): 0-based indices to mutate.
        n_decoys (int): Number of decoys to generate.
        
    Returns:
        list: Mutated sequence strings.
    """
    if not contact_positions:
        return [None] * n_decoys

    decoys = []
    for _ in range(n_decoys):
        seq_list = list(sequence)
        for pos in contact_positions:
            # Ensure the mutation actually changes the residue (optional)
            # current_res = seq_list[pos]
            # new_res = random.choice([aa for aa in STANDARD_AA if aa != current_res])
            
            new_res = random.choice(STANDARD_AA)
            seq_list[pos] = new_res
            
        decoys.append("".join(seq_list))
    
    return decoys

def parse_index_list(x):
    """Converts semicolon-separated strings or lists into a list of integers."""
    if isinstance(x, list):
        return x
    if pd.isna(x) or str(x).strip() == '':
        return []
    return [int(i) for i in str(x).split(';') if i.strip()]

# ─────────────────────────────────────────────
# MAIN EXECUTION
# ─────────────────────────────────────────────

def main():
    print("Loading data...")
    df_contacts = pd.read_csv(CONTACTS_CSV)
    df_sites = pd.read_csv(SITES_CSV)

    # Prepare output structure
    results = []

    print(f"Generating {N_DECOYS} random decoy(s) per cutoff...")
    
    for _, row in tqdm(df_sites.iterrows(), total=len(df_sites)):
        pdb_id = row['Entry']
        original_seq = row['sequence']
        
        # Filter contact data for this specific PDB
        pdb_contacts = df_contacts[df_contacts['Entry'] == pdb_id]
        
        row_results = row.to_dict()

        for cutoff in CUTOFFS:
            # Extract positions for the current cutoff
            col_name = f"contacts_{cutoff}A"
            if col_name in pdb_contacts.columns:
                raw_pos = pdb_contacts.iloc[0][col_name]
                positions = parse_index_list(raw_pos)
                
                # Generate mutations
                decoys = generate_random_decoys(original_seq, positions, N_DECOYS)
                
                # Store results (handling N_DECOYS > 1 if necessary)
                for i, d in enumerate(decoys):
                    suffix = f"_{i}" if N_DECOYS > 1 else ""
                    row_results[f"decoy_{cutoff}A{suffix}"] = d
            else:
                row_results[f"decoy_{cutoff}A"] = None

        results.append(row_results)

    # Save to CSV
    out_df = pd.DataFrame(results)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Done! Saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
