import os
import sys
import time
import csv
import argparse
import yaml
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pymol import cmd
import requests
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from kneed import KneeLocator
from scipy.cluster.hierarchy import linkage, leaves_list


ions = [
    # Alkaline metals
    "LI", "NA", "K", "RB", "CS",
    # Alkaline-earth metals
    "MG", "CA", "SR", "BA",
    # Transition metals
    "MN", "FE", "CO", "NI", "CU", "ZN",
    "CD", "HG", "PT", "AU", "AG",
    # Rare transition metals and lanthanides
    "Y", "ZR", "MO", "RU", "RH", "PD", "W", "OS", "IR",
    "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
    # Actinides (rare in PDB, but occur)
    "U", "TH",
    # Halides
    "F", "CL", "BR", "I",
    # Other anions
    "SO4", "PO4", "NO3", "CO3", "MO4", "WO4", "SE4"
]


modified_residues = [
    "CSO",  # S-hydroxycysteine – oxidized cysteine (SG-OH)
    "MSE",  # Selenomethionine – selenomethionine (used in structural biology)
    "SEP",  # Phosphoserine – phosphorylated serine
    "TPO",  # Phosphothreonine – phosphorylated threonine
    "PTR",  # Phosphotyrosine – phosphorylated tyrosine

    "HIC",  # Hydroxyisoleucine or Hydroxylysine – hydroxylated isoleucine/lysine
    "HYP",  # Hydroxyproline – hydroxyproline (common in collagen)
    "SEC",  # Selenocysteine – selenocysteine (21st amino acid)
    "PYL",  # Pyrrolysine – pyrrolysine (22nd amino acid, occurs in archaea)

    "CME",  # S-methylcysteine – methylated cysteine
    "FME",  # N-formylmethionine – initiating methionine in prokaryotes
    "CSS",  # Disulfide bond – formalized disulfide bridge record
    "CSD",  # Dithiothreitol-modified cysteine – reduced form of modified cysteine
    "CSX",  # Mixed disulfide – mixed forms of modified cysteines

    "KCX",  # Carboxyglutamate – γ-carboxyglutamate (occurs in calcium-binding proteins)
    "LLP",  # Lysylpyridinoline – crosslinked lysine (in collagen)
    "MLY",  # N-methyl-lysine – methylated lysine
    "MLZ",  # Di- or tri-methylated lysine
    "ALY",  # N-acetyl-lysine – acetylated lysine

    "PCA",  # Pyroglutamate – cyclized glutamine/glutamate at the N-terminus
]

# Buffer components
buffer_residues = [
    "TRS",  # Tris  
    "MES",  # 2-(N-Morpholino)ethanesulfonic acid  
    "HEP",  # HEPES  
    "ACT",  # Acetate  
    "MOP",  # MOPS  
    "PIP",  # PIPES  
    "CAP",  # CAPS  
    "BIS",  # Bis-Tris  
    "TES",  # TES  
    "ADA",  # ADA  
    "AEC",  # ACES  
]

# Cryoprotectants
cryo_residues = [
    "GOL",  # Glycerol  
    "EDO",  # Ethylene glycol  
    "MPD",  # 2-Methyl-2,4-pentanediol  
    "DMS",  # DMSO  
    "BME",  # β-mercaptoethanol  
    "PGO",  # Propylene glycol  
    "EGD",  # Ethylene glycol dimer  
    "PRP",  # n-Propanol  
]

# Sulfates / Phosphates
sulfate_phosphate_residues = [
    "SO4",  # Sulfate  
    "PO4",  # Phosphate  
    "POM",  # Pyrophosphate  
    "P5P",  # Polyphosphate  
]

# Reducing agents
reductant_residues = [
    "DTT",   # Dithiothreitol  
    "BME",   # β-mercaptoethanol  
    "TCEP",  # Tris(2-carboxyethyl)phosphine  
]


def fetch_pdbs_by_uniprot(uniprot_id, output_folder):
    print(f"\n Searching for PDB structures for UniProt ID: {uniprot_id}")
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
        response = requests.get(url)
        if response.status_code != 200:
            print("❌ Failed to retrieve data from PDBe API.")
            return []

        pdb_list = response.json().get(uniprot_id, [])
        pdb_ids = sorted(set(entry["pdb_id"].upper() for entry in pdb_list))  # Remove duplicates

        print(f"✅ Found {len(pdb_ids)} unique PDB structures: {', '.join(pdb_ids)}")

        downloaded = []
        for pdb_id in pdb_ids:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            dest_path = os.path.join(output_folder, f"{pdb_id}.pdb")

            if not os.path.exists(dest_path):
                r = requests.get(url)
                if r.status_code == 200:
                    with open(dest_path, 'w') as f:
                        f.write(r.text)
                    print(f"Downloaded: {pdb_id} to {output_folder}")
                    downloaded.append(pdb_id)
                else:
                    print(f"⚠️ Failed to download: {pdb_id}")
            else:
                print(f" Already exists: {pdb_id}")
                downloaded.append(pdb_id)

        return downloaded

    except Exception as e:
        print(f"❌ Error during download: {e}")
        return []

def fetch_uniprot_by_pdb(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    data = response.json()
    mappings = data.get(pdb_id.lower(), {}).get("UniProt", {})
    if mappings:
        return list(mappings.keys())[0]
    return None

def compute_rmsd(sel1, sel2, method="super"):
    n1 = cmd.count_atoms(sel1)
    n2 = cmd.count_atoms(sel2)

    if n1 == 0 or n2 == 0:
        if "_CA_" not in sel1:
            print(f"⚠️ Skipping RMSD between {sel1.replace('pocket_', '')} | {sel2.replace('pocket_', '')} — empty selection.")
        return None, f"Empty selection: {n1} vs {n2}"
    
    # allow_mismatch = int(n1 * 0.9)
    # if abs(n1 - n2) > allow_mismatch:
    #     if "_CA_" not in sel1:
    #         print(f"⚠️ Skipping RMSD between {sel1} and {sel2} — too different atom counts ({n1} vs {n2}).")
    #     return None, f"Too different atom counts: {n1} vs {n2}"

    try:
        if method in ["rms", "rms_cur"] and n1 != n2:
            msg = f"⚠️ Method {method} requires equal atom counts: {n1} vs {n2}"
            print(msg)
            return None, msg

        if method == "align":
            rmsd = cmd.align(sel1, sel2)[0]
        elif method == "cealign":
            result = cmd.cealign(sel1, sel2)
            rmsd = result['RMSD']
        elif method == "super":
            rmsd = cmd.super(sel1, sel2)[0]
        elif method == "rms":
            rmsd = cmd.rms(sel1, sel2)
        elif method == "rms_cur":
            rmsd = cmd.rms_cur(sel1, sel2)
        else:
            raise ValueError(f"Unknown RMSD method: {method}")

        return rmsd, f"Method: {method}, atoms: {n1} vs {n2}"

    except Exception as e:
        print(f"❌ Error calculating RMSD using method {method} between {sel1} and {sel2}: {e}")
        return None, f"{e} | Method: {method}, atoms: {n1} vs {n2}"

def print_directory_contents_pretty(path, columns=3):
    items = sorted(os.listdir(path))
    print(f"\n Working directory: {path}\n")

    # Add icons and sort
    pretty_items = []
    for name in items:
        full_path = os.path.join(path, name)
        if os.path.isdir(full_path):
            pretty_items.append(f"[📁] {name}")
        else:
            pretty_items.append(f"[📄] {name}")

    # Print in 3 columns
    print(" Directory contents:\n")
    for i in range(0, len(pretty_items), columns):
        row = pretty_items[i:i+columns]
        print("   ".join(f"{x:<40}" for x in row))  # <40 — left-align with width

def clean_structure(name, mode):
    if mode == "1":
        cmd.remove(f"{name} and solvent")
    elif mode == "2":
        ions_selection = " or ".join([f"{name} and resn {ion}" for ion in ions])
        cmd.remove(f"{name} and solvent")
        cmd.remove(ions_selection)
    elif mode == "3":
        cmd.remove(f"{name} and not polymer")
    elif mode == "4":
        pass
    else:
        print("⚠️ Invalid cleaning mode. Skipping.")

def load_yaml_config(config_path):
    print(f"\nAttempting to load YAML configuration from: {os.path.abspath(config_path)}")
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)  # Use safe_load instead of load
            print(f"✅ Configuration loaded:\n{config}\n")
            return config
    except Exception as e:
        print(f"❌ Error loading YAML file: {e}\n")
        return None
        
def parse_args():
    parser = argparse.ArgumentParser(description="Script for analyzing PDB files with RMSD and clustering.")
    parser.add_argument('--config', type=str, help="Path to configuration file (JSON).")
    parser.add_argument('--mode', type=str, choices=[1,2,3], help="Operation mode: 1 - local folder, 2 - UniProt ID, 3 - PDB ID")
    parser.add_argument('--folder_path', type=str, help="Path to folder with PDB files")
    parser.add_argument('--uniprot_id', type=str, help="UniProt ID to download PDB structures")
    parser.add_argument('--pdb_id', type=str, help="PDB ID to search for UniProt ID")
    parser.add_argument('--init_align', type=str, choices=['1', '2'], help="init_align method")
    parser.add_argument('--init_align_chain_id', type=str, help=" init_align_chain_id chain ID (e.g., 'A')")
    parser.add_argument('--do_preprocess', type=str, choices=['0', '1'], help="Perform preprocessing (1 - yes, 0 - no).")
    parser.add_argument('--clean_options', type=str, nargs='*', help="Cleaning options (1-6, 8 to keep).")
    parser.add_argument('--save_dir', type=str, help="Folder to save processed structures.")
    parser.add_argument('--ref_pocket', type=int, help="Index of reference pocket structure (1-based).")
    parser.add_argument('--ref_align', type=int, help="Index of reference alignment structure (1-based).")
    parser.add_argument('--pocket_method', type=str, choices=['1', '2', '3', '4'], help="Pocket definition method.")
    parser.add_argument('--ligand_resi', type=str, help="Residue number (e.g., '40').")
    parser.add_argument('--ligand_chain', type=str, help="chain ID (e.g., 'A').")
    parser.add_argument('--radius', type=float, help="Pocket radius (in Å).")
    parser.add_argument('--resi_chain', type=str, help="List of residues for method 2 (e.g., 'PRO 46 A, ASN 61 A').")
    parser.add_argument('--chain_id', type=str, help="Chain identifier for method 3.")
    parser.add_argument('--comparison_mode', type=str, choices=['1', '2'], help="Comparison mode (1 - all vs all, 2 - all vs ref).")
    parser.add_argument('--rmsd_method', type=str, choices=['1', '2', '3', '4', '5'], help="RMSD calculation method.")
    parser.add_argument('--linkage_method', type=str, choices=['1', '2', '3', '4', '5', '6', '7'], help="Clustering method.")
    parser.add_argument('--cl_choice', type=str, choices=['1', '2', '3'], help="Clustering approach.")
    parser.add_argument('--num_clusters', type=int, help="Number of clusters (for cl_choice 1).")
    parser.add_argument('--distance_threshold', type=float, help="Distance threshold (for cl_choice 2).")
    return parser.parse_args()


def save_config_to_yaml(config_dict, path="run_config.yaml"):
    # Remove None pairs
    cleaned_config = {k: v for k, v in config_dict.items() if v is not None}
    try:
        with open(path, "w") as f:
            yaml.dump(cleaned_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        print(f"\n✅ Configuration saved to file: {path}")
    except Exception as e:
        print(f"❌ Error saving configuration: {e}")


# PyMOL environment cleanup
cmd.reinitialize()

print("\n----------------------------------------------------------------------------------------------------")
print("                                        Loading Arguments                                         ")
print("----------------------------------------------------------------------------------------------------")

args = parse_args()
config = None
if args.config:
    config = load_yaml_config(args.config)
    if not config:
        print("Failed to load configuration (config). Switching to interactive mode.")
        config = {}
else:
    config = {}

# # Fill missing arguments from config file
# for key, value in config.items():
#     if getattr(args, key, None) is None:
#         setattr(args, key, value)

print(f"\n Arguments - args: {args}")
print(f"\n Arguments from Configuration: {config}")

ref_align_name = None
init_align = None
init_align_chain_id = None
failed_align = None
resi_chain = None
chain_id = None
num_clusters = None
distance_threshold = None


print("\n----------------------------------------------------------------------------------------------------")
print("                                        Loading PDB files                                         ")
print("----------------------------------------------------------------------------------------------------")

mode = args.mode or config.get('mode')
folder_path = args.folder_path or config.get('folder_path')
uniprot_id = args.uniprot_id or config.get('uniprot_id')
pdb_id = args.pdb_id or config.get('pdb_id')

# # Show current directory contents
# print_directory_contents_pretty(os.getcwd())

if not mode:
    # interactive input
    print("\nChoose operation mode:")
    print("1 — Local PDB files: Use local folder with PDB files")
    print("2 — UniProt ID: Download all corresponding PDB structures based on UniProt ID")
    print("3 — PDB ID: Determine UniProt ID from PDB ID and download all corresponding PDB structures")
    while True:
        mode = input("Enter mode number (1/2/3): ").strip()
        if mode not in {"1", "2", "3"}:
            print("⚠️ Invalid input. Try again.")
        else:
            mode = int(mode)
            break

mode = str(mode)

if mode == "1":
    if not folder_path:
        print_directory_contents_pretty(os.getcwd())
        folder_path = input("\nEnter path to PDB files folder (default current): ").strip()
        if not folder_path:
            folder_path = os.getcwd()
        print(f"PDB files folder - {folder_path}")

elif mode == "2":
    if not uniprot_id:
        uniprot_id = input("Enter UniProt ID: ").strip()
    folder_path = os.path.join(os.getcwd(), f"{uniprot_id}_pdbs")
    os.makedirs(folder_path, exist_ok=True)
    dwn=fetch_pdbs_by_uniprot(uniprot_id, folder_path)
    if dwn == []:
        print("❌ Failed to find UniProt ID")
        os.rmdir(folder_path)
        exit()
    print(f"✅ UniProt ID: {uniprot_id}")


elif mode == "3":
    if not pdb_id:
        pdb_id = input("🔹 Enter PDB ID: ").strip().lower()
    uniprot_id = fetch_uniprot_by_pdb(pdb_id)
    if not uniprot_id:
        print("❌ Failed to find UniProt ID for the given PDB ID.")
        exit()
    print(f"✅ UniProt ID: {uniprot_id}")
    folder_path = os.path.join(os.getcwd(), f"{uniprot_id}_pdbs")
    os.makedirs(folder_path, exist_ok=True)
    fetch_pdbs_by_uniprot(uniprot_id, folder_path)

# Check that folder exists and read PDB files
if not folder_path or not os.path.exists(folder_path):
    print(f"❌ Folder {folder_path} not found or not specified.")
    exit()

files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
files.sort()

# print(f"\n✅ Files successfully loaded from directory: {folder_path}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                       Preliminary cleaning                                      ")
print("----------------------------------------------------------------------------------------------------")

# --- Ask for preliminary cleaning ---
do_preprocess = args.do_preprocess or config.get('do_preprocess')
if do_preprocess is not None:
    do_preprocess = str(do_preprocess)
if not do_preprocess:
    while True:
        do_preprocess = input("\nPerform preliminary structure cleaning? (1 - yes, 0 - no): ").strip()
        if do_preprocess not in {"1", "0"}:
            print("⚠️ Invalid input. Try again.")
        else:
            break

    if do_preprocess == '1':
        print("\nPlease select appropriate preprocessing options:")
        print("1 — Remove water (solvent)")
        print("2 — Remove (Mg²⁺, Cl⁻, etc.)")
        print("3 — Remove sulfates and phosphates (SO4, PO4, etc.)")
        print("4 — Remove buffer components (TRS, MES, HEP, etc.)")
        print("5 — Remove cryoprotectants (GOL, EDO, MPD, etc.)")
        print("6 — Remove reducing agents (DTT, BME, TCEP)")
        print("7 — Remove all water, ions, buffers, cryoprotectants, phosphates, reducing agents")
        print("8 — Remove modified amino acid residues (CSO, MSE, SEP, TPO, PTR, etc.)")
        print("9 — Remove everything except protein (keep only polymer chain)")
        print("10 — Remove alternative conformations (altloc)")
        print("11 — Remove anisotropic parameters (ANISOU)")
        print("12 — Remove hydrogen atoms (H)")
        print("13 - Save processed structures to a specified folder")
        print("14 — Do not clean / Finish selection")

        save_dir = None
        to_remove = set()
        
        while True:
            choice = input("\nEnter option number (1–13) to select action, or 14 to finish: ").strip()
            if choice == "14":
                break
            elif choice in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]:
                to_remove.add(choice)
                print(f"Option {choice} will be applied.")
            elif choice == "13":
                save_dir = input("Enter path to folder for saving processed structures: ").strip()
                if not os.path.exists(save_dir):
                    try:
                        os.makedirs(save_dir)
                        print(f"Folder not found. Creating new directory: {os.path.abspath(save_dir)}")
                    except Exception as e:
                        print(f"Failed to create folder. Error: {e}")
                        save_dir = None
                else:
                    print(f"Structures will be saved to existing folder: {os.path.abspath(save_dir)}")
            else:
                print("Invalid input. Please enter a number from 1 to 14")

        print("\nSelected options will be applied: " + (", ".join(sorted(to_remove)) if to_remove else "nothing"))
    else:
        save_dir = None
        to_remove = None

else:
    save_dir = args.save_dir or config.get('save_dir')
    if save_dir:
        print(f"\nStructures after preliminary cleaning will be saved to existing folder: {os.path.abspath(save_dir)}")
    else:
        print("\n⚠️ Structures after preliminary cleaning will not be saved to a separate folder. Save folder (\"save_dir\") not specified.")

    to_remove = args.clean_options or config.get('clean_options') or []
    to_remove = set(to_remove)
    print("\nThe following structure processing options will be applied: " + (", ".join(sorted(to_remove)) if to_remove else "--> ⚠️ No preliminary processing options were selected"))


# --- Load and process all structures ---
names = []
for f in files:
    name = os.path.splitext(f)[0]
    names.append(name)
    cmd.load(os.path.join(folder_path, f), name)


    # --- Apply selected cleaning options ---
    if do_preprocess == '1':
        if "1" in to_remove:
            cmd.remove(f"{name} and solvent")
        if "2" in to_remove:
            for ion in ions:
                cmd.remove(f"{name} and resn {ion}")
        if "3" in to_remove:
            for resn in sulfate_phosphate_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "4" in to_remove:
            for resn in buffer_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "5" in to_remove:
            for resn in cryo_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "6" in to_remove:
            for resn in reductant_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "7" in to_remove:
            # Combined list of "unnecessary" ligands
            all_small_residues = set(
                buffer_residues +
                cryo_residues +
                sulfate_phosphate_residues +
                reductant_residues)

            # Remove water
            cmd.remove(f"{name} and solvent")
            
            # Remove ions
            for ion in ions:
                cmd.remove(f"{name} and resn {ion}")
            
            # Remove all other "small" molecules
            for resn in all_small_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "8" in to_remove:
            for resn in modified_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "9" in to_remove:
            cmd.remove(f"{name} and not polymer")
        if "10" in to_remove:
            cmd.remove(f"{name} and not (alt '' + A)")
        if "11" in to_remove:
            ani_state = 0
        if "12" in to_remove:
            cmd.remove(f"{name} and elem H")

        # Save processed structure if path specified
        if save_dir is not None:
            save_path = os.path.join(save_dir, f"{name}_processed.pdb")
            cmd.save(save_path, name, state=0)
            print(f"Saved: {save_path}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                      Selecting reference structures                                    ")
print("----------------------------------------------------------------------------------------------------")

# --- Reference structure selection ---
ref_name = args.ref_pocket or config.get('ref_pocket')

if not ref_name:
    print("\nAvailable PDB files:")
    for i, f in enumerate(files, 1):
        print(f"{i}. {f}")

    ref_index = int(input("\nSelect reference structure for alignment by entering its number from the provided list: ")) - 1
    ref_name = names[ref_index]
else:
    print(f"\nReference structure for alignment: {ref_name}")




# list of all HET groups in selected reference structure
het_atoms = cmd.get_model(f"{ref_name} and not polymer and not solvent").atom
het_residues = sorted({(a.resn, a.resi, a.chain) for a in het_atoms}, key=lambda x: (int(x[1]), x[2]))

# print("\n----------------------------------------------------------------------------------------------------")
# print("                   Preliminary alignment of all structures relative to reference                ")
# print("----------------------------------------------------------------------------------------------------")

# # --- Align all to reference ---
# init_align = args.init_align or config.get('init_align')
# if not init_align:    
#     print(f"\nChoose method for aligning all structures to reference")
#     print("1 - By chain (specify chain ID)")
#     print("2 - Automatic (no chain specification)")
#     init_align = input("\nEnter method (1 or 2): ").strip()
    
#     if init_align == "1":
#         init_align_chain_id = input("\nEnter chain ID, e.g.: A ").strip().upper()
#     elif init_align == "2":
#         print("\nAutomatic selected")
# else:
#     init_align =str(init_align)   
#     if init_align == "1":
#         init_align_chain_id = args.init_align_chain_id or config.get('init_align_chain_id')
#         print(f"\nAligning all structures to reference: 1 - By chain {init_align_chain_id}")
#     elif init_align == "2":
#         print("\nAutomatic selected")

print("\n----------------------------------------------------------------------------------------------------")
print("                         Selecting method for defining alignment region                         ")
print("----------------------------------------------------------------------------------------------------")


# --- Pocket definition method selection ---
ligand_resi = None
ligand_chain = None
pocket_method = args.pocket_method or config.get('pocket_method')
if not pocket_method:
    # Show list of all HET groups 
    if het_residues:
        print(f" \nHET groups detected in structure {ref_name}:")
        for resn, resi, chain in het_residues:
            print(f"  - {resn:>3} {resi:>4} {chain}")
    else:
        print(f"\nNo HET groups detected in structure {ref_name}.")
    
    print(f"\nChoose method for defining alignment region: region is defined ... ")
    print("1 – On ref. structure by specified residue ID and radius (Å), then this region is searched for and aligned in all structures")
    print("2 – After preliminary alignment of all structures, for each structure around the selected residue within a given radius (Å)")
    print("3 – For each structure around its HET groups within a specified radius (Å)")
    print("4 – By user-provided residue list, and then this region is searched for and aligned in all structures.")
    print("5 – On ref. structure by specified chain identifier, then this region is searched for and aligned in all structures.")
    
    while True:
        pocket_method = input("\nEnter method (1, 2, 3, 4 or 5): ").strip()
        if pocket_method == "1":
            while True:
                ligand_input = input("\nEnter residue number and chain ID separated by space, e.g. '40 A': ").strip()
                parts = ligand_input.split()
                if len(parts) != 2:
                    print("Error: enter residue number and chain ID separated by space, e.g. '40 A': Try again.")
                    continue
                ligand_resi, ligand_chain = parts 
                if not ligand_resi.isdigit():
                    print("Error: residue number must be a number. Try again.")
                    continue
                break
            try:
                radius = float(input("\nEnter pocket radius (in Å, default 7): ").strip())
            except ValueError:
                radius = 7.0
            break
        elif pocket_method == "4":
            resi_chain = input("\nEnter residue list in format (resn resi chain), e.g.: PRO 46 A, ASN 61 A: ").strip()
            break
        elif pocket_method == "5":
            chain_id = input("\nEnter chain ID, e.g.: A: ").strip().upper()
            break
        elif pocket_method == "3":
            try:
                radius = float(input("\nEnter radius for pockets around all HET groups (in Å, default 7): ").strip())
            except ValueError:
                radius = 7.0
            chain_id = input("\nEnter chain ID, e.g.: A: ").strip().upper()
            break
        elif pocket_method == "2":

            print("\n----------------------------------------------------------------------------------------------------")
            print("                   Preliminary alignment of all structures relative to reference                ")
            print("----------------------------------------------------------------------------------------------------")

            print("\nAvailable PDB files:")
            for i, f in enumerate(files, 1):
                print(f"{i}. {f}")
            ref_align_index = int(input("\nSelect reference structure for preliminary alignment by entering its number: ")) - 1
            ref_align_name = names[ref_align_index]
        
            # --- Align all to reference ---    
            print(f"\nChoose method for preliminary alignment of all structures to reference")
            print("1 - By chain (specify chain ID)")
            print("2 - By entire molecule")
            init_align = input("\nEnter method (1 or 2): ").strip()
            
            if init_align == "1":
                init_align_chain_id = input("\nEnter chain ID, e.g.: A ").strip().upper()
            elif init_align == "2":
                print("\nBy entire molecule")

            while True:
                ligand_input = input("\nTo define pockets, enter reference residue number and chain ID separated by space, e.g. '40 A': ").strip()
                parts = ligand_input.split()
                if len(parts) != 2:
                    print("Error: enter residue number and chain ID separated by space, e.g. '40 A'. Try again.")
                    continue
                ligand_resi, ligand_chain = parts
                if not ligand_resi.isdigit():
                    print("Error: residue number must be a number. Try again.")
                    continue
                break
            try:
                radius = float(input("\nEnter pocket radius (in Å, default 7): ").strip())
            except ValueError:
                radius = 7.0
            break
        else:
            print("Invalid input. Please enter a number from 1 to 5")


else:
    pocket_method =str(pocket_method)   
    if pocket_method == "1":
        print("1 – On reference structure by specified residue ID and radius (Å), then this region is searched for and aligned in all structures.")
        ligand_resi = args.ligand_resi or config.get('ligand_resi')
        print(f"ligand_resi - {ligand_resi}")
        ligand_chain = args.ligand_chain or config.get('ligand_chain')
        print(f"ligand_chain - {ligand_chain}")
        radius = args.radius or config.get('radius')
        print(f"radius - {radius}")
    elif pocket_method == "4":
        print("4 – By user-provided residue list, and then this region is searched for and aligned in all structures.")
        resi_chain = args.resi_chain or config.get('resi_chain')
        print(f"resi_chain - {resi_chain}")
    elif pocket_method == "5":
        print("5 – On reference structure by specified chain identifier, then this region is searched for and aligned in all structures.")
        chain_id = args.chain_id or config.get('chain_id')
        print(f"chain_id - {chain_id}")
    elif pocket_method == "3":
        print("3 – For each structure around its HET groups within a specified radius (Å).")
        radius = args.radius or config.get('radius')
        chain_id = args.chain_id or config.get('chain_id')
        print(f"chain_id - {chain_id}")
        print(f"radius - {radius}")
    if pocket_method == "2":
        print("2 – After preliminary alignment of all structures, for each structure around the selected residue within a given radius (Å).")
        ref_align_name = args.ref_align or config.get('ref_align')
        print(f"\nReference structure for preliminary alignment: {ref_align_name}")
        init_align = args.init_align or config.get('init_align')
        init_align = str(init_align)   
        if init_align == "1":
            init_align_chain_id = args.init_align_chain_id or config.get('init_align_chain_id')
            print(f"\n Preliminary alignment of all structures to reference: 1 - By chain {init_align_chain_id}")
        elif init_align == "2":
            print("\nBy entire molecule")
        ligand_resi = args.ligand_resi or config.get('ligand_resi')
        print(f"ligand_resi - {ligand_resi}")
        ligand_chain = args.ligand_chain or config.get('ligand_chain')
        print(f"ligand_chain - {ligand_chain}")
        radius = args.radius or config.get('radius')
        print(f"radius - {radius}")
    
print("\n----------------------------------------------------------------------------------------------------")
print("                                        Comparison mode                                             ")
print("----------------------------------------------------------------------------------------------------")

# --- Comparison mode ---
comparison_mode = args.comparison_mode or config.get('comparison_mode')
if not comparison_mode:
    print("\nChoose comparison mode:")
    print("1 - All vs all")
    print("2 - All vs reference")
    comparison_mode = input("\nEnter mode (1 or 2): ").strip()
else:
    comparison_mode = str(comparison_mode)
    if comparison_mode =="1":
            print("\nSelected comparison mode: 1 - All vs all")
    elif comparison_mode =="2":
            print("\nSelected comparison mode: 2 - All vs reference")

print("\n----------------------------------------------------------------------------------------------------")
print("                                          RMSD method                                                ")
print("----------------------------------------------------------------------------------------------------")

# --- RMSD method ---
rmsd_method_input = args.rmsd_method or config.get('rmsd_method')
if not rmsd_method_input:
    print("\nChoose RMSD calculation method:")
    print("1 - align   (Strict RMSD evaluation, used for precise alignment, excluding unmatched atoms)")
    print("2 - cealign (Geometric-based structural alignment, effective even with low similarity.)")
    print("3 - super   (Flexible option, allows mismatches and automatically matches corresponding atoms.)")
    print("4 - rms     (Accurate RMSD, works only with complete atom correspondence (no mismatches))")
    print("5 - rms_cur (Simplified and faster RMSD calculation, also requires full atom correspondence)")
    rmsd_method_input = input("\nEnter method (1, 2, 3, 4 or 5): ").strip()

rmsd_method_input =str(rmsd_method_input)
if rmsd_method_input == "1":
    rmsd_method = "align"
elif rmsd_method_input == "2":
    rmsd_method = "cealign"
elif rmsd_method_input == "3":
    rmsd_method = "super"
elif rmsd_method_input == "4":
    rmsd_method = "rms"
elif rmsd_method_input == "5":
    rmsd_method = "rms_cur"
else:
    raise ValueError("Invalid RMSD calculation method.")

print(f"\nSelected RMSD calculation method: {rmsd_method_input} - {rmsd_method}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                    Clustering method                                             ")
print("----------------------------------------------------------------------------------------------------")

# --- Clustering method ---

linkage_method_choice = args.linkage_method or config.get('linkage_method')
if not linkage_method_choice:
    print("\nChoose hierarchical clustering method:")
    print("1 - ward (minimizes intra-cluster variance, requires Euclidean distance)")
    print("2 - single (minimum distance between clusters)")
    print("3 - complete (maximum distance between clusters)")
    print("4 - average (average distance between clusters (UPGMA))")
    print("5 - centroid (distance between cluster centroids)")
    print("6 - median (median distance between clusters)")
    print("7 - weighted (weighted average distance (WPGMA))")

    help_table = """
    | Method       |  Compact clusters   | Elongated clusters | Outlier sensitivity |
    | ------------ | ------------------- | ------------------ | ------------------- |
    | 1. ward      | ✅ excellent        | ❌ poor            | ⚠️ medium            |
    | 2. single    | ❌ poor             | ✅ excellent       | ⚠️ high              |
    | 3. complete  | ✅ good             | ❌ poor            | ✅ robust           |
    | 4. average   | ✅ universal        | ✅ fair            | ⚠️ medium            |
    | 5. centroid  | ⚠️ unstable          | ⚠️ unstable         | ⚠️ unstable          |
    | 6. median    | ⚠️ unstable          | ⚠️ unstable         | ⚠️ unstable          |
    | 7. weighted  | ✅ good             | ✅ good            | ⚠️ medium            |
    """

    print(help_table)

    linkage_method_choice = input("\nEnter method number (1-7): ").strip()

linkage_method_choice = str(linkage_method_choice)
method_map = {
    "1": "ward",
    "2": "single",
    "3": "complete",
    "4": "average",
    "5": "centroid",
    "6": "median",
    "7": "weighted"
}

if linkage_method_choice in method_map:
    linkage_method = method_map[linkage_method_choice]
    print(f"\nSelected clustering method: {linkage_method_choice} - {linkage_method}")
else:
    linkage_method = "ward"  # default method
    print(f"\nInvalid input, default method selected: 1 - {linkage_method}")

print("\n----------------------------------------------------------------------------------------------------")
print("                         Clustering parameters: number or threshold                              ")
print("----------------------------------------------------------------------------------------------------")

# --- Clustering approach ---
cl_choice = args.cl_choice if args.cl_choice else config.get('cl_choice')
if not cl_choice:
    print("\nChoose clustering parameters:")
    print("1 — Specify number of clusters (maxclust)")
    print("2 — Specify distance threshold (distance)")
    print("3 — Automatic: Threshold based on 70% of maximum merge distance") 
    print("4 — Automatic: Elbow Method")

    cl_choice = input("\nEnter a number (1, 2, 3 or 4): ").strip()

    if cl_choice == '1':
        try:
            num_clusters = int(input("\nEnter desired number of clusters: ").strip())
        except ValueError:
            print("Invalid cluster count.")
    elif cl_choice == '2':
        try:
            distance_threshold = float(input("\nEnter distance threshold (e.g., 3.0): ").strip())
        except ValueError:
            print("Invalid distance threshold.")
    elif cl_choice == '3':
        print(f"Distance threshold will be automatically selected! ")
    elif cl_choice == '4':
        print(f"Distance threshold will be automatically selected! ")
    else:
        print("Invalid choice.")
else:
    cl_choice = str(cl_choice)
    if cl_choice == '1':
        print(f"\nSelected clustering approach: 1 — Specify number of clusters (maxclust)")
        num_clusters = args.num_clusters or config.get('num_clusters')
        print(f"num_clusters: {num_clusters}")
    elif cl_choice == '2':
        print(f"\nSelected clustering approach: 2 — Specify distance threshold (distance)")
        distance_threshold = args.distance_threshold or config.get('distance_threshold')
        print(f"distance_threshold: {distance_threshold}")
    elif cl_choice == "3":
        print(f"\nSelected clustering approach: 3 — Automatic: Threshold based on 70% of maximum merge distance") 
    elif cl_choice == "4":
        print(f"\nSelected clustering approach: 4 — Automatic: Elbow Method") 

#===================================================================================================================================================================
print("\n" + "-" * 100)
print( " " * 40 + "Script starting!" +  " " * 40)
print("-" * 100 + "\n")
#===================================================================================================================================================================

for name in names:
    print(f"✅ {name} Starting processing ")
    time.sleep(0.02)

#--------------------------------------------------------------------
# text_list = [f"✅ {name} Starting processing" for name in names]

# for text in text_list:
#     for char in text:
#         sys.stdout.write(char)
#         sys.stdout.flush()
#         time.sleep(0.005)  # delay between characters 
#     sys.stdout.write("\n")
#     sys.stdout.flush()
#     time.sleep(0.03)  # delay between lines
#--------------------------------------------------------------------

# --- Align all to reference ---

os.makedirs(f"{folder_path}/output", exist_ok=True)
os.makedirs(f"{folder_path}/output/aligned_structures", exist_ok=True)

# --- Get pocket ---

if pocket_method == "1":
    lig_sel = f"{ref_name} and resi {ligand_resi} and chain {ligand_chain}"
    cmd.select("ligand_ref", lig_sel)
    cmd.select("pocket_ref", f"(byres (ligand_ref around {radius})) and {ref_name} and polymer")
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model("pocket_ref").atom})
elif pocket_method == "4":
    resid_list = []
    for res in resi_chain.split(","):
        res = res.strip()
        if not res:
            continue
        try:
            resn, resi, chain = res.split()
            resid_list.append((resn, resi, chain))
        except ValueError:
            print(f"⚠️ Skipping invalid format: {res}")
    if not resid_list:
        raise ValueError("No valid residues specified.")

elif pocket_method == "5":
    chain_sel = f"{ref_name} and chain {chain_id} and polymer"
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model(chain_sel).atom})

elif pocket_method == "3":

    for name in names:
        het_atoms = cmd.get_model(f"{name} and not polymer and not solvent and chain {chain_id}").atom
        unique_residues = set((a.resn, a.resi, a.chain) for a in het_atoms)

        # Filter: keep only HET groups not in modified residues list
        filtered_residues = [
            (resn, resi, chain)
            for (resn, resi, chain) in unique_residues
            if resn not in modified_residues
        ]

        for resn, resi, chain in filtered_residues:
            lig_sel = f"{name} and resn {resn} and resi {resi} and chain {chain}"
            pocket_sel = f"(byres ({lig_sel} around {radius})) and {name} and polymer"
            ca_sel = f"{pocket_sel} and name CA"
            cmd.select(f"pocket_{name}", pocket_sel)
            cmd.select(f"pocket_CA_{name}", ca_sel)

    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model(f"pocket_{name}").atom})

    #         all_het_selections.append(pocket_sel)

    # cmd.select("pocket_ref", " or ".join(all_het_selections))

    # resid_list = list({(a.resn, a.resi, a.chain) for a in cmd.get_model("pocket_ref").atom})
    # resid_list.sort(key=lambda x: (x[2], int(x[1])))

    # print(f"\nCommon pocket formed from all HET groups (radius {radius} Å):")
    # for resn, resi, chain in resid_list:
    #     print(f"  - {resn:>3} {resi:>4} {chain}")

elif pocket_method == "2":
    failed_align = []  # list for proteins that failed to align

    if init_align == "1":
        ref_sel = f"{ref_align_name} and chain {init_align_chain_id}"
        for name in names:
            if name == ref_align_name:
                continue
            mobil_str = f"{name} and chain {init_align_chain_id}"
            try:
                cmd.align(mobil_str, ref_sel)
                # aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{name}_{init_align_chain_id}_aligned_to_{ref_align_name}_{init_align_chain_id}.pdb")
                # cmd.save(aligned_path, mobil_str)
                # print(f"✅ {name}_{init_align_chain_id} saved as {name}_{init_align_chain_id}_aligned_to_{ref_align_name}_{init_align_chain_id}.pdb")
            except:
                print(f"❌ Failed to align {name}_{init_align_chain_id} to {ref_align_name}_{init_align_chain_id}")
                failed_align.append(name)  # Add to failure list

    elif init_align == "2":
        for name in names:
            if name == ref_align_name:
                continue
            mobil_str = f"{name}"
            try:
                cmd.align(mobil_str, ref_align_name)
                # aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{name}_aligned_to_{ref_align_name}.pdb")
                # cmd.save(aligned_path, mobil_str)
                # print(f"✅ {name} saved as {name}_aligned_to_{ref_align_name}.pdb")
            except:
                print(f"❌ Failed to align {name} to {ref_align_name}")
                failed_align.append(name)  # Add to failure list

    lig_sel = f"{ref_name} and resi {ligand_resi} and chain {ligand_chain}"
    cmd.select("ligand_ref", lig_sel)
    cmd.select("pocket_ref", f"(byres (ligand_ref around {radius})) and {ref_name} and polymer")
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model("pocket_ref").atom})
else:
    raise ValueError("Invalid pocket definition method.")

resid_list.sort(key=lambda x: (int(x[1]), x[2]))


if pocket_method in ("1", "4", "5"):
    # --- Create pocket and Cα selections ---
    for name in names:
        parts = [f"(resi {resi} and chain {chain})" for resn, resi, chain in resid_list]
        pocket_sel = f"model {name} and polymer and (" + " or ".join(parts) + ")"
        ca_sel = f"{pocket_sel} and name CA"
        cmd.select(f"pocket_{name}", pocket_sel)
        cmd.select(f"pocket_CA_{name}", ca_sel)

elif pocket_method == "3":
    print("OK")

elif pocket_method == "2":
    for name in names:
        pocket_sel = f"(byres (ligand_ref around {radius})) and {name} and polymer"
        ca_sel = f"{pocket_sel} and name CA"
        cmd.select(f"pocket_{name}", pocket_sel)
        cmd.select(f"pocket_CA_{name}", ca_sel)

# --- Align all to reference ---
# ref_ca = f"pocket_CA_{ref_align_name}"
# os.makedirs(f"{folder_path}/output", exist_ok=True)
# os.makedirs(f"{folder_path}/output/aligned_structures", exist_ok=True)

# failed_align = []  # list for proteins that failed to align
# for name in names:
#     if name == ref_align_name:
#         continue
#     moving_ca = f"pocket_CA_{name}"
#     try:
#         cmd.align(moving_ca, ref_ca)
#         aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{name}_aligned_to_{ref_align_name}.pdb")
#         cmd.save(aligned_path, name)
#         print(f"✅ {name} saved as {name}_aligned_to_{ref_align_name}.pdb")
#     except:
#         print(f"❌ Failed to align {name} to {ref_align_name}")
#         failed_align.append(name)  # Add to failure list



print(f"{ref_name} - Pocket residues:", ', '.join([f"{resn} {resi} {chain}" for resn, resi, chain in resid_list]))
# Form PyMOL selection
pymol_sel_parts = [f"(resi {resi} and chain {chain})" for resn, resi, chain in resid_list]
pymol_sel = f"{ref_name} and polymer and (" + " or ".join(pymol_sel_parts) + ")"


# Commands for beautiful visualization
vis_code = f'''
# Highlight and color pocket
select pocket, {pymol_sel}
show sticks, pocket
color red, pocket

# Display Cα atoms as spheres
select pocket_CA, pocket and name CA
show spheres, pocket_CA
set sphere_scale, 0.3, pocket_CA
color tv_red, pocket_CA
set stick_radius, 0.15, selection=pocket


# Rest of protein — semi-transparent surface
select rest_protein, {ref_name} and polymer and not pocket
show surface, rest_protein
set transparency, 0.6, rest_protein
color gray80, rest_protein
'''
print("\nPyMOL commands for beautiful pocket visualization:\n")
print(vis_code.strip())
print(" \n ")


# --- RMSD calculation ---
rmsd_all = []
rmsd_ca = []

if comparison_mode == "1":
    # All vs all
    n = len(names)
    heatmap_matrix = np.full((n, n), np.nan)
    heatmap_labels = names

    for i in range(n):
        for j in range(n):
            if i == j:
                heatmap_matrix[i, j] = 0.0
                continue

            sel1_all = f"pocket_{names[i]}"
            sel2_all = f"pocket_{names[j]}"
            sel1_ca = f"pocket_CA_{names[i]}"
            sel2_ca = f"pocket_CA_{names[j]}"

            rms_all, reason_all = compute_rmsd(sel1_all, sel2_all, rmsd_method)
            rms_ca, reason_ca = compute_rmsd(sel1_ca, sel2_ca, rmsd_method)

            if sel2_all == f"pocket_{ref_name}":
                aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{names[i]}_aligned_to_{ref_name}.pdb")
                cmd.save(aligned_path, sel1_all)
                print(f"✅ {names[i]} saved as {names[i]}_aligned_to_{ref_name}.pdb")


            rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
            rms_ca_s = f"{rms_ca:.5f}" if rms_ca is not None else ""

            # print(f"{names[i]} vs {names[j]} | RMSD (all): {rms_all_s} | RMSD (Cα): {rms_ca_s}")
            
            rmsd_all.append((names[i], names[j], rms_all_s, reason_all))
            rmsd_ca.append((names[i], names[j], rms_ca_s, reason_ca))

            if rms_all is not None:
                heatmap_matrix[i, j] = rms_all
                heatmap_matrix[j, i] = rms_all  # symmetric

    # Remove rows and columns completely filled with NaN (except diagonal)
    def is_row_invalid(row, idx):
        return all(np.isnan(v) for j, v in enumerate(row) if j != idx)

    invalid_indices = [i for i, row in enumerate(heatmap_matrix) if is_row_invalid(row, i)]

    if invalid_indices:
        print("⚠️ Structures removed (alignment failed):")
        for i in invalid_indices:
            print(f"  - {heatmap_labels[i]}")

        heatmap_matrix = np.delete(heatmap_matrix, invalid_indices, axis=0)
        heatmap_matrix = np.delete(heatmap_matrix, invalid_indices, axis=1)
        heatmap_labels = [name for i, name in enumerate(heatmap_labels) if i not in invalid_indices]
        n = heatmap_matrix.shape[0]
#----------
    # 1) Create DataFrame from matrix for convenience
    df = pd.DataFrame(heatmap_matrix, index=heatmap_labels, columns=heatmap_labels)

    # 2) Fill NaN with maximum value (to avoid interfering with clustering)
    max_val = np.nanmax(df.values)
    df_filled = df.fillna(max_val)

    # 3) linkage for rows and columns
    row_Z = linkage(df_filled.values, method=linkage_method, metric='euclidean')
    col_Z = linkage(df_filled.values.T, method=linkage_method, metric='euclidean')

    # 4) Get leaf order
    row_order = leaves_list(row_Z)
    col_order = leaves_list(col_Z)

    # 5) Reorder DataFrame
    df_ord = df_filled.iloc[row_order, :].iloc[:, col_order]
    ordered_labels = df_ord.index.tolist()

    # --- Add mean values ---
    df_ord['Mean'] = df_ord.mean(axis=1)       # row means
    mean_row = df_ord.mean(axis=0)             # column means
    df_ord.loc['Mean'] = mean_row

    # --- Plot ordered heatmap ---
    plt.figure(figsize=(12,10))
    sns.heatmap(
        df_ord,
        cmap='coolwarm',
        annot=True, fmt=".2f",
        annot_kws={"size": 4},
        linewidths=0.5,
        cbar_kws={"label":"RMSD (Å)"}
    )
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title(f"Clustered RMSD heatmap ({linkage_method})")
    plt.tight_layout()

    heatmap_path = os.path.join(folder_path, "output", "rmsd_heatmap_clustered.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"Ordered heatmap saved as: rmsd_heatmap_clustered.png")

    # # Plot unordered heatmap
    # plt.figure(figsize=(12, 10))

    # # Calculate average RMSD per row (ignoring diagonal and nan)
    # mean_rmsd = []
    # for i in range(n):
    #     values = [heatmap_matrix[i, j] for j in range(n) if i != j and not np.isnan(heatmap_matrix[i, j])]
    #     avg = np.mean(values) if values else np.nan
    #     mean_rmsd.append(avg)

    # # Add mean as last column
    # heatmap_matrix_with_mean = np.hstack((heatmap_matrix, np.array(mean_rmsd).reshape(-1, 1)))
    # heatmap_labels_with_mean = heatmap_labels + ["Mean"]

    # sns.heatmap(heatmap_matrix_with_mean, xticklabels=heatmap_labels_with_mean, yticklabels=heatmap_labels,
    #             cmap='coolwarm', annot=True, fmt=".2f", linewidths=0.5, cbar_kws={"label": "RMSD (Å)"})

    # plt.title(f"RMSD (all atoms) Heatmap – Method: {rmsd_method}")
    # plt.xticks(rotation=45, ha='right')
    # plt.yticks(rotation=0)
    # plt.tight_layout()

    # heatmap_path = os.path.join(folder_path, "output", "rmsd_heatmap.png")
    # plt.savefig(heatmap_path, dpi=300)
    # print(f"\nRMSD heatmap saved as: rmsd_heatmap.png")

    # Prepare for dendrogram
    df = pd.DataFrame(heatmap_matrix, index=heatmap_labels, columns=heatmap_labels)


    # Remove rows and columns where more than 90% values are NaN
    threshold = 0.7  # 70%
    print(df)
    row_nan_fraction = df.isna().mean(axis=1)
    col_nan_fraction = df.isna().mean(axis=0)

    df_filtered = df.loc[row_nan_fraction < threshold, df.columns[col_nan_fraction < threshold]]

    print("Maximum RMSD value:", df_filtered.values.max())
    print("Minimum RMSD value:", df_filtered.values.min())


    valid_values = df_filtered.mask(df_filtered == 0).stack()  # exclude diagonal (0.0)
    max_rmsd = valid_values.max()

    df_filtered.replace([np.inf, -np.inf], np.nan, inplace=True)

    nan_locs = np.argwhere(np.isnan(df_filtered.values))
    for i, j in nan_locs:
        row_name = df_filtered.index[i]
        col_name = df_filtered.columns[j]
        print(f"🔍 NaN between structures: {row_name} ↔ {col_name}")

    df_filtered = df_filtered.fillna(max_rmsd)

    if df_filtered.shape[0] < 2:
        print("Insufficient data for clustering (less than two structures).")
    else:
        condensed = squareform(df_filtered.values, checks=False)
        Z = linkage(condensed, method=f'{linkage_method}')

        clusters = None
        clustering_description = ""

        if cl_choice == '1':
            try:
                clusters = fcluster(Z, t=num_clusters, criterion='maxclust')
                clustering_description = f"{num_clusters}_clusters"
            except ValueError:
                print("Invalid cluster count.")
        elif cl_choice == '2':
            try:
                clusters = fcluster(Z, t=distance_threshold, criterion='distance')
                clustering_description = f"threshold_{distance_threshold}"
            except ValueError:
                print("Invalid distance threshold.")

        elif cl_choice == '3':
            # Threshold based on 70% of maximum merge distance
            optimal_threshold = 0.7 * max(Z[:, 2])

            clusters = fcluster(Z, t=optimal_threshold, criterion='distance')
            clustering_description = f"knee_threshold_{optimal_threshold:.2f}"
            print(f"Threshold selected based on 70% of maximum merge distance: {optimal_threshold:.2f}")

        elif cl_choice == '4':
            # Automatic selection - Threshold using elbow method
            distances = Z[:, 2]
            kneedle = KneeLocator(range(1, len(distances)+1), distances, curve="convex", direction="increasing")
            optimal_threshold = distances[kneedle.knee]

            clusters = fcluster(Z, t=optimal_threshold, criterion='distance')
            clustering_description = f"knee_threshold_{optimal_threshold:.2f}"
            print(f"Threshold selected using Elbow method: {optimal_threshold:.2f}")
        else:
            print("Invalid choice.")

        plt.figure(figsize=(10, 7))
        dendrogram(Z, labels=df_filtered.index.tolist())
        plt.title(f"Hierarchical clustering based on RMSD ({linkage_method})")
        plt.tight_layout()
        dendro_path = os.path.join(folder_path, "output", f"RMSD_dendrogram_{linkage_method}.png")
        plt.savefig(dendro_path)
        plt.close()
        print(f"Dendrogram saved to: RMSD_dendrogram_{linkage_method}.png")

        if clusters is not None:
            cluster_df = pd.DataFrame({
                'Protein': df_filtered.index,
                'Cluster': clusters
            })
            output_csv_path = os.path.join(folder_path, "output", f"cluster_assignments_{clustering_description}_{linkage_method}.csv")
            cluster_df.to_csv(output_csv_path, index=False)
            print(f"Clustering results saved to: cluster_assignments_{clustering_description}_{linkage_method}.csv")

elif comparison_mode == "2":
    # All vs reference
    for name in names:
        if name == ref_align_name:
            continue
        sel1_all = f"pocket_{ref_align_name}"
        sel2_all = f"pocket_{name}"
        sel1_ca = f"pocket_CA_{ref_align_name}"
        sel2_ca = f"pocket_CA_{name}"

        rms_all, reason_all = compute_rmsd(sel1_all, sel2_all, rmsd_method)
        rms_ca, reason_ca = compute_rmsd(sel1_ca, sel2_ca, rmsd_method)
        
        # Format RMSD for nice output
        rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
        rms_ca_s = f"{rms_ca:.5f}" if rms_ca is not None else ""

        rmsd_all.append((ref_align_name, name, rms_all_s, reason_all))
        rmsd_ca.append((ref_align_name, name, rms_ca_s, reason_ca))

        # print(f"{ref_align_name} vs {name} | RMSD (all): {rms_all_s} | RMSD (Cα): {rms_ca_s}")

else:
    raise ValueError("Invalid comparison mode.")


# --- Save RMSD ---
with open(f"{folder_path}/output/rmsd_all_atoms.csv", "w", newline='', encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["Protein1", "Protein2", "RMSD_all_atoms", "Comment"])
    writer.writerows(rmsd_all)

with open(f"{folder_path}/output/rmsd_calpha.csv", "w", newline='', encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["Protein1", "Protein2", "RMSD_Calpha", "Comment"])
    writer.writerows(rmsd_ca)

# --- Save pocket information ---
with open(f"{folder_path}/output/info.txt", "w", encoding="utf-8-sig") as info_file:
    info_file.write(f"Alignment structure: {ref_name}\n")
    if pocket_method == "1":
        info_file.write(f"Alignment region definition method: On ref. structure (radius {radius} Å from residue {ligand_resi} {ligand_chain})\n")
    elif pocket_method == "2":
        info_file.write(f"Alignment region definition method: After preliminary alignment, around (radius {radius} Å from residue {ligand_resi} {ligand_chain})\n")
    elif pocket_method == "3":
        info_file.write(f"Alignment region definition method: For each structure around its HET groups (within {radius} Å radius)\n")
    elif pocket_method == "4":
        info_file.write(f"Alignment region definition method: By user-provided residue list - {resi_chain} \n")
    elif pocket_method == "5":
        info_file.write(f"Alignment region definition method: On ref. structure by specified chain ID {chain_id}\n")

    if pocket_method == '2':
        info_file.write(f"Reference structure for preliminary alignment: {ref_align_name}\n")
        if init_align == "1":
            info_file.write(f"Preliminary alignment method: by chain {init_align_chain_id}\n")
        elif init_align == "2":
            info_file.write(f"Preliminary alignment method: by entire molecule\n")

    info_file.write(f"RMSD calculation method: {rmsd_method}")

    info_file.write(f"\nSelected clustering method: {linkage_method}")

    if cl_choice == '1':
        info_file.write(f"\nSelected clustering approach: 1 — Specify number of clusters (maxclust)")
        info_file.write(f"\nnum_clusters: {num_clusters}")
    elif cl_choice == '2':
        info_file.write(f"\nSelected clustering approach: 2 — Specify distance threshold (distance)")
        info_file.write(f"\ndistance_threshold: {distance_threshold}")
    elif cl_choice == "3":
        info_file.write(f"\nSelected clustering approach: 3 — Automatic: Threshold based on 70% of maximum merge distance") 
    elif cl_choice == "4":
        info_file.write(f"\nSelected clustering approach: 4 — Automatic: Elbow or Knee Method") 

    if comparison_mode == "1":
        info_file.write("\nSelected comparison mode: 1 - All vs all")
    elif comparison_mode == "2":
        info_file.write("\nSelected comparison mode: 2 - All vs reference")


    if failed_align:
        info_file.write("\nProteins that failed to align:\n")
        for f_name in failed_align:
            info_file.write(f"- {f_name}\n")
    else:
        info_file.write("\n\nAll proteins successfully aligned.\n")
    
    info_file.write(f"\n\nPocket information (resn resi chain):\n\n")
    for name in names:
        pocket_atoms = cmd.get_model(f"pocket_{name}").atom
        residues = {(atom.resi, atom.resn, atom.chain) for atom in pocket_atoms}
        residues = sorted(residues, key=lambda x: (int(x[0]), x[2]))
        residue_lines = [f"{resn} {resi} {chain}" for resi, resn, chain in residues]
        info_file.write(f"{name}:\n")
        info_file.write("  " + ', '.join(residue_lines) + "\n\n")


    # --- Compare pockets to reference ---
    info_file.write("\nPocket comparison to reference structure:\n\n")
    
    # Form dictionary for reference: (resi, chain) → resn
    ref_residues = {
        (atom.resi, atom.chain): atom.resn
        for atom in cmd.get_model(f"pocket_{ref_name}").atom
    }

    results = []

    # Collect all diffs in list
    for name in names:
        if name == ref_name:
            continue
        target_atoms = cmd.get_model(f"pocket_{name}").atom
        target_residues = {
            (atom.resi, atom.chain): atom.resn
            for atom in target_atoms
        }

        diffs = []
        for key in ref_residues:
            ref_resn = ref_residues[key]
            tgt_resn = target_residues.get(key)
            if tgt_resn is None:
                diffs.append(f"{key[0]} {key[1]}: {ref_resn} → missing")
            elif tgt_resn != ref_resn:
                diffs.append(f"{key[0]} {key[1]}: {ref_resn} → {tgt_resn}")

        results.append((name, diffs))

    # Sort: first those with len(diffs)==0, then by len(diffs) ascending
    sorted_results = sorted(
        results,
        key=lambda x: (len(x[1]) != 0, len(x[1]))
    )

# Write to file in desired order
    for name, diffs in sorted_results:
        if not diffs:
            info_file.write(f"{name}: identical to reference pocket.\n\n")
        else:
            info_file.write(f"{name}: differs in {len(diffs)} residues:\n")
            for d in diffs:
                info_file.write(f"  {d}\n")
            info_file.write("\n")

# --- Plot RMSD histograms ---
try:
    import matplotlib.pyplot as plt
except ImportError:
    print("⚠️ matplotlib module not found. Install with: pip install matplotlib")
    plt = None

try:
    import seaborn as sns
except ImportError:
    print("⚠️ seaborn module not found. Install with: pip install seaborn")
    sns = None

# Set style
sns.set(style="whitegrid")

if plt and sns:
    def plot_rmsd_histogram(rmsd_values, title, filename):
        values = [float(val[2]) for val in rmsd_values if val[2]]
        if not values:
            print(f"⚠️ No data for histogram: {filename}")
            return

        plt.figure(figsize=(10, 6))
        sns.histplot(values, bins=30, kde=True, color='cornflowerblue', edgecolor='black')

        # Vertical line for mean RMSD
        mean_val = np.mean(values)
        plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean = {mean_val:.2f} Å')
        plt.legend(fontsize=12)
        
        plt.title(title, fontsize=16, fontweight='bold')
        plt.xlabel("RMSD (Å)", fontsize=14)
        plt.ylabel("Frequency", fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        # Add X-axis labels with rounded RMSD values
        min_val, max_val = min(values), max(values)
        plt.xticks([round(x, 2) for x in 
                    list(plt.xticks()[0]) 
                    if min_val <= x <= max_val])

        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f"{folder_path}/output/{filename}", dpi=300)
        plt.close()
        print(f"Histogram saved: {filename}")
    plot_rmsd_histogram(rmsd_all, "RMSD Histogram (all atoms)", "rmsd_all_atoms_hist.png")
    plot_rmsd_histogram(rmsd_ca, "RMSD Histogram (Cα)", "rmsd_calpha_hist.png")

# Collect configuration
config = {
    "mode": mode,
    "folder_path": folder_path,
    "uniprot_id": uniprot_id,
    "pdb_id": pdb_id,
    "save_dir": save_dir,
    "do_preprocess": do_preprocess,
    "clean_options": to_remove,
    "ref_pocket": ref_name,
    "ref_align": ref_align_name,
    "init_align": init_align,
    "init_align_chain_id": init_align_chain_id,
    "pocket_method": pocket_method,
    "ligand_resi": ligand_resi,
    "ligand_chain": ligand_chain,
    "radius": radius,
    "resi_chain": resi_chain,
    "chain_id": chain_id,
    "comparison_mode": comparison_mode,
    "rmsd_method": rmsd_method_input,
    "linkage_method": linkage_method_choice,
    "cl_choice": cl_choice,
    "num_clusters": num_clusters,
    "distance_threshold": distance_threshold,
}



# Save
if not args.config:
    save_config_to_yaml(config, f"{folder_path}/output/run_config.yaml")


print("\n✅ Done. Files saved:")
print("-----------------------------")
print(f"{folder_path}/output/")
print("- RMSD all atoms: rmsd_all_atoms.csv")
print("- RMSD Cα: rmsd_calpha.csv")
print("- Pocket residues (resn resi chain): info.txt")
print("- Aligned structures: aligned_structures/ folder")
print(f"- RMSD histogram (all atoms): rmsd_all_atoms_hist.png")
print(f"- RMSD histogram (Cα): rmsd_calpha_hist.png")
if comparison_mode == "1": 
    print(f"- Heatmap: rmsd_heatmap.png")
    print(f"- Dendrogram: RMSD_dendrogram_{linkage_method}.png ")
    print(f"- Clustering results: cluster_assignments_{clustering_description}_{linkage_method}.csv")
print("\n")