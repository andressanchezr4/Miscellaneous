import math
import os
import pandas as pd
import gzip 
from Bio.PDB.MMCIFParser import MMCIFParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import is_aa
from Bio.SeqUtils import molecular_weight
import gzip
import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

path2cif = '/home/andres/Desktop/mmCIF/'
cofactors =  ['TDP', 'FAD', 'FMN', 'NAD', 'PNS', 'COA', 'PLP', 'GSH', 'BTN', 'FFO', 'B12', 
                 'ASC', 'MQ7', 'UQ1', 'MGD', 'H4B', 'MDO', 'SAM', 'F43', 'COM', 'TP7', 'HEA', 
                 'DPM', 'PQQ', 'TPQ', 'TRQ', 'LPA', 'HEM', 'SAH', 'TPP', 'BIA']
import requests

# your starting list
user_list = ['TDP', 'FAD', 'FMN', 'NAD', 'PNS', 'COA', 'PLP', 'GSH', 'BTN', 'FFO', 'B12',
             'ASC', 'MQ7', 'UQ1', 'MGD', 'H4B', 'MDO', 'SAM', 'F43', 'COM', 'TP7', 'HEA',
             'DPM', 'PQQ', 'TPQ', 'TRQ', 'LPA', 'HEM', 'SAH', 'TPP', 'BIA']

api_url = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors"

resp = requests.get(api_url, timeout=30)
resp.raise_for_status()
data = resp.json()

all_cof = []
for entry in data:
    all_cof.extend(data[entry][0].get('cofactors', []))

cofactors = set(cofactors + all_cof)
metal_ions = ['BA', 'SR', 'CD', 'CR', 'MG', 'ZN', 'CU', 'CA', 'FE', 'MN', 'NI', 'CO', 'MO', 'K', 'NA', 'CL', 'V', 'W', 'SE'] 

substitutions = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', '5OW':'LYS', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MK8':'LEU', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'MLY': 'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR',
    'CAS': 'CYS', 'MSE': 'MET', 'SEP': 'SER', 'M3L': 'LYS', 'MLZ': 'LYS'
}


def uniprot_chain_dict(df, pdb_id):
    """
    Build dictionary {uniprot: [chains]} for a given PDB ID.
    """
    # Filter rows for this PDB
    pdb_df = df[df["PDB"] == pdb_id]

    # Group by UniProt, collect chains (unique)
    return (
        pdb_df.groupby("SP_PRIMARY")["CHAIN"]
        .apply(lambda x: sorted(set(x)))  # remove duplicates, sort for consistency
        .to_dict()
    )

def process_pdb(pdb):
    """'1n2t.cif.gz'
    Procesa un solo archivo PDB y devuelve filas con [pdb_id, cofactor, metals_combination]."""
    try:
        with gzip.open(path2cif + pdb.lower(), "rt") as f:
            mmcif_dict = MMCIF2Dict(f)
            
        with gzip.open(path2cif + pdb.lower(), 'rt') as f:
            parser = MMCIFParser(QUIET=True, auth_chains=True)
            struct = parser.get_structure(pdb, f)[0]
            
        all_residues = mmcif_dict.get("_atom_site.label_comp_id", [])
        aa_residues = [n for n,i in enumerate(all_residues) if i in cofactors]
        
        label_id = [mmcif_dict.get("_atom_site.label_comp_id", [])[aa_residues[n]] for n in range(len(aa_residues))]
        auth_chains = [mmcif_dict.get("_atom_site.auth_asym_id", [])[aa_residues[n]] for n in range(len(aa_residues))]
        label_chains = [mmcif_dict.get("_atom_site.label_asym_id", [])[aa_residues[n]] for n in range(len(aa_residues))]
        # set(mmcif_dict.get("_atom_site.label_asym_id", [])[:aa_residues])
        
        # Build a mapping: auth_asym_id -> set of label_asym_id(s)
        # auth_to_label = {}
        correspondence = set()
        for auth, label, id in zip(auth_chains, label_chains, label_id):
            correspondence.add((auth, label, id))
            # print(auth, label, id)
        
        chains_of_interest = uniprot_chain_dict(uniprot2pdb_raw, pdb.split('.')[0])

        conn_type = mmcif_dict.get("_struct_conn.conn_type_id", [])
        ptnr1_label_comp_id = mmcif_dict.get("_struct_conn.ptnr1_label_comp_id", [])
        ptnr2_label_comp_id = mmcif_dict.get("_struct_conn.ptnr2_label_comp_id", [])
        ptnr1_label_asym_id = mmcif_dict.get("_struct_conn.ptnr1_label_asym_id", [])
        ptnr2_label_asym_id = mmcif_dict.get("_struct_conn.ptnr2_label_asym_id", [])
        ptnr1_label_seq_id = mmcif_dict.get("_struct_conn.ptnr1_label_seq_id", [])
        ptnr2_label_seq_id = mmcif_dict.get("_struct_conn.ptnr2_label_seq_id", [])
        # con_df = pd.DataFrame([conn_type, ptnr1_label_comp_id,
        #                        ptnr2_label_comp_id, ptnr1_label_asym_id,
        #                        ptnr2_label_asym_id, ptnr1_label_seq_id,
        #                        ptnr2_label_seq_id])
        
        covalent_cofactors = []
        for ctype, comp1, comp2, chain1, chain2, id1, id2 in zip(
                conn_type, ptnr1_label_comp_id, ptnr2_label_comp_id,
                ptnr1_label_asym_id, ptnr2_label_asym_id,
                ptnr1_label_seq_id, ptnr2_label_seq_id):
        
            if ctype == "covale":  # mmCIF tag for covalent bond
                if comp1 in cofactors:
                    covalent_cofactors.append(f'{comp1}_{id1}_{chain1}')
                if comp2 in cofactors:
                    covalent_cofactors.append(f'{comp2}_{id2}_{chain2}')
        
        covalent_cofactors = set(covalent_cofactors) 
        
        het_res = [res for res in struct.get_residues() if res.id[0][2:].upper() in cofactors]
        met_het_res = [res for res in struct.get_residues() if res.id[0][2:].upper() in metal_ions]
        
        ns_prot = NeighborSearch([atom for res in struct.get_residues() if res.id[0] == ' ' for atom in res or res.id[0][2:] == substitutions.values()])
        
        corrected_res = set()
        for cc in covalent_cofactors:
            if '.' in cc:
                cc_id = cc.split('_')
                for h_r in het_res:
                    h_r_auth_chain = h_r.get_parent().id
                    h_r_label_chain = [c[1] for c in correspondence if h_r_auth_chain in c]
                    if h_r_label_chain:
                        h_r_label_chain = h_r_label_chain[0]
                    else:
                        continue
                    if h_r.id[0][2:].upper() == cc_id[0] and h_r_label_chain == cc_id[2]:
                        corrected_res.add(f'{h_r.id[0][2:].upper()}_{h_r.id[1]}_{h_r.get_parent().id}')
        
        covalent_cofactors = covalent_cofactors | corrected_res      
        
        for h_r in het_res:
            for atom in h_r:
                if len(ns_prot.search(atom.coord, 1.5, level="A")) > 0:
                    # print(h_r)
                    covalent_cofactors.add(f'{h_r.id[0][2:].upper()}_{h_r.id[1]}_{h_r.get_parent().id}')
                    break
                
        toreturn = []
        for key, value in chains_of_interest.items():
            
            # if not het_res + met_het_res:
            #     return None
            
            for chain in struct.get_chains():
                
                if chain.id not in value:
                    continue
                
                atoms_chain = [atom for residue in chain
                              if residue.id[0] == ' '
                              for atom in residue]
                ns = NeighborSearch(atoms_chain)
                for cofactor_res in het_res:
                    cofactor_name = f'{cofactor_res.id[0][2:].upper()}_{cofactor_res.id[1]}_{cofactor_res.get_parent().id}'
                    
                    if cofactor_name in covalent_cofactors:
                        cofactor_name = f'{cofactor_res.id[0][2:].upper()}*_{cofactor_res.id[1]}_{cofactor_res.get_parent().id}'
                    
                    atoms_cofactor = [atom for atom in cofactor_res]
                    
                    # Verificar si cofactor está cerca de la proteína
                    if not any(len(ns.search(atom.coord, 5.0, level="A")) > 0 for atom in atoms_cofactor):
                        continue
                    
                    # Buscar metales cercanos al cofactor
                    ns_cofactor = NeighborSearch(atoms_cofactor)
                    nearby_metals = []
                    
                    for metal_res in met_het_res:
                        metal_name = f'{metal_res.id[0][2:].upper()}_{metal_res.id[1]}_{metal_res.get_parent().id}'
                        atoms_metal = [atom for atom in metal_res]
                        if any(len(ns_cofactor.search(atom_m.coord, 4.1, level="A")) > 0 for atom_m in atoms_metal):
                            if metal_name not in nearby_metals:
                                nearby_metals.append(metal_name)
                                
                    if '*' in cofactor_name and cofactor_name.split('_')[-1] != chain.id:
                        continue
                    
                    toreturn.append([key, pdb.lower().split('.')[0], chain.id, cofactor_name, ",".join(nearby_metals)])
                              
            return toreturn
        
    except Exception as e:
        print(f"Error procesando {pdb}: {e}")
        return None

all_pdbs = os.listdir(path2cif)
n = len(all_pdbs)
chunk_size = math.ceil(n / 3) 

results = []
path2uniprot2pdb_raw = '/home/andres/Desktop/prueba_andres_1/uniprot2pdb_raw.csv'
uniprot2pdb_raw = pd.read_csv(path2uniprot2pdb_raw)
for i in range(0, n, chunk_size):
    pdb_chunk = all_pdbs[i:i+chunk_size]
    print(f"Procesando chunk {i//chunk_size + 1} de 3 ({len(pdb_chunk)} PDBs)")
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_pdb, pdb) for pdb in pdb_chunk]
        for future in as_completed(futures):
            rows = future.result()
            if rows:
                results.extend(rows)
              
from collections import defaultdict
cofactors =  ['TDP', 'FAD', 'FMN', 'NAD', 'PNS', 'COA', 'PLP', 'GSH', 'BTN', 'FFO', 'B12', 
                 'ASC', 'MQ7', 'UQ1', 'MGD', 'H4B', 'MDO', 'SAM', 'F43', 'COM', 'TP7', 
                 'DPM', 'PQQ', 'TPQ', 'TRQ', 'LPA', 'HEM', 'TPP']

cofactor_to_parent = dict()
for parent_name, entries in data.items():
    for entry in entries:
        for c in cofactors:
            if c in entry.get("cofactors", []):
                cofactor_to_parent[c] = cofactor_to_parent.get(c, []) + entry.get("cofactors", [])
                    
df = pd.DataFrame(results, columns=['uniprot_id', "pdb_id", 'chain', "cofactor", 'metal']).drop_duplicates()

df["cofactor"] = df["cofactor"].str.replace(r'_.*$', '', regex=True)
df["metal"] = df["metal"].apply(
    lambda x: ",".join(re.sub(r'_.*$', '', c.strip()) for c in x.split(","))
)

# Filtrar metales
allowed_metals = ['MG', 'ZN', 'CU', 'CA', 'FE', 'MN', 'NI', 'CO', 'MO', 'V', 'W', 'BA', 'SR', 'CD', 'CR']
def process_metals(s):
    metals = s.rstrip(',').split(',')
    filtered = [m for m in metals if m in allowed_metals]  # keep only allowed
    unique_metals = sorted(filtered)
    return ','.join(unique_metals)

df['metal'] = df['metal'].apply(process_metals)

def parent_adder(row):
    cof = str(row).strip()[:3]  # first 3 letters
    for parent, child in cofactor_to_parent.items():
        if cof in child and '*' in str(row).strip():  # compare first 3 letters
            print(f"{cof}* → {parent}*")  # see the mapping
            return parent +'*'
        elif cof in child:
            print(f"{cof} → {parent}")
            return parent
    return ''  # if no match

df["parent_cofactor"] = df["cofactor"].apply(parent_adder)

summary = (
    df.groupby(['parent_cofactor', 'metal'])['pdb_id']
    .nunique()
    .reset_index()
    .rename(columns={'pdb_id': 'count'})
)

print(summary)  
summary.to_csv('/home/andres/Desktop/pdbs_ALL_PARENT_cofactors_and_metals_covalent.csv')

from collections import defaultdict
cofactors =  ['TDP', 'FAD', 'FMN', 'NAD', 'PNS', 'COA', 'PLP', 'GSH', 'BTN', 'FFO', 'B12', 
                 'ASC', 'MQ7', 'UQ1', 'MGD', 'H4B', 'MDO', 'SAM', 'F43', 'COM', 'TP7', 
                 'DPM', 'PQQ', 'TPQ', 'TRQ', 'LPA', 'HEM', 'TPP']

cofactor_to_parent = dict()
for parent_name, entries in data.items():
    for entry in entries:
        for c in cofactors:
            if c in entry.get("cofactors", []):
                cofactor_to_parent[c] = cofactor_to_parent.get(c, []) + entry.get("cofactors", [])
                    
df = pd.DataFrame(results, columns=['uniprot_id', "pdb_id", 'chain', "cofactor", 'metal']).drop_duplicates()

# Normalizar cofactor y metal
df["cofactor"] = df["cofactor"].str.replace(r'_.*$', '', regex=True)
df["metal"] = df["metal"].apply(
    lambda x: ",".join(re.sub(r'_.*$', '', c.strip()) for c in x.split(","))
)

# Filtrar metales
allowed_metals = ['MG', 'ZN', 'CU', 'CA', 'FE', 'MN', 'NI', 'CO', 'MO', 'V', 'W', 'BA', 'SR', 'CD', 'CR']
def process_metals2(s):
    metals = s.rstrip(',').split(',')
    filtered = [m for m in metals if m in allowed_metals]  # keep only allowed
    unique_metals = sorted(set(filtered))
    return ','.join(unique_metals)

df['metal'] = df['metal'].apply(process_metals2)

def parent_adder(row):
    cof = str(row).strip()[:3]  # first 3 letters
    for parent, child in cofactor_to_parent.items():
        if cof in child and '*' in str(row).strip():  # compare first 3 letters
            print(f"{cof} → {parent}")  # see the mapping
            return parent +'*'
        elif cof in child:
            return parent
    return ''  # if no match

df["parent_cofactor"] = df["cofactor"].apply(parent_adder)

summary_unique = (
    df.groupby(['parent_cofactor', 'metal'])['pdb_id']
    .nunique()
    .reset_index()
    .rename(columns={'pdb_id': 'count'})
)

print(summary_unique)  

summary_unique.to_csv('/home/andres/Desktop/pdbs_ALL_PARENT_cofactors_and_metals_unique_covalent.csv')
df.to_csv('/home/andres/Desktop/pdbs_ALL_PARENT_cofactors_and_metals_covalent_withoutcounts.csv')



import re
df = pd.DataFrame(results, columns=['uniprot_id', "pdb_id", 'chain', "cofactor", 'metal']).drop_duplicates()
df["cofactor"] = df["cofactor"].str.replace(r'_.*$', '', regex=True)
df["metal"] = df["metal"].apply(
    lambda x: ",".join(re.sub(r'_.*$', '', c.strip()) for c in x.split(","))
)

# === LIST OF COFACTORS ===
cofactor_list = ['TDP', 'FAD', 'FMN', 'NAD', 'PNS', 'COA', 'PLP', 'GSH', 'BTN', 'FFO', 'B12',
                 'ASC', 'MQ7', 'UQ1', 'MGD', 'H4B', 'MDO', 'SAM', 'F43', 'COM', 'TP7', 'HEA',
                 'DPM', 'PQQ', 'TPQ', 'TRQ', 'LPA', 'HEM', 'SAH', 'TPP', 'BIA']

# === Allowed metals ===
allowed_metals = ['MG', 'ZN', 'CU', 'CA', 'FE', 'MN', 'NI', 'CO', 'MO', 'V', 'W', 'BA', 'SR', 'CD', 'CR'] 

def process_metals(s):
    metals = s.rstrip(',').split(',')
    filtered = [m for m in metals if m in allowed_metals]  # keep only allowed
    unique_metals = sorted(filtered)
    return ','.join(unique_metals)

df['metal'] = df['metal'].apply(process_metals)


# df_check = df_exploded[df_exploded['metal_combination'] != '']

# --- Step 3: Group by (cofactor, metal_combination) and count pdbs ---
summary = (
    df.groupby(['cofactor', 'metal'])['pdb_id']
    .nunique()
    .reset_index()
    .rename(columns={'pdb_id': 'count'})
)

print(summary)
summary.to_csv('/home/andres/Desktop/pdbs_ALL_cofactors_and_metals_covalent.csv')

df = pd.DataFrame(results, columns=['uniprot_id', "pdb_id", 'chain', "cofactor", 'metal']).drop_duplicates()
df["cofactor"] = df["cofactor"].str.replace(r'_.*$', '', regex=True)
df["metal"] = df["metal"].apply(
    lambda x: ",".join(re.sub(r'_.*$', '', c.strip()) for c in x.split(","))
)

def process_metals2(s):
    metals = s.rstrip(',').split(',')
    filtered = [m for m in metals if m in allowed_metals]  # keep only allowed
    unique_metals = sorted(set(filtered))
    return ','.join(unique_metals)

df['metal'] = df['metal'].apply(process_metals2)

summary_unique = (
    df.groupby(['cofactor', 'metal'])['pdb_id']
    .nunique()
    .reset_index()
    .rename(columns={'pdb_id': 'count'})
)

print(summary_unique)
summary_unique.to_csv('/home/andres/Desktop/pdbs_ALL_cofactors_and_metals_unique_covalent.csv')

df.to_csv('/home/andres/Desktop/pdbs_ALL_cofactors_and_metals_covalent_withoutcounts.csv')
