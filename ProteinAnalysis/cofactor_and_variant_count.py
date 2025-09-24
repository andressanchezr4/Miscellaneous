import math
import os
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
