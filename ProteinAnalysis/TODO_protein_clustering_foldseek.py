#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 17:10:27 2026

@author: andres
"""

import os
import subprocess
import requests
import pandas as pd
from pathlib import Path
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import MMCIF2Dict
from Bio.PDB import Structure, Model, Chain, Residue, Atom, ShrakeRupley
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from requests import get

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

def selected_residues(chains2keep, struct, res2sel = '', het_sel = False):
    
    # select the chain and its surrounding het res
    atoms_prot = [atom for residue in struct.get_residues() 
                  if (residue.id[0] == ' ' or residue.id[0][2:] in substitutions) and residue.get_parent().id.lower() == chains2keep.lower() 
                  for atom in residue]
    ns = NeighborSearch(atoms_prot)
    
    selected_residues = []
    for res in struct.get_residues():
        
        if het_sel and res.id[0].startswith('H_'):
            tot_atm = 0
            close_atm = 0
            for atom in res:
                if len(ns.search(atom.coord, 5, level="R")) > 0:
                    close_atm += 1
                tot_atm += 1
            if close_atm / tot_atm > 0.3:
                selected_residues.append(res)
        
        if (res.id[0] == ' ' or res.id[0][2:] in list(substitutions)) and res.get_parent().id.lower() == chains2keep.lower():
            selected_residues.append(res)
            
        elif res2sel == res:
            selected_residues.append(res)
            
    for chain in struct.get_chains():
        for residue in list(chain):
            if residue not in selected_residues:
                chain.detach_child(residue.id)

    return struct

# def download_pdb(pdb_id, uniprot, out_dir):
#     url = f"https://files.rcsb.org/download/{pdb_id}.cif"
#     out_file = out_dir + f"/{uniprot}_{pdb_id}.cif"
    
#     r = requests.get(url)
#     if r.status_code == 200:
#         with open(out_file, 'w') as out_file_write:
#             out_file_write.write(r.text)
#         print(f"⬇ Descargado {pdb_id}")
#     else:
#         print(f"❌ Error descargando {pdb_id}")
    
    
#     uni_chain = uniprot2pdb_raw[(uniprot2pdb_raw.SP_PRIMARY == uniprot) & (uniprot2pdb_raw.PDB == pdb_id.lower())]['CHAIN'].tolist()
#     for chain in uni_chain:
#         parser = MMCIFParser(QUIET=True, auth_chains=False)
#         struct = parser.get_structure(pdb_id, out_file)
#         final_path = out_dir + f"/{uniprot}_{pdb_id}_{chain}.cif"
#         struct2save = selected_residues(chain, struct, res2sel='', het_sel = True)
#         io = MMCIFIO()
#         io.set_structure(struct2save)
#         io.save(final_path,
#                 preserve_atom_numbering=True)
        
#     os.remove(out_file)

def download_pdb(pdb_id, uniprot, out_dir):
    out_file = out_dir + f"/{uniprot}_{pdb_id}.cif"
    url = "https://www.ebi.ac.uk/pdbe/entry-files/download/"

    download_link = url + pdb_id.lower() + '_updated.cif'

    # Download
    response = get(download_link, allow_redirects=True)
    open(out_file, "wb").write(response.content)
    
    if response.status_code == 200:
        with open(out_file, 'w') as out_file_write:
            out_file_write.write(response.text)
        print(f"⬇ Descargado {pdb_id}")
    else:
        print(f"❌ Error descargando {pdb_id}")
    
    
    uni_chain = uniprot2pdb_raw[(uniprot2pdb_raw.SP_PRIMARY == uniprot) & (uniprot2pdb_raw.PDB == pdb_id.lower())]['CHAIN'].tolist()
    
    for chain in uni_chain:
        # uni_len = df[(df.SP_PRIMARY == uniprot) & (df.PDB == pdb_id.lower()) & (df.CHAIN == chain)]
        try:
            parser = MMCIFParser(QUIET=True, auth_chains=True)
            struct = parser.get_structure(pdb_id, out_file)
            final_path = out_dir + f"/{uniprot}_{pdb_id}_{chain}.cif"
            struct2save = selected_residues(chain, struct, res2sel='', het_sel = True)
            io = MMCIFIO()
            io.set_structure(struct2save)
            io.save(final_path,
                    preserve_atom_numbering=True)
        except:
            # try:
            #     parser = MMCIFParser(QUIET=True, auth_chains=True)
            #     struct = parser.get_structure(pdb_id, out_file)
            #     final_path = out_dir + f"/{uniprot}_{pdb_id}_{chain}.cif"
            #     struct2save = selected_residues(chain, struct, res2sel='', het_sel = True)
            #     io = MMCIFIO()
            #     io.set_structure(struct2save)
            #     io.save(final_path,
            #             preserve_atom_numbering=True)
            # except:
                print(f'SE LANZO {pdb_id} {uniprot} {chain}', )
                os.remove(out_file)
                return 
            
    os.remove(out_file)
    
# =========================
# CONFIG
# =========================
analysis_folder = '/home/andres/Desktop/foldseekresults3'
all_dirs = Path(analysis_folder)
uniprot2pdb = pd.read_csv('/home/andres/Desktop/human_bac_short_done2/uniprot2pdb.csv')  # o .pkl si usas pickle
uniprot2pdb_raw = pd.read_csv('/home/andres/Desktop/human_bac_short_done2/uniprot2pdb_raw.csv')
OUTPUT_DIR = Path("/home/andres/Desktop/pdbs_selected")
FOLDSEEK_OUT = "/home/andres/Desktop/foldseek_clusters"
TMP_DIR = "/home/andres/Desktop/tmp_foldseek"

df_uniprots = pd.read_csv('/home/andres/Desktop/prioritized_prots_to_cluster.csv', sep = ';')

tar2uniprots = (
    df_uniprots
    .groupby('tar_name')['uniprot_id']
    .apply(list)
    .to_dict()
)

print(tar2uniprots)
df = pd.read_csv("/home/andres/Desktop/human_bac_short_done2/my_tsv.tsv.gz", sep="\t", compression="gzip", skiprows = 1)

all_dirs.mkdir(exist_ok=True, parents=True)

for tar_name, UNIPROT_TO_SELECT in tar2uniprots.items():
    tar_name = tar_name.replace(' ', '_')
    
    mer = 'mono'
    OUTPUT_DIR = Path(analysis_folder+f"/{mer}_{tar_name}/pdbs_selected")
    FOLDSEEK_OUT = analysis_folder+f"/{mer}_{tar_name}/foldseek_clusters"
    TMP_DIR = analysis_folder+f"/{mer}_{tar_name}/tmp_foldseek"

    all_pdbs = dict()
    for u_i in set(UNIPROT_TO_SELECT):
        all_pdbs[u_i] = uniprot2pdb[uniprot2pdb.SP_PRIMARY.apply(lambda x: u_i in x.split('-'))].PDB.tolist()[0].split(' ') 
    
    # all_pdbs_df = pd.concat(all_pdbs, ignore_index=True)
    # all_pdbs = all_pdbs_df.PDB.str.upper().str.split().values.flatten().tolist()
    # all_pdbs = [j for i in all_pdbs for j in i]

    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    for u_i, pdbs in all_pdbs.items():
        for p in pdbs:
            # try:
                download_pdb(p, u_i, str(OUTPUT_DIR))
            # except:
            #     print('FAILED TO DOWNLOAD: ', p, u_i)
            #     continue
    
    cmd = [
        "foldseek",
        "easy-cluster",
        str(OUTPUT_DIR),
        FOLDSEEK_OUT,
        TMP_DIR,
        "-c", "0.8"
    ]
    
    print("\n Ejecutando Foldseek...\n")
    subprocess.run(cmd, check=True)
    
    print("\n✅ Clustering terminado.")
    print(f" Resultados en: {FOLDSEEK_OUT}_clu.tsv")







rep_df = pd.DataFrame(columns = ['uniprots', 'tar_name', 'pdb'])
for folder in os.listdir(analysis_folder):
    for file in os.listdir(analysis_folder + f'/{folder}'):
        if file.endswith('rep_seq.fasta'):
            with open(analysis_folder+f'/{folder}/{file}') as fasta:
                for line in fasta:
                    if line.startswith('>'):
                        pdb = line.split(' ')[0][1:].split('_')
                        print(pdb[0], folder[5:], "_".join(pdb[1:]))
                        chain = "_".join(pdb[1:3])
                        rep_df.loc[len(rep_df)] = pdb[0], folder[5:], chain
    
                                                                             
                        # for tar_name, UNIPROT_TO_SELECT in tar2uniprots.items():
                        #     for u_i in set(UNIPROT_TO_SELECT):
                        #         all_pdbs = uniprot2pdb_raw[uniprot2pdb_raw.SP_PRIMARY.apply(lambda x: u_i in x.split('-'))]
                        #         # all_pdbs_df = pd.concat(all_pdbs, ignore_index=True)
                        #         all_pdbs = all_pdbs[['PDB','CHAIN']].astype(str).agg('_'.join, axis=1).tolist()
                        #         all_pdbs = [i.lower() for i in all_pdbs]
                        #         if pdb.lower() in all_pdbs:
                        #             rep_df.loc[len(rep_df)] = u_i, tar_name, pdb
                                
rep_df.columns =['uniprot_id', 'tar_name', 'pdb']                           
df_uniprots2 = df_uniprots.merge(rep_df, on=['uniprot_id', 'tar_name'], how = 'left')         
df_uniprot_pdb = (
    df_uniprots2
    .groupby(['uniprot_id', 'tar_name'], as_index=False)
    .agg({
        'pdb': lambda s: '-'.join(sorted(s.astype(str).unique()))
    })
)

df_uniprots3 = df_uniprots.merge(df_uniprot_pdb, on=['uniprot_id', 'tar_name'], how = 'left')

n_clus = 0
clus = open('/home/andres/Desktop/foldseekresults3/mono_Tyrosine_decarboxylase/foldseek_clusters_all_seqs.fasta', 'r').read()
import pandas as pd

def parse_cluster_fasta(folder, path):
    rows = []
    cluster_id = -1
    last_header = None
    folder = folder.split('_')[1][:6]
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                continue

            header = line[1:]  # remove '>'

            # If repeated header → new cluster, representative
            if line.startswith('>') and next(f).startswith('>'):
                cluster_id += 1

                # Parse representative
                parts = header.split("_")
                rows.append({
                    "uniprot": parts[0],
                    "pdb": parts[1],
                    "chain": parts[2],
                    "cluster_id": f'{folder}_cluster_{cluster_id}',
                    "is_representative": True
                })

            else:
                # Normal member line
                parts = header.split("_")
                rows.append({
                    "uniprot": parts[0],
                    "pdb": parts[1],
                    "chain": parts[2],
                    "cluster_id": f'{folder}_cluster_{cluster_id}',
                    "is_representative": False
                })

            last_header = header

    return pd.DataFrame(rows)

# ver = parse_cluster_fasta('/home/andres/Desktop/foldseekresults3/mono_Tyrosine_decarboxylase/foldseek_clusters_all_seqs.fasta')

contador = []
for folder in os.listdir(analysis_folder):
    for file in os.listdir(analysis_folder + f'/{folder}'):
        if file.endswith('rep_seq.fasta'):
            contador.append(parse_cluster_fasta(folder, analysis_folder + f'/{folder}/foldseek_clusters_all_seqs.fasta'))

contador = pd.concat(contador, axis = 0)
contador.columns = ['uniprot_id', 'pdb', 'chain', 'cluster_id', 'is_representative']
df2save = df_uniprots.merge(contador, on = 'uniprot_id', how = 'left')

ver = df2save[df2save.is_representative == True]
df2save[['uniprot_id', 'pdb', 'chain', 'cluster_id', 'is_representative', 'tar_name', 'organism', 'organism_type', 'tar_cl',
       'Gene_name', 'natural ligands', 'tier', 'Shekarabi']].to_csv('/home/andres/Desktop/prioritized_prots_to_cluster_pdb_multimer.csv', index = None, sep = ';')

contador = dict()
for folder in os.listdir(analysis_folder):
    contador[folder] = []
    for file in os.listdir(analysis_folder + f'/{folder}'):
        if file.endswith('rep_seq.fasta'):
            with open(analysis_folder+f'/{folder}/{file}') as fasta:
                for line in fasta:
                    if line.startswith('>'):
                        contador[folder] = contador.get(folder, []) + [line.split(' ')[0]]
                        

df_uniprots3[['uniprot_id', 'pdb', 'organism', 'organism_type', 'tar_cl', 'tar_name',
       'Gene_name', 'natural ligands', 'tier', 'Shekarabi']].to_csv('/home/andres/Desktop/prioritized_prots_to_cluster_pdb.csv', index = None, sep = ';')

def download_pdb(pdb_id, uniprot, out_dir):
    out_file = out_dir + f"/{uniprot}_{pdb_id}.cif"
    url = "https://www.ebi.ac.uk/pdbe/entry-files/download/"

    download_link = url + pdb_id.lower() + '_updated.cif'

    # Download
    response = get(download_link, allow_redirects=True)
    open(out_file, "wb").write(response.content)
    
    if response.status_code == 200:
        with open(out_file, 'w') as out_file_write:
            out_file_write.write(response.text)
        print(f"⬇ Descargado {pdb_id}")
    else:
        print(f"❌ Error descargando {pdb_id}")
    
    
# =========================
# CONFIG
# =========================
analysis_folder = '/home/andres/Desktop/foldseekresults4'
all_dirs = Path(analysis_folder)
uniprot2pdb = pd.read_csv('/home/andres/Desktop/human_bac_short_done2/uniprot2pdb.csv')  # o .pkl si usas pickle
uniprot2pdb_raw = pd.read_csv('/home/andres/Desktop/human_bac_short_done2/uniprot2pdb_raw.csv')
OUTPUT_DIR = Path("/home/andres/Desktop/pdbs_selected")
FOLDSEEK_OUT = "/home/andres/Desktop/foldseek_clusters"
TMP_DIR = "/home/andres/Desktop/tmp_foldseek"

df_uniprots = pd.read_csv('/home/andres/Desktop/prioritized_prots_to_cluster.csv', sep = ';')

uniprot2pdb[uniprot2pdb.SP_PRIMARY == 'Q97S36']

tar2uniprots = (
    df_uniprots
    .groupby('tar_name')['uniprot_id']
    .apply(list)
    .to_dict()
)

print(tar2uniprots)
df = pd.read_csv("/home/andres/Desktop/human_bac_short_done2/my_tsv.tsv.gz", sep="\t", compression="gzip", skiprows = 1)

all_dirs.mkdir(exist_ok=True, parents=True)

for tar_name, UNIPROT_TO_SELECT in tar2uniprots.items():
    tar_name = tar_name.replace(' ', '_')
    
    mer = 'mono'
    OUTPUT_DIR = Path(analysis_folder+f"/{mer}_{tar_name}/pdbs_selected")
    FOLDSEEK_OUT = analysis_folder+f"/{mer}_{tar_name}/foldseek_clusters"
    TMP_DIR = analysis_folder+f"/{mer}_{tar_name}/tmp_foldseek"

    all_pdbs = dict()
    for u_i in set(UNIPROT_TO_SELECT):
        all_pdbs[u_i] = uniprot2pdb[uniprot2pdb.SP_PRIMARY.apply(lambda x: u_i in x.split('-'))].PDB.tolist()[0].split(' ') 
    
    # all_pdbs_df = pd.concat(all_pdbs, ignore_index=True)
    # all_pdbs = all_pdbs_df.PDB.str.upper().str.split().values.flatten().tolist()
    # all_pdbs = [j for i in all_pdbs for j in i]

    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    for u_i, pdbs in all_pdbs.items():
        for p in pdbs:
            # try:
                download_pdb(p, u_i, str(OUTPUT_DIR))
            # except:
            #     print('FAILED TO DOWNLOAD: ', p, u_i)
            #     continue
    
    cmd = [
        "foldseek",
        "easy-multimercluster",
        str(OUTPUT_DIR),
        FOLDSEEK_OUT,
        TMP_DIR,
        "-c", "0.8"
    ]
    
    print("\n Ejecutando Foldseek...\n")
    subprocess.run(cmd, check=True)
    
    print("\n✅ Clustering terminado.")
    print(f" Resultados en: {FOLDSEEK_OUT}_clu.tsv")

def parse_cluster_fasta(folder, path):
    rows = []
    cluster_id = -1
    last_header = None
    folder = folder.split('_')[1][:6]
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                continue

            header = line[1:]  # remove '>'

            # If repeated header → new cluster, representative
            if line.startswith('>') and next(f).startswith('>'):
                cluster_id += 1

                # Parse representative
                parts = header.split("_")
                rows.append({
                    "uniprot": parts[0],
                    "pdb": parts[1],
                    "chain": parts[2],
                    "cluster_id": f'{folder}_cluster_{cluster_id}',
                    "is_representative": True
                })

            else:
                # Normal member line
                parts = header.split("_")
                rows.append({
                    "uniprot": parts[0],
                    "pdb": parts[1],
                    "chain": parts[2],
                    "cluster_id": f'{folder}_cluster_{cluster_id}',
                    "is_representative": False
                })

            last_header = header

    return pd.DataFrame(rows)

# ver = parse_cluster_fasta('/home/andres/Desktop/foldseekresults3/mono_Tyrosine_decarboxylase/foldseek_clusters_all_seqs.fasta')

contador = []
for folder in os.listdir(analysis_folder):
    for file in os.listdir(analysis_folder + f'/{folder}'):
        if file.endswith('rep_seq.fasta'):
            contador.append(parse_cluster_fasta(folder, analysis_folder + f'/{folder}/foldseek_clusters_all_seqs.fasta'))

contador = pd.concat(contador, axis = 0)
contador.columns = ['uniprot_id', 'pdb', 'chain', 'cluster_id', 'is_representative']
df2save = df_uniprots.merge(contador, on = 'uniprot_id', how = 'left')

ver = df2save[df2save.is_representative == True]
df2save[['uniprot_id', 'pdb', 'chain', 'cluster_id', 'is_representative', 'tar_name', 'organism', 'organism_type', 'tar_cl',
       'Gene_name', 'natural ligands', 'tier', 'Shekarabi']].to_csv('/home/andres/Desktop/prioritized_prots_to_cluster_pdb_multimer.csv', index = None, sep = ';')
