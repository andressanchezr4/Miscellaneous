# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 11:32:30 2024

@author: andres.sanchez
"""

# AttributeDetailsPDB: https://www.rcsb.org/docs/search-and-browse/advanced-search/attribute-details
# More details on attributes: https://search.rcsb.org/structure-search-attributes.html
# RCSB PDB api documentation: https://www.rcsb.org/docs/search-and-browse/advanced-search/attribute-details
# Usar también la busqueda avanzada del rcsb pdb (conretamente en el BOTON SearchAPI, en la parte superior derecha)
import requests
import rcsbsearchapi.search as rcsbsearch
from rcsbsearchapi.search import AttributeQuery
from rcsbsearchapi.search import ChemSimilarityQuery
import pandas as pd
from rcsbsearchapi.search import AttributeQuery
import time

def FromInchiToPDBids_noncovalent(inchi):
    # Takes a little bit longer than the request done
    # through requests library
                   
    # WITH THIS CALL WE RETRIEVE THE LIGANDS WITH THE EXACT 
    # INCHI THAT ARE NOT COVALENTLY ATTACHED
    
    # start = time.time()
    q = ChemSimilarityQuery(value=inchi,
                             query_type="descriptor",
                             descriptor_type="InChI",
                             match_type="graph-exact")
    
    q2 =  AttributeQuery("rcsb_nonpolymer_instance_annotation.type", "exact_match",
                         "HAS_NO_COVALENT_LINKAGE")
    q_final = q & q2
    
    try:
        n_pdbs = list(q_final(return_type='entry'))
        ligand_id = list(q_final(return_type='mol_definition'))
    
    except:
        print('SALTO EL LIMITE')
        time.sleep(5)
        n_pdbs = list(q_final(return_type='entry'))
        ligand_id = list(q_final(return_type='mol_definition'))
    # finish = time.time() - start
    return n_pdbs, ligand_id

df = pd.read_csv('C:/Users/andres.sanchez/Desktop/comp_sets_Review 1.csv', sep = ';')
# df[['pdb_id', 'ligand_id']] = df['inchi'].apply(lambda inchi: pd.Series(FromInchiToPDBids(inchi)))

n_pdbs = []
n_ligands = []
for n, inch in enumerate(df['inchi']): #1984
    n_p, l_i = FromInchiToPDBids_noncovalent(inch)
    print(n_p, l_i)
    n_pdbs.append(n_p)
    n_ligands.append(l_i)
    print(f'{n}/{len(df)}')
    

################
### REQUESTS ###
################

def noAPIFromInchiToPDBids_noncovalent(inchi):
    
    # WITH THIS CALL WE RETRIEVE THE LIGANDS WITH THE EXACT 
    # INCHI THAT ARE NOT COVALENTLY ATTACHED
    
    results_api_call = []
    
    for r_t in ['mol_definition', 'entry']:
        
        api_params = {
          "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
              {
            "type": "terminal",
            "service": "chemical",
            "parameters": {
              "value": inchi,
              "type": "descriptor",
              "descriptor_type": "InChI",
              "match_type": "graph-exact"
            }
              },
              {
                "type": "terminal",
                "service": "text",
                "parameters": {
                  "attribute": "rcsb_nonpolymer_instance_annotation.type",
                  "operator": "exact_match",
                  "value": "HAS_NO_COVALENT_LINKAGE"
                }
              }
            ]
          },
          "return_type": r_t,
          "request_options": {
            "results_verbosity": "compact"
          }
        }
        
        url = 'https://search.rcsb.org/rcsbsearch/v2/query'
        
        try:
            response = requests.post(url, json=api_params)
        except:
            time.sleep(5)
            response = requests.post(url, json=api_params)
        
        if response.status_code == 200:
            results_api_call.append(response.json()['result_set'])
        else:
            results_api_call.append([None])
    
    return results_api_call[0], results_api_call[1]


def FromInchiToPDBids_noncovalent(inchi):
    # Takes a little bit longer than the request done
    # through requests library
    
    # WITH THIS CALL WE RETRIEVE ALL THE LIGANDS WITH THE EXACT 
    # INCHI THAT ARE BOTH COVALENTLY AND NON COVALENTLY ATTACHED
    
    # start = time.time()
    q_final = ChemSimilarityQuery(value=inchi,
                             query_type="descriptor",
                             descriptor_type="InChI",
                             match_type="graph-exact")
    
    
    try:
        n_pdbs = list(q_final(return_type='entry'))
        ligand_id = list(q_final(return_type='mol_definition'))
    
    except:
        print('SALTO EL LIMITE')
        time.sleep(5)
        n_pdbs = list(q_final(return_type='entry'))
        ligand_id = list(q_final(return_type='mol_definition'))
    # finish = time.time() - start
    return n_pdbs, ligand_id

df = pd.read_csv('./comp_sets_Review 1.csv', sep = ';')
# df[['pdb_id', 'ligand_id']] = df['inchi'].apply(lambda inchi: pd.Series(FromInchiToPDBids(inchi)))

n_pdbs = []
n_ligands = []
for n, inch in enumerate(df['inchi']): #1984
    n_p, l_i = FromInchiToPDBids_noncovalent(inch)
    print(n_p, l_i)
    n_pdbs.append(n_p)
    n_ligands.append(l_i)
    print(f'{n}/{len(df)}')
    
df2 = pd.read_csv('./inchi_to_ligand_id.csv')  
df_final = df.merge(df2, on = 'inchi', how = 'inner')
df_txt = df_final[['hmdb_id', 'ligand_id', 'pdb_id', 'inchi']]

df['pdb_id'] = n_pdbs
df['ligand_id'] = ligand_id

df_api = df[['hmdb_id', 'pdb_id', 'ligand_id', 'inchi']]
df_api2 = df_api.copy()
df_api2['ligand_id'] = df_api2['ligand_id'].apply(lambda x: [i for i in x if len(i) == 3] if len(x) > 0 else x)

df['pdb_id'] = n_ligands
df['ligand_id'] = n_pdbs
df_final = df[['hmdb_id', 'ligand_id', 'pdb_id', 'inchi']].copy()

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
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}

solvents = 'EOH, MOH, GOL, 1GP, 2DP, 3PH, 6A8, DGA, DGD, DGG, DR9, DVG, G3P, HGP, \
            HGX, IGP, INB, L1P, L2P, L3P, L4P, LHG, LI1, LIO, LPC, PGM, SGL, SGM, SQD, TGL, \
            12P, 15P, 1PE, 2PE, CE9, CP4, DIO, P4C, P6G, PG4, PGE, VNY, DMS, EDO, PEG, TRS, \
            IPA, TBU, ACT, EEE, ACY, BME, MBN, NAG, SIA, FUK, IVA, STA, BMA, SO4, MAN, GAL, \
            DOD, SO3, IOD, PO4, PO3, TLA, PSA, MES, PG4, FUC, SOG, GLC, SF4'
solvents = [i.strip() for i in solvents.split(',')]

df_final_onlyTP = df_final[~df_final['ligand_id'].apply(lambda x: any(i in x for i in list(substitutions.keys()) + solvents))]
df_final_onlyTP = df_final_onlyTP[~((df_final_onlyTP['ligand_id'].apply(lambda x: x[0] is None)) & (df_final_onlyTP['pdb_id'].apply(lambda x: x[0] is None)))]

