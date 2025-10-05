# -*- coding: utf-8 -*-
"""
Created on 2021

@author: Andres Sanchez
"""

import time

import pandas as pd
#import random

import numpy as np

from rdkit import Chem

from chembl_structure_pipeline import standardizer as sdz
from chembl_structure_pipeline import checker

from chembl_webresource_client.new_client import new_client
molecule = new_client.molecule
activities = new_client.activity
molecule.set_format('json')
activities.set_format('json')


start = time.time()
def molblock_estandarizadotosmile(molblock_estandarizado):
    try:
        return Chem.MolToSmiles(Chem.MolFromMolBlock(molblock_estandarizado))
    except:
        return 'NaN'
    
# Generamos los molblock legibles
def moltomolblock(x):
    try:
        return Chem.MolToMolBlock(x)
    except:

        return 'NaN'

# Estandarizamos los molblock
def molblockestandarizer(molblock):
    try:
        return sdz.standardize_molblock(molblock)
    except:
        return 'NaN'   
    
# Buscamos flags 
def flags(molblock_estandarizado):
    try:
        return sdz.exclude_flag(molblock_estandarizado)
    except:
        return 'NaN'
    
# Generamos los smiles de los molblock para ver si son compuestos simples o comlejos
def moltosmile(molblock):
    try:
        return Chem.MolToSmiles(molblock)
    except:
        return 'NaN'
    
# Determinamos qué compuestos son simples y cúales complejos
def busquedaparent(smile):
    if "." in smile:
        return "Compuesto_complejo"
    return "Compuesto_simple"

# Buscamos problemas en la anotación de las moleculas        
def checkeoissues(molblock_estandarizado):
    issues = checker.check_molblock(molblock_estandarizado)
    if len(issues) > 0:
        return issues[0][0]
    if len(issues) == 0:
        return 0
# Hacemos la funcion para la búsqueda por el método FLEX
def buscadorflexCHEMBL(smile, contador_funcion=[1], 
                     contador_ausentes=[0], contador_encontrados=[0]):    
    print("Busqueda FLEX. He corrido la funcion", contador_funcion, "veces.", "He encontrado el ID del", round((contador_encontrados[0]/contador_funcion[0])*100,2),"%")  
    contador_funcion[0] += 1
    if smile == 'NaN' or smile == None:
        contador_ausentes[0] += 1
        return 'NaN'
    try:
        res_list = ""
        res = molecule.filter(molecule_structures__canonical_smiles__flexmatch=smile)
        res_list = ["".join(x['molecule_chembl_id']) for x in res]
        if len(res) == 0:
            contador_ausentes[0] += 1
            return 'Not_found'
        contador_encontrados[0] += 1
        return res_list
    except:
        contador_ausentes[0] += 1
        res = "Not_found"
        return res

#Buscador actividades
def buscadoractividadesCHEMBL(chemblid, contador_funcion=[1], 
                     contador_ausentes=[0], contador_encontrados=[0]):  
    print("Busqueda CHEMBL ids. He corrido la funcion", contador_funcion, "veces.", "He encontrado el ID del", round((contador_encontrados[0]/contador_funcion[0])*100,2),"%")  
    try:        
        res = activities.filter(molecule_chembl_id__in=chemblid).only(['molecule_chembl_id', 
                                                                       'target_chembl_id', 
                                                                       'molecule_pref_name',
                                                                       'pchembl_value', 
                                                                       'target_organism'])
        res_chemblid = [x['molecule_chembl_id'] for x in res]     
        res_target = [x['target_chembl_id'] for x in res]
        res_name = [x['molecule_pref_name'] for x in res]
        res_pchembl = [x['pchembl_value'] for x in res]
        res_org = [x['target_organism'] for x in res]
        return pd.Series([res_chemblid, res_target, res_name, res_pchembl, res_org])
    except:
        contador_ausentes[0] += 1
        res = "Not_found"
        return res
    
# Hacemos la búsqueda por GET
def buscadorenCHEMBL(inchikey, contador_funcion=[1], 
                     contador_ausentes=[0], contador_encontrados=[0]):    
    print("Busqueda GET. He corrido la funcion", contador_funcion, "veces.", "He encontrado el ID del", round((contador_encontrados[0]/contador_funcion[0])*100,2),"%")  
    contador_funcion[0] += 1
    if inchikey == 'NaN' or inchikey == None or inchikey == 'Compuesto_simple' or inchikey == 'No_apto':
        contador_ausentes[0] += 1
        return 'NaN'
    try: 
        res = molecule.get(inchikey)
        contador_encontrados[0] += 1
        return res['molecule_chembl_id']
    except:
        contador_ausentes[0] += 1
        res = "Molecule not found"
        return res
   
# Generamos el parent de aquellos compuestos complejos
def getparent(molblock_estandarizado):
    try: 
        parent, _ = sdz.get_parent_molblock(molblock_estandarizado)
        return parent
    except:
        parent = "No localizable"
        return parent

# Generamos los inchikey para el parent
def inchikeygenerator(molblock_estandarizado):
    if molblock_estandarizado == 'Compuesto_simple':
        return 'NaN'
    try:
        inchikey = Chem.InchiToInchiKey(Chem.MolBlockToInchi(molblock_estandarizado))
        return inchikey
    except:
        return 'NaN'
    return 'NaN'

def molfrommolblock(molblock):
    try:
        return Chem.MolFromMolBlock(molblock)
    except:
        return 'Compuesto Simple'

# Guardamos los foodb IDs
filename = 'C:/Users/Andres Sanchez/OneDrive - FUNDACION IMDEA-ALIMENTACION/Escritorio/foodb_sdf.sdf'     
foodb_ids = []
with open(filename, encoding = 'utf-8') as f:
    for line in f:
        if line.startswith('> <DATABASE_ID>'):
            nextLine = next(f)
            foodb_ids.append(nextLine)

#Leemos el archivo sdf directamente
pattern_start = '$'
pattern_end = ">"
molblock = """"""
lista_molblocks = []
flist = open(filename, encoding = 'utf-8').readlines()
parsing = False
for line in flist:
    if line.startswith(pattern_start):
        parsing = True
    if line.startswith(pattern_end):
        parsing = False
        if len(molblock) > 0: 
            if molblock[14] == "M":
                lista_molblocks.append(molblock[5:])
                molblock = """"""
                continue
            if molblock[9] == "M":
                lista_molblocks.append(molblock[5:])
                molblock = """"""
                continue
            if molblock[5] == " ":
                lista_molblocks.append(molblock[4:])
                molblock = """"""
                continue
            lista_molblocks.append(molblock[5:])
            molblock = """"""                
    if parsing:
        molblock += line

#random.seed(1234)
prueba = pd.DataFrame(foodb_ids, columns=['foodb_ID'])
prueba['Molblock'] = lista_molblocks
prueba['Mol'] = prueba['Molblock'].apply(lambda x: Chem.MolFromMolBlock(x))
#muestra = 5000
#random_rows = random.sample(range(len(lista_molblocks)), muestra)
#prueba = df.iloc[random_rows]

# COMENZAMOS PIPELINE
prueba['Exclude_flag'] = prueba['Molblock'].apply(lambda x: flags(x))
prueba['Issues'] = prueba['Molblock'].apply(lambda x: checkeoissues(x))
prueba['Molblock_estandarizado'] = np.where((prueba['Issues']>5.9) | (prueba['Exclude_flag'] == True), 
                                            'No_apto',
                                            prueba['Molblock'].apply(lambda x: molblockestandarizer(x)))
prueba['Mol_estandarizado'] = prueba['Molblock_estandarizado'].apply(lambda x: molfrommolblock(x))
prueba['Inchikeys'] = prueba['Molblock_estandarizado'].apply(lambda x: inchikeygenerator(x))
# Hasta aquí
prueba['Smiles'] = prueba['Mol'].apply(lambda x: moltosmile(x))
prueba['Smiles_canonizados'] = prueba['Mol_estandarizado'].apply(lambda x: moltosmile(x))
prueba['Tipo_compuesto'] = prueba['Smiles_canonizados'].apply(lambda x: busquedaparent(x))
prueba['ID_por_FLEX'] = prueba['Smiles_canonizados'].apply(lambda x: buscadorflexCHEMBL(x))
prueba['ID_por_GET'] = prueba.loc[prueba['ID_por_FLEX'] != 'Not_found', 'Inchikeys'].apply(lambda x: buscadorenCHEMBL(x))
prueba['Parent_molblock'] = np.where(prueba['Tipo_compuesto'] == 'Compuesto_complejo',
                            prueba['Molblock_estandarizado'].apply(lambda x: getparent(x)),
                            'Compuesto_simple')
prueba['Parent_mol'] = prueba['Parent_molblock'].apply(lambda x: Chem.MolFromMolBlock(x))
prueba['Parent_smiles'] = prueba['Parent_molblock'].apply(lambda x: molblock_estandarizadotosmile(x))
prueba['Parent_ID_por_FLEX'] = prueba['Parent_smiles'].apply(lambda x: buscadorflexCHEMBL(x))
prueba['Parent_ID_por_FLEX'].fillna("Not_found",inplace=True)
prueba['Parent_inchikeys'] = prueba['Parent_molblock'].apply(lambda x: inchikeygenerator(x))
prueba['Parent_ID_por_GET'] = prueba.loc[prueba['Parent_ID_por_FLEX'] != 'Not_found', 'Parent_inchikeys'].apply(lambda x: buscadorenCHEMBL(x))

prueba.to_csv('C:/Users/Andres Sanchez/OneDrive - FUNDACION IMDEA-ALIMENTACION/Escritorio/TODOS(15-05-2021).csv', index = False)


# Observamos la distribución de los issues 
prueba["Issues"].plot(kind="hist", bins = 7, figsize=(4,6))

# Fraccion de flags positivos y negativos
contador_true = prueba.loc[prueba["Exclude_flag"] == True].shape[0]
contador_false = prueba.loc[prueba["Exclude_flag"] == False].shape[0]
total = contador_false + contador_true
fraccion_exclude = "El", round((contador_false/total)*100,2), "% no tienen ninguna flag", round((contador_true/total)*100,2), "% tienen flag"

# Fraccion de compuestos simples y complejos
contador_simple = prueba.loc[prueba["Tipo_compuesto"] == "Compuesto_simple"].shape[0]
contador_complejo = prueba.loc[prueba["Tipo_compuesto"] == "Compuesto_complejo"].shape[0]
total = contador_simple + contador_complejo
fraccion_tipo_compuesto = "El", round((contador_simple/total)*100,2), "% de los compuestos son simples y el", round((contador_complejo/total)*100,2), "% son complejos"

end = time.time()
print("He tardado", ((end - start)/60), "min")  
# ~33 horas para los 70K







# ESTA PARTE ESTA BIEN HECHA EN EL OTRO SCRIPT DE BUSQUEDA ACTIVIDADES POR CHUNKS
# NO LO BORRO POR SI ACASO 

prueba = pd.read_csv('C:/Users/Andres Sanchez/OneDrive - FUNDACION IMDEA-ALIMENTACION/Escritorio/TODOS(15-05-2021).csv')

ids = prueba[['foodb_ID', 'ID_por_GET', 'Parent_ID_por_GET']]
ids['ID_por_GET'].fillna("Molecule not found",inplace=True)
encontrados = ids.loc[ids['ID_por_GET'] != 'Molecule not found']
#Cuantas moleculas simples se han encontrado
encontrados = encontrados.loc[encontrados['ID_por_GET'] != 'NaN']

#Cuantas moleculas complejas se han encontrado 
complejos = ids.loc[(ids['Parent_ID_por_GET'] != 'NaN') & (ids['Parent_ID_por_GET'] != 'Molecule not found')]
complejos['Parent_ID_por_GET'].fillna("Molecule not found",inplace=True)
complejos = complejos.loc[complejos['Parent_ID_por_GET'] != 'Molecule not found']

#Buscamos las actividades de los compuestos
#encontrados = encontrados[:20]

encontrados[['ChEMBL_ID',
             'Target_ChEMBL_ID', 
             'ChEMBL_name', #name no se corresponde con el target si no con el name del ChEMBL ID
             'Target_pChEMBL', 
             'Target_organism']] = encontrados['ID_por_GET'].apply(lambda x: buscadoractividadesCHEMBL(x))

#Generamos el vector con los foodb ids en funcion del numero de hits
longitud = [len(x) for x in encontrados['ChEMBL_ID']]
vector_foodb_IDs = [[x] * l for x, l in zip(encontrados['foodb_ID'], longitud)]
vector_foodb_IDs = [compuesto for lista in vector_foodb_IDs for compuesto in lista]
             
#Deshacemos las listas generadas en la búsqueda de actividades
lista_chemblid = [compuesto for lista in encontrados['ChEMBL_ID'] for compuesto in lista]
lista_target = [compuesto for lista in encontrados['Target_ChEMBL_ID'] for compuesto in lista]
lista_pchembl = [compuesto for lista in encontrados['Target_pChEMBL'] for compuesto in lista]
lista_name = [compuesto for lista in encontrados['ChEMBL_name'] for compuesto in lista]
lista_organism = [compuesto for lista in encontrados['Target_organism'] for compuesto in lista]

#Construimos finalmente el dataframe
actividades = pd.DataFrame(vector_foodb_IDs, columns=['Foodb_ID'])
actividades['ChEMBL_ID'] = lista_chemblid
actividades['ChEMBL_name'] =  lista_name
actividades['Target_ChEMBL_ID'] = lista_target
actividades['Target_pChEMBL'] =  lista_pchembl
actividades['Target_organism'] = lista_organism

#complejos = complejos[:20]
complejos[['Parent_ChEMBL_ID',
             'Parent_Target_ChEMBL_ID', 
             'Parent_ChEMBL_name', #name no se corresponde con el target si no con el name del ChEMBL ID
             'Parent_Target_pChEMBL', 
             'Parent_Target_organism']] = complejos['Parent_ID_por_GET'].apply(lambda x: buscadoractividadesCHEMBL(x))

#Generamos el vector con los foodb ids en funcion del numero de hits
longitud = [len(x) for x in complejos['Parent_ChEMBL_ID']]
vector_foodb_IDs = [[x] * l for x, l in zip(complejos['foodb_ID'], longitud)]
vector_foodb_IDs = [compuesto for lista in vector_foodb_IDs for compuesto in lista]
             
#Deshacemos las listas generadas en la búsqueda de actividades
lista_chemblid = [compuesto for lista in complejos['Parent_ChEMBL_ID'] for compuesto in lista]
lista_target = [compuesto for lista in complejos['Parent_Target_ChEMBL_ID'] for compuesto in lista]
lista_pchembl = [compuesto for lista in complejos['Parent_Target_pChEMBL'] for compuesto in lista]
lista_name = [compuesto for lista in complejos['Parent_ChEMBL_name'] for compuesto in lista]
lista_organism = [compuesto for lista in complejos['Parent_Target_organism'] for compuesto in lista]

actividades_p = pd.DataFrame(vector_foodb_IDs, columns=['Foodb_ID'])
actividades_p['ChEMBL_ID'] = lista_chemblid
actividades_p['ChEMBL_name'] =  lista_name
actividades_p['Target_ChEMBL_ID'] = lista_target
actividades_p['Target_pChEMBL'] =  lista_pchembl
actividades_p['Target_organism'] = lista_organism

#actividades.to_csv('C:/Users/Andres Sanchez/OneDrive - FUNDACION IMDEA-ALIMENTACION/Escritorio/pruebasactividades(15-05-2021).csv', index = False)
#actividades_p.to_csv('C:/Users/Andres Sanchez/OneDrive - FUNDACION IMDEA-ALIMENTACION/Escritorio/pruebasactividades_p(15-05-2021).csv', index = False)

