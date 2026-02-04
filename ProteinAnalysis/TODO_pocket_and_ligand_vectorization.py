#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 11:45:10 2026

@author: andres
"""

import os
import pandas as pd
from Bio.PDB import PDBParser
import numpy as np
import numpy as np
from scipy.stats import skew 
import random
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from torch.utils.data import DataLoader
import torch
from torch_geometric.data import Batch
from rdkit import Chem
import Smiles2Graphs as mg
from rdkit import Chem
import torch
from rdkit.Chem import AllChem
from torch_geometric.data import Data
torch.manual_seed(42)
random.seed(42)
    
# def count_atom_classes(df, mol_column='mol'):
#     """
#     Count atom classes in RDKit Mol objects in a DataFrame.

#     Args:
#         df: pandas DataFrame containing RDKit Mol objects
#         mol_column: name of the column containing the Mol objects

#     Returns:
#         List of dictionaries, one per molecule, with counts per atomic number
#         Example: [{'C': 6, 'O': 1, 'N': 2}, {...}, ...]
#     """
#     atom_counts_list = []

#     for mol in df[mol_column]:
#         counts = Counter()
#         for atom in mol.GetAtoms():
#             # Use symbol instead of atomic number for readability
#             counts[atom.GetSymbol()] += 1
#         atom_counts_list.append(dict(counts))

#     return atom_counts_list
# ver = count_atom_classes(df_all_mol)
# total_counts = Counter()
# for d in ver:
#     total_counts.update(d)

def mol_from_sdf(sdf_path):
    """
    Lee la primera molécula de un archivo SDF y devuelve
    el número de enlaces rotables según la definición de RDKit.
    """
    # Supplier para leer el SDF
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    
    # Extraer primera molécula válida
    mol = next((m for m in suppl if m is not None), None)
    return mol
    
def compute_usr_descriptor(coords):
    """
    Compute a 12-element USR descriptor from N×3 coordinates.

    Args:
        coords: NumPy array of shape (N, 3)
    Returns:
        1D NumPy array of length 12
    """

    # 1) Centroid of all atoms
    ctd = coords.mean(axis=0)

    # 2) Closest atom to centroid (cst)
    dists_to_ctd = np.linalg.norm(coords - ctd, axis=1)
    cst = coords[dists_to_ctd.argmin()]

    # 3) Farthest atom to centroid (fct)
    fct = coords[dists_to_ctd.argmax()]

    # 4) Atom farthest from fct (ftf)
    dists_to_fct = np.linalg.norm(coords - fct, axis=1)
    ftf = coords[dists_to_fct.argmax()]

    # Reference points list
    refs = [ctd, cst, fct, ftf]

    usr_vector = []

    for ref in refs:
        distances = np.linalg.norm(coords - ref, axis=1)
        usr_vector.append(distances.mean())
        usr_vector.append(distances.var())
        usr_vector.append(skew(distances))

    return np.array(usr_vector)

pathpockets = '/home/andres/Desktop/uni_prio/p2rank_batches/'
all_configs = [i for i in os.listdir(pathpockets) if i.endswith('.txt') and 'P2RANK' in i]

df_configs = pd.DataFrame(columns = ['uniprot', 'pdb', 'pdb_file', 'pocket_center', 'pocket_size'])
for config in all_configs:
    split = config.split('_')
    uniprot = split[0]
    pdb_file = config.split('.')[0] + '.pdb'
    pdb = "_".join(split[:5])
    c = open(pathpockets+config).read().split()
    size_pocket = float(c[-1]), float(c[-4]), float(c[-7]) 
    center_pocket = float(c[2]), float(c[5]), float(c[8])
    df_configs.loc[len(df_configs)] = uniprot, pdb, pdb_file, center_pocket, size_pocket
    
path2uniprots= '/home/andres/Desktop/uni_prio/uniprot/'
usr = []
for _, row in df_configs.iterrows():
    cx, cy, cz = row['pocket_center']    # center of the box
    sx, sy, sz = row['pocket_size']       # total length in each dimension
    
    xmin = cx - sx/2
    xmax = cx + sx/2
    ymin = cy - sy/2
    ymax = cy + sy/2
    zmin = cz - sz/2
    zmax = cz + sz/2
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", f"{path2uniprots}{row['uniprot']}/selected_cif_chains/{row['pdb_file']}")
    
    coords = []
    for atom in structure.get_atoms():
        x, y, z = atom.coord
        if (xmin <= x <= xmax and
            ymin <= y <= ymax and
            zmin <= z <= zmax):
            coords.append([x, y, z])
    
    coords = np.array(coords)
    usr.append(compute_usr_descriptor(coords))

df_configs['usr'] = usr
df_configs['Protein'] = df_configs['pdb_file'].apply(lambda x: "_".join(x.split('.')[0].split('_')[1:8]))

df_all = pd.read_csv('/home/andres/Desktop/protligeff.csv')
path2ligands = '/home/andres/Desktop/PRUEBAS_PROTEINSUPERIMP/analisis_box_size_fixed/raw_ligands2/'
all_mol = [AllChem.ComputeGasteigerCharges(mol_from_sdf(path2ligands+p)) for p in os.listdir(path2ligands)]

df_all_mol = pd.DataFrame([i.split('.')[0] for i in os.listdir(path2ligands)], columns = ['Ligand'])
df_all_mol['mol'] = all_mol
df_all = df_all.merge(df_all_mol, on = 'Ligand', how = 'left')
df_all_final = df_all.merge(df_configs[['Protein', 'usr']], on = 'Protein', how = 'left').dropna()
# unir al df_all la columna de cada ligando con su correspondiente mol
# unir al df_all la columna de cada proteina con su correspondiente prot
# splitearlos con los grupos train/test que ya tengo
# hacer los loader


test_prot = pd.read_csv('/home/andres/Desktop/prot2test.csv')['Protein'].tolist()
test_lig = pd.read_csv('/home/andres/Desktop/ligands2test2.csv')['Ligand'].tolist()
df_train = df_all_final[~(df_all_final.Protein.isin(test_prot) | df_all_final.Ligand.isin(test_lig))]
df_test = df_all_final[(df_all_final.Protein.isin(test_prot) & df_all_final.Ligand.isin(test_lig))]

def has_bromine(mol):
    """
    Returns True if the molecule contains at least one bromine atom (Br), False otherwise.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Br':
            return True
    return False

# df_train['is_br'] = df_train.mol.apply(has_bromine)
# df_test['is_br'] = df_test.mol.apply(has_bromine)

# ver = count_atom_classes(df_train)
# total_counts = Counter()
# for d in ver:
#     total_counts.update(d)

df_br = df_test[df_test.mol.apply(has_bromine)]
df_test = df_test[~df_test.mol.apply(has_bromine)]
df_train = pd.concat([df_train, df_br], axis=0, ignore_index=True)

# df_train.to_csv('/home/andres/train_df_timedock.csv', index = False)
# df_test.to_csv('/home/andres/test_df_timedock.csv', index = False)
# 0. Imports
# ===============================
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import random_split
from rdkit import Chem
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.nn import GINEConv
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# ===============================
# 1. RDKit Mol → Molecular Graph
# ===============================
def mol_to_graph(mol):
    node_feats = []
    atom_types=['C','N','O','S','F','Cl','Br','I','P']
    for atom in mol.GetAtoms():
        hyb = atom.GetHybridization()
        # q = 0.0 if np.isnan(float(atom.GetProp('_GasteigerCharge'))) else float(atom.GetProp('_GasteigerCharge'))
        symbol = atom.GetSymbol()
        ohe = [int(symbol == t) for t in atom_types]  # thi
        node_feats.append([
            atom.GetAtomicNum(),
            atom.GetTotalDegree(),
            atom.GetTotalNumHs(),
            int(atom.GetIsAromatic()),
            int(atom.IsInRing()),
            int(hyb == Chem.rdchem.HybridizationType.SP),
            int(hyb == Chem.rdchem.HybridizationType.SP2),
            int(hyb == Chem.rdchem.HybridizationType.SP3),
            int(hyb == Chem.rdchem.HybridizationType.SP3 and not atom.IsInRing()),
            # q,
            *ohe
            
        ])
    x = torch.tensor(node_feats, dtype=torch.float)

    edge_index = []
    edge_attr = []

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        bt = bond.GetBondType()
        bond_feat = [
            int(bt == Chem.rdchem.BondType.SINGLE),
            int(bt == Chem.rdchem.BondType.DOUBLE),
            int(bt == Chem.rdchem.BondType.TRIPLE),
            int(bt == Chem.rdchem.BondType.AROMATIC),
            int(bond.IsInRing()),
            int(bond.GetIsConjugated()),
            int(bt == Chem.rdchem.BondType.SINGLE and not bond.IsInRing())
        ]

        edge_index += [[i, j], [j, i]]
        edge_attr += [bond_feat, bond_feat]

    edge_index = torch.tensor(edge_index, dtype=torch.long).t()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    return x, edge_index, edge_attr

# ===============================
# 2. Dataset (Ligando + USR + y)
# ===============================
class DockingDataset(Dataset):
    def __init__(self, mols, usr_vectors, y_values):
        super().__init__()
        self.mols = mols
        self.usr = usr_vectors
        self.y = y_values

    def len(self):
        return len(self.mols)

    def get(self, idx):
        x, edge_index, edge_attr = mol_to_graph(self.mols[idx])

        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
            usr=torch.tensor(self.usr[idx], dtype=torch.float).unsqueeze(0),
            y=torch.tensor(self.y[idx], dtype=torch.float)
        )
        return data

# ===============================
# 3. GNN Encoder del Ligando
# ===============================
# class LigandEncoder(nn.Module):
#     def __init__(self, node_dim, emb_dim=128):
#         super().__init__()
#         self.conv1 = GINEConv(node_dim, 64)
#         self.conv2 = GINEConv(64, emb_dim)

#     def forward(self, data):
#         x, edge_index, batch = data.x, data.edge_index, data.batch
#         x = F.relu(self.conv1(x, edge_index))
#         x = self.conv2(x, edge_index)
#         x = global_mean_pool(x, batch)
#         return x

class LigandEncoder(nn.Module):
    def __init__(self, node_dim, edge_dim, emb_dim=128):
        super().__init__()

        nn1 = nn.Sequential(
            nn.Linear(node_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 64)
        )

        nn2 = nn.Sequential(
            nn.Linear(64, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )

        self.conv1 = GINEConv(nn1, edge_dim=edge_dim)
        self.conv2 = GINEConv(nn2, edge_dim=edge_dim)

    def forward(self, data):
        x, edge_index, edge_attr, batch = \
            data.x, data.edge_index, data.edge_attr, data.batch

        x = F.relu(self.conv1(x, edge_index, edge_attr))
        x = self.conv2(x, edge_index, edge_attr)
        return global_mean_pool(x, batch)

# ===============================
# 4. Modelo Final (Ligando + USR)
# ===============================
class DockingTimeModel(nn.Module):
    def __init__(self, node_dim, edge_dim, usr_dim=12, emb_dim=128):
        super().__init__()
        self.encoder = LigandEncoder(node_dim, edge_dim)

        self.mlp = nn.Sequential(
            nn.Linear(emb_dim + usr_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 16),
            nn.ReLU(),
            nn.Linear(16, 1)
        )

    def forward(self, data):
        ligand_emb = self.encoder(data)
        usr = data.usr
        combined = torch.cat([ligand_emb, usr], dim=1)
        return self.mlp(combined)

# ===============================
# 5. Entrenamiento
# ===============================
def train(model, loader, optimizer, loss_fn, device):
    model.train()
    total_loss = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred = model(batch)
        loss = loss_fn(pred, batch.y.view(-1, 1))
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)

# ===============================
# 6. Evaluación con métricas
# ===============================
def evaluate(model, loader):
    model.eval()
    y_true, y_pred = [], []

    with torch.no_grad():
        for batch in loader:
            pred = model(batch)
            y_true.append(batch.y.cpu().numpy())
            y_pred.append(pred.cpu().numpy())

    y_true = np.vstack(y_true).ravel()
    y_pred = np.vstack(y_pred).ravel()

    return {
        "RMSE": mean_squared_error(y_true, y_pred, squared=False),
        "MAE": mean_absolute_error(y_true, y_pred),
        "R2": r2_score(y_true, y_pred)
    }


def validate(model, loader, loss_fn, device):
    model.eval()
    total_loss = 0
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            pred = model(batch)
            loss = loss_fn(pred, batch.y.view(-1, 1))
            total_loss += loss.item()
    return total_loss / len(loader)

# ===============================
# 7. Ejecución completa
# ===============================
# mols        → lista de RDKit Mol
# usr_vectors → array shape (N, 12)
# y_values    → tiempos de docking

# train_dataset = DockingDataset(df_train['mol'].tolist(), 
#                          df_train['usr'].tolist(),
#                          df_train['Value'].tolist())

# test_dataset = DockingDataset(df_test['mol'].tolist(), 
#                          df_test['usr'].tolist(),
#                          df_test['Value'].tolist())

# train_loader = DataLoader(train_dataset, batch_size=256, shuffle = True)
# test_loader  = DataLoader(test_dataset, batch_size=256)

# model = DockingTimeModel(node_dim=4).to('cuda')
# optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
# loss_fn = nn.MSELoss()

# for epoch in range(50):
#     loss = train(model, train_loader, optimizer, loss_fn)
#     print(f"Epoch {epoch:03d} | Train Loss: {loss:.4f}")

# metrics = evaluate(model, test_loader)
# print("Test metrics:", metrics)
# ======================================================
# 6. Data preparation (train / val / test)
# ======================================================
device = torch.device("cuda")

full_train_dataset = DockingDataset(
    df_train['mol'].tolist(),
    df_train['usr'].tolist(),
    df_train['Value'].tolist()
)

n_total = len(full_train_dataset)
n_val = int(0.15 * n_total)
n_train = n_total - n_val

train_dataset, val_dataset = random_split(
    full_train_dataset,
    [n_train, n_val],
    generator=torch.Generator().manual_seed(42)
)

test_dataset = DockingDataset(
    df_test['mol'].tolist(),
    df_test['usr'].tolist(),
    df_test['Value'].tolist()
)

train_loader = DataLoader(train_dataset, batch_size=256, shuffle=True, pin_memory=True)
val_loader   = DataLoader(val_dataset, batch_size=256, shuffle=False, pin_memory=True)
test_loader  = DataLoader(test_dataset, batch_size=256, shuffle=False, pin_memory=True)

# ======================================================
# 7. Model, optimizer, step decay
# ======================================================
model = DockingTimeModel(node_dim=18, edge_dim = 7).to(device)

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

scheduler = torch.optim.lr_scheduler.StepLR(
    optimizer,
    step_size=10,   # decay every 10 epochs
    gamma=0.5
)

loss_fn = nn.MSELoss()


# ======================================================
# 8. Training loop
# ======================================================
n_epochs = 50

for epoch in range(n_epochs):
    train_loss = train(model, train_loader, optimizer, loss_fn, device)
    val_loss = validate(model, val_loader, loss_fn, device)
    scheduler.step()

    lr = optimizer.param_groups[0]["lr"]

    print(
        f"Epoch {epoch:03d} | "
        f"Train Loss: {train_loss:.4f} | "
        f"Val Loss: {val_loss:.4f} | "
        f"LR: {lr:.2e}"
    )


# ======================================================
# 9. Final test metrics
# ======================================================
metrics = evaluate(model, test_loader, device)
print("Test metrics:", metrics)






