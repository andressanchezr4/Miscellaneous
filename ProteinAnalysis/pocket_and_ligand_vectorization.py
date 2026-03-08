#!/usr/bin/env python3
"""
Pocket + Ligand vectorization pipeline (clean version).

Este script reemplaza `TODO_pocket_and_ligand_vectorization.py` con una versión
modular y ejecutable por argumentos.

Incluye dos subcomandos:
1) `prepare`: calcula descriptores USR del pocket, construye moléculas RDKit de
   ligandos, fusiona datasets y genera split train/test.
2) `train`: entrena un modelo GNN (Torch Geometric) con features de ligando + USR.

Dependencias para `prepare`:
- pandas, numpy, scipy, biopython, rdkit

Dependencias adicionales para `train`:
- torch, torch_geometric, scikit-learn

Ejemplo (solo preparación):
    python pocket_ligand_pipeline_clean.py prepare \
      --p2rank-config-dir /data/p2rank_batches \
      --structures-root /data/uniprot \
      --ligands-dir /data/raw_ligands2 \
      --protlig-csv /data/protligeff.csv \
      --test-proteins-csv /data/prot2test.csv \
      --test-ligands-csv /data/ligands2test2.csv \
      --out-dir /data/processed

Ejemplo (entrenamiento):
    python pocket_ligand_pipeline_clean.py train \
      --train-pkl /data/processed/train.pkl \
      --test-pkl /data/processed/test.pkl \
      --target-col Value \
      --epochs 50
"""

from __future__ import annotations

import argparse
import logging
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import skew


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def compute_usr_descriptor(coords: np.ndarray) -> np.ndarray:
    """Calcula descriptor USR de 12 dimensiones desde coordenadas N x 3."""
    if coords.size == 0:
        return np.full(12, np.nan, dtype=float)

    ctd = coords.mean(axis=0)

    dists_to_ctd = np.linalg.norm(coords - ctd, axis=1)
    cst = coords[dists_to_ctd.argmin()]
    fct = coords[dists_to_ctd.argmax()]

    dists_to_fct = np.linalg.norm(coords - fct, axis=1)
    ftf = coords[dists_to_fct.argmax()]

    refs = [ctd, cst, fct, ftf]
    usr_vector: list[float] = []

    for ref in refs:
        distances = np.linalg.norm(coords - ref, axis=1)
        usr_vector.append(float(distances.mean()))
        usr_vector.append(float(distances.var()))
        usr_vector.append(float(skew(distances)))

    return np.array(usr_vector, dtype=float)


def _extract_floats(text: str) -> list[float]:
    floats = []
    for tok in text.replace(",", " ").split():
        try:
            floats.append(float(tok))
        except ValueError:
            continue
    return floats


def parse_p2rank_config(config_path: Path) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """
    Devuelve (center_xyz, size_xyz).

    Conserva la lógica del script original usando índices esperados en tokens,
    con fallback por posición relativa si el formato varía levemente.
    """
    tokens = config_path.read_text(encoding="utf-8", errors="ignore").split()
    floats = _extract_floats(" ".join(tokens))

    # Camino principal: replicar intención del script original.
    try:
        center = (float(tokens[2]), float(tokens[5]), float(tokens[8]))
        size = (float(tokens[-7]), float(tokens[-4]), float(tokens[-1]))
        return center, size
    except Exception:
        pass

    # Fallback robusto: usar primeras 3 y últimas 3 cifras detectadas.
    if len(floats) >= 6:
        center = tuple(floats[:3])
        size = tuple(floats[-3:])
        return center, size  # type: ignore[return-value]

    raise ValueError(f"No se pudo parsear center/size en {config_path}")


def build_pocket_config_table(p2rank_config_dir: Path) -> pd.DataFrame:
    rows = []
    for cfg in sorted(p2rank_config_dir.glob("*.txt")):
        if "P2RANK" not in cfg.name:
            continue

        split = cfg.stem.split("_")
        if len(split) < 5:
            logging.warning("Formato inesperado de nombre de config: %s", cfg.name)
            continue

        uniprot = split[0]
        pdb = "_".join(split[:5])
        pdb_file = f"{cfg.stem}.pdb"

        try:
            center, size = parse_p2rank_config(cfg)
        except Exception as exc:
            logging.warning("Config inválida %s: %s", cfg.name, exc)
            continue

        rows.append(
            {
                "uniprot": uniprot,
                "pdb": pdb,
                "pdb_file": pdb_file,
                "pocket_center": center,
                "pocket_size": size,
            }
        )

    return pd.DataFrame(rows)


def extract_pocket_coords(structure_path: Path, center: tuple[float, float, float], size: tuple[float, float, float]) -> np.ndarray:
    cx, cy, cz = center
    sx, sy, sz = size

    xmin, xmax = cx - sx / 2, cx + sx / 2
    ymin, ymax = cy - sy / 2, cy + sy / 2
    zmin, zmax = cz - sz / 2, cz + sz / 2

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", str(structure_path))

    coords = []
    for atom in structure.get_atoms():
        x, y, z = atom.coord
        if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
            coords.append([x, y, z])

    return np.array(coords, dtype=float)


def add_usr_to_configs(df_configs: pd.DataFrame, structures_root: Path) -> pd.DataFrame:
    usr_vectors = []
    protein_ids = []

    for _, row in df_configs.iterrows():
        structure_path = structures_root / row["uniprot"] / "selected_cif_chains" / row["pdb_file"]
        if not structure_path.exists():
            usr_vectors.append(np.full(12, np.nan, dtype=float))
            protein_ids.append("_".join(Path(row["pdb_file"]).stem.split("_")[1:8]))
            logging.warning("No existe estructura: %s", structure_path)
            continue

        coords = extract_pocket_coords(structure_path, row["pocket_center"], row["pocket_size"])
        usr_vectors.append(compute_usr_descriptor(coords))
        protein_ids.append("_".join(Path(row["pdb_file"]).stem.split("_")[1:8]))

    out = df_configs.copy()
    out["usr"] = usr_vectors
    out["Protein"] = protein_ids
    return out


def mol_from_sdf(sdf_path: Path):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is not None:
        AllChem.ComputeGasteigerCharges(mol)
    return mol


def build_ligand_table(ligands_dir: Path) -> pd.DataFrame:
    rows = []
    for sdf in sorted(ligands_dir.glob("*.sdf")):
        mol = mol_from_sdf(sdf)
        if mol is None:
            logging.warning("Ligando inválido: %s", sdf.name)
            continue
        rows.append({"Ligand": sdf.stem, "mol": mol})
    return pd.DataFrame(rows)


def has_bromine(mol) -> bool:
    return any(atom.GetSymbol() == "Br" for atom in mol.GetAtoms())


def split_train_test(df: pd.DataFrame, test_proteins: set[str], test_ligands: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_train = df[~(df["Protein"].isin(test_proteins) | df["Ligand"].isin(test_ligands))].copy()
    df_test = df[(df["Protein"].isin(test_proteins) & df["Ligand"].isin(test_ligands))].copy()

    # Mantener la lógica original: moléculas con Br pasan a train.
    df_br = df_test[df_test["mol"].apply(has_bromine)].copy()
    df_test = df_test[~df_test["mol"].apply(has_bromine)].copy()
    df_train = pd.concat([df_train, df_br], ignore_index=True)

    return df_train, df_test


def ensure_columns(df: pd.DataFrame, required: Iterable[str], df_name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} no contiene columnas requeridas: {missing}")


def cmd_prepare(args: argparse.Namespace) -> None:
    set_seed(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df_configs = build_pocket_config_table(args.p2rank_config_dir)
    if df_configs.empty:
        raise RuntimeError("No se encontraron configs P2RANK válidas.")

    df_configs = add_usr_to_configs(df_configs, args.structures_root)

    df_all = pd.read_csv(args.protlig_csv)
    ensure_columns(df_all, ["Protein", "Ligand", args.target_col], "protlig_csv")

    df_lig = build_ligand_table(args.ligands_dir)
    df_merged = df_all.merge(df_lig, on="Ligand", how="left")
    df_merged = df_merged.merge(df_configs[["Protein", "usr"]], on="Protein", how="left")

    # Filtrar filas sin mol o sin usr válido.
    df_merged = df_merged[df_merged["mol"].notna()].copy()
    df_merged = df_merged[df_merged["usr"].apply(lambda x: isinstance(x, np.ndarray) and np.isfinite(x).all())].copy()

    test_proteins = set(pd.read_csv(args.test_proteins_csv)["Protein"].astype(str).tolist())
    test_ligands = set(pd.read_csv(args.test_ligands_csv)["Ligand"].astype(str).tolist())

    df_train, df_test = split_train_test(df_merged, test_proteins, test_ligands)

    # Guardados: pickle (completo con objetos RDKit) + csv (sin columna mol)
    merged_pkl = args.out_dir / "merged.pkl"
    train_pkl = args.out_dir / "train.pkl"
    test_pkl = args.out_dir / "test.pkl"

    df_merged.to_pickle(merged_pkl)
    df_train.to_pickle(train_pkl)
    df_test.to_pickle(test_pkl)

    df_merged.drop(columns=["mol"]).to_csv(args.out_dir / "merged.csv", index=False)
    df_train.drop(columns=["mol"]).to_csv(args.out_dir / "train.csv", index=False)
    df_test.drop(columns=["mol"]).to_csv(args.out_dir / "test.csv", index=False)

    logging.info("Preparación completada.")
    logging.info("merged: %s", merged_pkl)
    logging.info("train:  %s", train_pkl)
    logging.info("test:   %s", test_pkl)


def cmd_train(args: argparse.Namespace) -> None:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from torch.utils.data import random_split
    from torch_geometric.data import Data, Dataset
    from torch_geometric.loader import DataLoader
    from torch_geometric.nn import GINEConv, global_mean_pool

    set_seed(args.seed)
    torch.manual_seed(args.seed)

    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))

    df_train = pd.read_pickle(args.train_pkl)
    df_test = pd.read_pickle(args.test_pkl)
    ensure_columns(df_train, ["mol", "usr", args.target_col], "train_pkl")
    ensure_columns(df_test, ["mol", "usr", args.target_col], "test_pkl")

    def mol_to_graph(mol):
        atom_types = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P"]
        node_feats = []
        for atom in mol.GetAtoms():
            hyb = atom.GetHybridization()
            symbol = atom.GetSymbol()
            ohe = [int(symbol == t) for t in atom_types]
            node_feats.append(
                [
                    atom.GetAtomicNum(),
                    atom.GetTotalDegree(),
                    atom.GetTotalNumHs(),
                    int(atom.GetIsAromatic()),
                    int(atom.IsInRing()),
                    int(hyb == Chem.rdchem.HybridizationType.SP),
                    int(hyb == Chem.rdchem.HybridizationType.SP2),
                    int(hyb == Chem.rdchem.HybridizationType.SP3),
                    int(hyb == Chem.rdchem.HybridizationType.SP3 and not atom.IsInRing()),
                    *ohe,
                ]
            )
        x = torch.tensor(node_feats, dtype=torch.float)

        edge_index, edge_attr = [], []
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bt = bond.GetBondType()
            bf = [
                int(bt == Chem.rdchem.BondType.SINGLE),
                int(bt == Chem.rdchem.BondType.DOUBLE),
                int(bt == Chem.rdchem.BondType.TRIPLE),
                int(bt == Chem.rdchem.BondType.AROMATIC),
                int(bond.IsInRing()),
                int(bond.GetIsConjugated()),
                int(bt == Chem.rdchem.BondType.SINGLE and not bond.IsInRing()),
            ]
            edge_index += [[i, j], [j, i]]
            edge_attr += [bf, bf]

        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)
        return x, edge_index, edge_attr

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
            return Data(
                x=x,
                edge_index=edge_index,
                edge_attr=edge_attr,
                usr=torch.tensor(self.usr[idx], dtype=torch.float).unsqueeze(0),
                y=torch.tensor(self.y[idx], dtype=torch.float),
            )

    class LigandEncoder(nn.Module):
        def __init__(self, node_dim, edge_dim, emb_dim=128):
            super().__init__()
            nn1 = nn.Sequential(nn.Linear(node_dim, 64), nn.ReLU(), nn.Linear(64, 64))
            nn2 = nn.Sequential(nn.Linear(64, emb_dim), nn.ReLU(), nn.Linear(emb_dim, emb_dim))
            self.conv1 = GINEConv(nn1, edge_dim=edge_dim)
            self.conv2 = GINEConv(nn2, edge_dim=edge_dim)

        def forward(self, data):
            x = F.relu(self.conv1(data.x, data.edge_index, data.edge_attr))
            x = self.conv2(x, data.edge_index, data.edge_attr)
            return global_mean_pool(x, data.batch)

    class DockingTimeModel(nn.Module):
        def __init__(self, node_dim, edge_dim, usr_dim=12, emb_dim=128):
            super().__init__()
            self.encoder = LigandEncoder(node_dim, edge_dim, emb_dim)
            self.mlp = nn.Sequential(
                nn.Linear(emb_dim + usr_dim, 128),
                nn.ReLU(),
                nn.Linear(128, 64),
                nn.ReLU(),
                nn.Linear(64, 32),
                nn.ReLU(),
                nn.Linear(32, 16),
                nn.ReLU(),
                nn.Linear(16, 1),
            )

        def forward(self, data):
            ligand_emb = self.encoder(data)
            combined = torch.cat([ligand_emb, data.usr], dim=1)
            return self.mlp(combined)

    def train_epoch(model, loader, optimizer, loss_fn):
        model.train()
        total = 0.0
        for batch in loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            pred = model(batch)
            loss = loss_fn(pred, batch.y.view(-1, 1))
            loss.backward()
            optimizer.step()
            total += loss.item()
        return total / max(len(loader), 1)

    def validate_epoch(model, loader, loss_fn):
        model.eval()
        total = 0.0
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(device)
                pred = model(batch)
                total += loss_fn(pred, batch.y.view(-1, 1)).item()
        return total / max(len(loader), 1)

    def evaluate(model, loader):
        model.eval()
        y_true, y_pred = [], []
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(device)
                pred = model(batch)
                y_true.append(batch.y.cpu().numpy())
                y_pred.append(pred.cpu().numpy())

        y_true = np.vstack(y_true).ravel()
        y_pred = np.vstack(y_pred).ravel()
        return {
            "RMSE": float(mean_squared_error(y_true, y_pred, squared=False)),
            "MAE": float(mean_absolute_error(y_true, y_pred)),
            "R2": float(r2_score(y_true, y_pred)),
        }

    full_train_dataset = DockingDataset(
        df_train["mol"].tolist(),
        df_train["usr"].tolist(),
        df_train[args.target_col].tolist(),
    )
    test_dataset = DockingDataset(
        df_test["mol"].tolist(),
        df_test["usr"].tolist(),
        df_test[args.target_col].tolist(),
    )

    n_total = len(full_train_dataset)
    n_val = int(args.val_fraction * n_total)
    n_train = n_total - n_val

    train_dataset, val_dataset = random_split(
        full_train_dataset,
        [n_train, n_val],
        generator=torch.Generator().manual_seed(args.seed),
    )

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False, pin_memory=True)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False, pin_memory=True)

    model = DockingTimeModel(node_dim=18, edge_dim=7, usr_dim=12, emb_dim=args.emb_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=args.lr_step_size, gamma=args.lr_gamma)
    loss_fn = torch.nn.MSELoss()

    for epoch in range(args.epochs):
        train_loss = train_epoch(model, train_loader, optimizer, loss_fn)
        val_loss = validate_epoch(model, val_loader, loss_fn)
        scheduler.step()
        lr = optimizer.param_groups[0]["lr"]
        logging.info("Epoch %03d | Train %.4f | Val %.4f | LR %.2e", epoch, train_loss, val_loss, lr)

    metrics = evaluate(model, test_loader)
    logging.info("Test metrics: %s", metrics)

    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": model.state_dict(), "metrics": metrics}, args.model_out)
    logging.info("Modelo guardado en: %s", args.model_out)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pocket + ligand vectorization and optional GNN training.")
    parser.add_argument("--verbose", action="store_true", help="Activa logs detallados.")
    parser.add_argument("--seed", type=int, default=42, help="Semilla global reproducible.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    p_prepare = subparsers.add_parser("prepare", help="Preparar features de pocket y ligando + split train/test.")
    p_prepare.add_argument("--p2rank-config-dir", type=Path, required=True, help="Directorio con configs .txt de P2RANK.")
    p_prepare.add_argument("--structures-root", type=Path, required=True, help="Raíz con /<uniprot>/selected_cif_chains/<pdb_file>.")
    p_prepare.add_argument("--ligands-dir", type=Path, required=True, help="Directorio con ligandos .sdf.")
    p_prepare.add_argument("--protlig-csv", type=Path, required=True, help="CSV maestro con columnas Protein, Ligand y target.")
    p_prepare.add_argument("--test-proteins-csv", type=Path, required=True, help="CSV con columna Protein para test.")
    p_prepare.add_argument("--test-ligands-csv", type=Path, required=True, help="CSV con columna Ligand para test.")
    p_prepare.add_argument("--target-col", default="Value", help="Nombre de columna objetivo en protlig_csv.")
    p_prepare.add_argument("--out-dir", type=Path, required=True, help="Directorio de salida para merged/train/test.")
    p_prepare.set_defaults(func=cmd_prepare)

    p_train = subparsers.add_parser("train", help="Entrenar GNN con datos preparados.")
    p_train.add_argument("--train-pkl", type=Path, required=True, help="Archivo train.pkl generado por prepare.")
    p_train.add_argument("--test-pkl", type=Path, required=True, help="Archivo test.pkl generado por prepare.")
    p_train.add_argument("--target-col", default="Value", help="Nombre de columna objetivo.")
    p_train.add_argument("--device", default="", help="Dispositivo PyTorch, ej: cuda o cpu. Vacío = auto.")
    p_train.add_argument("--epochs", type=int, default=50)
    p_train.add_argument("--batch-size", type=int, default=256)
    p_train.add_argument("--lr", type=float, default=1e-3)
    p_train.add_argument("--lr-step-size", type=int, default=10)
    p_train.add_argument("--lr-gamma", type=float, default=0.5)
    p_train.add_argument("--val-fraction", type=float, default=0.15)
    p_train.add_argument("--emb-dim", type=int, default=128)
    p_train.add_argument("--model-out", type=Path, default=Path("./docking_time_model.pt"))
    p_train.set_defaults(func=cmd_train)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    configure_logging(args.verbose)
    args.func(args)


if __name__ == "__main__":
    main()
