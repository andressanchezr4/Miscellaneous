#!/usr/bin/env python3
"""
Pipeline limpio para clustering estructural con Foldseek.

Este script reemplaza el flujo de trabajo del archivo TODO original:
1) Lee un listado de proteínas priorizadas por target.
2) Mapea UniProt -> PDB.
3) Descarga estructuras CIF desde PDBe.
4) En modo mono: recorta cadena de interés + hetero-residuos cercanos.
5) Ejecuta Foldseek (easy-cluster o easy-multimercluster).
6) Parsea los clústeres desde *_clu.tsv.
7) Exporta un CSV final con metadatos + asignación de clúster.

Uso mínimo:
    python foldseek_pipeline_clean.py \
      --df-uniprots prioritized_prots_to_cluster.csv \
      --uniprot2pdb uniprot2pdb.csv \
      --uniprot2pdb-raw uniprot2pdb_raw.csv \
      --analysis-folder ./foldseek_results \
      --output-csv ./prioritized_prots_to_cluster_pdb.csv

Notas:
- Requiere: pandas, requests, biopython y foldseek en PATH.
- El CSV de entrada debe contener al menos: tar_name, uniprot_id.
- En modo mono se usa uniprot2pdb_raw para resolver cadenas por PDB.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests
from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure


# Sustituciones de residuos no estándar a aminoácidos estándar
SUBSTITUTIONS = {
    "2AS": "ASP", "3AH": "HIS", "5HP": "GLU", "5OW": "LYS", "ACL": "ARG", "AGM": "ARG", "AIB": "ALA", "ALM": "ALA",
    "ALO": "THR", "ALY": "LYS", "ARM": "ARG", "ASA": "ASP", "ASB": "ASP", "ASK": "ASP", "ASL": "ASP", "ASQ": "ASP",
    "AYA": "ALA", "BCS": "CYS", "BHD": "ASP", "BMT": "THR", "BNN": "ALA", "BUC": "CYS", "BUG": "LEU", "C5C": "CYS",
    "C6C": "CYS", "CAS": "CYS", "CCS": "CYS", "CEA": "CYS", "CGU": "GLU", "CHG": "ALA", "CLE": "LEU", "CME": "CYS",
    "CSD": "ALA", "CSO": "CYS", "CSP": "CYS", "CSS": "CYS", "CSW": "CYS", "CSX": "CYS", "CXM": "MET", "CY1": "CYS",
    "CY3": "CYS", "CYG": "CYS", "CYM": "CYS", "CYQ": "CYS", "DAH": "PHE", "DAL": "ALA", "DAR": "ARG", "DAS": "ASP",
    "DCY": "CYS", "DGL": "GLU", "DGN": "GLN", "DHA": "ALA", "DHI": "HIS", "DIL": "ILE", "DIV": "VAL", "DLE": "LEU",
    "DLY": "LYS", "DNP": "ALA", "DPN": "PHE", "DPR": "PRO", "DSN": "SER", "DSP": "ASP", "DTH": "THR", "DTR": "TRP",
    "DTY": "TYR", "DVA": "VAL", "EFC": "CYS", "FLA": "ALA", "FME": "MET", "GGL": "GLU", "GL3": "GLY", "GLZ": "GLY",
    "GMA": "GLU", "GSC": "GLY", "HAC": "ALA", "HAR": "ARG", "HIC": "HIS", "HIP": "HIS", "HMR": "ARG", "HPQ": "PHE",
    "HTR": "TRP", "HYP": "PRO", "IAS": "ASP", "IIL": "ILE", "IYR": "TYR", "KCX": "LYS", "LLP": "LYS", "LLY": "LYS",
    "LTR": "TRP", "LYM": "LYS", "LYZ": "LYS", "MAA": "ALA", "MEN": "ASN", "MHS": "HIS", "MIS": "SER", "MK8": "LEU",
    "MLE": "LEU", "MPQ": "GLY", "MSA": "GLY", "MSE": "MET", "MVA": "VAL", "NEM": "HIS", "NEP": "HIS", "NLE": "LEU",
    "NLN": "LEU", "NLP": "LEU", "NMC": "GLY", "OAS": "SER", "OCS": "CYS", "OMT": "MET", "PAQ": "TYR", "PCA": "GLU",
    "PEC": "CYS", "PHI": "PHE", "PHL": "PHE", "PR3": "CYS", "PRR": "ALA", "PTR": "TYR", "PYX": "CYS", "SAC": "SER",
    "SAR": "GLY", "SCH": "CYS", "SCS": "CYS", "SCY": "CYS", "SEL": "SER", "SEP": "SER", "SET": "SER", "SHC": "CYS",
    "SHR": "LYS", "SMC": "CYS", "SOC": "CYS", "STY": "TYR", "SVA": "SER", "TIH": "ALA", "TPL": "TRP", "TPO": "THR",
    "TPQ": "ALA", "TRG": "LYS", "MLY": "LYS", "TRO": "TRP", "TYB": "TYR", "TYI": "TYR", "TYQ": "TYR", "TYS": "TYR",
    "TYY": "TYR", "M3L": "LYS", "MLZ": "LYS"
}


@dataclass
class TargetPaths:
    target_name: str
    target_slug: str
    output_dir: Path
    foldseek_out_prefix: Path
    tmp_dir: Path


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pipeline de clustering estructural (Foldseek) para proteínas priorizadas."
    )
    parser.add_argument("--df-uniprots", type=Path, required=True, help="CSV con columnas tar_name y uniprot_id.")
    parser.add_argument("--uniprot2pdb", type=Path, required=True, help="CSV de mapeo UniProt->PDB.")
    parser.add_argument("--uniprot2pdb-raw", type=Path, required=True, help="CSV con columnas SP_PRIMARY, PDB, CHAIN.")
    parser.add_argument("--analysis-folder", type=Path, required=True, help="Carpeta raíz de resultados.")
    parser.add_argument("--output-csv", type=Path, required=True, help="CSV final con asignación de clústeres.")
    parser.add_argument("--mode", choices=["mono", "multimer"], default="mono", help="Tipo de clustering Foldseek.")
    parser.add_argument("--foldseek-bin", default="foldseek", help="Binario de Foldseek en PATH.")
    parser.add_argument("--coverage", type=float, default=0.8, help="Umbral -c para Foldseek.")
    parser.add_argument("--input-sep", default=";", help="Separador del CSV --df-uniprots.")
    parser.add_argument("--output-sep", default=";", help="Separador del CSV de salida.")
    parser.add_argument("--request-timeout", type=int, default=60, help="Timeout HTTP en segundos.")
    parser.add_argument("--keep-raw-cif", action="store_true", help="No borrar CIF temporal descargado.")
    parser.add_argument("--skip-existing", action="store_true", help="No re-descargar/re-generar CIF si ya existe.")
    parser.add_argument("--verbose", action="store_true", help="Logs detallados.")
    return parser.parse_args()


def sanitize_target_name(name: str) -> str:
    return "_".join(name.split())


def is_standard_or_substituted_residue(residue) -> bool:
    hetfield = residue.id[0]
    return hetfield == " " or hetfield[2:] in SUBSTITUTIONS


def select_residues_for_chain(chain_id: str, struct: Structure, include_het: bool = True) -> Structure:
    chain_id_lower = chain_id.lower()
    protein_atoms = [
        atom
        for residue in struct.get_residues()
        if is_standard_or_substituted_residue(residue) and residue.get_parent().id.lower() == chain_id_lower
        for atom in residue
    ]

    if not protein_atoms:
        return struct

    neighbor_search = NeighborSearch(protein_atoms)
    selected_residues = []

    for residue in struct.get_residues():
        if include_het and residue.id[0].startswith("H_"):
            total_atoms = 0
            close_atoms = 0
            for atom in residue:
                total_atoms += 1
                if neighbor_search.search(atom.coord, 5, level="R"):
                    close_atoms += 1
            if total_atoms > 0 and (close_atoms / total_atoms) > 0.3:
                selected_residues.append(residue)

        if is_standard_or_substituted_residue(residue) and residue.get_parent().id.lower() == chain_id_lower:
            selected_residues.append(residue)

    selected_set = set(selected_residues)
    for chain in struct.get_chains():
        for residue in list(chain):
            if residue not in selected_set:
                chain.detach_child(residue.id)

    return struct


def download_cif_from_pdbe(pdb_id: str, out_file: Path, timeout: int) -> bool:
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id.lower()}_updated.cif"
    response = requests.get(url, allow_redirects=True, timeout=timeout)
    if response.status_code != 200:
        logging.warning("No se pudo descargar %s (HTTP %s)", pdb_id, response.status_code)
        return False
    out_file.write_bytes(response.content)
    logging.info("Descargado %s -> %s", pdb_id, out_file.name)
    return True


def get_chains_for_uniprot_pdb(uniprot2pdb_raw: pd.DataFrame, uniprot: str, pdb_id: str) -> list[str]:
    mask = (
        (uniprot2pdb_raw["SP_PRIMARY"] == uniprot)
        & (uniprot2pdb_raw["PDB"].str.lower() == pdb_id.lower())
    )
    chains = uniprot2pdb_raw.loc[mask, "CHAIN"].dropna().astype(str).tolist()
    return list(dict.fromkeys(chains))


def download_and_prepare_mono(
    pdb_id: str,
    uniprot: str,
    output_dir: Path,
    uniprot2pdb_raw: pd.DataFrame,
    timeout: int,
    keep_raw_cif: bool,
    skip_existing: bool,
) -> None:
    raw_cif = output_dir / f"{uniprot}_{pdb_id}.cif"
    if not raw_cif.exists() or not skip_existing:
        ok = download_cif_from_pdbe(pdb_id, raw_cif, timeout)
        if not ok:
            return

    chains = get_chains_for_uniprot_pdb(uniprot2pdb_raw, uniprot, pdb_id)
    if not chains:
        logging.warning("Sin cadenas para %s %s en uniprot2pdb_raw", uniprot, pdb_id)
        if raw_cif.exists() and not keep_raw_cif:
            raw_cif.unlink(missing_ok=True)
        return

    for chain in chains:
        final_path = output_dir / f"{uniprot}_{pdb_id}_{chain}.cif"
        if final_path.exists() and skip_existing:
            continue
        try:
            parser = MMCIFParser(QUIET=True, auth_chains=True)
            struct = parser.get_structure(pdb_id, str(raw_cif))
            selected_struct = select_residues_for_chain(chain, struct, include_het=True)
            io = MMCIFIO()
            io.set_structure(selected_struct)
            io.save(str(final_path), preserve_atom_numbering=True)
        except Exception as exc:
            logging.warning("Error procesando %s %s cadena %s: %s", uniprot, pdb_id, chain, exc)

    if raw_cif.exists() and not keep_raw_cif:
        raw_cif.unlink(missing_ok=True)


def download_for_multimer(
    pdb_id: str,
    uniprot: str,
    output_dir: Path,
    timeout: int,
    skip_existing: bool,
) -> None:
    final_path = output_dir / f"{uniprot}_{pdb_id}.cif"
    if final_path.exists() and skip_existing:
        return
    download_cif_from_pdbe(pdb_id, final_path, timeout)


def foldseek_command(mode: str, foldseek_bin: str, output_dir: Path, out_prefix: Path, tmp_dir: Path, coverage: float) -> list[str]:
    subcommand = "easy-cluster" if mode == "mono" else "easy-multimercluster"
    return [
        foldseek_bin,
        subcommand,
        str(output_dir),
        str(out_prefix),
        str(tmp_dir),
        "-c",
        str(coverage),
    ]


def run_foldseek(cmd: list[str]) -> None:
    logging.info("Ejecutando: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def map_uniprot_to_pdbs(uniprot2pdb: pd.DataFrame, uniprot_id: str) -> list[str]:
    matches = uniprot2pdb[
        uniprot2pdb["SP_PRIMARY"].astype(str).apply(lambda x: uniprot_id in x.split("-"))
    ]
    if matches.empty:
        return []
    pdb_field = str(matches.iloc[0]["PDB"]).strip()
    if not pdb_field:
        return []
    return [p.strip() for p in pdb_field.split() if p.strip()]


def parse_foldseek_id(seq_id: str) -> tuple[str, str, str]:
    header = seq_id.strip().lstrip(">")
    id_part = header.split()[0]
    parts = id_part.split("_")
    uniprot_id = parts[0] if len(parts) > 0 else ""
    pdb = parts[1] if len(parts) > 1 else ""
    chain = "_".join(parts[2:]) if len(parts) > 2 else ""
    return uniprot_id, pdb, chain


def parse_clusters_from_tsv(cluster_tsv: Path, target_name: str) -> pd.DataFrame:
    if not cluster_tsv.exists():
        logging.warning("No existe el archivo de clusters: %s", cluster_tsv)
        return pd.DataFrame(columns=["uniprot_id", "pdb", "chain", "cluster_id", "is_representative", "tar_name"])

    pairs = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["representative", "member"]) 
    if pairs.empty:
        return pd.DataFrame(columns=["uniprot_id", "pdb", "chain", "cluster_id", "is_representative", "tar_name"])

    rep_order = list(dict.fromkeys(pairs["representative"].tolist()))
    cluster_index = {rep: idx for idx, rep in enumerate(rep_order)}
    target_prefix = sanitize_target_name(target_name)[:6]

    rows = []
    for _, row in pairs.iterrows():
        rep = str(row["representative"])
        member = str(row["member"])
        uniprot_id, pdb, chain = parse_foldseek_id(member)
        rows.append(
            {
                "uniprot_id": uniprot_id,
                "pdb": pdb,
                "chain": chain,
                "cluster_id": f"{target_prefix}_cluster_{cluster_index[rep]}",
                "is_representative": member == rep,
                "tar_name": target_name,
            }
        )

    return pd.DataFrame(rows)


def build_target_paths(analysis_folder: Path, target_name: str, mode: str) -> TargetPaths:
    target_slug = sanitize_target_name(target_name)
    root = analysis_folder / f"{mode}_{target_slug}"
    return TargetPaths(
        target_name=target_name,
        target_slug=target_slug,
        output_dir=root / "pdbs_selected",
        foldseek_out_prefix=root / "foldseek_clusters",
        tmp_dir=root / "tmp_foldseek",
    )


def ensure_required_columns(df: pd.DataFrame, required: Iterable[str], df_name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} no contiene columnas requeridas: {missing}")


def main() -> None:
    args = parse_args()
    configure_logging(args.verbose)

    args.analysis_folder.mkdir(parents=True, exist_ok=True)

    df_uniprots = pd.read_csv(args.df_uniprots, sep=args.input_sep)
    uniprot2pdb = pd.read_csv(args.uniprot2pdb)
    uniprot2pdb_raw = pd.read_csv(args.uniprot2pdb_raw)

    ensure_required_columns(df_uniprots, ["tar_name", "uniprot_id"], "df_uniprots")
    ensure_required_columns(uniprot2pdb, ["SP_PRIMARY", "PDB"], "uniprot2pdb")
    ensure_required_columns(uniprot2pdb_raw, ["SP_PRIMARY", "PDB", "CHAIN"], "uniprot2pdb_raw")

    tar2uniprots = (
        df_uniprots.groupby("tar_name")["uniprot_id"].apply(list).to_dict()
    )

    all_cluster_tables = []

    for tar_name, uniprots in tar2uniprots.items():
        paths = build_target_paths(args.analysis_folder, tar_name, args.mode)
        paths.output_dir.mkdir(parents=True, exist_ok=True)
        paths.tmp_dir.mkdir(parents=True, exist_ok=True)

        logging.info("Procesando target: %s", tar_name)

        for uniprot_id in sorted(set(uniprots)):
            pdb_ids = map_uniprot_to_pdbs(uniprot2pdb, uniprot_id)
            if not pdb_ids:
                logging.warning("Sin PDBs para %s", uniprot_id)
                continue

            for pdb_id in pdb_ids:
                if args.mode == "mono":
                    download_and_prepare_mono(
                        pdb_id=pdb_id,
                        uniprot=uniprot_id,
                        output_dir=paths.output_dir,
                        uniprot2pdb_raw=uniprot2pdb_raw,
                        timeout=args.request_timeout,
                        keep_raw_cif=args.keep_raw_cif,
                        skip_existing=args.skip_existing,
                    )
                else:
                    download_for_multimer(
                        pdb_id=pdb_id,
                        uniprot=uniprot_id,
                        output_dir=paths.output_dir,
                        timeout=args.request_timeout,
                        skip_existing=args.skip_existing,
                    )

        cmd = foldseek_command(
            mode=args.mode,
            foldseek_bin=args.foldseek_bin,
            output_dir=paths.output_dir,
            out_prefix=paths.foldseek_out_prefix,
            tmp_dir=paths.tmp_dir,
            coverage=args.coverage,
        )
        run_foldseek(cmd)

        cluster_tsv = Path(f"{paths.foldseek_out_prefix}_clu.tsv")
        cluster_df = parse_clusters_from_tsv(cluster_tsv, tar_name)
        all_cluster_tables.append(cluster_df)

    if all_cluster_tables:
        cluster_all = pd.concat(all_cluster_tables, ignore_index=True)
    else:
        cluster_all = pd.DataFrame(columns=["uniprot_id", "pdb", "chain", "cluster_id", "is_representative", "tar_name"])

    result = df_uniprots.merge(cluster_all, on=["uniprot_id", "tar_name"], how="left")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_csv, index=False, sep=args.output_sep)

    logging.info("Pipeline completado. CSV final: %s", args.output_csv)


if __name__ == "__main__":
    main()
