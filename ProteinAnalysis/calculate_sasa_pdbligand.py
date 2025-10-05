# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 12:16:02 2025

@author: andres.sanchez
"""

from Bio.PDB import ShrakeRupley, Structure, Model, Chain, Residue

def compute_relative_sasa(structure, chain_id, ligand_resname, ligand_resid):
    """
    Compute the relative SASA of a ligand in a protein complex.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        Input structure object.
    chain_id : str
        Chain ID containing the ligand.
    ligand_resname : str
        Residue name of the ligand (e.g., 'ATP').
    ligand_resid : int
        Residue number of the ligand.

    Returns
    -------
    float
        Relative SASA (SASA_in / SASA_out).
    """

    # --- 1. Get the ligand residue ---
    ligand = None
    for res in structure[0][chain_id]:
        hetflag, resseq, icode = res.id
        if hetflag.strip() != "" and res.get_resname() == ligand_resname and resseq == ligand_resid:
            ligand = res
            break
    if ligand is None:
        raise ValueError("Ligand not found in structure")

    # --- 2. Remove water from the structure ---
    struct_copy = structure.copy()
    for chain in list(struct_copy.get_chains()):
        for residue in list(chain):
            hetflag = residue.id[0]
            if residue.get_resname() == "HOH":  # water
                chain.detach_child(residue.id)

    # --- 3. Compute SASA of ligand in protein context ---
    sr = ShrakeRupley()
    sr.compute(struct_copy, level="R")
    sasa_in = [res.sasa for res in struct_copy.get_residues() if res == ligand][0]

    # --- 4. Compute SASA of ligand in isolation ---
    empty_struct = Structure.Structure("ligand_only")
    empty_model = Model.Model(0)
    empty_chain = Chain.Chain(chain_id)
    empty_model.add(empty_chain)
    empty_chain.add(ligand.copy())
    empty_struct.add(empty_model)

    sr = ShrakeRupley()
    sr.compute(empty_struct, level="R")
    sasa_out = [res.sasa for res in empty_struct.get_residues()][0]

    # --- 5. Return relative SASA ---
    return round(sasa_in / sasa_out, 3)

from Bio.PDB import PDBParser, MMCIFParser

parser = MMCIFParser(QUIET=True)
struct = parser.get_structure("test", "C:/Users/andres.sanchez/Downloads/3USN.cif")

rel_sasa = compute_relative_sasa(struct, chain_id="A", ligand_resname="ATT", ligand_resid=174)
print("Relative SASA:", rel_sasa)
