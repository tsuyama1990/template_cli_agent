import re

import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.data.models import DFTCompute


def create_qe_input_from_atoms(atoms: Atoms, config: DFTCompute) -> str:
    """Creates a Quantum Espresso input file from an ASE Atoms object."""
    symbols = sorted(list(set(atoms.get_chemical_symbols())))
    nat = len(atoms)
    ntyp = len(symbols)

    input_lines = [
        "&CONTROL",
        "    calculation = 'scf'",
        "    restart_mode = 'from_scratch'",
        "    prefix = 'pwscf'",
        "    outdir = './out'",
        "    pseudo_dir = './pseudos'",
        "/",
        "&SYSTEM",
        f"    ibrav = 0, nat = {nat}, ntyp = {ntyp}",
        f"    ecutwfc = {config.ecutwfc}",
        f"    ecutrho = {config.ecutrho}",
        "/",
        "&ELECTRONS",
        "    conv_thr = 1.0e-8",
        "    mixing_beta = 0.7",
        "/",
        "ATOMIC_SPECIES",
    ]

    for symbol in symbols:
        input_lines.append(f"    {symbol}  1.0  {symbol}.{config.pseudopotentials}.UPF")

    input_lines.append("ATOMIC_POSITIONS {angstrom}")
    for atom in atoms:
        input_lines.append(
            f"    {atom.symbol}  {atom.position[0]:.8f}  {atom.position[1]:.8f}  {atom.position[2]:.8f}"
        )

    input_lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.cell:
        input_lines.append(f"    {vector[0]:.8f}  {vector[1]:.8f}  {vector[2]:.8f}")

    k_points = np.ceil(
        np.linalg.norm(atoms.cell, axis=1) * config.kpoints_density
    ).astype(int)
    input_lines.append("K_POINTS {automatic}")
    input_lines.append(f"    {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0")

    return "\n".join(input_lines)


def parse_qe_output(output: str) -> dict:
    """Parses a Quantum Espresso output file to extract energy, forces, and stress."""
    results = {}

    energy_match = re.search(r"!    total energy\s+=\s+(-?[\d\.]+)", output)
    if energy_match:
        # QE output is in Ry, convert to eV for ASE
        results["energy"] = float(energy_match.group(1)) * 13.605693122994

    forces_block_match = re.search(
        r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n\s*\n(.*?)(?:\n\n|Total force|add F_BERRY|\Z)",
        output,
        re.DOTALL,
    )
    if forces_block_match:
        forces_block = forces_block_match.group(1)
        # QE forces are in Ry/au, convert to eV/Angstrom for ASE
        ry_au_to_ev_ang = 13.605693122994 / 0.529177210903
        forces_match = re.findall(
            r"atom\s+\d+\s+type\s+\d+\s+force =\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",
            forces_block,
        )
        if forces_match:
            results["forces"] = (
                np.array(forces_match, dtype=float) * ry_au_to_ev_ang
            ).tolist()

    stress_match = re.search(
        r"total stress .*? p\(kbar\) =\s*\n\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",
        output,
        re.DOTALL,
    )
    if stress_match:
        # Convert from kbar to eV/A^3
        kbar_to_ev_ang3 = 1.0 / 1602.1766208
        stress_voigt = (
            np.array([float(s) for s in stress_match.groups()]) * kbar_to_ev_ang3
        )
        stress_matrix = np.array(
            [
                [stress_voigt[0], stress_voigt[3], stress_voigt[4]],
                [stress_voigt[3], stress_voigt[1], stress_voigt[5]],
                [stress_voigt[4], stress_voigt[5], stress_voigt[2]],
            ]
        ).tolist()
        results["stress"] = stress_matrix

    return results
