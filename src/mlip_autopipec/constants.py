"""
This module defines physical and chemical constants used throughout the application.

Storing these values in a central location avoids magic numbers and makes the
code more maintainable and readable.
"""

# Recommended plane-wave cutoffs (in Ry) for SSSP Precision pseudopotentials
# These values are critical for ensuring the accuracy of DFT calculations.
SSSP_PRECISION_CUTOFFS: dict[str, tuple[int, int]] = {
    "H": (50, 400),
    "He": (60, 480),
    "Li": (40, 320),
    "Be": (50, 400),
    "B": (60, 480),
    "C": (60, 480),
    "N": (70, 560),
    "O": (70, 560),
    "F": (80, 640),
    "Ne": (90, 720),
    "Na": (40, 320),
    "Mg": (40, 320),
    "Al": (40, 320),
    "Si": (40, 320),
    "P": (50, 400),
    "S": (50, 400),
    "Cl": (60, 480),
    "Ar": (70, 560),
    "K": (40, 320),
    "Ca": (40, 320),
    "Sc": (60, 480),
    "Ti": (70, 560),
    "V": (70, 560),
    "Cr": (70, 560),
    "Mn": (70, 560),
    "Fe": (80, 640),
    "Co": (80, 640),
    "Ni": (80, 640),
    "Cu": (70, 560),
    "Zn": (70, 560),
    "Ga": (60, 480),
    "Ge": (60, 480),
    "As": (60, 480),
    "Se": (60, 480),
    "Br": (60, 480),
    "Kr": (70, 560),
    "Rb": (40, 320),
    "Sr": (40, 320),
    "Y": (60, 480),
    "Zr": (60, 480),
    "Nb": (70, 560),
    "Mo": (70, 560),
    "Tc": (70, 560),
    "Ru": (70, 560),
    "Rh": (70, 560),
    "Pd": (70, 560),
    "Ag": (60, 480),
    "Cd": (60, 480),
    "In": (60, 480),
    "Sn": (60, 480),
    "Sb": (60, 480),
    "Te": (60, 480),
    "I": (60, 480),
    "Xe": (70, 560),
    "Cs": (40, 320),
    "Ba": (40, 320),
    "La": (60, 480),
    "Ce": (60, 480),
    "Pr": (60, 480),
    "Nd": (60, 480),
    "Pm": (60, 480),
    "Sm": (60, 480),
    "Eu": (60, 480),
    "Gd": (70, 560),
    "Tb": (70, 560),
    "Dy": (70, 560),
    "Ho": (70, 560),
    "Er": (70, 560),
    "Tm": (70, 560),
    "Yb": (70, 560),
    "Lu": (70, 560),
    "Hf": (70, 560),
    "Ta": (70, 560),
    "W": (80, 640),
    "Re": (80, 640),
    "Os": (80, 640),
    "Ir": (80, 640),
    "Pt": (90, 720),
    "Au": (80, 640),
    "Hg": (70, 560),
    "Tl": (60, 480),
    "Pb": (60, 480),
    "Bi": (60, 480),
    "Po": (60, 480),
    "At": (60, 480),
    "Rn": (70, 560),
}

# Set of elements for which ferromagnetic calculations are enabled by default
MAGNETIC_ELEMENTS: set[str] = {"Fe", "Co", "Ni"}

# Dictionary of elemental melting points (in Kelvin) used for heuristic guesses
ELEMENT_MELTING_POINTS: dict[str, int] = {
    "Fe": 1811,
    "Pt": 2041,
    "Si": 1687,
    "Al": 933,
}

# --- Structure Generation Constants ---
ALLOY_TARGET_ATOMS: int = 64
ALLOY_STRAIN_MAGNITUDE: float = 0.05
