from __future__ import annotations

from ase.build import bulk, fcc111, molecule

from mlip_autopipec.core.physics import detect_vacuum


def test_detect_vacuum_bulk():
    """Test that a bulk FCC crystal is correctly classified as 'bulk'."""
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    classification = detect_vacuum(atoms)
    assert classification == "bulk"


def test_detect_vacuum_slab():
    """Test that a slab with a vacuum layer is correctly classified as 'slab'."""
    atoms = fcc111("Au", size=(2, 2, 3), vacuum=10.0)
    classification = detect_vacuum(atoms)
    assert classification == "slab"


def test_detect_vacuum_molecule():
    """Test that a single molecule is correctly classified as 'molecule'."""
    atoms = molecule("H2O")
    atoms.set_cell([10, 10, 10])
    atoms.set_pbc(True)
    classification = detect_vacuum(atoms)
    assert classification == "molecule"


def test_detect_vacuum_no_pbc():
    """Test that a system with no PBC is classified as 'molecule'."""
    atoms = molecule("H2O")
    atoms.set_pbc(False)
    classification = detect_vacuum(atoms)
    assert classification == "molecule"


def test_detect_vacuum_with_small_voids():
    """Test a porous-like material is still 'bulk' if vacuum is not connected."""
    atoms = bulk("Si", "diamond", a=5.43)
    # This structure has natural voids, but they are not connected
    # across the periodic boundary, so it should be classified as bulk.
    classification = detect_vacuum(atoms)
    assert classification == "bulk"
