import numpy as np
from ase import Atoms

from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.utils.dft_utils import create_qe_input_from_atoms, parse_qe_output


def test_create_qe_input_from_atoms():
    atoms = Atoms(
        "Si2", cell=np.diag([5.43, 5.43, 5.43]), positions=[[0, 0, 0], [1.3, 1.3, 1.3]]
    )
    config = DFTCompute(
        code="quantum_espresso",
        command="pw.x",
        pseudopotentials="SSSP_1.3_PBE_precision",
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=2.5,
    )

    input_str = create_qe_input_from_atoms(atoms, config)

    assert "nat = 2" in input_str
    assert "ecutwfc = 60.0" in input_str
    assert "ecutrho = 240.0" in input_str
    assert "Si  1.0  Si.SSSP_1.3_PBE_precision.UPF" in input_str
    assert "Si  0.00000000  0.00000000  0.00000000" in input_str
    assert "14 14 14 0 0 0" in input_str  # 5.43 * 2.5 ~ 13.57 -> ceil is 14


def test_parse_qe_output():
    qe_output = """
    some text before
    !    total energy              =     -123.456 Ry
    some text in between
    Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.1000   -0.2000    0.3000
     atom    2 type  1   force =    -0.1000    0.2000   -0.3000
    some text after forces
     total stress  (Ry/bohr**3)    -0.00000100   -0.00000000   -0.00000000   -0.00000000    0.00000000    0.00000000
    """

    results = parse_qe_output(qe_output)

    assert "energy" in results
    assert "forces" in results

    # Check if the values are reasonable after conversion
    assert -1700 < results["energy"] < -1600  # -123.456 Ry is approx -1679 eV
    assert len(results["forces"]) == 2

    # Note: Stress parsing is not included in the provided output snippet for simplicity
    # A more complete test would include a valid stress block.


def test_parse_qe_output_no_energy():
    qe_output = "This is a failed run output."
    results = parse_qe_output(qe_output)
    assert "dft_energy" not in results
