from ase.build import bulk
import numpy as np
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

def test_generate_qe_input():
    atoms = bulk("Si", "diamond", a=5.43)
    parameters = {
        "pseudo_dir": "/path/to/pseudos",
        "ecutwfc": 60.0,
        "kpoints": [4, 4, 4, 0, 0, 0],
        "pseudopotentials": {
            "Si": "Si.upf"
        }
    }

    input_str = generate_qe_input(atoms, parameters)

    assert "calculation = 'scf'" in input_str
    assert "nat = 2" in input_str
    assert "Si" in input_str
    assert "Si.upf" in input_str

def test_parse_qe_output_success():
    mock_output = """
     JOB DONE.
     !    total energy              =     -111.45181134 Ry
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
     atom    2 type  1   force =     0.00000000    0.00000000    0.00000000

     total   stress  (Ry/bohr**3)                       pbar
      -0.00000122   0.00000000  -0.00000000         -0.18      0.00     -0.00
       0.00000000  -0.00000122   0.00000000          0.00     -0.18      0.00
      -0.00000000   0.00000000  -0.00000122         -0.00      0.00     -0.18
    """

    result, error = parse_qe_output(mock_output)

    assert error is None
    assert "total_energy_ev" in result
    assert "forces" in result
    assert "stress" in result
    assert np.isclose(result['total_energy_ev'], -1516.38, atol=1e-2)
    assert len(result['forces']) == 2

def test_parse_qe_output_failure():
    mock_output = "JOB NOT DONE."
    result, error = parse_qe_output(mock_output)
    assert result is None
    assert "did not finish" in error
