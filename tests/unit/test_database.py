import os
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.db import connect

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult


class TestAseDB(unittest.TestCase):
    def setUp(self):
        self.db_path = "test.db"
        if os.path.exists(self.db_path):
            os.remove(self.db_path)
        self.db = AseDB(self.db_path)
        self.atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])

    def tearDown(self):
        if os.path.exists(self.db_path):
            os.remove(self.db_path)

    def test_add_atoms(self):
        atoms_id = self.db.add_atoms(self.atoms, state="test_state")
        self.assertIsInstance(atoms_id, int)
        with connect(self.db_path) as con:
            row = con.get(id=atoms_id)
        self.assertEqual(row.key_value_pairs["state"], "test_state")

    def test_get_atoms_to_label(self):
        self.db.add_atoms(self.atoms, state="initial")
        self.db.add_atoms(self.atoms, state="labelled")
        atoms_to_label = self.db.get_atoms_to_label()
        self.assertEqual(len(atoms_to_label), 1)
        self.assertEqual(atoms_to_label[0].get_chemical_symbols(), ["H", "H"])

    def test_write_dft_result(self):
        atoms_id = self.db.add_atoms(self.atoms)
        dft_result = DFTResult(
            energy=-1.0,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
            status="success",
            metadata={"test_key": "test_value"},
        )
        self.db.write_dft_result(atoms_id, dft_result)
        with connect(self.db_path) as con:
            row = con.get(id=atoms_id)
        self.assertEqual(row.key_value_pairs["state"], "labelled")
        self.assertIn("dft_result", row.key_value_pairs)
        retrieved_result = DFTResult.model_validate_json(
            row.key_value_pairs["dft_result"]
        )
        self.assertEqual(retrieved_result.energy, -1.0)
        self.assertEqual(retrieved_result.metadata["test_key"], "test_value")

    def test_update_state(self):
        atoms_id = self.db.add_atoms(self.atoms)
        self.db.update_state(atoms_id, "new_state")
        with connect(self.db_path) as con:
            row = con.get(id=atoms_id)
        self.assertEqual(row.key_value_pairs["state"], "new_state")

    def test_get_training_data(self):
        atoms_id = self.db.add_atoms(self.atoms, state="initial")
        dft_result = DFTResult(
            energy=-1.0,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
            status="success",
        )
        self.db.write_dft_result(atoms_id, dft_result)

        atoms_id_2 = self.db.add_atoms(self.atoms, state="initial")
        self.db.update_state(atoms_id_2, "failed")

        atoms_list, results_list = self.db.get_training_data()
        self.assertEqual(len(atoms_list), 1)
        self.assertEqual(len(results_list), 1)
        self.assertEqual(results_list[0].energy, -1.0)


if __name__ == "__main__":
    unittest.main()
