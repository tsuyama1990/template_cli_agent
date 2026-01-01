import pytest
import os
from tempfile import TemporaryDirectory
from mlip_autopipec.database import AseDBWrapper
from ase.build import molecule

@pytest.fixture(scope="function")
def temp_db():
    """Provides a temporary database for a single test function."""
    with TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test.db")
        db_wrapper = AseDBWrapper(db_path)

        # Add a sample atom to the database for tests that need it
        h2o = molecule("H2O")
        h2o.set_cell([10, 10, 10])
        h2o.center()
        db_wrapper.add_atoms(h2o)

        yield db_wrapper
