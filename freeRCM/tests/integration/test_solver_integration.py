"""
Integration tests for the solver functionality.
"""

import sys
import os
import numpy as np
import pytest

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'python'))

from core.data_structures import dict2struct
from core.thermodynamics import orgProps


def create_mock_simulation_data():
    """Create mock simulation data for testing."""
    # Mock component data
    comps = np.array(['Ethanol', 'Water', 'Acetone'])
    selected_comps = np.array(['Ethanol', 'Water', 'Acetone'])

    # Mock thermodynamic properties
    allProps = dict2struct({
        'antoine': np.array([
            [8.20417, 1642.89, 230.3],  # Ethanol
            [8.07131, 1730.63, 233.426],  # Water
            [7.632, 1321.0, 192.8]  # Acetone
        ]),
        'NRTL_aij': np.random.rand(3, 3),
        'NRTL_bij': np.random.rand(3, 3),
        'NRTL_cij': np.random.rand(3, 3),
        'PLXANT': np.random.rand(3, 7),
        'TcCel': np.array([243.1, 374.1, 235.0]),
        'Pc': np.array([63.8, 220.6, 47.0]),
        'omega': np.array([0.645, 0.344, 0.307])
    })

    # Mock solver options
    opts = dict2struct({
        'dxi': 0.02,
        'n_it': 50,  # Reduced for testing
        'lines': 5,   # Reduced for testing
        'antMethod': 1,
        'activity': 3,
        'lmopts': {'maxiter': 100, 'ftol': 1e-12, 'xtol': 1e-12}
    })

    return comps, selected_comps, allProps, opts


def test_thermodynamics_data_processing():
    """Test that thermodynamics data is processed correctly."""
    comps, selected_comps, allProps, opts = create_mock_simulation_data()

    # Test Antoine parameter extraction
    antoine_props = orgProps(1, comps, selected_comps, allProps)
    assert 'antoine' in antoine_props
    assert antoine_props.antoine.shape == (3, 3)

    # Test NRTL parameter extraction
    nrtl_props = orgProps(2, comps, selected_comps, allProps)
    assert 'NRTL_aij' in nrtl_props
    assert nrtl_props.NRTL_aij.shape == (3, 3)

    # Test critical property extraction
    crit_props = orgProps(4, comps, selected_comps, allProps)
    assert 'TcCel' in crit_props
    assert len(crit_props.TcCel) == 3


def test_solver_initialization():
    """Test that the solver can be initialized without errors."""
    try:
        from core.solver import RCM
        # Just test that the function exists and can be called
        # (we won't actually run it due to complexity of full simulation)
        assert callable(RCM)
        print("Solver initialization test passed")
    except ImportError as e:
        pytest.skip(f"Solver import failed (expected on some systems): {e}")
    except Exception as e:
        pytest.fail(f"Unexpected error during solver initialization: {e}")


def test_data_structures():
    """Test dict2struct functionality."""
    from core.data_structures import dict2struct

    # Test basic functionality
    data = {'a': 1, 'b': {'nested': True}, 'c': [1, 2, 3]}
    ds = dict2struct(data)

    assert ds.a == 1
    assert ds.b.nested == True
    assert ds.c == [1, 2, 3]

    # Test that it still works as a dict
    assert ds['a'] == 1
    ds['new_key'] = 'value'
    assert ds.new_key == 'value'


def test_module_imports():
    """Test that all modules can be imported."""
    try:
        from core import data_structures, thermodynamics, solver
        from gui import plot_widget
        print("All module imports successful")
    except ImportError as e:
        pytest.fail(f"Module import failed: {e}")


if __name__ == "__main__":
    # Run basic tests
    test_data_structures()
    test_thermodynamics_data_processing()
    test_module_imports()
    test_solver_initialization()
    print("All integration tests passed!")