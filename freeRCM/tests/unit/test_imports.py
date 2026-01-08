"""
Unit tests for module imports and basic functionality.
"""

import sys
import os
import pytest

# Add src to path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'python'))


def test_core_imports():
    """Test that core modules can be imported."""
    try:
        from core.data_structures import dict2struct
        from core.thermodynamics import orgProps
        from core.solver import RCM
        assert True
    except ImportError as e:
        pytest.fail(f"Failed to import core modules: {e}")


def test_gui_imports():
    """Test that GUI modules can be imported."""
    try:
        from gui.plot_widget import RCMplot
        # Note: main_window import may fail due to PySide6/Qt dependencies
        # but we can test the plot functionality
        assert True
    except ImportError as e:
        pytest.fail(f"Failed to import GUI modules: {e}")


def test_data_structures():
    """Test dict2struct functionality."""
    from core.data_structures import dict2struct

    # Test basic functionality
    d = {'a': 1, 'b': 2}
    ds = dict2struct(d)

    assert ds.a == 1
    assert ds.b == 2
    assert ds['a'] == 1  # Still works as dict


def test_thermodynamics_basic():
    """Test basic thermodynamics functionality."""
    from core.thermodynamics import orgProps
    import numpy as np

    # Create mock data
    comps = np.array(['Component1', 'Component2', 'Component3'])
    selected_comps = np.array(['Component1', 'Component2', 'Component3'])
    allProps = {
        'antoine': np.random.rand(3, 3),
        'NRTL_aij': np.random.rand(3, 3),
        'NRTL_bij': np.random.rand(3, 3),
        'NRTL_cij': np.random.rand(3, 3),
        'PLXANT': np.random.rand(3, 7),
        'TcCel': np.random.rand(3),
        'Pc': np.random.rand(3),
        'omega': np.random.rand(3)
    }

    # Test parameter fetching
    props = orgProps(1, comps, selected_comps, allProps)
    assert 'antoine' in props

    props = orgProps(2, comps, selected_comps, allProps)
    assert 'NRTL_aij' in props