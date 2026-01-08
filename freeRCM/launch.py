#!/usr/bin/env python3
"""
FreeRCM - Simple Launcher

A simple launcher that avoids relative import issues.
"""

import sys
import os

# Add the src/python directory to the path
project_root = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(project_root, 'src', 'python')
sys.path.insert(0, src_path)

# Add core and gui paths for direct imports
core_path = os.path.join(src_path, 'core')
gui_path = os.path.join(src_path, 'gui')
if core_path not in sys.path:
    sys.path.insert(0, core_path)
if gui_path not in sys.path:
    sys.path.insert(0, gui_path)

def main():
    """Launch the FreeRCM application."""
    try:
        # Import PySide6
        from PySide6.QtWidgets import QApplication
        print("✓ PyQt6 imported")

        # Import our modules step by step
        import data_structures
        print("✓ data_structures imported")

        import thermodynamics
        print("✓ thermodynamics imported")

        import solver
        print("✓ solver imported")

        import plot_widget
        print("✓ plot_widget imported")

        import main_window
        print("✓ main_window imported")

        # Create and run the application
        app = QApplication(sys.argv)
        window = main_window.GetStartedWindow()
        window.show()
        sys.exit(app.exec())

    except ImportError as e:
        print(f"Import error: {e}")
        import traceback
        traceback.print_exc()
        print("Make sure all dependencies are installed:")
        print("  pip install -r build/requirements/pyReqs.txt")
        sys.exit(1)

if __name__ == "__main__":
    main()