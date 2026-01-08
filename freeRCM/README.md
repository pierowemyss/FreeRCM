# FreeRCM

FreeRCM is a professional, graphical-user-interfaced residue curve mapping tool for chemical engineering distillation design. It specializes in analyzing distillation columns, particularly for entrainer screening in extractive distillation.

## Features

- **Interactive Residue Curve Mapping**: Click anywhere on ternary diagrams to generate distillation trajectories
- **Modern GUI**: Built with PyQt6 for cross-platform compatibility
- **Thermodynamic Models**: NRTL activity coefficients, SRK equations of state, Antoine vapor pressure correlations
- **High-Performance Computing**: C/Fortran numerical solvers for fast, accurate calculations
- **Professional Architecture**: Modular design with proper separation of concerns
- **Cross-Platform**: macOS, Linux, and Windows support
- **Data Persistence**: Save and load simulation files (.rcm format)

## Installation

### Prerequisites

- Python 3.13 (PySide6 compatibility - not available for 3.14+)
- GCC and GFortran compilers
- GNU Scientific Library (GSL)
- MINPACK library

### Quick Start

1. **Clone and setup**:
   ```bash
   git clone <repository-url>
   cd freeRCM
   pip install -r build/requirements/pyReqs.txt
   ```

2. **Build native libraries**:
   ```bash
   cd build
   make
   cd ..
   ```

3. **Run the application**:
   ```bash
   python3.13 launch.py
   ```

## Project Structure

```
freeRCM/
├── src/
│   ├── python/
│   │   ├── core/          # Core functionality
│   │   ├── gui/           # User interface
│   │   └── utils/         # Utilities
│   └── native/            # C/Fortran sources
├── lib/                   # Compiled libraries
├── tests/                 # Test suite
├── build/                 # Build system
├── examples/              # Sample data
├── docs/                  # Documentation
└── packaging/             # Distribution files
```

## Usage

### Basic Workflow

1. **Component Selection**: Choose 3 components for ternary system analysis
2. **Parameter Input**: Enter thermodynamic properties (Antoine coefficients, NRTL parameters, etc.)
3. **Interactive Plotting**: Click on the ternary diagram to generate residue curves
4. **Analysis**: Study distillation trajectories and phase behavior

### Advanced Features

- **Auto-generation**: Batch curve generation for complete mapping
- **Parameter Validation**: Built-in checks for thermodynamic consistency
- **Export Capabilities**: Save plots and data for reports

## Development

### Running Tests

```bash
# Unit tests
python -m pytest tests/unit/

# Integration tests
python tests/integration/test_solver_integration.py

# All tests
python -m pytest tests/
```

### Building from Source

The project uses a Makefile-based build system:

```bash
# Clean build
cd build && make clean && make

# Development build with debug symbols
make CFLAGS="-g -Wall -fPIC"
```

### Code Style

- **Imports**: Standard library → third-party → local modules
- **Naming**: `snake_case` for functions/variables, `PascalCase` for classes
- **Formatting**: 4-space indentation, no trailing whitespace
- **Documentation**: NumPy-style docstrings

## Technical Details

### Architecture

FreeRCM uses a hybrid performance architecture:
- **Python**: GUI (PySide6), data management, high-level logic
- **C/Fortran**: Numerical solving, thermodynamic calculations
- **Library Loading**: Proper dependency resolution with rpath

### Numerical Methods

- **ODE Integration**: Residue curves via `dx/dξ = x - y`
- **Nonlinear Solving**: MINPACK for phase equilibrium calculations
- **Thermodynamic Models**: NRTL, SRK, extended Antoine equations

### Performance

The native C/Fortran backend provides significant performance improvements over pure Python implementations, enabling interactive analysis of complex distillation systems.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## License

See LICENSE file for details.

## References

- Doherty, M. F., & Malone, M. F. (2001). Conceptual design of Distillation Systems. Boston: McGraw-Hill.
- [MINPACK](https://github.com/fortran-lang/minpack)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)