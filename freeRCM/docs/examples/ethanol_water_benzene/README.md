# Ethanol-Water-Benzene Ternary System Example

This example demonstrates residue curve mapping for the ethanol-water-benzene system, a classic ternary mixture with azeotropic behavior used in chemical engineering education and industrial applications.

## System Description

**Components:**
- **Ethanol (C₂H₅OH)**: Polar alcohol with hydrogen bonding
- **Water (H₂O)**: Universal solvent with strong hydrogen bonding
- **Benzene (C₆H₆)**: Non-polar aromatic hydrocarbon

**Key Features:**
- **Ethanol-Water Minimum-Boiling Azeotrope**: Forms at ~78°C with composition ~82% ethanol, 18% water
- **Complex Phase Behavior**: Three distillation regions separated by distillation boundaries
- **Industrial Relevance**: Solvent recovery, fuel blending, extractive distillation design

## Thermodynamic Models

- **Vapor Pressure**: Extended Antoine equation (PLXANT)
- **Activity Coefficients**: NRTL model with temperature dependence
- **Critical Properties**: Used for SRK equation of state (if needed)

## Usage

1. **Load the Simulation:**
   ```bash
   # Launch FreeRCM and use "Open Simulation" to load:
   docs/examples/ethanol_water_benzene/simulation.rcm
   ```

2. **Run Analysis:**
   ```bash
   cd docs/examples/ethanol_water_benzene
   python3 analysis.py
   ```

3. **Expected Results:**
   - Residue curves showing distillation boundary
   - Azeotropic point at ~82 mol% ethanol, 18 mol% water
   - Three distinct distillation regions
   - Univolatility line connecting pure benzene to azeotrope

## Educational Value

This example teaches:
- Setting up ternary systems with azeotropes
- Interpreting residue curve maps
- Understanding distillation boundaries
- Practical distillation column design considerations

## References

- Gmehling, J., et al. "Azeotropic Data." VCH Publishers, 1994
- Poling, B.E., et al. "The Properties of Gases and Liquids." McGraw-Hill, 2001
- DECHEMA Chemistry Data Series, Dortmund Data Bank

## Files in This Example

- `simulation.rcm`: Pre-configured FreeRCM simulation file
- `parameters.py`: Thermodynamic parameters and data
- `analysis.py`: Analysis script demonstrating key concepts
- `README.md`: This documentation file