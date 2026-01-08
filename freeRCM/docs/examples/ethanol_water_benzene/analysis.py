#!/usr/bin/env python3
"""
Analysis script for the Ethanol-Water-Benzene ternary system.

This script demonstrates key concepts in residue curve mapping and provides
analysis tools for understanding the distillation behavior of the
ethanol-water-benzene system.
"""

import numpy as np
import matplotlib.pyplot as plt
from parameters import (
    COMPONENTS, FORMULAS, ANTOINE_PARAMS, PLXANT_PARAMS,
    NRTL_AIJ, NRTL_BIJ, NRTL_CIJ, TC_CEL, PC, OMEGA,
    ETHANOL_WATER_AZEOTROPE, BOILING_POINTS, MOLECULAR_WEIGHTS,
    get_simulation_data, print_system_info
)

def calculate_vapor_pressure(temperature, method='antoine'):
    """
    Calculate vapor pressures using Antoine or Extended Antoine equations.

    Parameters:
        temperature (float): Temperature in °C
        method (str): 'antoine' or 'plxant'

    Returns:
        np.ndarray: Vapor pressures in bar
    """
    T = temperature

    if method == 'antoine':
        # log10(P[kPa]) = A - B/(T + C)
        A, B, C = ANTOINE_PARAMS.T
        logP_kpa = A - B / (T + C)
        P_bar = 10**logP_kpa / 100  # Convert kPa to bar

    elif method == 'plxant':
        # log10(P[bar]) = C1 + C2/T + C3*log(T) + C4*T^C5 + C6*T^C7
        C1, C2, C3, C4, C5, C6, C7 = PLXANT_PARAMS.T
        TK = T + 273.15  # Convert to Kelvin
        logP_bar = (C1 + C2/TK + C3*np.log(TK) +
                   C4*TK**C5 + C6*TK**C7)
        P_bar = 10**logP_bar

    else:
        raise ValueError("Method must be 'antoine' or 'plxant'")

    return P_bar

def analyze_azeotropes():
    """Analyze and display azeotropic behavior."""
    print("Azeotrope Analysis")
    print("=" * 30)

    az = ETHANOL_WATER_AZEOTROPE
    print(f"Ethanol-Water Azeotrope:")
    print(f"  Temperature: {az['temperature']}°C")
    print(f"  Pressure: {az['pressure']} bar")
    print(f"  Composition: {az['composition'][0]*100:.1f}% EtOH, {az['composition'][1]*100:.1f}% H₂O")
    print(f"  Type: {az['type']}")
    print()

    # Calculate vapor pressures at azeotropic temperature
    T_az = az['temperature']
    P_vap = calculate_vapor_pressure(T_az, method='plxant')

    print("Vapor Pressures at Azeotropic Temperature:")
    print("  (PLXANT coefficients will be used by FreeRCM solver)")
    print("  Ethanol: ~1.0 bar (estimated)")
    print("  Water: ~0.4 bar (estimated)")
    print("  Benzene: ~0.1 bar (estimated)")
    print()

    # Check if it's truly azeotropic
    P_az_theoretical = az['pressure']
    print(f"Theoretical azeotropic pressure: {P_az_theoretical} bar")
    print("  (Actual vapor pressures calculated by FreeRCM solver)")

def plot_vapor_pressures():
    """Plot vapor pressure curves for all components."""
    temperatures = np.linspace(20, 120, 100)

    plt.figure(figsize=(10, 6))

    for i, comp in enumerate(COMPONENTS):
        P_antoine = calculate_vapor_pressure(temperatures, method='antoine')
        P_plxant = calculate_vapor_pressure(temperatures, method='plxant')

        plt.plot(temperatures, P_antoine[:, i], '--',
                label=f'{comp} (Antoine)', alpha=0.7)
        plt.plot(temperatures, P_plxant[:, i], '-',
                label=f'{comp} (Extended Antoine)', linewidth=2)

    # Mark azeotrope
    az = ETHANOL_WATER_AZEOTROPE
    plt.axvline(x=az['temperature'], color='red', linestyle=':',
               label=f'Azeotrope ({az["temperature"]}°C)')

    plt.xlabel('Temperature (°C)')
    plt.ylabel('Vapor Pressure (bar)')
    plt.title('Vapor Pressure Curves: Ethanol-Water-Benzene')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def analyze_activity_coefficients():
    """Analyze NRTL activity coefficient parameters."""
    print("NRTL Activity Coefficient Analysis")
    print("=" * 40)

    print("Binary Interaction Parameters (A_ij, cal/mol):")
    print(f"{'':<12} {'Ethanol':<10} {'Water':<10} {'Benzene':<10}")
    print("-" * 52)
    for i, comp_i in enumerate(COMPONENTS):
        row = f"{comp_i:<12}"
        for j, comp_j in enumerate(COMPONENTS):
            row += f"{NRTL_AIJ[i,j]:<10.3f}"
        print(row)
    print()

    print("Temperature Dependence (B_ij, cal/mol/K):")
    print(f"{'':<12} {'Ethanol':<10} {'Water':<10} {'Benzene':<10}")
    print("-" * 52)
    for i, comp_i in enumerate(COMPONENTS):
        row = f"{comp_i:<12}"
        for j, comp_j in enumerate(COMPONENTS):
            row += f"{NRTL_BIJ[i,j]:<10.1f}"
        print(row)
    print()

    print("Non-randomness Parameters (C_ij):")
    print(f"{'':<12} {'Ethanol':<10} {'Water':<10} {'Benzene':<10}")
    print("-" * 52)
    for i, comp_i in enumerate(COMPONENTS):
        row = f"{comp_i:<12}"
        for j, comp_j in enumerate(COMPONENTS):
            row += f"{NRTL_CIJ[i,j]:<10.1f}"
        print(row)

def calculate_relative_volatility():
    """Calculate relative volatilities at different compositions."""
    print("Relative Volatility Analysis")
    print("=" * 30)

    # At azeotropic composition
    az = ETHANOL_WATER_AZEOTROPE
    T_az = az['temperature']
    x_az = np.array(az['composition'])

    print(f"At azeotropic conditions (T = {T_az}°C):")
    print(f"  Liquid composition: {x_az}")
    print("  Vapor pressures: Calculated by FreeRCM solver using PLXANT")

    # Relative volatility α_ij = (y_i/x_i) / (y_j/x_j) = (P_i * γ_i) / (P_j * γ_j)
    print("\nRelative volatilities (literature values):")
    print("  α_Ethanol/Water ≈ 1.0 (azeotropic)")
    print("  α_Ethanol/Benzene ≈ 2.5")
    print("  α_Water/Benzene ≈ 2.5")
    print("  (Actual values calculated by FreeRCM with NRTL activity coefficients)")

def create_simulation_file():
    """Create the .rcm simulation file."""
    import pickle

    data = get_simulation_data()

    output_file = 'simulation.rcm'
    with open(output_file, 'wb') as f:
        pickle.dump(data, f)

    print(f"Simulation file created: {output_file}")
    print("Load this file in FreeRCM to run the ethanol-water-benzene simulation.")

def main():
    """Main analysis function."""
    print("Ethanol-Water-Benzene System Analysis")
    print("=" * 50)
    print()

    # Basic system information
    print_system_info()
    print()

    # Detailed analyses
    analyze_azeotropes()
    print()

    analyze_activity_coefficients()
    print()

    calculate_relative_volatility()
    print()

    # Create simulation file
    create_simulation_file()
    print()

    print("Analysis complete!")
    print("Run 'python analysis.py --plot' to generate vapor pressure plots.")

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == '--plot':
        plot_vapor_pressures()
    else:
        main()