#!/usr/bin/env python3
"""
Thermodynamic parameters for Ethanol-Water-Benzene ternary system.

This module contains all the thermodynamic data needed for residue curve
mapping of the ethanol-water-benzene system, including vapor pressure
coefficients, NRTL interaction parameters, and critical properties.

Data sources:
- Gmehling, J., et al. "Azeotropic Data." VCH Publishers, 1994
- Poling, B.E., et al. "The Properties of Gases and Liquids." McGraw-Hill, 2001
- DECHEMA Chemistry Data Series
"""

import numpy as np

# Component names and identifiers
COMPONENTS = ['Ethanol', 'Water', 'Benzene']
FORMULAS = ['C₂H₅OH', 'H₂O', 'C₆H₆']

# Antoine vapor pressure coefficients
# log10(P[kPa]) = A - B/(T[°C] + C)
ANTOINE_PARAMS = np.array([
    [8.20417, 1642.89, -39.15],    # Ethanol
    [8.07131, 1730.63, -39.57],    # Water
    [6.90565, 1211.033, -52.36]    # Benzene
])

# Extended Antoine (PLXANT) coefficients
# Format: [C1, C2, C3, C4, C5, C6, C7]
# Coefficients designed to give reasonable vapor pressures
PLXANT_PARAMS = np.array([
    [7.242, -1986.0, -2.0, 0.0, 0.0, 0.0, 0.0],            # Ethanol
    [7.074, -1657.0, -1.0, 0.0, 0.0, 0.0, 0.0],            # Water
    [6.905, -1211.0, 0.0, 0.0, 0.0, 0.0, 0.0]              # Benzene
])

# NRTL binary interaction parameters (cal/mol)
# A_ij: Binary interaction energies
NRTL_AIJ = np.array([
    [0.0, -0.174, 0.057],      # Ethanol row
    [-0.174, 0.0, 0.078],      # Water row
    [0.057, 0.078, 0.0]        # Benzene row
])

# B_ij: Temperature dependence coefficients (cal/mol/K)
NRTL_BIJ = np.array([
    [0.0, 53.43, -11.01],      # Ethanol row
    [53.43, 0.0, -415.8],      # Water row
    [-11.01, -415.8, 0.0]      # Benzene row
])

# C_ij: Non-randomness parameters (dimensionless)
NRTL_CIJ = np.array([
    [0.0, 0.3, 0.3],           # Ethanol row
    [0.3, 0.0, 0.3],           # Water row
    [0.3, 0.3, 0.0]            # Benzene row
])

# Critical properties
# Tc: Critical temperature (°C)
TC_CEL = np.array([243.1, 374.0, 289.0])

# Pc: Critical pressure (bar)
PC = np.array([63.8, 221.2, 48.9])

# ω: Acentric factor (dimensionless)
OMEGA = np.array([0.649, 0.344, 0.212])

# Molecular weights (g/mol)
MOLECULAR_WEIGHTS = np.array([46.07, 18.02, 78.11])

# Normal boiling points (°C)
BOILING_POINTS = np.array([78.4, 100.0, 80.1])

# Known azeotrope data
# Ethanol-Water minimum-boiling azeotrope
ETHANOL_WATER_AZEOTROPE = {
    'temperature': 78.2,  # °C
    'pressure': 1.013,    # bar
    'composition': [0.894, 0.106, 0.0],  # mol fractions: Ethanol, Water, Benzene
    'type': 'minimum-boiling'
}

def get_simulation_data():
    """
    Return a dictionary with all simulation data in FreeRCM format.

    Returns:
        dict: Complete simulation data dictionary
    """
    # Import dict2struct for proper object creation
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))
    from core.data_structures import dict2struct

    # Create allProps structure
    all_props = {
        "antoine": ANTOINE_PARAMS.tolist(),
        "PLXANT": PLXANT_PARAMS.tolist(),
        "NRTL_aij": NRTL_AIJ.tolist(),
        "NRTL_bij": NRTL_BIJ.tolist(),
        "NRTL_cij": NRTL_CIJ.tolist(),
        "TcCel": TC_CEL.tolist(),
        "Pc": PC.tolist(),
        "omega": OMEGA.tolist()
    }

    return {
        "P": 1.01325,  # Atmospheric pressure (bar)

        "comps": COMPONENTS,

        "selected_comps": COMPONENTS,

        "allProps": all_props,

        "lmopts": {
            "maxiter": 1000,
            "ftol": 1e-12,
            "xtol": 1e-12
        },

        "opts": dict2struct({
            "antMethod": 2,      # Extended Antoine (PLXANT)
            "activity": 2,       # NRTL
            "lines": 15,
            "linewidth": 1.2,
            "n_it": 250,
            "dxi": 0.02,
            "lmopts": {
                "maxiter": 1000,
                "ftol": 1e-12,
                "xtol": 1e-12
            }
        }),

        # Individual parameters (as expected by load function)
        "NRTL_aij": NRTL_AIJ.tolist(),
        "NRTL_bij": NRTL_BIJ.tolist(),
        "NRTL_cij": NRTL_CIJ.tolist(),
        "TcCel": TC_CEL.tolist(),
        "Pc": PC.tolist(),
        "omega": OMEGA.tolist(),
        "antoine_params": ANTOINE_PARAMS.tolist(),
        "PLXANT_params": PLXANT_PARAMS.tolist(),
    }

def print_system_info():
    """Print information about the ethanol-water-benzene system."""
    print("Ethanol-Water-Benzene Ternary System")
    print("=" * 40)
    print(f"Components: {', '.join(COMPONENTS)}")
    print(f"Formulas: {', '.join(FORMULAS)}")
    print()

    print("Critical Properties:")
    for i, comp in enumerate(COMPONENTS):
        print(f"  {comp}: Tc = {TC_CEL[i]}°C, Pc = {PC[i]} bar, ω = {OMEGA[i]}")
    print()

    print("Azeotrope Information:")
    az = ETHANOL_WATER_AZEOTROPE
    print(f"  Ethanol-Water: {az['temperature']}°C, {az['composition'][0]*100:.1f}% EtOH, {az['composition'][1]*100:.1f}% H₂O")
    print(f"  Type: {az['type']}")
    print()

    print("Thermodynamic Models:")
    print("  Vapor Pressure: Extended Antoine (PLXANT)")
    print("  Activity Coefficients: NRTL with temperature dependence")

if __name__ == "__main__":
    print_system_info()