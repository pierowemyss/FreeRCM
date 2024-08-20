# FreeRCM

FreeRCM is a free, graphical-user-interfaced residue curve mapping tool that is helpful for designing distillation columns, especially for entrainer screening for extractive distillation columns. This specific program also allows for simulations to be saved as "*.rcm" files and opened later.

## How to use:

Download the project and run with Python (shown below for MacOS/Linux) then follow the steps enumerated below.

```bash
#For MacOS/Linux
cd path/to/project/folder
python3 freeRCM.py
```

NOTE: I've been experimenting with building a standalone app using Nuitka; however, I'm still resolving some startup time issues with that so this is the best way for now.



### 1. Create a new or open an existing simulation:

The "Open Simulation" button will prompt a "*.rcm" file to be selected. An example simulation, "example.rcm," has been provided for you to play around with the program.

<p align="center">
<img width="827" alt="StartScreen" src="https://github.com/user-attachments/assets/487fb408-4c9a-444b-b062-a93d824a7744">
</p>

### 2. Set up or modify a simulation:

Add components by clicking the "Add Components" button and delete them by selecting a component and clicking "Delete Component." Make sure three components are in the "Selected Components" box before continuing. Components can be added to or removed from here by selecting a component and clicking ">>" or "<<."

<p align="center">
<img width="798" alt="SimSetUp" src="https://github.com/user-attachments/assets/bd2dbf04-9efb-44ab-b188-7b4fd0517e48">
</p>

Additionally, you must also ensure that for your selected thermodynamic model, all parameters have been input. Click on "Input Parameters" and copy/paste parameters from excel (or enter manually) to do this.

IMPORTANT: All methods deal with temperature in units of degrees Celsius and pressure in units of bar.

<p align="center">
<img width="1103" alt="InpParams" src="https://github.com/user-attachments/assets/941ed499-2ac9-4f01-ac6d-6fcfaea4d149">
</p>

### 3. Create plots, have fun!

Either auto-generate curves or click on the plot to generate a curve. Click on "Save Simulation" to save a "*.rcm" to return to your simulation later.

<p align="center">
<img width="1008" alt="Screenshot 2024-08-19 at 11 03 56 PM" src="https://github.com/user-attachments/assets/75fc6f0c-5b51-42c8-8171-c855e66a066c">
</p>

## How it works:

Residue curve mapping is used to determine, from any starting composition, the only possible path forward and backward in a distillation. From any point, there is only one possible path because residue curves do not intersect. Residue curves are determined from the the differential equation:

$$\frac{dx}{d\xi} = x - y$$

This can be solved forward and backwards in warped time, $\xi$, through the difference equation:

$$x_{i+1} = x_i + \Delta\xi(x_i - y_i)$$

During each step in warped time, $y_i$ is found via vapor-liquid equilibria relations with $x_i$. The result of this is displayed below:

<p align="center">
<img src="https://github.com/user-attachments/assets/e8499168-d793-439f-af49-38e0e3130ce9" alt="firstExample" style="width:50%;">
</p>

NOTE: The NRTL activity coefficient model sometimes comes in different forms. The form used in this program is as follows:

$$\ln{\gamma_i} = \frac{\sum_j x_j\tau_{ji}G_{ji}}{\sum_k x_kG_{ki}} + \sum_j \frac{x_jG_{ij}}{\sum_k x_kG_{{kj}}}\left(\tau_{ij} - \frac{\sum_m x_m\tau_{mj}G_{mj}}{\sum_k x_kG_{kj}}\right)$$

where:
$G_{ij} = \text{exp}(-c_{ij}\tau_{ij})$
and
$\tau_{ij} = a_{ij} + b_{ij}/T$

ADDITIONAL NOTE: The Extended Antoine equation used is as follows:

$$\ln{P_i^{SAT}} = C_{1,i} + \frac{C_{2,i}}{T+C_{3,i}} + C_{4,i}T + C_{5,i}\ln{T} + C_{6,i}T^{C_{7,i}}$$

Other methods (like regular Antoine and SRK) use their most standard forms, so there should not be much ambiguity when inputting parameters. Please remember, as mentioned before, that all methods deal with temperature in units of degrees Celsius and pressure in units of bar.

## Features to be added in the future:

1. Automatic acentricity prediction for SRK using Pitzer correlation.
2. Add useful information to "Help" buttons.
3. Build standalone app (with Nuitka or something). I've already done this, but startup time is horrendous (probably because of MatPlotLib or NumPy, if anyone has any experience with this PLEASE help).
4. Parameter auto-fill from databank?

## References:
- Doherty, M. F., & Malone, M. F. (2001). Conceptual design of Distillation Systems. Boston: McGraw-Hill. 
- [MINPACK](https://github.com/fortran-lang/minpack/tree/main)
