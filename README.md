# Codes to perform the analyses presented in the paper "Global Sensitivity Analysis of State Estimation for Power Distribution Systems"

MAIN ROUTINE: 
MAIN.m, which allows to perform sensitivity analysis (SA) of state estimation (SE) and SE-based applications. 

LIST OF SUBROUTINES:

1) Sensitivity analysis of SE performed on IEEE 6-bus system: provides Sobol' indexes of the estimates (voltage magnitudes and phase angles)

2) Sensitivity analysis of SE performed on IEEE 15-bus system: provides Sobol' indexes of the estimates (voltage magnitudes and phase angles)

3) Sensitivity analysis of SE performed on IEEE 69-bus system: provides Sobol' indexes of the estimates (voltage magnitudes and phase angles). The user can set the desired DG penetration level.

4) Sensitivity analysis of total loss, performed on IEEE 69-bus system: provides Sobol' indexes of total losses inn the system. 13 different DG penetration levels are considered.

BE ALWAYS SURE TO BE IN THE CORRECT PATH.
MATPOWER, UQLAB AND PARALLEL COMPUTING TOOLBOX™ REQUIRED. PLACE THEIR FOLDERS IN THE "Mains" FOLDER AND MAKE SURE TO HAVE THEM INSTALLED.

CREDITS:
Riccardo Scalabrin and Gianmarco Cocchi
riccardo.scalabrin@mail.polimi.it
gianmarco.cocchi@mail.polimi.it
Mirko Ginocchi 
mirko.ginocchi@eonerc.rwth-aachen.de

MATPOWER: 
R. D. Zimmerman, C. E. Murillo-Sanchez. MATPOWER User-s Manual, Version 7.1. 2020. [Online]. Available: https://matpower.org/docs/MATPOWER-manual-7.1.pdf // Link to the repository: https://github.com/MATPOWER/matpower

UQLAB: 
UQLab: A Framework for Uncertainty Quantification in MATLAB, Stefano Marelli and Bruno Sudret, In The 2nd International Conference on Vulnerability and Risk Analysis and Management (ICVRAM 2014), University of Liverpool, United Kingdom, July 13-16, 2014, pp. 2554–2563. DOI: 10.1061/9780784413609.257

The adopted WLS State Estimation algorithm is based on a function originally provided by: Praviraj PG (2022). Power System State Estimation using WLS (https://www.mathworks.com/matlabcentral/fileexchange/23052-power-system-state-estimation-using-wls), MATLAB Central File Exchange.
