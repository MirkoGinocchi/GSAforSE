# Codes to perform the analyses presented in the paper "Global Sensitivity Analysis of State Estimation for Power Distribution Systems"

MAIN ROUTINE: 

File "MAIN.m", which allows to perform sensitivity analysis (SA) of state estimation (SE) and SE-based applications. 

LIST OF SUBROUTINES contained in "MAIN.m":

1) Sensitivity analysis of SE performed on IEEE 6-bus system: provides Sobol' indexes of the estimates (voltage magnitudes and phase angles)

2) Sensitivity analysis of SE performed on IEEE 15-bus system: provides Sobol' indexes of the estimates (voltage magnitudes and phase angles)

3) Sensitivity analysis of SE performed on IEEE 69-bus system: provides Sobol' indexes of the estimates (voltage magnitudes and phase angles). The user can set the desired DG penetration level.

4) Sensitivity analysis of total power losses, performed on IEEE 69-bus system: provides Sobol' indexes of total losses in the system. 13 different DG penetration levels are considered.

Be always sure to be in the correct path. MATPOWER, UQLab and Parallel Computing Toolbox™ required.

CREDITS:
Riccardo Scalabrin (riccardo.scalabrin@mail.polimi.it), Gianmarco Cocchi (gianmarco.cocchi@mail.polimi.it), Mirko Ginocchi (mirko.ginocchi@eonerc.rwth-aachen.de)

MATPOWER: 
R. D. Zimmerman, C. E. Murillo-Sanchez. MATPOWER User-s Manual, Version 7.1. 2020. [Online]. Available: https://matpower.org/docs/MATPOWER-manual-7.1.pdf // Link to the repository: https://github.com/MATPOWER/matpower

UQLAB: 
UQLab: A Framework for Uncertainty Quantification in MATLAB, Stefano Marelli and Bruno Sudret, In The 2nd International Conference on Vulnerability and Risk Analysis and Management (ICVRAM 2014), University of Liverpool, United Kingdom, July 13-16, 2014, pp. 2554–2563. DOI: https://doi.org/10.1061/9780784413609.257. Instruction for UQLab installations: https://www.uqlab.com/  

The adopted WLS State Estimation algorithm is based on a function originally provided by: Praviraj PG (2022). Power System State Estimation using WLS (https://www.mathworks.com/matlabcentral/fileexchange/23052-power-system-state-estimation-using-wls), MATLAB Central File Exchange.
