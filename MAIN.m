%% STATE ESTITMATOR GLOBAL SENSITIVITY ANALYSIS TOOLBOX: MAIN ROUTINE %%
% Software to perform sensitivity analysis (SA) of state estimation (SE) and SE-based
% applications. 
% *********************************************************************** %
%% LIST OF SUBROUTINES:
% 1) Sensitivity analysis of SE performed on IEEE 6-bus system: provides Sobol' indexes
% of the estimates (voltage magnitudes and phase angles)
% *********************************************************************** %

% 2) Sensitivity analysis of SE performed on IEEE 15-bus system: provides Sobol' indexes
% of the estimates (voltage magnitudes and phase angles)
% *********************************************************************** %

% 3) Sensitivity analysis of SE performed on IEEE 69-bus system: provides Sobol' indexes
% of the estimates (voltage magnitudes and phase angles). The user can set
% the desired DG penetration level.
% *********************************************************************** %

% 4) Sensitivity analysis of total loss, performed on IEEE 69-bus system: provides Sobol' indexes
% of total losses inn the system. 13 different DG penetration levels are
% considered.
% *********************************************************************** %

% BE ALWAYS SURE TO BE IN THE CORRECT PATH.
% MATPOWER AND UQLAB REQUIRED. PLACE THEIR FOLDERS IN THE "Mains" FOLDER
% AND MAKE SURE TO HAVE THEM INSTALLED.
% 
% *********************************************************************** %
%% CREDITS:
%% Riccardo Scalabrin - Gianmarco Cocchi - Mirko Ginocchi
%    riccardo.scalabrin@polimi.it
%    gianmarco.cocchi@mail.polimi.it 
%    mirko.ginocchi@eonerc.rwth-aachen.de
%% MATPOWER: 
%   R. D. Zimmerman, C. E. Murillo-Sanchez. MATPOWER User-s Manual, Version 7.1. 2020.
%   [Online]. Available: https://matpower.org/docs/MATPOWER-manual-7.1.pdf
%% UQLAB: 
%   UQLab: A Framework for Uncertainty Quantification in MATLAB,
%   Stefano Marelli and Bruno Sudret,
%   In The 2nd International Conference on Vulnerability and Risk Analysis and Management (ICVRAM 2014),
%   University of Liverpool, United Kingdom, July 13-16, 2014, pp. 2554–2563.
%   DOI: 10.1061/9780784413609.257
%% Our WLS State Estimation algorithm is based on a function originally provided by: 
%   Praviraj PG (2022). Power System State Estimation using WLS 
%   (https://www.mathworks.com/matlabcentral/fileexchange/23052-power-system-state-estimation-using-wls), MATLAB Central File Exchange.
% *********************************************************************** %
%% SUBROUTINE 1: SA of IEEE 6-bus system
clc
clear
run("ieee_6_bus_SA.m")

%% SUBROUTINE 2: SA of IEEE 15-bus system
clc
clear
run("ieee_15_bus_SA.m")

%% SUBROUTINE 3: SA of IEEE 69-bus system
clc
clear
pen_level = 30; % Select the desired DG penetration level [%] (0 % -> passive grid)
run("ieee_69_bus_SA.m")

%% SUBROUTINE 4: SA of total loss - IEEE 69-bus system
clc
clear
run("total_loss_SA_ieee_69_bus.m")
