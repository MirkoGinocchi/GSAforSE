%% IEEE 15-bus system SA routine
%% GENERAL DESCRIPTION: 
% This script is used to compute the Sensitivity Analysis of the State Estimator
% on IEEE 15-bus system
%
%% ********************************************************************* %%
%% Measurement configuration and meter placement

V_meas_bus = [1]; % voltage meters locations
P_meas_bus = 2:15; % active power inj. meters and pseudo-measurements locations
Q_meas_bus = 2:15; % reactive power inj. meters and pseudos-measurements locations
nvi = length(V_meas_bus); % # of voltage meters and pseudo-meters
npi = length(P_meas_bus); % # of active power inj meters and pseudo-meters
nqi = length(Q_meas_bus); % # of reactive power inj meters and pseudo-meters

%% Parameters
network_selection = 'case15da';
mpc = loadcase(network_selection); % load network information
pf = runpf(mpc); % run power flow
Nsamples = 1000; % UA sample size
states_true = [pf.bus(:,9), pf.bus(:,8)]; % nominal values of angles and voltage magnitudes
Sbus = makeSbus(pf.baseMVA, pf.bus, pf.gen); % vector of complex power inj. @ all buses
nbus = length(mpc.bus(:,1)); % # of buses
nbranches = length(mpc.branch(:,1)); % # of branches

%% Inputs uncertainty characterization
    
voltmeter_class = 1; % voltage meters error [%]
PM_error = 50; % Pseudo Measurement Error [%]
parameters_tol = 0.1; % tolerance of line parameters

% --------| V_mag | P_inj | Qinj |--------
Z_true = [pf.bus(V_meas_bus,8); real(Sbus(P_meas_bus)); imag(Sbus(Q_meas_bus))]; % true values of measurements

% std. devs. values for inputs statistical characterization & UA
sigma_voltmeter = voltmeter_class/100.*Z_true(1:nvi,1);
sigma_Pinj = PM_error/3/100.*abs(real(Sbus(P_meas_bus)));
sigma_Qinj = PM_error/3/100.*abs(imag(Sbus(P_meas_bus)));

Z_true_sigma = abs([sigma_voltmeter; sigma_Pinj; sigma_Qinj]); 

%----------| R | X |------------%
RXB_true = [mpc.branch(:,3), mpc.branch(:,4)]; % true values of line parameters

%----------| R | X |------------%
RXB_tol = abs(parameters_tol*RXB_true);  % line parameters std dev

%% WLS SE parameters

Parameters.Network = network_selection;
Parameters.Nsamples = Nsamples;
Parameters.Nbus = nbus;
Parameters.fbus = mpc.branch(:,1);
Parameters.tbus = mpc.branch(:,2);
Parameters.Nbranches = nbranches;
Parameters.V_meas_bus = V_meas_bus;
Parameters.P_meas_bus = P_meas_bus;
Parameters.Q_meas_bus = Q_meas_bus;
Parameters.PQf_meas_branch = 0;
Parameters.nvi = nvi;
Parameters.npi = npi;
Parameters.nqi = nqi;
Parameters.npf = 0;
Parameters.nqf = 0;
Parameters.sigma_meas = Z_true_sigma;
Parameters.states_true_phasors = states_true(:,2).*exp(1j.*(states_true(:,1)/180*pi)); %[states_true(:,1)*pi/180; states_true(:,2)];
   
    
%% Uncertainty Analysis (UA)

uqlab;
modelopts.mFile = 'WLS_ieee_15';
modelopts.Parameters = Parameters; 
model = uq_createModel(modelopts);

    for i = 1:nbranches % branch resistances distributions 
        PDFs.Marginals(i).Name = 'R';
        PDFs.Marginals(i).Type = 'Uniform';
        PDFs.Marginals(i).Parameters = [RXB_true(i,1) - RXB_tol(i,1), RXB_true(i,1) + RXB_tol(i,1)];
    end    

    for i = 1:nbranches % branch reactances distributions 
        PDFs.Marginals(i + nbranches).Name = 'X';
        PDFs.Marginals(i + nbranches).Type = 'Uniform';
        PDFs.Marginals(i + nbranches).Parameters = [RXB_true(i,2) - RXB_tol(i,2), RXB_true(i,2) + RXB_tol(i,2)];
    end    
            
    for i = 1:nvi
        PDFs.Marginals(2*nbranches + i).Name = 'V_m';
        PDFs.Marginals(2*nbranches + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + i).Parameters = [Z_true(i), Z_true_sigma(i)];
    end % voltage meters distributions

    for i = 1:npi
        PDFs.Marginals(2*nbranches + nvi + i).Name = 'P_in';
        PDFs.Marginals(2*nbranches + nvi + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + nvi + i).Parameters = [Z_true(i + nvi), Z_true_sigma(i + nvi)];
    end % P_inj pseuso-m. distributions  

    for i = 1:nqi
        PDFs.Marginals(2*nbranches + nvi + npi + i).Name = 'Q_in';
        PDFs.Marginals(2*nbranches + nvi + npi + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + nvi + npi + i).Parameters = [Z_true(i + nvi + npi), Z_true_sigma(i + nvi + npi)];
    end % Q_inj pseuso-m. distributions

PDFs.Copula.Type = 'Independent'; % independent inputs (NO CORRELATION)   
Z = uq_createInput(PDFs);
X = uq_getSample(Z, Nsamples); % input space MC sampling 
estimates = uq_evalModel(X); % outputs distribution matrix

%% Sensitivity Analysis
% PCE SURROGATE MODEL

clear model modelopts; % clears the previous model 

PCEOpts.Type = 'Metamodel'; % select the metamodel tool
PCEOpts.MetaType = 'PCE';
PCEOpts.Input = Z; % inputs statistical properties
PCEOpts.ExpDesign.X = X; % input samples
PCEOpts.ExpDesign.Y = estimates; % output samples
PCEOpts.Method = 'LARS';
PCEOpts.TruncOptions.qNorm = 0.7; % q-norm hyperbolic truncation parameter
PCEOpts.Degree = 3; % polynomials degree
myPCE = uq_createModel(PCEOpts); % PCE surrogate model

%% GSA: Sobol' variance-based sensitivity analysis

PCESobol.Type = 'Sensitivity';
PCESobol.Method = 'Sobol';
PCESobol.Sobol.Order = 2; % maximum order of Sobol' indexes
PCESobolAnalysis = uq_createAnalysis(PCESobol);


