%% IEEE 6-bus system SA routine
%% GENERAL DESCRIPTION: 
% This script is used to compute the Sensitivity Analysis of the State Estimator
% on IEEE 6-bus system
%
%% ********************************************************************* %%
%% Measurement configuration and meter placement

V_meas_bus = 1:6;    % voltage meters locations
P_meas_bus = [1,2,3]; % active power inj. meters locations
Pf_meas_branch = [1;2;3;4;5;6;7;8;9]; % active power flow meters locations (branches)
Qf_meas_branch = [1;4]; % reactive power flow meters locations (branches)
nvi = length(V_meas_bus); % # of voltage meters and pseudo-meters
npi = length(P_meas_bus); % # of active power inj meters and pseudo-meters
nqi = 0; % # of reactive power inj meters and pseudo-meters
npf = length(Pf_meas_branch); % # of active power flow meters and pseudo-meters
nqf = length(Qf_meas_branch); % # of reactive power flow meters and pseudo-meters

%% Parameters

network_selection = 'case6ww'; %IEEE 6-bus system
mpc = loadcase(network_selection); % load network information
pf = runpf(mpc); % run power flow
Nsamples = 1000; % UA sample size
states_true = [pf.bus(:,9), pf.bus(:,8)]; % nominal values of angles and voltage magnitudes
Sbus = makeSbus(pf.baseMVA, pf.bus, pf.gen); % vector of complex power inj. @ all buses
nbus = length(mpc.bus(:,1)); % # of buses
nbranches = length(mpc.branch(:,1)); % # of branches

%% Inputs uncertainty characterization

voltmeter_class = 1; % voltage meters error [%]
P_inj_error = 2; % active power inj. meters error [%]
PQA_P_class = 2; % active power flow meters error [%]
PQA_Q_class = 2; % reactive power flow meters error [%]
parameters_tol = 0; % tolerance of line parameters

Z_true = [1.1 1.1 1.098 1.018 1.006 1.034 1.325 1.625 0.585 0.231 0.584 0.509 0.205 0.746 0.352 0.529 0.256 0.547 -0.131 -0.064]; % true values of measurements

% std. devs. values for inputs statistical characterization & UA
sigma_voltmeter = voltmeter_class/100.*Z_true(1:nvi)'; % voltage meters std. dev.
sigma_Pinj = P_inj_error/100.*Z_true(nvi+1:nvi+npi)'; % active power inj. meters std. dev.
sigma_Pf = PQA_P_class/100.*Z_true(nvi+npi+1:nvi+npi+npf)'; % active power flow meters std. dev.
sigma_Qf = PQA_Q_class/100.*Z_true(nvi+npi+npf+1:nvi+npi+npf+nqf)'; % reactive power flow meters std. dev.

Z_true_sigma = abs([sigma_voltmeter; sigma_Pinj; sigma_Pf; sigma_Qf]);

%----------| R | X | B |------------%
RXB_true = [mpc.branch(:,3), mpc.branch(:,4), mpc.branch(:,5)]; % true values of line parameters

%----------| R | X | B |------------%
RXB_sigma = abs(parameters_tol*RXB_true);  % line parameters std dev

%% WLS SE parameters

Parameters.Network = network_selection;
Parameters.Nsamples = Nsamples;
Parameters.Nbus = nbus;
Parameters.fbus = mpc.branch(:,1);
Parameters.tbus = mpc.branch(:,2);
Parameters.Nbranches = nbranches;
Parameters.V_meas_bus = V_meas_bus;
Parameters.P_meas_bus = P_meas_bus;
Parameters.PQf_meas_branch = [Pf_meas_branch; Qf_meas_branch];
Parameters.nvi = nvi;
Parameters.npi = npi;
Parameters.nqi = nqi; 
Parameters.npf = npf;
Parameters.nqf = nqf;
Parameters.sigma_meas = Z_true_sigma;

%% Uncertainty Analysis (UA)

uqlab;
modelopts.mFile = 'WLS_ieee_6';
modelopts.Parameters = Parameters; 
model = uq_createModel(modelopts);

    for i = 1:nbranches % branch resistances distributions 
        PDFs.Marginals(i).Name = 'R';
        PDFs.Marginals(i).Type = 'Gaussian';
        PDFs.Marginals(i).Parameters = [RXB_true(i,1), RXB_sigma(i,1)];
    end    

    for i = 1:nbranches % branch reactances distributions 
        PDFs.Marginals(i + nbranches).Name = 'X';
        PDFs.Marginals(i + nbranches).Type = 'Gaussian';
        PDFs.Marginals(i + nbranches).Parameters = [RXB_true(i,2), RXB_sigma(i,2)];
    end    

    for i = 1:nbranches % branch suscteptances distributions 
        PDFs.Marginals(i + 2*nbranches).Name = 'B';
        PDFs.Marginals(i + 2*nbranches).Type = 'Gaussian';
        PDFs.Marginals(i + 2*nbranches).Parameters = [RXB_true(i,3), RXB_sigma(i,3)];
    end
        
    for i = 1:nvi % voltage meters distributions
        PDFs.Marginals(3*nbranches + i).Name = 'V_m';
        PDFs.Marginals(3*nbranches + i).Type = 'Gaussian';
        PDFs.Marginals(3*nbranches + i).Parameters = [Z_true(i), Z_true_sigma(i)];
    end 

    for i = 1:npi % P_inj distributions 
        PDFs.Marginals(3*nbranches + nvi + i).Name = 'P_in';
        PDFs.Marginals(3*nbranches + nvi + i).Type = 'Gaussian';
        PDFs.Marginals(3*nbranches + nvi + i).Parameters = [Z_true(i + nvi), Z_true_sigma(i + nvi)];
    end  

    for i = 1:npf % P_f distributions
        PDFs.Marginals(3*nbranches + nvi + npi + nqi + i).Name = 'P_f';
        PDFs.Marginals(3*nbranches + nvi + npi + nqi + i).Type = 'Gaussian';
        PDFs.Marginals(3*nbranches + nvi + npi + nqi + i).Parameters = [Z_true(i + nvi + npi + nqi), Z_true_sigma(i + nvi + npi + nqi)];
    end 

    for i = 1:nqf % Q_f distributions
        PDFs.Marginals(3*nbranches + nvi + npi + npf + nqi + i).Name = 'Q_f';
        PDFs.Marginals(3*nbranches + nvi + npi + npf + nqi + i).Type = 'Gaussian';
        PDFs.Marginals(3*nbranches + nvi + npi + npf + nqi + i).Parameters = [Z_true(i + nvi + npi + nqi + npf), Z_true_sigma(i + nvi + npi + nqi + npf)];
    end  

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


