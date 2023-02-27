%% GSA study of total loss - IEEE 69-bus system
%% GENERAL DESCRIPTION: 
% This script is used to compute the Sensitivity Analysis of the total
% losses on IEEE 69-bus system. 13 different DG penetration levels are
% considered, starting from a passive grid (0% DG penetration) and increasing the DG
% penetration by 10% for each level.
%
%% ********************************************************************* %%
%% Parameters

network_selection = 'case69';

for time_sector = 1:13
    
mpc = loadcase(network_selection); % load network information
pf = runpf(mpc); % run power flow with initial loads (necessary to get P_load)
Sbus = makeSbus(pf.baseMVA, pf.bus, pf.gen); % vector of complex power inj. @ all buses
P_meas_bus = [find(mpc.bus(:,3) ~= 0)];
sigma_loads = 50/3/100.*real(Sbus(P_meas_bus)); % loads active power std. dev. 

%% Distributed Generators - placement and settings

Pd_tot = sum(mpc.bus(:,3)); % Total active power required by all the loads
dg_penetration = [0:10:120]./100; % DG penetration level
Pg_to_install = dg_penetration(time_sector)*Pd_tot; % Total active power injected by DGs
dg_position(:,1) = [7;18;24;28;33;39;45;48]; % DGs locations (buses)
% dg_position(:,2) = [0.3;0.15;0.15;0.1;0.1;0.1;0.05;0.05]; % Total active P. inj. by DGs is shared among the single units 
dg_position(:,2) = ones(1,8)*0.125; % Total active P. inj. by DGs is shared among the single units 
mpc.bus(dg_position(1,1),3) = mpc.bus(dg_position(1,1),3) - dg_position(1,2)*Pg_to_install; % DGs provide power to loads connected where DGs are located, 
mpc.bus(dg_position(2,1),3) = mpc.bus(dg_position(2,1),3) - dg_position(2,2)*Pg_to_install; % decreasing their demand
mpc.bus(dg_position(3,1),3) = mpc.bus(dg_position(3,1),3) - dg_position(3,2)*Pg_to_install;
mpc.bus(dg_position(4,1),3) = mpc.bus(dg_position(4,1),3) - dg_position(4,2)*Pg_to_install;
mpc.bus(dg_position(5,1),3) = mpc.bus(dg_position(5,1),3) - dg_position(5,2)*Pg_to_install;
mpc.bus(dg_position(6,1),3) = mpc.bus(dg_position(6,1),3) - dg_position(6,2)*Pg_to_install;
mpc.bus(dg_position(7,1),3) = mpc.bus(dg_position(7,1),3) - dg_position(7,2)*Pg_to_install;
mpc.bus(dg_position(8,1),3) = mpc.bus(dg_position(8,1),3) - dg_position(8,2)*Pg_to_install;
sigma_gen = 50/3/100.*dg_position(:,2)*Pg_to_install; % gens. active power std. dev.

pf = runpf(mpc); % run new power flow with new load values
Nsamples = 2000; % UA sample size
states_true = [pf.bus(:,9), pf.bus(:,8)]; % nominal values of angles and voltage modules 
nbus = length(mpc.bus(:,1)); % # of buses
nbranches = length(mpc.branch(:,1)); % # of branches

%% Measurement configuration and meter placement

V_meas_bus = [1]'; % voltage meters locations
P_meas_bus = [find(mpc.bus(:,3) ~= 0)]; % active power inj. meters and pseudo-measurements locations
Q_meas_bus = [find(mpc.bus(:,4) ~= 0)]; % reactive power inj. meters and pseudo-measurements locations
PQf_meas_branch = [1]'; % active power flow meters and pseudo-measurements locations
PQ_transit = find(abs(Sbus) == 0); % transit buses
nvi = length(V_meas_bus); % # of voltage meters and pseudo-meters
npi = length(P_meas_bus); % # of active power inj meters and pseudo-meters
nqi = length(Q_meas_bus); % # of reactive power inj meters and pseudo-meters
npf = length(PQf_meas_branch); % # of active power flow meters and pseudo-meters
nqf = length(PQf_meas_branch); % # of reactive power flow meters and pseudo-meters
npt = length(PQ_transit); % # of transit buses
nqt = length(PQ_transit); % # of transit buses

%% Inputs uncertainty characterization

voltmeter_class = 1; % voltage meters error [%]
PM_error = 50/3; % Pseudo Measurement Error [%]
PQA_P_class = 1; % PQA active power error [%]
PQA_Q_class = 1; % PQA reactive power error [%]
parameters_tol = 0.1; % tolerance of line parameters

%--------| V_mag | P_inj | Qinj | P_transit | Q_transit | P_from | Q_from |--------
Z_true = [pf.bus(V_meas_bus,8); real(Sbus(P_meas_bus)); imag(Sbus(Q_meas_bus)); real(Sbus(PQ_transit)); imag(Sbus(PQ_transit));...
            pf.branch(PQf_meas_branch,14)./pf.baseMVA; pf.branch(PQf_meas_branch,15)./pf.baseMVA]; % true values of measurements

% std. devs. values for inputs statistical characterization & UA
sigma_voltmeter = voltmeter_class/100.*Z_true(1:nvi,1);
sigma_Pinj = PM_error/100.*real(Sbus(P_meas_bus));

for ii = 1:length(P_meas_bus) % active power uncertainty combination @ DG buses
    for kk = 1:length(dg_position(:,1))
        if P_meas_bus(ii) == dg_position(kk,1)
            sigma_Pinj(ii) =  sqrt(sigma_gen(kk)^2 + sigma_loads(ii)^2);
        end
    end
end

sigma_Qinj = PM_error/100.*imag(Sbus(P_meas_bus));
sigma_Pf = PQA_P_class/100.*pf.branch(PQf_meas_branch,14)./pf.baseMVA;
sigma_Qf = PQA_Q_class/100.*pf.branch(PQf_meas_branch,15)./pf.baseMVA;
sigma_P_transit = 1e-4*ones(npt,1);
sigma_Q_transit = 1e-4*ones(nqt,1);

Z_true_sigma = abs([sigma_voltmeter; sigma_Pinj; sigma_Qinj; sigma_P_transit; sigma_Q_transit; sigma_Pf; sigma_Qf]);

% std. devs. values for WLS weights
sigma_voltmeter_wls = voltmeter_class/100.*ones(nvi,1);
sigma_Pinj_wls = PM_error/100.*ones(npi,1);
sigma_Qinj_wls = PM_error/100.*ones(nqi,1);
sigma_Pf_wls = PQA_P_class/100.*ones(npf,1);
sigma_Qf_wls = PQA_Q_class/100.*ones(nqf,1);
Z_true_sigma_WLS = abs([sigma_voltmeter_wls; sigma_Pinj_wls; sigma_Qinj_wls; sigma_P_transit; sigma_Q_transit; sigma_Pf_wls; sigma_Qf_wls]);

%----------| R | X |------------%
RXB_true = [mpc.branch(:,3), mpc.branch(:,4)]; % true values of line parameters

%----------| R | X |------------%
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
Parameters.Q_meas_bus = Q_meas_bus;
Parameters.PQf_meas_branch = PQf_meas_branch;
Parameters.PQ_transit = PQ_transit;
Parameters.nvi = nvi;
Parameters.npi = npi;
Parameters.nqi = nqi;
Parameters.npf = npf;
Parameters.nqf = nqf;
Parameters.npt = npt;
Parameters.nqt = nqt;
Parameters.sigma_meas = Z_true_sigma_WLS;
Parameters.states_true_phasors = states_true(:,2).*exp(1j.*(states_true(:,1)/180*pi));
% Parameters.states_true_phasors = [states_true(:,1)*pi/180, states_true(:,2)]; %states_true(:,2).*exp(1j.*(states_true(:,1)/180*pi));

%% Uncertainty Analysis (UA)

uqlab;
modelopts.mFile = 'WLS_plus_losses_ieee_69';
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
    
% _________________________________________________________________________
        
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

    for i = 1:npt
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + i).Name = 'P_t';
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + i).Parameters = [Z_true(i + nvi + npi + nqi), 0*Z_true_sigma(i + nvi + npi + nqi)];
    end % virtual measurements (P) distributions
    
    for i = 1:nqt
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + npt + i).Name = 'Q_t';
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + npt + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + npt + i).Parameters = [Z_true(i + nvi + npi + nqi + npt), 0*Z_true_sigma(i + nvi + npi + nqi + npt)];
    end % virtual measurements (Q) distributions        
    
    for i = 1:npf
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + npt + nqt + i).Name = 'P_f';
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + npt + nqt + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + nvi + npi + nqi + npt + nqt + i).Parameters = [Z_true(i + nvi + npi + nqi + npt + nqt), Z_true_sigma(i + nvi + npi + nqi + npt + nqt)];
    end % P_f pseuso-m. distributions

    for i = 1:nqf
        PDFs.Marginals(2*nbranches + nvi + npi + npf + nqi + npt + nqt + i).Name = 'Q_f';
        PDFs.Marginals(2*nbranches + nvi + npi + npf + nqi + npt + nqt + i).Type = 'Gaussian';
        PDFs.Marginals(2*nbranches + nvi + npi + npf + nqi + npt + nqt + i).Parameters = [Z_true(i + nvi + npi + nqi + npf + npt + nqt), Z_true_sigma(i + nvi + npi + nqi + npf + npt + nqt)];
    end % Q_f pseuso-m. distributions 

PDFs.Copula.Type = 'Independent';    
Z = uq_createInput(PDFs);
X = uq_getSample(Z, Nsamples); % input space MC sampling 
[Y, mean_I_flowing(time_sector,:)] = uq_evalModel(X);
loss_dist(:,time_sector) = transpose(Y);
mean_loss(time_sector) = mean(Y);

%% Sensitivity Analysis
% PCE SURROGATE MODEL

clear model modelopts; % clears the previous model 

PCEOpts.Type = 'Metamodel'; % select the metamodel tool
PCEOpts.MetaType = 'PCE';
PCEOpts.Input = Z; % inputs statistical properties
PCEOpts.ExpDesign.X = X; % input samples
PCEOpts.ExpDesign.Y = transpose(Y); % output samples
PCEOpts.Method = 'LARS';
PCEOpts.TruncOptions.qNorm = 0.7; % q-norm hyperbolic truncation parameter
PCEOpts.Degree = 3; % polynomials degree
myPCE = uq_createModel(PCEOpts); % PCE surrogate model

%% GSA: Sobol' variance-based sensitivity analysis

PCESobol.Type = 'Sensitivity';
PCESobol.Method = 'Sobol';
PCESobol.Sobol.Order = 2; % maximum order of Sobol' indexes
PCESobolAnalysis = uq_createAnalysis(PCESobol);

total(:,time_sector) = PCESobolAnalysis.Results.Total;
first(:,time_sector) = PCESobolAnalysis.Results.FirstOrder;
second(:,time_sector) = PCESobolAnalysis.Results.AllOrders{1, 2};
second_ind = PCESobolAnalysis.Results.VarIdx{2, 1};
error(time_sector) = myPCE.Error.ModifiedLOO;
second(:,14:15) = second_ind;
clear myPCE PCEOpts PCESobolAnalysis PCESobol mpc pf;

end


