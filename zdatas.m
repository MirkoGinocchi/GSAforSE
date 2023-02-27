%% HELP ZDATAS
% This function provides the correct measurement data arrangement to
% solve the WLS State Estimation. To add other networks, just add a new
% case structure. Keep in mind that the slack bus must be labelled as bus
% #1, so please verify that the bus numbering is correct.
%% ********************************************************************* %%

function zdt = zdatas(Z_real, sigma_meas, nvi, npi, nqi, npt, nqt, npf, nqf,...
                                V_index, P_index, Q_index, PQf_index, PQ_transit, bus_indices_branch)

nbranches = length(bus_indices_branch(:,1));

switch nbranches

    case 68 % IEEE 69-bus system 
        meas_index = 1:length(Z_real);
        meas_type = [ones(1,nvi), 2*ones(1,npi), 2*ones(1,npt), 3*ones(1,nqi), 3*ones(1,nqt), 4*ones(1,npf), 5*ones(1, nqf)]';
        from_bus = [V_index; P_index; PQ_transit; Q_index; PQ_transit; bus_indices_branch(PQf_index,1); bus_indices_branch(PQf_index,1)];
        to_bus = [zeros(1, nvi + npi + nqi + npt + nqt), bus_indices_branch(PQf_index,2)', bus_indices_branch(PQf_index,2)']';

    case 14 % IEEE 15-bus system
        meas_index = 1:length(Z_real);
        meas_type = [ones(1,nvi), 2*ones(1,npi), 3*ones(1,nqi)]';
        from_bus = [V_index, P_index, Q_index]';
        to_bus = [zeros(1, nvi + npi + nqi)]';

    case 11 % IEEE 6-bus system
        meas_index = 1:length(Z_real);
        meas_type = [ones(1,nvi), 2*ones(1,npi), 4*ones(1,npf), 5*ones(1, nqf)]';
        from_bus = [V_index, P_index, bus_indices_branch(PQf_index,1)']';
        to_bus = [zeros(1, nvi + npi + nqi), bus_indices_branch(PQf_index,2)']';        

    otherwise
end
    
    
%       |Msnt |     Type |      Value     |     From  |     To |    Rii (sigma) |
    zdt = [meas_index'   meas_type   Z_real'     from_bus    to_bus  sigma_meas];


end
