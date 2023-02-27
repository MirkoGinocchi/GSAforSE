function [losses, mean_I_flowing] = WLS_plus_losses_ieee_69(X, P)
% Power System State Estimation using Weighted Least Square Method

num = P.Nbus;
nbus = num; % Get number of buses
nvi = P.nvi; % Number of Voltage measurements
npi = P.npi; % Number of Real Power Injection measurements
nqi = P.nqi; % Number of Reactive Power Injection measurements
npf = P.npf; % Number of Real Power Flow measurements
nqf = P.nqf; % Number of Reactive Power Flow measurements
npt = P.npt; % Number of Real virtual measurements
nqt = P.nqt; % Number of Reactive virtual measurements

bus_indices_branch = [P.fbus, P.tbus];
RXB_samples = X(:,1:2*P.Nbranches);
Z_real = X(:,2*P.Nbranches + 1:end);

V_e_mod = zeros(P.Nsamples, num);
V_e_ang = V_e_mod;

parfor MCit = 1:length(X(:,1))

    ybus = ybusppg(RXB_samples(MCit,:), bus_indices_branch, num); % Get YBus
    bpq = bbusppg(RXB_samples(MCit,:), bus_indices_branch, num); % Get B data
    G = real(ybus);
    B = imag(ybus);
          
    zdata = zdatas(Z_real(MCit,:), P.sigma_meas, P.nvi, P.npi, P.nqi, P.npt, P.nqt, P.npf, P.nqf,...
                            P.V_meas_bus, P.P_meas_bus, P.Q_meas_bus, P.PQf_meas_branch, P.PQ_transit, bus_indices_branch); % Get Measurement data..  
    fbus = zdata(:,4); % From bus
    tbus = zdata(:,5); % To bus
    type = zdata(:,2); % Type of measurement, Vi - 1, Pi - 2, Qi - 3, Pij - 4, Qij - 5, Iij - 6
    z = zdata(:,3); % Measuement values

    Ri = diag(zdata(:,6).^2); % Measurement Error
    W_onehalf = chol(Ri\eye(size(Ri,1))); % W^1/2
    V = ones(nbus,1);
    del = zeros(nbus,1);
    E = [del(1:nbus-1); V];   % State Vector

    vi = find(type == 1); % Index of voltage magnitude measurements
    ppi = find(type == 2); % Index of real power injection measurements
    qi = find(type == 3); % Index of reactive power injection measurements
    pf = find(type == 4); % Index of real powerflow measurements
    qf = find(type == 5); % Index of reactive powerflow measurements
    
    iter = 1;
    tol_abs = 5;
    
    while(tol_abs > 1e-4)

        %Measurement Function, h
        h1 = V(fbus(vi),1);
        h2 = zeros(npi + npt,1);
        h3 = zeros(nqi + nqt,1);
        h4 = zeros(npf,1);
        h5 = zeros(nqf,1);
        
        for i = 1:npi + npt
            m = fbus(ppi(i));
            for k = 1:nbus
                h2(i) = h2(i) + V(m)*V(k)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
            end
        end

        for i = 1:nqi + nqt
            m = fbus(qi(i));
            for k = 1:nbus
                h3(i) = h3(i) + V(m)*V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
            end
        end

        for i = 1:npf
            m = fbus(pf(i));
            p = tbus(pf(i));
            h4(i) = -V(m)^2*G(m,p) - V(m)*V(p)*(-G(m,p)*cos(del(m)-del(p)) - B(m,p)*sin(del(m)-del(p)));
        end

        for i = 1:nqf
            m = fbus(qf(i));
            p = tbus(qf(i));
            h5(i) = -V(m)^2*(-B(m,p)+bpq(m,p)) - V(m)*V(p)*(-G(m,p)*sin(del(m)-del(p)) + B(m,p)*cos(del(m)-del(p)));
        end

        h = [h1; h2; h3; h4; h5];

        % Residue
        r = z - h;

        % Jacobian..
        % H11 - Derivative of V with respect to angles All Zeros
        H11 = zeros(nvi,nbus-1);

        % H12 - Derivative of V with respect to V.. 
        H12 = zeros(nvi,nbus);
        for k = 1:nvi
            for n = 1:nbus
                if n == k
                    H12(k,n) = 1;
                end
            end
        end

        % H21 - Derivative of Real Power Injections with Angles
        H21 = zeros(npi + npt,nbus-1);
        for i = 1:npi + npt
            m = fbus(ppi(i));
            for k = 1:(nbus-1)
                if k+1 == m
                    for n = 1:nbus
                        H21(i,k) = H21(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                    end
                    H21(i,k) = H21(i,k) - V(m)^2*B(m,m);
                else
                    H21(i,k) = V(m)* V(k+1)*(G(m,k+1)*sin(del(m)-del(k+1)) - B(m,k+1)*cos(del(m)-del(k+1)));
                end
            end
        end

        % H22 - Derivative of Real Power Injections with V
        H22 = zeros(npi + npt,nbus);
        for i = 1:npi + npt
            m = fbus(ppi(i));
            for k = 1:(nbus)
                if k == m
                    for n = 1:nbus
                        H22(i,k) = H22(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                    end
                    H22(i,k) = H22(i,k) + V(m)*G(m,m);
                else
                    H22(i,k) = V(m)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
                end
            end
        end

        % H31 - Derivative of Reactive Power Injections with Angles
        H31 = zeros(nqi + nqt,nbus-1);
        for i = 1:nqi + nqt
            m = fbus(qi(i));
            for k = 1:(nbus-1)
                if k+1 == m
                    for n = 1:nbus
                        H31(i,k) = H31(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                    end
                    H31(i,k) = H31(i,k) - V(m)^2*G(m,m);
                else
                    H31(i,k) = V(m)* V(k+1)*(-G(m,k+1)*cos(del(m)-del(k+1)) - B(m,k+1)*sin(del(m)-del(k+1)));
                end
            end
        end

        % H32 - Derivative of Reactive Power Injections with V
        H32 = zeros(nqi + nqt,nbus);
        for i = 1:nqi + nqt
            m = fbus(qi(i));
            for k = 1:(nbus)
                if k == m
                    for n = 1:nbus
                        H32(i,k) = H32(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                    end
                    H32(i,k) = H32(i,k) - V(m)*B(m,m);
                else
                    H32(i,k) = V(m)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
                end
            end
        end

        % H41 - Derivative of Real Power Flows with Angles
        H41 = zeros(npf,nbus-1);
        for i = 1:npf
            m = fbus(pf(i));
            p = tbus(pf(i));
            for k = 1:(nbus-1)
                if k+1 == m
                    H41(i,k) = V(m)* V(p)*(-G(m,p)*sin(del(m)-del(p)) + B(m,p)*cos(del(m)-del(p)));
                else 
                    if k+1 == p
                    H41(i,k) = -V(m)* V(p)*(-G(m,p)*sin(del(m)-del(p)) + B(m,p)*cos(del(m)-del(p)));
                    else
                        H41(i,k) = 0;
                    end
                end
            end
        end

        % H42 - Derivative of Real Power Flows with V
        H42 = zeros(npf,nbus);
        for i = 1:npf
            m = fbus(pf(i));
            p = tbus(pf(i));
            for k = 1:nbus
                if k == m
                    H42(i,k) = -V(p)*(-G(m,p)*cos(del(m)-del(p)) - B(m,p)*sin(del(m)-del(p))) - 2*G(m,p)*V(m);
                else 
                    if k == p
                    H42(i,k) = -V(m)*(-G(m,p)*cos(del(m)-del(p)) - B(m,p)*sin(del(m)-del(p)));
                    else
                        H42(i,k) = 0;
                    end
                end
            end
        end

        % H51 - Derivative of Reactive Power Flows with Angles
        H51 = zeros(nqf,nbus-1);
        for i = 1:nqf
            m = fbus(qf(i));
            p = tbus(qf(i));
            for k = 1:(nbus-1)
                if k+1 == m
                    H51(i,k) = -V(m)* V(p)*(-G(m,p)*cos(del(m)-del(p)) - B(m,p)*sin(del(m)-del(p)));
                else 
                    if k+1 == p
                    H51(i,k) = V(m)* V(p)*(-G(m,p)*cos(del(m)-del(p)) - B(m,p)*sin(del(m)-del(p)));
                    else
                        H51(i,k) = 0;
                    end
                end
            end
        end

        % H52 - Derivative of Reactive Power Flows with V
        H52 = zeros(nqf,nbus);
        for i = 1:nqf
            m = fbus(qf(i));
            p = tbus(qf(i));
            for k = 1:nbus
                if k == m
                    H52(i,k) = -V(p)*(-G(m,p)*sin(del(m)-del(p)) + B(m,p)*cos(del(m)-del(p))) - 2*V(m)*(-B(m,p)+ bpq(m,p));
                else 
                    if k == p
                    H52(i,k) = -V(m)*(-G(m,p)*sin(del(m)-del(p)) + B(m,p)*cos(del(m)-del(p)));
                    else
                        H52(i,k) = 0;
                    end
                end
            end
        end

        % Measurement Jacobian, H
        H = [H11 H12; H21 H22; H31 H32; H41 H42; H51 H52];
        H_tilde = W_onehalf * H;
        r_tilde = W_onehalf * r;
        [Q_n,U] = qr(H_tilde,0);
        r_q = Q_n' * r_tilde;

        % solve by backward substitution
        m = length(r_q);
        dE = zeros(m,1);
        dE(m) = r_q(m)/U(m,m);
            for k = m-1:-1:1
                dE(k) = (r_q(k)-U(k,k+1:m)*dE(k+1:m))/U(k,k);
            end

        E = E + dE;
        del(2:end) = E(1:nbus-1);
        V = E(nbus:end);
        iter = iter + 1;
        tol_abs = max(abs(dE));

       
    end

%     CvE(:,MCit) = diag(inv(H'*inv(Ri)*H)); % Covariance matrix..

    Del = del;

    V_e_mod(MCit,:) = V;
    V_e_ang(MCit,:) = Del;

    V_r = V.*cos(Del);
    V_i = V.*sin(Del);

    V_r_from = V_r(bus_indices_branch(:,1));
    V_i_from = V_i(bus_indices_branch(:,1));
    V_r_to = V_r(bus_indices_branch(:,2));
    V_i_to = V_i(bus_indices_branch(:,2));

    I_flowing_r_1 = (((V_r_from - V_r_to))'./...
        (RXB_samples(MCit,1:P.Nbranches).^2 + RXB_samples(MCit,P.Nbranches+1:end).^2)).*RXB_samples(MCit,1:P.Nbranches);
    I_flowing_r_2 = (((V_i_from - V_i_to))'./...
        (RXB_samples(MCit,1:P.Nbranches).^2 + RXB_samples(MCit,P.Nbranches+1:end).^2)).*RXB_samples(MCit,P.Nbranches+1:end);
    
    I_flowing_i_1 = (((V_i_from - V_i_to))'./...
        (RXB_samples(MCit,1:P.Nbranches).^2 + RXB_samples(MCit,P.Nbranches+1:end).^2)).*RXB_samples(MCit,1:P.Nbranches);
    I_flowing_i_2 = -(((V_r_from - V_r_to))'./...
        (RXB_samples(MCit,1:P.Nbranches).^2 + RXB_samples(MCit,P.Nbranches+1:end).^2)).*RXB_samples(MCit,P.Nbranches+1:end);    
 
    I_flowing_r = I_flowing_r_1 + I_flowing_r_2;
    I_flowing_i = I_flowing_i_1 + I_flowing_i_2;
    
    I_flowing_square = I_flowing_r.^2 + I_flowing_i.^2;
    loss(MCit) = RXB_samples(MCit,1:P.Nbranches)*(I_flowing_square)';
    
    fprintf('  %d  \n', MCit);    
    I_flowing_real(MCit,:) = I_flowing_r;
    I_flowing_imag(MCit,:) = I_flowing_i;
end
    mean_I_flowing = mean(I_flowing_real,1) + 1j.*mean(I_flowing_real,1);
    losses = loss;
end