function V_est = WLS_ieee_15(X, P)
% Power System State Estimation using Weighted Least Square Method..

num = P.Nbus;
nbus = num; % Get number of buses..
nvi = P.nvi; % Number of Voltage measurements
npi = P.npi; % Number of Real Power Injection measurements
nqi = P.nqi; % Number of Reactive Power Injection measurements
npf = P.npf; % Number of Real Power Flow measurements
nqf = P.nqf; % Number of Reactive Power Flow measurements
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
          
    zdata = zdatas(Z_real(MCit,:), P.sigma_meas, P.nvi, P.npi, P.nqi, [], [], P.npf, P.nqf,...
                            P.V_meas_bus, P.P_meas_bus, P.Q_meas_bus, P.PQf_meas_branch, [], bus_indices_branch); % Get Measurement data
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

    iter = 1;
    tol_abs = 5;
    
    while(tol_abs > 1e-4)

        %Measurement Function, h
        h1 = V(fbus(vi),1);
        h2 = zeros(npi,1);
        h3 = zeros(nqi,1);

        for i = 1:npi
            m = fbus(ppi(i));
            for k = 1:nbus
                h2(i) = h2(i) + V(m)*V(k)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
            end
        end

        for i = 1:nqi
            m = fbus(qi(i));
            for k = 1:nbus
                h3(i) = h3(i) + V(m)*V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
            end
        end


        h = [h1; h2; h3];

        % Residue
        r = z - h;

        % Jacobian..
        % H11 - Derivative of V with respect to angles.. All Zeros
        H11 = zeros(nvi,nbus-1);

        % H12 - Derivative of V with respect to V
        H12 = zeros(nvi,nbus);
        for k = 1:nvi
            for n = 1:nbus
                if n == k
                    H12(k,n) = 1;
                end
            end
        end

        % H21 - Derivative of Real Power Injections with Angles
        H21 = zeros(npi,nbus-1);
        for i = 1:npi
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
        H22 = zeros(npi,nbus);
        for i = 1:npi
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
        H31 = zeros(nqi,nbus-1);
        for i = 1:nqi
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
        H32 = zeros(nqi,nbus);
        for i = 1:nqi
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


        % Measurement Jacobian, H
        H = [H11 H12; H21 H22; H31 H32];
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

%     CvE(:,MCit) = diag(inv(H'*inv(Ri)*H)); % Covariance matrix

    Del = 180/pi*del;

    V_e_mod(MCit,:) = V;
    V_e_ang(MCit,:) = Del;
    
    fprintf('  %d  \n', MCit);
    
    
end
    
    V_est = [V_e_ang, V_e_mod];
end