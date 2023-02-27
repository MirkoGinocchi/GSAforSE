%% HELP YBUSPPG
% This function computes the Ybus matrix of the network.
%% ********************************************************************* %%

function ybus = ybusppg(branch, bus_indices_branch, num)  % Returns ybus

nbus = num;    % no. of buses

fb = bus_indices_branch(:,1);     % From bus number
tb = bus_indices_branch(:,2);     % To bus number
nbranch = length(fb);           % no. of branches
b = zeros(1,nbranch);
r = branch(1:nbranch);      % Resistance, R
x = branch(nbranch + 1:2*nbranch);      % Reactance, X
if nbus == 49 || nbus == 6 || nbus == 81 || nbus == 40 % These grids have ground admittances specified. 
                                                       % If lines are modelled using R-X only, it is not necessary to include bus admittances
    b = branch(2*nbranch + 1:3*nbranch)./2;   % Ground Admittance, B/2
    b = 1i*b;                % Make B imaginary
end
a = ones(length(fb),1);        % Tap setting value
z = r + 1i*x;            % Z matrix
y = 1./z;               % To get inverse of each element



ybus = zeros(nbus,nbus);        % Initialise YBus
 
 % Formation of the Off Diagonal Elements
 for k=1:nbranch
     ybus(fb(k),tb(k)) = ybus(fb(k),tb(k))-y(k)/a(k);
     ybus(tb(k),fb(k)) = ybus(fb(k),tb(k));
 end
 
 % Formation of Diagonal Elements
 for m =1:nbus
     for n =1:nbranch
         if fb(n) == m
             ybus(m,m) = ybus(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             ybus(m,m) = ybus(m,m) + y(n) + b(n);
         end
     end
 end
 %ybus;                  % Bus Admittance Matrix
 %zbus = inv(ybus);      % Bus Impedance Matrix