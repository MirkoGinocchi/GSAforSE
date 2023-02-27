%% HELP BBUSPPG
% This function returns shunt admittances.
%% ********************************************************************* %%
function bbus = bbusppg(branch, bus_indices_branch, num)     % Returns B-bus


fb = bus_indices_branch(:,1);
tb = bus_indices_branch(:,2);
nbranch = length(fb);           % no. of branches
b = zeros(1,nbranch);
if num == 49 || num == 6 || num == 81 || num == 40 
    b = branch(2*nbranch + 1:end)./2;
end
nbus = num;   % no. of buses
bbus = zeros(nbus,nbus);

 for k=1:nbranch
     bbus(fb(k),tb(k)) = b(k);
     bbus(tb(k),fb(k)) = bbus(fb(k),tb(k));
 end