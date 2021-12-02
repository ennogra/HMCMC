function [xp,pp,dPotE,PotE] = leapfrog(xp,pp,dPotE,epsilon,L,dU,dK)
% Function: leap frog function (no need to flip sign of pp because of the symmetry of the algorithm)
%[dPotE,~] = dU(xp);
pp = pp - epsilon/2*dPotE;     % inital half update for p
for jj=1:L-1                    % L-1 update for x and p
    xp = xp + epsilon*dK(pp);
    [dPotE,~] = dU(xp);    
    pp = pp - epsilon*dPotE;
end
xp = xp + epsilon*dK(pp);       % last update for x
[dPotE,PotE] = dU(xp);    
pp = pp - epsilon/2*dPotE;     % last half update for p
end

