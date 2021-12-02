function [xp,pp,dPotE,PotE] = leapfrogonce(xp,pp,dPotE,epsilon,dU,dK)
% leap frog function (Do NOT flip the sign of pp, the algorithm already takes care of the reversibility)
pp = pp - epsilon/2*dPotE;
xp = xp + epsilon*dK(pp);
[dPotE,PotE] = dU(xp);
pp = pp - epsilon/2*dPotE;
