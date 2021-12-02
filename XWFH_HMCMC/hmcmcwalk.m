function [xp,logp_xp,dPotE,logr] = hmcmcwalk(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,lambda)
% --- for simple HMCMC (no NUTS)

% --- propose pp
pp = mvnrnd(zeros(size(xp)),M);

logp_pp = -K(pp);

% --- move one step along MC
Lm = max(1,round(lambda/epsilon));  % # of hamiltonian walks
[xp1,pp1,dPotE1,PotE1] = leapfrog(xp,pp,dPotE,epsilon,Lm,dU,dK);
%logp_xp1 = -U(xp1);
logp_xp1 = -PotE1;
logp_pp1 = -K(pp1);
logr = logp_xp1 + logp_pp1 - (logp_xp+logp_pp);
if log(rand)<logr   % update; otherwise, return the same xp and logp_xp values
    xp = xp1;
    logp_xp = logp_xp1;     % update logp_xp so no need to recalculate it in the next iteration
    dPotE = dPotE1;
end