function [epsilon,counter] = heuristic4epsilon(xp_init,pp_init,dPotE_init,epsilon_init,U,K,dU,dK,logp_init)
% Function: Heuristic for initial epsilon 
% this method works best if only one probability mode exists.
epsilon = epsilon_init;
% [xp1,pp1] = leapfrogonce(xp_init,pp_init,epsilon,dU,dK);  % leap frog once
% logp1 = -(U(xp1)+K(pp1));
[~,pp1,~,PotE1] = leapfrogonce(xp_init,pp_init,dPotE_init,epsilon,dU,dK);  % leap frog once
logp1 = -(PotE1+K(pp1));
logr = logp1 - logp_init;
a = 2*(exp(logr)>0.5)-1;
counter = 0;
while a*logr>-a*log(2) % equivalent to exp(logr)^a > 2^(-a)
    epsilon = 2^a*epsilon;
    %     [xp1,pp1] = leapfrogonce(xp_init,pp_init,epsilon,dU,dK);
    %     logp1 = -(U(xp1)+K(pp1));
    [~,pp1,~,PotE1] = leapfrogonce(xp_init,pp_init,dPotE_init,epsilon,dU,dK);
    logp1 = -(PotE1+K(pp1));
    logr = logp1 - logp_init;
    counter  = counter+1;
end



