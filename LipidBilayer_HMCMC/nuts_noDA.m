function [xp,logp_xp,dPotE] = nuts_noDA(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max)
% one-step NUTS without dual averaging
pp = mvnrnd(zeros(size(xp)),M);
logp_pp = -K(pp);
logp = logp_xp + logp_pp;
u = unifrnd(0,exp(logp));
xp_neg = xp;
xp_pos = xp;
pp_neg = pp;
pp_pos = pp;
jj = 0;
n = 1;
s = true;
while s==1
    vj = randi([0,1])*2-1;  % randomly choose -1 or 1
    if vj==-1  % grow to negative side
        [xp_neg,pp_neg,~,~,xp1,n1,s1,logp_xp1,dPotE1] = buildtree(xp_neg,pp_neg,dPotE,u,vj,jj,epsilon,U,K,dU,dK,delta_max);
    else % grow tree to positive side
        [~,~,xp_pos,pp_pos,xp1,n1,s1,logp_xp1,dPotE1] = buildtree(xp_pos,pp_pos,dPotE,u,vj,jj,epsilon,U,K,dU,dK,delta_max);
    end
    if s1==1 && rand<n1/n
        xp = xp1;
        logp_xp = logp_xp1;
        dPotE = dPotE1;
    end
    n = n + n1;
    s = s1 & ((xp_pos-xp_neg)*pp_neg'>=0) & ((xp_pos-xp_neg)*pp_pos'>=0);
    jj = jj + 1;
end
end % EOF


function [xp_neg,pp_neg,xp_pos,pp_pos,xp1,n1,s1,logp_xp1,dPotE1] = buildtree(xp,pp,dPotE,u,v,jj,epsilon,U,K,dU,dK,delta_max)
if jj==0
    [xp1,pp1,dPotE1,PotE1] = leapfrogonce(xp,pp,dPotE,v*epsilon,dU,dK);  % one-step leapfrog and direction given by v
    logp_xp1 = -PotE1;
    logp_pp1 = -K(pp1);
    logp1 = logp_xp1 + logp_pp1;    
    n1 = double(log(u)<=logp1);       % numeric
    s1 = log(u)<(logp1+delta_max);    % logical
    xp_neg = xp1;
    xp_pos = xp1;
    pp_neg = pp1;
    pp_pos = pp1;
else
    [xp_neg,pp_neg,xp_pos,pp_pos,xp1,n1,s1,logp_xp1,dPotE1] = buildtree(xp,pp,dPotE,u,v,jj-1,epsilon,U,K,dU,dK,delta_max);
    if s1==1
        if v == -1
            [xp_neg,pp_neg,~,~,xp2,n2,s2,logp_xp2,dPotE2] = buildtree(xp_neg,pp_neg,dPotE1,u,v,jj-1,epsilon,U,K,dU,dK,delta_max);
        else
            [~,~,xp_pos,pp_pos,xp2,n2,s2,logp_xp2,dPotE2] = buildtree(xp_pos,pp_pos,dPotE1,u,v,jj-1,epsilon,U,K,dU,dK,delta_max);
        end
        if rand<n2/(n1+n2)
            xp1 = xp2;
            logp_xp1 = logp_xp2;
            dPotE1 = dPotE2;
        end
        s1 = s2 & ((xp_pos-xp_neg)*pp_neg'>=0) & ((xp_pos-xp_neg)*pp_pos'>=0);
        n1 = n1 + n2;
    end
end
end    % EOF