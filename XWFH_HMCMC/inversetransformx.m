function x=inversetransformx(xtrans,lbub)
% inverse transform to x
% transx: p x n (p is the dimension of parameter space, n is the number of
% data points)

% (-Inf,Inf)
x = xtrans;             
% (a,b)
idx1 = lbub(:,3) == 1;  
x(idx1,:) = ( lbub(idx1,2) + lbub(idx1,1).*exp(-xtrans(idx1,:)) ) ./ (1+exp(-xtrans(idx1,:)));
% (a,Inf)
idx2 = lbub(:,3) == 2;   
x(idx2,:) = lbub(idx2,1) + exp(xtrans(idx2,:));
% (-Inf,b)
idx3 = lbub(:,3) == 3;
x(idx3,:) = lbub(idx3,2) - exp(xtrans(idx3,:));
% no transformation
idx5 = lbub(:,3) == 5;
x(idx5,:) = xtrans(idx5,:);
