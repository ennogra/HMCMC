function xtrans = transformx(x,lbub)
% Transform x so it becomes unbounded
% x: p x n (p is the dimension of parameter space, n is the number of
% data points)

% (-Inf,Inf)
xtrans = x;             
% (a,b)
idx1 = lbub(:,3) == 1;  
xtrans(idx1,:) = log( (x(idx1,:)-lbub(idx1,1)) ./ (lbub(idx1,2)-x(idx1,:)) );
% (a,Inf)
idx2 = lbub(:,3) == 2;   
xtrans(idx2,:) = log( x(idx2,:)-lbub(idx2,1) );
% (-Inf,b)
idx3 = lbub(:,3) == 3;
xtrans(idx3,:) = log( lbub(idx3,2) - x(idx3,:) );
% no transformation
idx5 = lbub(:,3) == 5;
xtrans(idx5,:) = x(idx5,:);
