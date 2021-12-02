function U = potentialE(X,vz,I,weight,sspFlag,XFix)
XAll = nan(1,length(sspFlag));
XAll(sspFlag) = X;
XAll(~sspFlag) = XFix;
Sigmavi = exp(XAll(end));
logp = -( log10(fcn_flu_spline_hmcmc(XAll(1:end-1),vz))-log10(I) ).^2 / (2*Sigmavi) - log(sqrt(2*pi*Sigmavi));
U = -sum(weight.*logp);  
U = U - cbs_regularization(XAll);

function reg = cbs_regularization(X)
global beta1 beta2
global fitFlag x2
global lbub
x=zeros(length(fitFlag),1);
x(fitFlag==1) = X(1:end-1);
x(fitFlag==0) = x2;
x = inversetransformx(x,lbub);
a = x(17:end);
a = a(:);

A1 = sum( abs(a) );
a_diff = diff(a);       
A2 = sum( abs(a_diff).^2 );        % LASSO
% Nc1 = (1-w)*A1;         % due to smoothness
% Nc2 = w*A2;             % due to deviation to intial amean
% Nc_a = Nc1+Nc2;
reg = beta1*A1 + beta2*A2;