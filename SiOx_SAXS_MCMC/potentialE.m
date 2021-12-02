function U = potentialE(X,qdata,Iqdata,weight)
Sigmavi = exp(X(end));
logp = -( log10(fcn_saxs_sphere_hmcmc(X(1:end-1),qdata))-log10(Iqdata) ).^2 / (2*Sigmavi) - log(sqrt(2*pi*Sigmavi));
U = -sum(weight.*logp);  
