function U = potentialE(X,nsets,alphai_set,Iset,Iweight_set)
Sigmavi = exp(X(end-nsets+1:end));
I_cal_set = fcn_xray_ref_hmcmc(X(1:end-nsets),alphai_set);
U = 0;
for iset = 1:nsets
    I_cal = I_cal_set{iset};
    I = Iset{iset};
    Iweight = Iweight_set{iset};
    logp = -( log10(I_cal)-log10(I) ).^2 / (2*Sigmavi(iset)) - log(sqrt(2*pi*Sigmavi(iset)));
    U = U - sum(Iweight.*logp);
end
