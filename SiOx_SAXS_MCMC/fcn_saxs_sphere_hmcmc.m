function I_cal = fcn_saxs_sphere_hmcmc(x1,qdata)
global fitFlag x2
global lbub

% --- Generate all the fitting parameters
x=zeros(length(fitFlag),1);
x(fitFlag==1) = x1;
x(fitFlag==0) = x2;
x = inversetransformx(x,lbub);

Inorm       = x(1);
Ibk          = x(2);
R           = x(3);
sigma_R     = x(4);
sigma_q     = x(5);

% if sigma_R ~= 0
%     if polymodel_flag == 1      % gaussain
%         nPoly = 20;         % # of points
%         fPoly = 3;          % factor (number of sigmas to account)
%         R_span = linspace(max(eps,R-sigma_R*fPoly),R+sigma_R*fPoly,nPoly);
%         gPoly = 1/(sqrt(2*pi)*sigma_R)*exp(-(R_span-R).^2/(2*sigma_R^2));  % Gaussian function
%     elseif polymodel_flag == 2  % schultz        
%         nPoly = 20;         % # of points
%         fPoly = 3;          % factor (number of sigmas to account)
%         sigma_z_core = sigma_R;
%         z_core = R;
%         schultz_rho = sigma_z_core/z_core;
%         schultz_z = (1/schultz_rho)^2-1;
%         R_span = linspace(max(eps,z_core-sigma_z_core*fPoly),z_core+sigma_z_core*fPoly*1.5,nPoly);
%         %gPoly = ((schultz_z+1)/z_core)^(schultz_z+1)*R_span.^schultz_z.*exp(-(schultz_z+1)*R_span/z_core) / gamma(schultz_z+1);
%         gPoly = ((schultz_z+1)/z_core)^(schultz_z+1).*exp(schultz_z.*log(R_span)-(schultz_z+1)*R_span/z_core) / gamma(schultz_z+1);
%         
%         
%     end
% else
%     R_span = R;
%     gPoly = 1;
% end
%  
% --- gaussian polydispersity
if sigma_R ~= 0
    nPoly = 11;         % # of points
    fPoly = 3;          % factor (number of sigmas to account)
    R_span = linspace(max(eps,R-sigma_R*fPoly),R+sigma_R*fPoly,nPoly);
    gPoly = 1/(sqrt(2*pi)*sigma_R)*exp(-(R_span-R).^2/(2*sigma_R^2));  % Gaussian function
else
    R_span = R;
    gPoly = 1;
end
% figure
% plot(R_span,gPoly);
% trapz(R_span,gPoly)

% --- loop through q's
nq = length(qdata);
I_cal = nan(nq,1);
for ii = 1:nq
    iq = qdata(ii);
    if sigma_q >1e-9
        nRes = 10;         % # of points
        fRes = 3;          % factor (number of sigmas to account)
        q_span = linspace(max(eps,iq-sigma_q*fRes),iq+sigma_q*fRes,nRes);
        q_span = q_span(:);
        gRes = 1/(sqrt(2*pi)*sigma_q)*exp(-(q_span-iq).^2/(2*sigma_q^2));  % Gaussian function
    else
        nRes = 1;
        q_span = iq;
        gRes = 1;
    end
%     figure
%     plot(q_span,gRes);
%     trapz(q_span,gRes)    
    % --- calcualte F matrix for q_span and R_span
    [R_span_2D,q_span_2D] = meshgrid(R_span,q_span);
    gPoly_2D = repmat(gPoly,nRes,1);
    f = formfactor_sphere(q_span_2D,R_span_2D);
    I = trapz(R_span,abs(f).^2.*gPoly_2D,2);
    
    % --- volume normalization
    V = 4/3*pi*R_span_2D.^3;
    VV = trapz(R_span,abs(V).^2.*gPoly_2D,2); 
    I = I./VV;

    % --- integration over resolution
    if sigma_q == 0
        I_cal(ii) = I;
    else
        I_cal(ii) = trapz(q_span,I.*gRes);
    end
end
I_cal = I_cal*Inorm+Ibk;    

function f = formfactor_sphere(q,R)
qR = q.*R;
f = 4*pi*R.^3.*(sin(qR)-qR.*cos(qR))./(qR).^3;