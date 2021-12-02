function [Y,Y1,Y2,Y_elastic_norm,edp0,p0] = fcn_flu_spline_hmcmc(x1,vz)
global fitFlag x2
global lbub
global sdd xenergy nSlice alpha0 plotEDPFlag
global N nx B
global si cr pd ptba ps au  % properties of materials
global rho

% --- Generate all the fitting parameters
x=zeros(length(fitFlag),1);
x(fitFlag==1) = x1;
x(fitFlag==0) = x2;
x = inversetransformx(x,lbub);

Inorm       = x(1);
Ibk         = x(2);
alphai_res  = x(3);         % exit angle resolution
vz0         = x(4);
sigma_si    = x(5);
d_cr        = x(6);
sigma_cr    = x(7);
d_pd        = x(8);
sigma_pd    = x(9);
d_ptba      = x(10);
sigma_ptba  = x(11);
d_ps        = x(12);
sigma_ps    = x(13);
d_au        = x(14);
xe_ratio    = x(15);
elastic_ratio = x(16);
a = x(17:end);
a = a(:);

%% exit angle and resolution span
alphai = atand((vz-vz0)/sdd);       % exit angle 

% - list of alphai due to resolution 
if alphai_res > eps
    nspan_alphai = 11;
    fspan_alphai = 3;
    alphai_span_2D = nan(length(alphai),nspan_alphai);
    gspan = nan(length(alphai),nspan_alphai);
    for ii=1:length(alphai)
        alphai_span_2D(ii,:) = linspace(alphai(ii)- fspan_alphai*alphai_res,alphai(ii)+ fspan_alphai*alphai_res,nspan_alphai);
        gspan(ii,:) = 1/sqrt(2*pi)/alphai_res * exp(-(alphai_span_2D(ii,:)  -  alphai(ii) ).^2/2/alphai_res^2);
    end    
    alphai_span = alphai_span_2D(:);
    
else
    alphai_span = alphai;
end


%% spline profile for gold
p_au0 = B*a;        % un-normalized spline profile of gold layer


%% calcualte intensity
% alpha0 is the incident angle
[EF1,edp1,p1] = flu(alphai_span,1,...
    d_au,p_au0,...
    d_ps,d_ptba,d_pd,d_cr,...
    sigma_ps,sigma_ptba,sigma_pd,sigma_cr,sigma_si);

[EF2,edp2,p2] = flu(alphai_span,2,...
    d_au,p_au0,...
    d_ps,d_ptba,d_pd,d_cr,...
    sigma_ps,sigma_ptba,sigma_pd,sigma_cr,sigma_si);


% --- elastic scattering for Pd and Au
% % incident angle side
% [EF0i,edp0i,p0i] = flu(alpha0,3,...
%     d_au,p_au0,...
%     d_ps,d_ptba,d_pd,d_cr,...
%     sigma_ps,sigma_ptba,sigma_pd,sigma_cr,sigma_si);
% % exit angle side
% [EF0f,edp0,p0] = flu(alphai_span,3,...
%     d_au,p_au0,...
%     d_ps,d_ptba,d_pd,d_cr,...
%     sigma_ps,sigma_ptba,sigma_pd,sigma_cr,sigma_si);

% calculate both incident and exit at the same time to save time
[EF0if,edp0,p0] = flu([alpha0;alphai_span],3,...
    d_au,p_au0,...
    d_ps,d_ptba,d_pd,d_cr,...
    sigma_ps,sigma_ptba,sigma_pd,sigma_cr,sigma_si);
EF0i = EF0if(:,1);
EF0f = EF0if(:,2:end);

%% %%%%%%%%%%%%%%%%%%%%%%
% --- fluoresnce signals from au
% xe1
Y_tmp_1 = abs(EF0i).^2.*p1(:,1)*ones(1,length(alphai_span)) .*abs(EF1).^2;
Y1 = NaN*ones(length(alphai_span),1);
% xe2
Y_tmp_2 = abs(EF0i).^2.*p2(:,1)*ones(1,length(alphai_span)) .*abs(EF2).^2;
Y2 = NaN*ones(length(alphai_span),1);
% loop angles
for ii=1:length(alphai_span)
    Y1(ii) = trapz(edp1(:,1),Y_tmp_1(:,ii));
    Y2(ii) = trapz(edp2(:,1),Y_tmp_2(:,ii));
end
Y2 = Y2*xe_ratio;
% get resolution convolutioned result
if alphai_res > eps
    Y1_tmp = NaN*ones(length(alphai),1);
    Y2_tmp = NaN*ones(length(alphai),1);
    Y1 = reshape(Y1,length(alphai),[]);
    Y2 = reshape(Y2,length(alphai),[]);
    for ii=1:length(alphai)
        Y1_tmp(ii) = trapz(alphai_span_2D(ii,:),Y1(ii,:).*gspan(ii,:));
        Y2_tmp(ii) = trapz(alphai_span_2D(ii,:),Y2(ii,:).*gspan(ii,:));
    end
    Y1 = Y1_tmp;
    Y2 = Y2_tmp;
end
Y_flu = Y1+Y2;

% --- elastic scattering from au,PS, ptba, Pd, Cr, Si
Y_tmp_au     = abs(EF0i).^2.*p0(:,1)*ones(1,length(alphai_span)) .*abs(EF0f).^2;
Y_tmp_ps     = abs(EF0i).^2.*p0(:,2)*ones(1,length(alphai_span)) .*abs(EF0f).^2;
Y_tmp_ptba   = abs(EF0i).^2.*p0(:,3)*ones(1,length(alphai_span)) .*abs(EF0f).^2;
Y_tmp_pd     = abs(EF0i).^2.*p0(:,4)*ones(1,length(alphai_span)) .*abs(EF0f).^2;
Y_tmp_cr     = abs(EF0i).^2.*p0(:,5)*ones(1,length(alphai_span)) .*abs(EF0f).^2;
Y_tmp_si     = abs(EF0i).^2.*p0(:,6)*ones(1,length(alphai_span)) .*abs(EF0f).^2;
Y_au = NaN*ones(length(alphai_span),1);
Y_ps = NaN*ones(length(alphai_span),1);
Y_ptba = NaN*ones(length(alphai_span),1);
Y_pd = NaN*ones(length(alphai_span),1);
Y_cr = NaN*ones(length(alphai_span),1);
Y_si = NaN*ones(length(alphai_span),1);
for ii=1:length(alphai_span)
    Y_au(ii)    = trapz(edp0(:,1),Y_tmp_au(:,ii));  
    Y_ps(ii)    = trapz(edp0(:,1),Y_tmp_ps(:,ii));
    Y_ptba(ii)  = trapz(edp0(:,1),Y_tmp_ptba(:,ii));
    Y_pd(ii)    = trapz(edp0(:,1),Y_tmp_pd(:,ii));
    Y_cr(ii)    = trapz(edp0(:,1),Y_tmp_cr(:,ii));
    Y_si(ii)    = trapz(edp0(:,1),Y_tmp_si(:,ii));
end
% get resolution convolutioned result
if alphai_res > eps
    Y_au_tmp = NaN*ones(length(alphai),1);
    Y_ps_tmp = NaN*ones(length(alphai),1);
    Y_ptba_tmp = NaN*ones(length(alphai),1);
    Y_pd_tmp = NaN*ones(length(alphai),1);
    Y_cr_tmp = NaN*ones(length(alphai),1);
    Y_si_tmp = NaN*ones(length(alphai),1);
    
    Y_au = reshape(Y_au,length(alphai),[]);
    Y_ps = reshape(Y_ps,length(alphai),[]);
    Y_ptba = reshape(Y_ptba,length(alphai),[]);
    Y_pd = reshape(Y_pd,length(alphai),[]);
    Y_cr = reshape(Y_cr,length(alphai),[]);
    Y_si = reshape(Y_si,length(alphai),[]);
    
    for ii=1:length(alphai)
        Y_au_tmp(ii) = trapz(alphai_span_2D(ii,:),Y_au(ii,:).*gspan(ii,:));
        Y_ps_tmp(ii) = trapz(alphai_span_2D(ii,:),Y_ps(ii,:).*gspan(ii,:));
        Y_ptba_tmp(ii) = trapz(alphai_span_2D(ii,:),Y_ptba(ii,:).*gspan(ii,:));
        Y_pd_tmp(ii) = trapz(alphai_span_2D(ii,:),Y_pd(ii,:).*gspan(ii,:));
        Y_cr_tmp(ii) = trapz(alphai_span_2D(ii,:),Y_cr(ii,:).*gspan(ii,:));
        Y_si_tmp(ii) = trapz(alphai_span_2D(ii,:),Y_si(ii,:).*gspan(ii,:));
    end
    Y_au = Y_au_tmp;
    Y_ps = Y_ps_tmp;
    Y_ptba = Y_ptba_tmp;
    Y_pd = Y_pd_tmp;
    Y_cr = Y_cr_tmp;  
    Y_si = Y_si_tmp;
end

% --- molar mass  of pd and au
M_ps = 104.14912;
M_ptba = 128.16898;
M_pd = 106.42;
M_cr = 7.19;
M_si = 28.0855;
M_au = 196.97;

% --- number of atoms of pd and au of a unit volume (or unit thickness)
N_ps = ps.massDensity/M_ps;
N_ptba = ptba.massDensity/M_ptba;
N_pd = pd.massDensity/M_pd;
N_cr = cr.massDensity/M_cr;
N_si = si.massDensity/M_si;
N_au = au.massDensity/M_au;

% --- scattering cross section of each atom (at 12.1kev)
cs_h = 4.33248907910636e-08;
cs_c = 6.10032192574906e-06;
cs_o = 1.31496033880765e-05;
cs_si =  5.91047383957732e-05 ;
cs_cr = 0.000230778062489777 ;
cs_pd = 0.00121999632238972;
cs_au = 0.00501882819931888;

% --- scattering cross section of each molecule (at 12.1 Kev)
cs_ps = 8*cs_c + 8*cs_h;
cs_ptba = 7*cs_c + 12*cs_h + 2*cs_o;

% --- elastic is normalized with respect to that of gold
Y_elastic_ps    = Y_ps*(N_ps*cs_ps); 
Y_elastic_ptba  = Y_ptba*(N_ptba*cs_ptba); 
Y_elastic_pd    = Y_pd*(N_pd*cs_pd); 
Y_elastic_cr    = Y_cr*(N_cr*cs_cr); 
Y_elastic_si    = Y_si*(N_si*cs_si); 
Y_elastic_au    = Y_au*(N_au*cs_au); 

Y_elastic_total = Y_elastic_ps + Y_elastic_ptba + Y_elastic_pd + Y_elastic_cr + Y_elastic_si + Y_elastic_au;
%Y_elastic_total = Y_elastic_pd +Y_elastic_au;


% % ---- plot elastic
% figure
% hold on;
% plot(alphai,Y_elastic_ps);
% plot(alphai,Y_elastic_ptba);
% plot(alphai,Y_elastic_pd);
% plot(alphai,Y_elastic_cr);
% plot(alphai,Y_elastic_si);
% plot(alphai,Y_elastic_au);
% plot(alphai,Y_elastic_total,'linewidth',2);
% legend('PS','PtBA','Pd','Cr','Si','Au','Total');
% title('elastic')

Y_elastic_norm = Y_elastic_total/(N_au*cs_au)*elastic_ratio; % elastic is normalized by gold amount

Y = (Y_flu + Y_elastic_norm)*Inorm + Ibk;



function [EF,edp,p] = flu(alphai,ixe,...
    d_au,p_au0,...
    d_ps,d_ptba,d_pd,d_cr,...
    sigma_ps,sigma_ptba,sigma_pd,sigma_cr,sigma_si)
global nSlice plotEDPFlag
global nx
global si cr pd ptba ps au      % properties of materials
global xenergy
global rho

% --- slicing without gold
rawedp          = [...
    0                       0               0               sigma_ps    0;
    d_ps                    ps.dispersion(ixe)   ps.absorption(ixe)  sigma_ptba  0;                    
    d_ptba+d_ps             ptba.dispersion(ixe) ptba.absorption(ixe) sigma_pd    0;            
    d_pd+d_ptba+d_ps        pd.dispersion(ixe)   pd.absorption(ixe)  sigma_cr    0;        
    d_cr+d_pd+d_ptba+d_ps   cr.dispersion(ixe)   cr.absorption(ixe)   sigma_si    0;    
    NaN                     si.dispersion(ixe)  si.absorption(ixe)   NaN         0];
[edp_0,dispersion,~,~] = edprofile(rawedp,nSlice);         % electron density profile without gold

% --- profile without gold
edp = edp_0;


% --- calculate volume fraction of each element (the intergration over z
% must be the effective total thickness)
% For gold:
ind_z = find(edp_0(:,1)>=0 & edp_0(:,1)<=d_ps+d_ptba);
z_mapped_from_index = ((1:nx)' - 1)/(nx-1)*(d_ps+d_ptba);
p_au = zeros(size(edp_0,1),1);
p_au(ind_z) = interp1(z_mapped_from_index,p_au0,edp_0(ind_z,1));
p_au = p_au/trapz(edp_0(:,1),p_au) * d_au;
% For other elements
p_ps = dispersion(:,2)/ps.dispersion(ixe);
p_ptba = dispersion(:,3)/ptba.dispersion(ixe);
p_pd = dispersion(:,4)/pd.dispersion(ixe);
p_cr = dispersion(:,5)/cr.dispersion(ixe);
p_si = dispersion(:,6)/si.dispersion(ixe);

p = [p_au,p_ps,p_ptba,p_pd,p_cr,p_si];      


% --- caculate EDP
edp(:,2) = edp_0(:,2) + (au.dispersion(ixe)-edp_0(:,2)).*p_au;
edp(:,3) = edp_0(:,3) + (au.absorption(ixe)-edp_0(:,3)).*p_au;

lambda = 12.3984171668278/xenergy(ixe);
if ixe == 3     % calculate for elastic energy only
    re = 2.814e-15;
    rho = [
        d_ps+d_ptba-edp(:,1), ...
        edp(:,2)*2*pi/(lambda*1e-10)^2/re/1e30, ...
        p_au/d_au];
end
if plotEDPFlag == 1 & ixe == 3 % for elastic energy only
  %  re = 2.814e-15;
  %  rho = edp(:,2)*2*pi/(lambda*1e-10)^2/re/1e30;
    figure
    subplot(2,1,1);
   % plot(d_ps+d_ptba-edp(:,1),rho,'-','linewidth',2);
   plot(rho(:,1),rho(:,2),'-','LineWidth',2);
    xlabel('Heigth above mirror (A)');
    ylabel('Electron density (A^{-3})');
    set(gca,'xminorTick','on');
    subplot(2,1,2);
   % plot(d_ps+d_ptba-edp(:,1),p_au/d_au,'-','linewidth',2);
   plot(rho(:,1),rho(:,3),'-','LineWidth',2);   
    xlabel('Heigth above mirror (A)');
    ylabel('\phi');
    set(gca,'xminorTick','on');
    set(gca,'xlim',[0,600]);
   
    set(gcf,'Position',[397   515   400   250]);

    plotEDPFlag = 0;
   
    
%     ddd = d_ps+d_ptba-edp(:,1);
%     ppp = p_au/d_au;
%     assignin('base','ddd',ddd);
%     assignin('base','ppp',ppp);
%     assignin('base','rho',rho);
%     assignin('base','edp',edp);
end

% --- Electic field
[EF,~,~,~] = efield(alphai',edp,lambda);
%  figure
%  imagesc(alphai,edp(:,1),abs(EF).^2);
%keyboard
% Y_tmp = p_au*ones(1,length(alphai)) .*abs(EF).^2;
% Y = NaN*ones(length(alphai),1);
% for ii=1:length(alphai)
%     Y(ii) = trapz(edp(:,1),Y_tmp(:,ii));
% end

