clear all

global fitFlag x2
global lbub
global sdd xenergy nSlice alpha0 plotEDPFlag
global N nx B
global rho
global si cr pd ptba ps au  % properties of materials
global beta1 beta2

nSlice = 650;

% --- needed for the computation of refraction index
addpath(fullfile(pwd,'xrayrefraction'),fullfile(pwd,'xrayrefraction','AtomicScatteringFactor'));

%% some flags ans settings
% --- plot EDP (electron density profile during iterations)
plotEDPFlag = 0;

% --- about data
weightFlag = 0;       % 0/1 use weight for intensity

% --- about HMCMC
updatePlotFlag = 1;    % 0/1 update plot during MCMC
dualAvgFlag = 0;       % 0/1 use dual averaging
nutsFlag = 0;          % 0/1 use NUTS 

%% Load image data
% --- flat field
FF = imread('FF_p1m0108_E11p9_T9p9_vrf_m0p3.tif');  % flat field
mask = imread('pilatus1mf_mask_20110110.tif');      % mask

fpath = pwd;
fname = 'img.tif';

f = fullfile(fpath,fname);
img = double(imread(f));

img = img.*FF;
img(~mask) = NaN;

% --- rotate image to account for tilt of the camera
angle = -atand(5/981);    % title angle of the camera 
img2 = double(imrotate(img,angle,'nearest','crop'));

% figure
% imagesc(img,[0,20]);

%% performce linecut (collapse horizontally)
d_ypixel = 0.172/cosd(abs(angle));
sumRange = [10:480,510:980];
Is = mean(img2(430:950,sumRange),2);
Is = flipud(Is);
vzs = (0:length(Is)-1)';
vzs = vzs*d_ypixel;

ind_nan = isnan(Is);
Is(ind_nan) = [];
vzs(ind_nan) = [];

% --- get background level from low (negative) angle data
[~,vzs_max_ind] = max(Is);
ind_bk = vzs_max_ind-60:vzs_max_ind-40;
Ibk_mean = mean(Is(ind_bk));

% --- cut off low angle data
fitRange = [80:1:160,161:2:210,213:4:350]; %length(Is);

I = Is(fitRange);
vz = vzs(fitRange);

% figure
% hold on;
% plot(vzs,Is,'*-');
% plot(vz,I,'o-');
% 
% xlabel('height (mm)');
% ylabel('fluoresence intensity');

% return

%% add weight to higher order peaks
if weightFlag == 0
    weight = ones(size(vz));
else
    error('no I_error column');
end

%% initialize and prepare parameters
sdd = 3180;             % sample detector distance
xe0 = 12.1;             % inicident elastic energy
xe1 = 9.70484684684685;      % weighted mean of L_alpha(1,2)
xe2 = 11.5847;          % weighted L_beta2,15 (no L_beta1 because the incident energy is not high enough to excite L2 shell)
alpha0 = 0.13;      % incident angle
xenergy = [xe1,xe2,xe0];    % for low to high energy; last one is elastic energy


sdd_air = sdd-1388;     % air path (vacuum length is 1388)
air = refrac('N78O21Ar',[xe1,xe2],0.0011839);   
abs_ratio = exp(-sdd_air*0.1*diff(1./air.attLength) );  % absortion ratio of air
quantum_ratio = 0.92/0.8;   % account for quantum efficienc for the two energies; from Pilatus spec sheet
yield_ratio = 23/(100+11);  % relative emission yield ration of L_beta(2,15)/L_alpha(1,2); from emission table
xe_ratio = yield_ratio*abs_ratio*quantum_ratio;    % total ratio of with air absorption and detector quantum efficiency included
elastic_ratio =    0.04;   % contribution of elastic scattering with respect to fluorescent of xe1

Inorm   =   0.19;
Ibk =     Ibk_mean;  
alphai_res = eps(eps);          % resolution for exit angle alpha;  (no fitting)
vz0  =    15.7;     % offset for exit angle (to find sample horizon) in terms of pixel

% --- parameter about samples
sigma_si  = 10.7426539747676;   % silicon substrate roughness
rho_si = 2.33;                  % electron density of silicon

% Cr layer
d_cr =     59.9955682709922;   
sigma_cr = 11.2238440521098;
rho_cr = 7.19;

% about Pd layer
d_pd =     633;
sigma_pd =  16.5211447363043 ;
rho_pd = 12.56; 

% about PTBA layer
d_ptba =       630;
sigma_ptba = 1.8;  
rho_ptba = 1.109;

% about PS capping
d_ps =         180;
sigma_ps =       35;
rho_ps = 1.04;

% about gold layer
d_au =     5.75;

% --- cubic b-spline settings
N = 30;         % # of linear splines
nx = 800;       % # of points in each spline

B = cbsbase(N+1,nx);          % create N+2 splines
B = B(:,2:end-1);           % use the middle N splines (i.e. discard the 1st and last one because both ends are D2O)

% - uniform shape 
a = eps+ones(N,1)/N*0.01;

% regularization
beta1 = 1e2;
beta2 = 1e4;

rho_au = 19.3;  % gold density

si      = refrac('Si',xenergy,rho_si);
cr      = refrac('Cr',xenergy,rho_cr);
pd      = refrac('Pd',xenergy,rho_pd);
ptba    = refrac('C7H12O2',xenergy,rho_ptba);
ps      = refrac('C8H8',xenergy,rho_ps);
au      = refrac('Au',xenergy,rho_au);

% --- construct fitting parameters
% [initial, lower bound, upper bound, fit(1/0)?]
fitparam = [...
    Inorm           0.05        1           1
    Ibk             0           1           0
    alphai_res      0           0.01        0
    vz0             10          20          1
    
    sigma_si        1           20          0   
        
    d_cr            1           100         0
    sigma_cr        1           50          0
    
    d_pd            100         1000        0
    sigma_pd        1           50          0        
   
    d_ptba          400         800         1
    sigma_ptba      0.1         100         1         

    d_ps            100         400         1 
    sigma_ps        0.1         100         1         

    d_au            1           10          1
    
    xe_ratio        0           1           0        
    elastic_ratio   0           0.1         1
    
    a               0*ones(N,1) 1e-2*ones(N,1) ones(N,1)
];

x       = fitparam(:,1);
lb      = fitparam(:,2);
ub      = fitparam(:,3);
fitFlag = fitparam(:,4);

% --- transform parameter so it become unbounded
% for any unbounded parameters, better rescale its initial value to 1. This
% needs to be done manually.
if any(lb-ub>=0)  % ub must be larger than lb
    error('incorrect bound');
end
if any(x<=lb | x>=ub)
    error('start value cannot be at or outside the bounds')
end
lbub = [lb,ub];
lbub(:,3) = NaN;  % Last column indicate the transformation method; 
lbub(isfinite(lb) & isfinite(ub),3)     = 1;    % (a, b)
lbub(isfinite(lb) & ~isfinite(ub),3)    = 2;    % (a, Inf)
lbub(~isfinite(lb) & isfinite(ub),3)    = 3;    % (-Inf,b);
lbub(~isfinite(lb) & ~isfinite(ub),3)   = 4;    % (-Inf,Inf);
lbub(fitFlag==0,3) = 5;                         % no transformation
xtrans = transformx(x,lbub);

% --- assign to-be and not-to-be fitted parameters
x1 = xtrans(fitFlag==1);
x2 = xtrans(fitFlag==0);

%% preview data
tic
[I_cal,Y1,Y2,Y_elastic] = fcn_flu_spline_hmcmc(x1,vz);     % preview initial conditions before fitting
toc

figure
subplot(5,1,1:3);
hold on; box on;
plot(atand((vz-vz0)/sdd),(I-Ibk)/Inorm,'o','markersize',6);
plot(atand((vz-vz0)/sdd),(I_cal-Ibk)/Inorm,'-','linewidth',2);
plot(atand((vz-vz0)/sdd),Y1,'-','linewidth',1);
plot(atand((vz-vz0)/sdd),Y2,'-','linewidth',1);
plot(atand((vz-vz0)/sdd),Y_elastic,'-','linewidth',1);
ylabel('Au L\alpha and L\beta_1 Fluorescence');
xlabel('Exit Angle (deg)');
set(gca,'xminortick','on')
set(gca,'xlim',[0,1]);
legend('data','fit','9.7kev','11.58kev','Elastic')
subplot(5,1,4:5);
plot(rho(:,1),rho(:,3),'-','LineWidth',2);
xlabel('Heigth above mirror (A)');
ylabel('\phi');
set(gca,'xminorTick','on');
set(gca,'xlim',[0,800]);


% return


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Start HMCMC setup 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize variance of residual to guessed values
Sigmavi = var(log10(I_cal) - log10(I));

% --- initial U and dU values
xp_init = [x1(:)',log(Sigmavi)];

%% Define potential, kinetic energies and their differentials
sspFlag = true(1,length(xp_init));
XFix = [];

% --- define potential and differential
% U = @(X)potentialE(X,vz,I,weight,sspFlag,XFix);
opts = admOptions('i', 1);  % only differentiate the 1st input argument (not param)
% dU = @(X)admDiffFor(@potentialE,1,X,vz,I,weight,opts);    % Forward mode
% dU = @(X)admDiffVFor(@potentialE,1,X,vz,I,weight,opts);     % Vectorized forward mode
% dU = @(X)admDiffRev(@potentialE,1,X,vz,I,weight,opts);    % Reverse mode
dU = @(X)admDiffFD(@potentialE,1,X,vz,I,weight,sspFlag,XFix,opts);     % Numerical finite difference method (not an auto diff method)

%% Initialize parameters 
plotEDPFlag = 0;    % force EDP update plot to stop
tic
[dPotE_init,PotE_init] = dU(xp_init);
fprintf('One differetial takes %.5f sec\n',toc);

% --- plot potential differentials
figure
plot(log10(abs(dPotE_init)+eps),'o-')
xlabel('index');
ylabel('log10(|dU/dx|)')

%% determine subspaces and epsilon values 
 nssp = 1;   % use only one subspaces
 sspFlag = true(nssp,length(dPotE_init));
 epsilon = 0.05;
 %  sspFlag(1,[1:3,5,7,end]) = true;          % subspace for Inorm, vz0, d_ptba, d_au, and variance
%  sspFlag(2,:) = ~sspFlag(1,:);       % subspace for other parameters
%  epsilon = [0.005,0.05];   % small epsilon for large differentials; and large epsilon for small differentials

% --- define Mass matrix values
iM_list = nan(1,length(dPotE_init));
iM_flag = false(nssp,length(dPotE_init));        % 0 for parameters with small gradients, 1 for parameters with large parameters
iM_flag([1:3,5,7,end]) = true;
iM_list(iM_flag) = 400;         % 
iM_list(~iM_flag) = 1;          % use large variance values for paremters with large parameters, i.e. larger momentum values
iM = diag(iM_list);


%% --- dual averaging settings
if dualAvgFlag == 1
    delta = 0.8;   % optimal acceptance probability
    gamma = 0.05;
    t0 = 10;
    kappa = 0.75;
    Madap = Inf;    % Inf is preferred especially in presence of multiple modes in the probablility distribution
    H_init = 0;
    Hm = ones(1,nssp)*H_init;
    logepsilon_init = log(10*epsilon);
    logepsilon_mean = log(epsilon);
end

%% --- NUTS or simple HMCMC settings
if nutsFlag == 1  
    delta_max = 500;        % stopping rule (equation 3), large enough so not to interfere
else
    nstep_target = 10;
    lambda_target = nstep_target*epsilon;    % target total number simulation length
end

%% final preparation for HMCMC
% --- assign initial values
xp = xp_init;
dPotE = dPotE_init;
logp_xp = -PotE_init;

% --- start MCMC 
nmcmc = 2000;       % total # of iterations
param_mcmc = nan(nmcmc,length(xp));     % to sava param, epsilon, and potential gradient
epsilon_mcmc = nan(nmcmc,nssp);
dPotE_mcmc = nan(nmcmc,length(xp));

param_mcmc0 = xp;  
epslion_mcmc0 = epsilon;
dPotE_mcmc0 = dPotE_init;

%% start iteration
warning('off','MATLAB:Figure:FigureSavedToMATFile');
tic
for ii=1:nmcmc    
    
    fprintf('%d of %d iterations:\n', ii,nmcmc);
    
    for issp = 1:nssp
        sspFlag_active = sspFlag(issp,:);        
        xp_active = xp(sspFlag_active);
        xp_fixed = xp(~sspFlag_active);
        
        
        % --- define U, dU, K, dK
        iU  = @(X)potentialE(X,vz,I,weight,sspFlag_active,xp_fixed);
        idU = @(X)admDiffFD(@potentialE,1,X,vz,I,weight,sspFlag_active,xp_fixed,opts);     
     %   iM  = eye(nnz( sspFlag_active ));  % default
        iK  = @(p)p*inv(iM)*p'/2 + log(det(iM));    % Euclidean-Gaussian kinetic energy
        idK = @(p)p*inv(iM);
        
        % --- update dPotE
        if nssp ~= 1     % no need to update if only one subspace
            [dPotE,PotE] = idU(xp_active);
            logp_xp = -PotE;
        end

        % --- move one step in HMCMC
        trial_hmcmc = 0;
        done_hmcmc = 0;
        while done_hmcmc == 0  % to prevent ill random pp values
            try
                if nutsFlag == 0    % use simple HMCMC
                    [xp_active,logp_xp,dPotE,logr(issp)] = hmcmcwalk(xp_active,logp_xp,dPotE,iU,iK,idU,idK,iM,epsilon(issp),lambda_target(issp));
                elseif nutsFlag == 1 && dualAvgFlag == 1    % use NUTS with dual averaging for epsilon
                    [xp_active,logp_xp,dPotE,alpha(issp),nalpha(issp)] = nuts(xp_active,logp_xp,dPotE,iU,iK,idU,idK,iM,epsilon(issp),delta_max);
                elseif nutsFlag == 1 && dualAvgFlag == 0    % use NUTS without dual averaging
                    [xp_active,logp_xp,dPotE] = nuts_noDA(xp_active,logp_xp,dPotE,iU,iK,idU,idK,iM,epsilon(issp),delta_max);
                end
                done_hmcmc = 1;
                fprintf('      %d of %d subspace (%.1f sec) with epsilon=%.5f, xp: ',issp,nssp,toc,epsilon(issp));
                fprintf(['[',repmat('%.10f ',1,length(xp_active)),'\b]\n'],xp_active);
                fprintf(['                                                      dPotE: [',repmat('%.10f ',1,length(dPotE)),'\b]\n'],dPotE);                
            catch
                trial_hmcmc = trial_hmcmc + 1;
                fprintf('   For the %dth iteration and %dth subspace: re-try %d time(s)\n', ii,issp,trial_hmcmc);
            end
        end

        % --- adapt epsilon with dual averaging
        if dualAvgFlag == 1
            if ii<=Madap
                if nutsFlag == 1    % update Hm for nuts
                    Hm(issp) = (1-1/(ii+t0))*Hm(issp) + 1/(ii+t0)*(delta-alpha(issp)/nalpha(issp));
                else                % update Hm for simple HMCMC
                    Hm(issp) = (1-1/(ii+t0))*Hm(issp) + 1/(ii+t0)*(delta-min(1,exp(logr(issp))));
                end
                logepsilon(issp) = logepsilon_init(issp) - sqrt(ii)/gamma*Hm(issp);
                logepsilon_mean(issp) = ii^(-kappa)*logepsilon(issp) +(1-ii^(-kappa))*logepsilon_mean(issp);  % moving average
                epsilon(issp) = exp(logepsilon(issp));
            else
                epsilon(issp) = exp(logepsilon_mean(issp));
            end
        end
        
        
        % assmelbe and collect xp, and dPotE_mcmc
        xp(sspFlag_active) = xp_active;
        param_mcmc(ii,:) = xp;
        dPotE_mcmc(ii,sspFlag_active) = dPotE;
        epsilon_mcmc(ii,issp) = epsilon(issp);
        
    end
    
    % --- adapt epsilon with dual averaging
    if dualAvgFlag == 1
        if ii<=Madap
            if nutsFlag == 1    % update Hm for nuts
                Hm = (1-1/(ii+t0))*Hm + 1/(ii+t0)*(delta-alpha/nalpha);
            else                % update Hm for simple HMCMC
                Hm = (1-1/(ii+t0))*Hm + 1/(ii+t0)*(delta-min(1,exp(logr)));
            end 
            logepsilon = logepsilon_init - sqrt(ii)/gamma*Hm;
            logepsilon_mean = ii^(-kappa)*logepsilon +(1-ii^(-kappa))*logepsilon_mean;  % moving average
            epsilon = exp(logepsilon);
        else
            epsilon = exp(logepsilon_mean);
        end
    end
    
    % --- update plot
    if updatePlotFlag == 1
        [I_cal,Y1,Y2,Y_elastic] = fcn_flu_spline_hmcmc(xp(1:end-1),vz);     % preview initial conditions before fitting
        
        if ii==1
            hfig = figure;
            hax1 = subplot(6,1,1:3);
            hold on; box on;
            plot(atand((vz-vz0)/sdd),(I-Ibk)/Inorm,'o','markersize',6);
            hline1a = plot(atand((vz-vz0)/sdd),(I_cal-Ibk)/Inorm,'-','linewidth',2);
            hline1b = plot(atand((vz-vz0)/sdd),Y1,'-','linewidth',1);
            hline1c = plot(atand((vz-vz0)/sdd),Y2,'-','linewidth',1);
            hline1d = plot(atand((vz-vz0)/sdd),Y_elastic,'-','linewidth',1);

            ylabel('Au L\alpha and L\beta_1 Fluorescence');
            xlabel('Exit Angle (deg)');
            set(gca,'xminortick','on')
            set(gca,'xlim',[0,1.2]);
            legend({'data','fit','9.7kev','11.58kev','Elastic'},'location','best')            
            hold off;
            
            hax2 = subplot(6,1,4:5);
            hold on; box on;
            hline2 = plot(rho(:,1),rho(:,3),'r-','LineWidth',1);
            xlabel('Heigth above mirror (A)');
            ylabel('\phi');
            set(gca,'xminorTick','on');
            set(gca,'xlim',[0,1000]);
            
            hax3 = subplot(6,1,6);
            hline3 = plot(epsilon_mcmc);
            set(gca,'xlim',[0,nmcmc]);
            xlabel('iteration');
            ylabel('\epsilon');
            set(gca,'yscale','log');
            htitle = sgtitle(sprintf('yan12 RT %d of %d: %.1f sec elapsed',ii,nmcmc,toc));            
        else
            hline1a.YData = (I_cal-Ibk)/Inorm;
            hline1b.YData = Y1;
            hline1c.YData = Y2;
            hline1d.YData = Y_elastic;
            set(hline2,'Color','b','linewidth',0.5); % change previous color to blue and linewidth to thin
            hline2 = plot(rho(:,1),rho(:,3),'r-','linewidth',1,'Parent',hax2);  % add new line            
            for ihline3 = 1:length(hline3)
                set(hline3(ihline3),'XData',1:nmcmc,'YData',epsilon_mcmc(:,ihline3));
            end            
            htitle.String = sprintf('yan12 RT %d of %d: %.1f sec elapsed',ii,nmcmc,toc);
        end
        pause(0.01);        
    end
    
    save tmp_randoma.mat        % iteratively save result
end
total_time =toc;

% --- save result
if nutsFlag == 1 && dualAvgFlag == 1
    save fresult_NUTS_DA_randoma_nontrivalmass_v1.mat
elseif nutsFlag == 1 && dualAvgFlag == 0
    save fresult_NUTS_NoDA_randoma_nontrivalmass_v1.mat    
elseif nutsFlag == 0 && dualAvgFlag == 1    
    save fresult_Simple_DA_randoma_nontrivalmass_v1.mat        
elseif nutsFlag == 0 && dualAvgFlag == 0
    save fresult_Simple_NoDA_randoma_nontrivalmass_v1.mat            
end


