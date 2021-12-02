% X-ray Reflectivity
% HMCMC with Non-trivial diagonal mass matrix M

%%%%%%%%%%%%%%%%%%%
% do NOT fit nsets of multiple energies. Must be same energy
% fit one set only
%%%%%%%%%%%%%%%%%%%%%

clear all
feature('numcores');

% --- needed for the computation of refraction index
addpath(fullfile(pwd,'xrayrefraction'),fullfile(pwd,'xrayrefraction','AtomicScatteringFactor'));


%% define global parameters 
global fitFlag x2 lbub
global lambda
global plotEDPFlag
global nsets nSlice nlayer
global edp_0 rawedp
global interface_profile
global conv_method
global binedp_delta_cutoff      % delta cutoff for edp binning

%% some flags ans settings
% --- about density profile
binedp_delta_cutoff = 1e-15;
conv_method = 1;            % 0/1/2/3: no resolution/vectorization/forloop/conv: 3 is the fastest (only 1 is implemented for now)
interface_profile   = 0;        % 0/1/2: error / hyperbolic tangent/cumulative skewed norm
nSlice = 150;
plotEDPFlag = 1;                % 0/1: plot edp in the call function
                                % must disable during MCMC
% --- about beam
xe = 20.017;                % @@@@@@ x-ray energy in KeV     
lambda      = 12.3984171668278/xe;     % wavelength (A)
                               
% --- about data
weightFlag = 1;       % 0/1 use weight for intensity
typeRes = 1;          % 1/2: gaussian/box

% --- about HMCMC
updatePlotFlag = 1;    % 0/1 update plot during MCMC
dualAvgFlag = 1;        % 0/1 use dual averaging
nutsFlag = 0;         % 0/1 use NUTS 

% --- plot settings for multiple curves
plotshift = 0.1;    % shift factor for multiple data sets when plotting

%% Load Data
fpath = pwd;
flist = {
   'A2DOPC.dat'; 
    };

% flist = dir('*.ABS');
% flist = {flist.name}';

nsets = length(flist); 
alphai_set = cell(1,nsets);  % collect data sets
Iset = cell(1,nsets); 
Ierr_set = cell(1,nsets); 
Iweight_set = cell(1,nsets);
sigma_q_set = cell(1,nsets); 

% --- points to delete
ind_bad = {
    []         % uncomment if no bad point to remove for this set
%     [30 34 49 54]
   };
% ind_bad = cell(length(flist),1);

for iset=1:nsets
    % --- for DAT format
    if strcmpi(flist{iset}(end-2:end),'ABS')
        s = importdata(fullfile(fpath,flist{iset}),' ',5);
        s = s.data;
    elseif strcmpi(flist{iset}(end-2:end),'DAT')
        s = load(fullfile(fpath,flist{iset}));     % for DAT format
    end

    % -- remove bad sigma q
    s(ind_bad{iset},:) = [];

    [~,ia,~] = unique(s(:,1),'rows','first');       % remove points of identical qz
    s = s(ia,:);
    ind = s(:,2)<=0 | s(:,3)< 0 | s(:,1)>2.0;   % % remove negative I and too noise data
%     ind = s(:,2)<=0 | s(:,3)< 0;
    s(ind,:) = [];
    
    % --- reduce data density
   fitRange = 1:1:size(s,1);
%     fitRange = [1:60,61:2:size(s,1)];
    s = s(fitRange,:);
    
    % convert from qz to alpha
    alphai00 = asind(s(:,1)*lambda/4/pi);     % inicident angle values (deg)
    
    I00 = s(:,2);        % intensity
    I00_err = s(:,3);       % weight
    
    % --- normalize intensity
    I00_max = max(I00);
    I00 = I00/I00_max;
    I00_err = I00_err/I00_max;
    
    % --- add weight
    if weightFlag == 0
        I00_weight = ones(length(q00),1);
    elseif weightFlag == 1
        I00_weight = I00./I00_err;
        I00_weight = I00_weight/sum(I00_weight)*length(alphai00);     % normalize so the sum is the number of points
    end
        
    % --- collect
    alphai_set{iset} = alphai00;
    Iset{iset} = I00;
    Ierr_set{iset} = I00_err;
    Iweight_set{iset} = I00_weight;
end

% % --- view raw data
figure
hold on; box on;
for iset=1:nsets
    errorbar(alphai_set{iset},Iset{iset}*plotshift^(iset-1),Ierr_set{iset}*plotshift^(iset-1),'o-');
end
set(gca,'yscale','log');
legend(flist,'Location','best','Interpreter','none');

%% initialize and prepare parameters
Inorm = 1*ones(nsets,1);      % intensity normalization
Ibk =   0*ones(nsets,1);     % constant background    

% --- resolution 
resdawave = 0.0;                        % relative dispersion of wavelength; keep it zero for monochromatic beam
% resconstant = 4.05869172548923e-05;      % constant resoltuion for qz in unit of A^-1; 
        % fitting will be much fast if brushh values are set zero and they are
        % not fitted (fit flag is zero), i.e., no resolution convolution
resconstant = 0;      % constant resoltuion for qz in unit of A^-1; 

        
foot_angle =  0.0001;                         % footprint angle, angle where the full beam matches the sample size
alphai_offset = 0;      % angle offset to account for incident angle misalignment. You need to fit or change this value in
            % order to match the critical angle of the substrate; otherwise
            % keep it zero.
              
massdensity_Si = 2.33;                  % g/cm^3 mass density for initial calculation of delta and beta
massdensity_SiOx = 2.2;
massdensity_H2O = 0.9973;
si = refrac('Si',xe,massdensity_Si);
siox = refrac('SiO2',xe,massdensity_SiOx);
h2o = refrac('H2O',xe,massdensity_H2O);   % buffer

massdensity_DOPC = 0.95;    % guessed value
DOPC = refrac('C44HO83P1N1',xe,massdensity_DOPC);       
           
delta_Si = si.dispersion;           % get delta; you can define your own initial value 
beta_Si  = si.absorption;            % get beta;  you can define your own initial value

sigma_Si = 2.5;
asymmetry_Si = 0;
CDF_Si = 0;
% %sigma_Si_SiOx       =  4.34513427878109;  % roughness of Si/SiOx  
% asymmetry_Si_SiOx   = 0;              % asymmetry factor (+ toward z+)
% CDF_Si_SiOx         = 0;

delta_buffer = h2o.dispersion;
beta_buffer =  h2o.absorption;

% --- layered edp from top (buffer but excluding buffer) to bottom (Si substrate, but exclusing substrate)
% [thickness, delta, beta, sigma, asymmetry] (sigma is the roughness of the
% top interface of the layer)
layeredp = [
12          6.5e-07         DOPC.absorption    4.5
24          4.5e-07         DOPC.absorption    4
9           9e-07           DOPC.absorption    8.5
];

nlayer = size(layeredp,1);
fitFlag_layeredp = [...
    ones(nlayer,1) ones(nlayer,1) zeros(nlayer,1) ones(nlayer,1)
    ];
lb_layeredp = [ ...
    eps*ones(nlayer,1)        1e-9*ones(nlayer,1)   1e-11*ones(nlayer,1)        eps*ones(nlayer,1)        %-10*ones(nlayer,1)          -eps*ones(nlayer,1)
    ];     % lower boundary 
ub_layeredp = [ ...
    40*ones(nlayer,1)       5e-6*ones(nlayer,1)     1e-8*ones(nlayer,1)         50*ones(nlayer,1)       %eps*ones(nlayer,1)          50*ones(nlayer,1)
    ];     % lower boundary 


% --- [start, lb, ub, fitflag]
fitparam = [...
    Inorm               ones(nsets,1)*0.0001    ones(nsets,1)*10       ones(nsets,1)*1        
    Ibk                 -ones(nsets,1)*1e-16     ones(nsets,1)*1e-6     ones(nsets,1)*0
    resdawave           -eps                    1e-3                    0
    resconstant         -eps                    0.01                    0
    foot_angle          0                       5                       0
    alphai_offset       -0.01                   0.01                    0
    
	delta_Si            0                       1e-5                    0
    beta_Si             0                       1e-5                    0
    sigma_Si            0                       20                      1
    asymmetry_Si   -10                     eps                     0
   CDF_Si          -eps                    50                      0
    delta_buffer        1e-10                   1e-5                    0
    beta_buffer         1e-10                   1e-9                   0
    
    layeredp(:)         lb_layeredp(:)          ub_layeredp(:)          fitFlag_layeredp(:)
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
    error('start value cannot be at the bounds')
end
lbub = [lb,ub];
lbub(:,3) = NaN;  % Last column indicate the transformation method; 
lbub(isfinite(lb) & isfinite(ub),3)     = 1;    % (a, b)
lbub(isfinite(lb) & ~isfinite(ub),3)    = 2;    % (a, Inf)
lbub(~isfinite(lb) & isfinite(ub),3)    = 3;    % (-Inf,b);
lbub(~isfinite(lb) & ~isfinite(ub),3)   = 4;    % (-Inf,Inf);
lbub(fitFlag==0,3) = 5;                         % no transformation
xtrans = transformx(x,lbub);    % transform 

% --- assign to-be and not-to-be fitted parameters
x1 = xtrans(fitFlag==1);
x2 = xtrans(fitFlag==0);


%% prevew data
plotEDPFlag = 1;
tic
[I_cal_set,qset] = fcn_xray_ref_hmcmc(x1,alphai_set);     
toc
plotEDPFlag = 0;

% % --- plot
figure
hold on; box on;
% plot data
for iset = 1:nsets
    errorbar(alphai_set{iset},Iset{iset}.*qset{iset}.^4*plotshift^(iset-1),Ierr_set{iset}.*qset{iset}.^4*plotshift^(iset-1),'o');
end
% plot calculation
for iset = 1:nsets
    plot(alphai_set{iset},I_cal_set{iset}.*qset{iset}.^4*plotshift^(iset-1),'k-','linewidth',1.5);
end
set(gca,'yscale','log');
legend(flist,'Location','best','Interpreter','none')

%% Initialize parameters 
plotEDPFlag = 0;    % disable online plot of EDP

% ---  initialize variance of residual to guessed values
Sigmavi = nan(1,nsets);
for iset=1:nsets   
    Sigmavi(iset) = var(log10(I_cal_set{iset}) - log10(Iset{iset}));
end

% --- initial U and dU values
xp_init = [x1(:)',log(Sigmavi)];
% xp_init = [-4.12776380423215,-17.6106277100161,-2.12284409489526,-1.96824145580311,-0.184713254832867,-1.26989898621000,-3.43606201206975,0.610705031822836,-1.85020262250683,-3.95000471741201];

%% Define potential, kinetic energies and their differentials
% --- define potention  and differential
U = @(X)potentialE(X,nsets,alphai_set,Iset,Iweight_set);
opts = admOptions('i', 1);  % only differentiate the 1st input argument (not param)
% dU = @(X)admDiffFor(@potentialE,1,X,nsets,alphai_set,Iset,Iweight_set,opts);   % Forward mode
% dU = @(X)admDiffVFor(@potentialE,1,X,nsets,alphai_set,Iset,Iweight_set,opts);     % Vectorized forward mode
% dU = @(X)admDiffRev(@potentialE,1,X,nsets,alphai_set,Iset,Iweight_set,opts);    % Reverse mode
dU = @(X)admDiffFD(@potentialE,1,X,nsets,alphai_set,Iset,Iweight_set,opts);     % Numerical finite difference method (not an auto diff method)

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

%% epsilon, and define mass matrix
epsilon = 0.01;  % Use as larger epsilon as possible

% Identity mass matrix
% M = eye(length(x1)+nsets);  % default

% non-trivial mass matrix determined by dPotE
M = (abs(dPotE_init)); 
M = diag(M);

%% define kinetic energy and differential
K = @(p)p*inv(M)*p'/2;
dK = @(p)p*inv(M);          % differation of kinetic energy w.r.t momentum


%% --- find good initial value for epsilon (may skip this part)
% pp_init = mvnrnd(zeros(size(xp_init)),M);
% logp_init = -PotE_init - K(pp_init);
% epsilon_init = 0.1;   % large start value is preferred than small start value
% tic
% [epsilon,counter] = heuristic4epsilon(xp_init,pp_init,dPotE_init,epsilon_init,U,K,dU,dK,logp_init);
% fprintf('Inialized epsilon=%f, after %d iterations of heuristic done in %f sec\n',epsilon,counter,toc);
% 
% % --- force epsilon if a different inital value is desired
% epsilon = 0.01;

%% --- dual averaging settings
if dualAvgFlag == 1
    delta = 0.8;   % optimal acceptance probability
    gamma = 0.05;
    t0 = 10;
    kappa = 0.75;
    Madap = Inf;    % Inf is preferred especially in presence of multiple modes in the probablility distribution
    H_init = 0;
    Hm = H_init;
    logepsilon_init = log(10*epsilon);
    logepsilon_mean = log(epsilon);
end

%% --- NUTS or simple HMCMC settings
if nutsFlag == 1  
    delta_max = 500;        % stopping rule (equation 3), large enough so not to interfere
else
%     lambda_target = 0.1;    % target total number simulation length
    nstep_target = 50;      %  number of steps
    lambda_target = nstep_target*epsilon;    % target total number simulation length
end


%% prepare HMCMC
% --- assign initial values
xp = xp_init;
dPotE = dPotE_init;
logp_xp = -PotE_init;

% --- start MCMC
nmcmc = 3000;   % # of iterations
param_mcmc = nan(nmcmc,length(xp));    
epsilon_mcmc = nan(nmcmc,1);
dPotE_mcmc = nan(nmcmc,length(xp));

param_mcmc0 = xp;
epslion_mcmc0 = epsilon;
dPotE_mcmc0 = dPotE_init;

%% Start HMCMC
warning('off','MATLAB:Figure:FigureSavedToMATFile');
tic
for ii=1:nmcmc    
    % --- move one step in HMCMC
    trial_hmcmc = 0;
    done_hmcmc = 0;
    while done_hmcmc == 0  % to prevent ill random pp values
        try
            if nutsFlag == 0    % use simple HMCMC
                [xp,logp_xp,dPotE,logr] = hmcmcwalk(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,lambda_target);                
            elseif nutsFlag == 1 && dualAvgFlag == 1    % use NUTS with dual averaging for epsilon
                [xp,logp_xp,dPotE,alpha,nalpha] = nuts(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max);
            elseif nutsFlag == 1 && dualAvgFlag == 0    % use NUTS without dual averaging
                [xp,logp_xp,dPotE] = nuts_noDA(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max);
            end
            done_hmcmc = 1;
            fprintf('%d of %d iterations (%.1f sec) with epsilon=%.8f: ',ii,nmcmc,toc,epsilon);
            fprintf(['[',repmat('%.3f ',1,6),'\b]\n'],xp([1:5,end]));     % only display 6 parameters
        catch
            trial_hmcmc = trial_hmcmc + 1;                        
            fprintf('For the %dth iteration: re-try %d time(s)\n', ii,trial_hmcmc);
        end
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
    
    % --- collect result
    param_mcmc(ii,:) = xp;
    epsilon_mcmc(ii) = epsilon;
    dPotE_mcmc(ii,:) = dPotE;
    
    % --- update plot
    if updatePlotFlag == 1
        [I_cal_set,qset] = fcn_xray_ref_hmcmc(xp(1:end-nsets),alphai_set);     
        if ii==1
            qset1 = qset;
            hfig = figure('position',[200,10,600,700]);
            % --- plot curve
            hax1 = subplot(11,1,1:4);
            hold on; box on;
            for iset = 1:nsets
                errorbar(alphai_set{iset},Iset{iset}.*qset1{iset}.^4*plotshift^(iset-1),Ierr_set{iset}.*qset1{iset}.^4*plotshift^(iset-1),'o');
            end
            % plot calculation
           %  hline1 = [];
            for iset = 1:nsets
                hline1(iset) = plot(alphai_set{iset},I_cal_set{iset}.*qset1{iset}.^4*plotshift^(iset-1),'k-','linewidth',1.5);
            end
            legend(flist,'Location','best','Interpreter','none')
            set(gca,'yscale','log');
            set(gca,'XMinorTick','on','YMinorTick','on');            
            xlabel('angle (deg)');
            ylabel('R\timesq^4');
            % --- plot EDP
            hax2 = subplot(11,1,6:7);            
            hold on; box on;
            hline2 = plot(edp_0(:,1)-rawedp(end-1,1),edp_0(:,2),'r-','linewidth',1);
            ylabel('\delta');      
            set(gca,'xticklabel',[]);
            set(gca,'XMinorTick','on','YMinorTick','on');
            % --- plot EDP beta
            hax3 = subplot(11,1,8:9);            
            hold on; box on;            
            hline3 = plot(edp_0(:,1)-rawedp(end-1,1),edp_0(:,3),'r-','linewidth',1);             
            ylabel('\beta')    
            xlabel(['depth (',char(197),')']);    
            set(gca,'XMinorTick','on','YMinorTick','on'); 
            linkaxes([hax2,hax3],'x');
%            set(gca,'xlim',[-500,1500]);
            % --- plot epsilon
            hax4 =subplot(11,1,11);
            hline4 = plot(param_mcmc(:,end));
            set(gca,'xlim',[0,nmcmc]);
            xlabel('iteration');
            ylabel('log(\sigma^2)');
%             set(gca,'yscale','log');
            htitle = sgtitle(sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc));
        else
            for iset = 1:nsets
                hline1(iset).YData = I_cal_set{iset}.*qset1{iset}.^4*plotshift^(iset-1);
            end
            set(hline2,'Color','b','linewidth',0.5); % change previous color to blue and linewidth to thin            
            set(hline3,'Color','b','linewidth',0.5); % change previous color to blue and linewidth to thin                        
            hline2 = plot(edp_0(:,1)-rawedp(end-1,1),edp_0(:,2),'r-','linewidth',1,'Parent',hax2);  % add new line
            hline3 = plot(edp_0(:,1)-rawedp(end-1,1),edp_0(:,3),'r-','linewidth',1,'Parent',hax3);  % add new line            
            set(hline4,'XData',1:nmcmc,'YData',param_mcmc(:,end));             
            htitle.String = sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc);
        end
        pause(0.01);        
    end
    
    save tmp_LB_HMCMC.mat        % iteratively save result
end
total_time =toc;

% --- save result
if nutsFlag == 1 && dualAvgFlag == 1
    save fresult_LB_NUTS_DA_v1.mat
elseif nutsFlag == 1 && dualAvgFlag == 0
    save fresult_LB_NUTS_NoDA_v1.mat    
elseif nutsFlag == 0 && dualAvgFlag == 1    
    save fresult_LB_Simple_DA_v1.mat        
elseif nutsFlag == 0 && dualAvgFlag == 0
    save fresult_LB_Simple_NoDA_v1.mat            
end

