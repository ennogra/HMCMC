clear

global fitFlag x2
global lbub

%% some flags ans settings
% --- about data
weightFlag = 1;       % 0/1 use weight for intensity

% --- about HMCMC
updatePlotFlag = 1;    % 0/1 update plot during MCMC
dualAvgFlag = 0;        % 0/1 use dual averaging
nutsFlag = 1;         % 0/1 use NUTS (must use unit Mass and cannot update Mass)

% force mass update to stop if NUTS is used
if nutsFlag == 1
    updateMassType = 0;
end

%% Load data
s = load('mrgSAXS_1percent_150nm.dat');

% reduce # of data points in log space
% [~,fitIndex] = ismember(...
%     interp1(s(:,1),s(:,1), logspace(log10(s(1,1)),log10(s(end,1)),133),'nearest'),...
%     s(:,1));
fitIndex = 1:5:size(s,1);
qdata =     s(fitIndex,1);
Idata = s(fitIndex,2);
Idata_err = s(fitIndex,3);

% normalize 
Idata_max = max(Idata);
Idata = Idata/Idata_max;
Idata_err = Idata_err/Idata_max;

if weightFlag == 0
    weight = ones(size(qdata));
elseif weightFlag == 1
    weight = Idata./Idata_err;
    weight = weight/sum(weight)*length(qdata);
end
% 
% figure
% errorbar(qdata,Idata,Idata_err,'o-');
% set(gca,'yscale','log');
% set(gca,'xscale','log');
% return
%% initialize and prepare parameters
Inorm = 1;
Ibk = 1e-5;
R = 600;
sigma_R = 20;
sigma_q = 1e-4;

% --- [start, lb, ub, fitflag]
fitparam = [...
    Inorm           0.1                 10                  1
    Ibk             0                   1e-3                1     
    R               100                 1000                 1
    sigma_R         1                   100                 1
    sigma_q         0                   1e-3                1
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
I_cal = fcn_saxs_sphere_hmcmc(x1,qdata);     % preview initial conditions before fitting
toc

% --- plot
figure
hold on;
errorbar(qdata,Idata,Idata_err,'o');
plot(qdata,I_cal,'-','linewidth',1.5);
hold off; box on;
set(gca,'xscale','log');
set(gca,'yscale','log');



%% Define potential, kinetic energies and their differentials
% --- define potential and differential
U = @(X)potentialE(X,qdata,Idata,weight);
opts = admOptions('i', 1);  % only differentiate the 1st input argument (not param)
% dU = @(X)admDiffFor(@potentialE,1,X,qdata,Idata,weight,opts);    % Forward mode
% dU = @(X)admDiffVFor(@potentialE,1,X,qdata,Idata,weight,opts);     % Vectorized forward mode
% dU = @(X)admDiffRev(@potentialE,1,X,qdata,Idata,weight,opts);    % Reverse mode
dU = @(X)admDiffFD(@potentialE,1,X,qdata,Idata,weight,opts);     % Numerical finite difference method (not an auto diff method)

%% Initialize parameters 
% ---  initialize variance of residual to guessed values
Sigmavi = var(log10(I_cal) - log10(Idata));

% --- initial U and dU values
xp_init = [x1(:)',log(Sigmavi)];
[dPotE_init,PotE_init] = dU(xp_init);

%% define M, K dK
% --- force epsilon if a different inital value is desired
epsilon = 0.05;

% M is estimated from previous run 
xp_cov = [0.000164021708211037,-0.000179455922662253,7.38702533447674e-05,-4.30634836927171e-06,0.000888335438219122,3.87183695892452e-05;-0.000179455922662253,0.0122934571185113,0.000134535806980337,-0.000107277461955319,0.00555197820936063,-0.000330074179482761;7.38702533447674e-05,0.000134535806980337,0.000129380499721589,-0.000320048414513545,0.00131327954980660,3.81042392878774e-05;-4.30634836927171e-06,-0.000107277461955319,-0.000320048414513545,0.00481127319900645,-0.00230447419487296,0.000197792440075726;0.000888335438219122,0.00555197820936063,0.00131327954980660,-0.00230447419487296,3.70161857895349,0.00649332775510657;3.87183695892452e-05,-0.000330074179482761,3.81042392878774e-05,0.000197792440075726,0.00649332775510657,0.0154936454118666];
M = inv(diag(diag(xp_cov)));

% % --- define kinetic energy and differential
K = @(p)p*inv(M)*p'/2;% + log(det(M));    % Euclidean-Gaussian kinetic energy
dK = @(p)p*inv(M);          % differation of kinetic energy w.r.t momentum

%% --- find good initial value for epsilon (may skip this part)
% pp_init = mvnrnd(zeros(size(xp_init)),M);
% logp_init = -PotE_init - K(pp_init);
% epsilon_init = 0.1;   % large start value is preferred than small start value
% tic
% [epsilon,counter] = heuristic4epsilon(xp_init,pp_init,dPotE_init,epsilon_init,U,K,dU,dK,logp_init);
% fprintf('Inialized epsilon=%f, after %d iterations of heuristic done in %f sec\n',epsilon,counter,toc);

% --- force epsilon if a different inital value is desired
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
%     lambda_target = 10*epsilon;%0.1;    % target total number simulation length
  %  lambda_target = 20*epsilon;%0.1;    % target total number simulation length    
    lambda_target = 1;
end

%% start HMCMC
% --- assign initial values
xp = xp_init;
dPotE = dPotE_init;
logp_xp = -PotE_init;

% --- start MCMC
nmcmc = 5000;
param_mcmc = nan(nmcmc,length(xp));
epsilon_mcmc = nan(nmcmc,1);
dPotE_mcmc = nan(nmcmc,length(xp));
pp_mcmc = nan(nmcmc,length(xp));

param_mcmc0 = xp_init;
epsilon_mcmc0 = epsilon;
dPotE_mcmc0 = dPotE;

%%
warning('off','MATLAB:Figure:FigureSavedToMATFile');
tic
for ii=1:nmcmc

   % --- move one step in HMCMC
    trial_hmcmc = 0;
    done_hmcmc = 0;
    while done_hmcmc == 0   % to prevent ill random pp values
        try
            if nutsFlag == 0    % use simple HMCMC
                [xp,logp_xp,dPotE,logr] = hmcmcwalk(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,lambda_target);                
            elseif nutsFlag == 1 && dualAvgFlag == 1    % use NUTS with dual averaging for epsilon
                [xp,logp_xp,dPotE,alpha,nalpha] = nuts(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max);
            elseif nutsFlag == 1 && dualAvgFlag == 0    % use NUTS without dual averaging
                [xp,logp_xp,dPotE] = nuts_noDA(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max);
            end
            done_hmcmc = 1;
            fprintf('%d of %d iterations (%.1f sec) with epsilon=%f: ',ii,nmcmc,toc,epsilon);
            fprintf(['[',repmat('%f ',1,length(xp)),'\b]\n'],xp);
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
        I_cal = fcn_saxs_sphere_hmcmc(xp(1:end-1),qdata);     % preview initial conditions before fitting
        if ii==1
            hfig = figure;
            subplot(3,1,1:2);
            hold on;
            errorbar(qdata,Idata,Idata_err,'o');
            hline = plot(qdata,I_cal,'-','linewidth',1.5);
            hold off; box on;
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            subplot(3,1,3);
            hline2 = plot(epsilon_mcmc);
            set(gca,'xlim',[0,nmcmc]);
            xlabel('iteration');
            ylabel('\epsilon');
            set(gca,'yscale','log');
            htitle = sgtitle(sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc));            
        end
        hline.YData = I_cal;
        htitle.String = sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc);
        set(hline2,'XData',1:nmcmc,'YData',epsilon_mcmc);
        pause(0.01);        
    end
    
    save tmp.mat        % iteratively save result
end
total_time =toc;




% add the inital values
param_mcmc = [param_mcmc0; param_mcmc];
epsilon_mcmc = [epsilon_mcmc0; epsilon_mcmc];

% --- save result
if nutsFlag == 1 && dualAvgFlag == 1
    save fresult_NUTS_DA_v1.mat
elseif nutsFlag == 1 && dualAvgFlag == 0
    save fresult_NUTS_NoDA_v1.mat    
elseif nutsFlag == 0 && dualAvgFlag == 1    
    save fresult_Simple_DA_v1.mat        
elseif nutsFlag == 0 && dualAvgFlag == 0
    save fresult_Simple_NoDA_v1.mat            
end