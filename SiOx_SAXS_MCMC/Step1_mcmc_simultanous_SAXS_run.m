% use Cov of warmup
clear

global fitFlag x2
global lbub

%% some flags ans settings
% --- about data
weightFlag = 1;       % 0/1 use weight for intensity

% --- about HMCMC
updatePlotFlag = 1;    % 0/1 update plot during MCMC


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
    R               300                 1000                 1
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

%% Define potential
% U = -log(P);
U = @(X)potentialE(X,qdata,Idata,weight); % X(end) is the Sigmavi


%% Initialize parameters 
% ---  initialize variance of residual to guessed values
Sigmavi = var(log10(I_cal) - log10(Idata));

% --- initial U and dU values
xp_init = [x1(:)',log(Sigmavi)];
Ui = U(xp_init);

%% proposal distribution (assume normal)
% epsilon = 0.02;
% Sigma_prop = epsilon*eye(length(xp_init));

M = [0.000165062153113839,-0.000169036430577784,7.50411945895435e-05,-2.54289456127242e-08,0.000308522967838566,1.72989067618196e-05;-0.000169036430577784,0.0123588122294425,0.000140910739148036,-0.000128428124146926,0.00492068467629250,-0.000365124905626878;7.50411945895435e-05,0.000140910739148036,0.000130920153942672,-0.000317503113541892,0.000968925650428361,1.89130710100616e-05;-2.54289456127242e-08,-0.000128428124146926,-0.000317503113541892,0.00481101233495571,-0.00399635628278936,0.000161046220063323;0.000308522967838566,0.00492068467629250,0.000968925650428361,-0.00399635628278936,3.92880764407997,0.00879487237066152;1.72989067618196e-05,-0.000365124905626878,1.89130710100616e-05,0.000161046220063323,0.00879487237066152,0.0155888150165789];
Sigma_prop = diag(diag(M));
Sigma_prop(5,5) = 0.01;

%% start MCMC
% --- assign initial values
xp = xp_init;

% --- start MCMC
nmcmc = 100000;
param_mcmc = nan(nmcmc,length(xp));
U_mcmc = nan(nmcmc,length(xp));
param_mcmc0 = xp_init;

counter = 0;

%%
warning('off','MATLAB:Figure:FigureSavedToMATFile');
tic
for ii=1:nmcmc
    
    % propose a move from proposal joint-distribution
    xp_tmp = mvnrnd(xp,Sigma_prop);
            
    % calculate log ratio (Metropolis)
    %logr = log(mvnpdf(xp_tmp,mu,Sigma)) - log(mvnpdf(xp,mu,Sigma));
    logr = U(xp) - U(xp_tmp);
 
    % take it or not
    if log(rand)<logr
        xp = xp_tmp;
        counter = counter+1;
    end

    fprintf('%d of %d iterations (%.1f sec) with counter=%d: ',ii,nmcmc,toc,counter);
    fprintf(['[',repmat('%f ',1,length(xp)),'\b]\n'],xp);
    
    % --- collect
    param_mcmc(ii,:) = xp;

    
        
    % --- update plot
    if updatePlotFlag == 1
        I_cal = fcn_saxs_sphere_hmcmc(xp(1:end-1),qdata);     % preview initial conditions before fitting
        if ii==1
            hfig = figure;
       %     subplot(3,1,1:2);
            hold on;
            errorbar(qdata,Idata,Idata_err,'o');
            hline = plot(qdata,I_cal,'-','linewidth',1.5);
            hold off; box on;
            set(gca,'xscale','log');
            set(gca,'yscale','log');
%             subplot(3,1,3);
%             hline2 = plot(epsilon_mcmc);
%             set(gca,'xlim',[0,nmcmc]);
            xlabel('iteration');
%            ylabel('\epsilon');
            set(gca,'yscale','log');
            htitle = sgtitle(sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc));            
        end
        hline.YData = I_cal;
        htitle.String = sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc);
%        set(hline2,'XData',1:nmcmc,'YData',epsilon_mcmc);
        pause(0.001);        
    end
    
    save tmp.mat        % iteratively save result
end
total_time =toc;

% add the inital values
param_mcmc = [param_mcmc0; param_mcmc];
save fresult_MCMC_simultaneous_v1.mat

