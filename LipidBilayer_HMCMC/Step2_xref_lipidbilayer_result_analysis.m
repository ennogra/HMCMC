% addpath(fullfile(pwd,'xrayrefraction'),fullfile(pwd,'xrayrefraction','AtomicScatteringFactor'));

%% ========================== post MCMC analysis ==========================
load fresult_LB_Simple_DA_v1.mat

%% add initial values
param_mcmc = [param_mcmc0; param_mcmc];
epsilon_mcmc = [epslion_mcmc0; epsilon_mcmc];
dPotE_mcmc = [dPotE_init; dPotE_mcmc];

%% trace plot
% % --- clear NaN (in case of incomplete iteration loops);
nmcmc = nnz(all(~isnan(param_mcmc),2));
param_mcmc(nmcmc+1:end,:) = [];
epsilon_mcmc(nmcmc+1:end) = [];
dPotE_mcmc(nmcmc+1:end,:) = [];

% --- inverse transform to normal scale
x1_mcmc = inversetransformx(param_mcmc(:,1:end-1)',lbub(logical(fitFlag),:))';
Sigmavi_mcmc = exp(param_mcmc(:,end));
param_mcmc_inverted = [x1_mcmc,Sigmavi_mcmc];

% --- select parameters for trace plot
idx_plot = 1:size(param_mcmc,2);
% idx_plot = [1:8,30,size(param_mcmc,2)];

ylabel_str = {'I_0','\sigma_{Si}','h_o','h_t','h_i','\delta_o','\delta_t','\delta_i','\sigma_o','\sigma_t','\sigma_i','\sigma^2'};

nsubplots = length(idx_plot);
 figure('position',[488,430,425,280]);
 ha = tight_subplot(size(param_mcmc,2),1,[0.005 0.01],[.15 .01],[.1 .025]);

kk=1;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel(ylabel_str{kk});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[1.01,1.1]);

% si roughness
kk = 2;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[2.3,2.9]);

kk = 3;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[9,15.8]);

kk = 4;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[21,25]);

kk = 5;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[9,13]);

kk = 6;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk)*1e7);
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[5.8,7.3]);

kk = 7;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk)*1e7);
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
set(gca,'ylim',[3.9,4.6]);

kk = 8;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
%  set(gca,'ylim',[9,9.5]);

kk = 9;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
% set(gca,'ylim',[3,7]);

kk = 10;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
% set(gca,'ylim',[2.8,5]);

kk = 11;
axes(ha(kk));
plot(param_mcmc_inverted(:,kk));
ylabel({ylabel_str{kk};});
set(gca,'XTickLabel',[]);
% set(gca,'ylim',[7.8,9.4]);

kk = 12;
axes(ha(kk));
plot(log10(param_mcmc_inverted(:,kk)));
ylabel(['log',ylabel_str{kk}]);
set(gca,'ylim',[-3.8,-1.8]);
xlabel('iteration');

linkaxes(ha,'x');
set(ha,'xlim',[0,3000]);
set(ha,'XMinorTick','on');
set(ha,'yticklabel',[]);



%% statistics of parameters
% --- burnin
nburnin= 20;
nstable = nmcmc-nburnin;
param_postburnin = param_mcmc(nburnin+1:end,:);
param_mcmc_inverted_postburnin = param_mcmc_inverted(nburnin+1:end,:);

% --- mean and variance of parameters
E_param = mean(param_mcmc_inverted_postburnin);
V_param = var(param_mcmc_inverted_postburnin);
tstat = E_param./sqrt(V_param); 
pvalue = (1-tcdf(E_param./sqrt(V_param),nmcmc))*2;

% --- MCMC standard error by effective sample size
[ESS,Sigma] = multiESS(param_mcmc_inverted_postburnin);    % book assume no correlation between thetas
%[ESS,Sigma] = multiESS(param_postburnin);    % book assume no correlation between thetas
ESS = round(ESS);
% MCMCSE = sqrt(diag(cov(param_mcmc))'/ESS);    % useing covariance of mcmc
MCMCSE = sqrt(diag(Sigma)'/ESS);    % useing covariance from mutliESS

%array2table([E_param; sqrt(V_param);MCMCSE],'RowNames',{'Mean','SD','MCMCSE'})


fprintf('Effective sample size is %d\n',ESS);



%table_var = {'I0','z_offset','d_ptba','sigma_dtba','d_ps','sigma_ps','d_au','f_e','sigma2','d_film'};
array2table([E_param; sqrt(V_param); tstat; pvalue;MCMCSE ],'RowNames',{'Mean','SD','t-stat','p-value','MCMCSE'})%,'VariableNames',table_var)


%% % --- matrix plot
 figure('position',[570,373,600,500]);
[S,AX,BigAx,H,HAx] = plotmatrix(param_mcmc_inverted_postburnin);

fontsize = 15;
for kk=1:size(param_mcmc_inverted,2)
    ylabel(AX(kk,1),ylabel_str{kk},'FontSize',fontsize);
    %xlabel(AX(end,kk),ylabel_str{kk});
    set(AX(1,kk),'XTickMode','auto','XAxisLocation','top');
    xlabel(AX(1,kk),ylabel_str{kk},'FontSize',fontsize);
end
set(AX(:,1),'YTicklabel',[]);


%% autocorrelation
% --- calculate
acf_param_mcmc = nan(size(param_postburnin));
for ii=1:size(param_postburnin,2)
    %     acf_param_mcmc(:,ii) = corrFFT(param_postburnin(:,ii));
    acf_param_mcmc(:,ii) = corrFFT(param_mcmc_inverted_postburnin(:,ii));
    
end

 figure('position',[773,125,260,320]);
 ha = tight_subplot(size(param_mcmc,2),1,[0.005 0.01],[.1 .01],[.15 .05]);

for ii=1:length(idx_plot)
    axes(ha(ii));
    %ha(ii) = subplot(length(idx_plot),1,ii);
    hold on; box on;
    stem(0:size(acf_param_mcmc,1)-1, acf_param_mcmc(1:end,idx_plot(ii)),'marker','none');
    %    legend(ylabel_list{ii},'Location','best','box','off','LineWidth',0.01);
    plot([0,size(acf_param_mcmc,1)],1.96/sqrt(size(acf_param_mcmc,1))*ones(1,2),'r--','linewidth',0.1);
    plot([0,size(acf_param_mcmc,1)],-1.96/sqrt(size(acf_param_mcmc,1))*ones(1,2),'r--','linewidth',0.1);
    
    hold off;
    if ii==round(length(idx_plot)/2)
        ylabel('ACF');
    end
    text(150,0.5,ylabel_str{ii},'FontSize',8);
    %     ylabel(['ACF (',ylabel_list{ii},')']);
    % ylabel(ylabel_str{ii});
     set(ha,'ylim',[-0.2,1]);
    set(ha(ii),'ytick',[0,.8]);
    set(ha(ii),'yticklabel',[0,.8]);

end

xlabel('Lags');
set(ha,'XMinorTick','on');
% set(ha,'YMinorTick','on');
%set(ha(end),'XTickMode','auto');
%linkaxes(ha,'xy');
set(ha,'xlim',[-1,200]);
set(ha(end),'xticklabel',get(ha(end),'xtick'));
%set(ha(1:end-1),'XTickLabel','');
% sgtitle('Autocorrelation');


%% effective sample size
ESS_by_acf = nan(1,size(acf_param_mcmc,2));
ESS_idx_by_acf = nan(1,size(acf_param_mcmc,2));
for ii=1:size(acf_param_mcmc,2)
    idx = find(abs(acf_param_mcmc(:,ii))<0.05,1,'first');
    ESS_idx_by_acf(ii) = idx;
    ESS_by_acf(ii) = round(size(acf_param_mcmc,1)/(1+2*sum((acf_param_mcmc(2:idx,ii) ) )));
end

ESS_by_multiESS = nan(1,size(acf_param_mcmc,2));
for ii=1:size(acf_param_mcmc,2)
    %    ESS_by_multiESS(ii) = multiESS(param_postburnin(:,ii));
    ESS_by_multiESS(ii) = round(multiESS(param_mcmc_inverted_postburnin(:,ii)));
end

array2table(round([ESS_by_acf',ESS_by_multiESS']),'VariableNames',{'ESS_by_ACF','ESS_by_multiESS'})
%fprintf('Effective sample size is %d\n',ESS);


%% collect result
counter = nnz(diff(param_postburnin(:,end))~=0);
nparam=size(param_mcmc_inverted_postburnin,2);

[repmat(size(param_mcmc_inverted_postburnin,1),nparam,1),...
    ESS_by_acf(:),ESS_by_multiESS(:),...
    ESS_by_acf(:)/size(param_mcmc_inverted_postburnin,1), ...
    ESS_by_multiESS(:)/size(param_mcmc_inverted_postburnin,1), ...
    repmat(counter,nparam,1),...
    repmat(counter/size(param_postburnin,1),nparam,1)];

%% Plot performance vs iterations: collect x-ray results
% --- collect result
I_cal_tmp = cell(nstable,nsets);
fresnel_tmp = cell(nstable,nsets);
q_tmp = cell(nstable,nsets);
edp_all = cell(nstable,1);
rawedp_all = cell(nstable,1);
for ii=1:nstable
    fprintf('%.5d',ii);
    %         if mod(ii,40) == 0
    %             fprintf('\n');
    %         end
    [I_cal_tmp(ii,:),q_tmp(ii,:),~,fresnel_tmp(ii,:)] = fcn_xray_ref_hmcmc(param_postburnin(ii,1:end-nsets),alphai_set);
    edp_all{ii} = edp_0;
    rawedp_all{ii} = rawedp;
    fprintf('\b\b\b\b\b');
end
fprintf('\n');
% --- convert to mat
I_cal_all = cell(1,nsets);
q_all = cell(1,nsets);
fresnel_all = cell(1,nsets);
for iset = 1:nsets
    I_cal_all{iset} = cell2mat(I_cal_tmp(:,iset)');
    q_all{iset} = cell2mat(q_tmp(:,iset)');
    fresnel_all{iset} = cell2mat(fresnel_tmp(:,iset)');
end

%save(fullfile(pwd,['reduced_',fresult_list{iresult}(8:end)]));
save msplot.mat;

%% Load back result and plot reflectivity result
load('msplot.mat');

normFlag = 2; %1/2/3: fresnel/q^4/log
figure('Position',[644,380,250,280]);
hold on; box on;

for iset = 1:nsets
    alphai = alphai_set{iset};
    qdata = q_all{iset};
    Idata = Iset{iset};
    Ierr = Ierr_set{iset};
    I_cal = I_cal_all{iset};
    I_fresnel = fresnel_all{iset};
    
    
    if normFlag == 1 % scale with fresnel
        plotshift = -0.8;
        % --- get fresnel of the readback alphai
        Ifresnel0 = fresnel(alphai,delta_Si,beta_Si,lambda);
        Ifresnel0 = Ifresnel0(:,2);
        
        errorbar(alphai,Idata./Ifresnel0 + (iset-1)*plotshift,Ierr./Ifresnel0,'o','markersize',6,'CapSize',4,'LineWidth',0.5);  % data
        I_cal_by_fresnel = I_cal./I_fresnel;
        lb_I_cal = quantile(I_cal_by_fresnel,0.025,2);   % 2.5% quantile
        ub_I_cal = quantile(I_cal_by_fresnel,0.975,2);   % 97.5% quantile
        patch_x = [alphai; flip(alphai)];
        patch_y = [lb_I_cal;flip(ub_I_cal)];
        patch(patch_x,patch_y + (iset-1)*plotshift ,'m','EdgeColor','m','FaceAlpha',.4);
        
    elseif normFlag ==2  % --- scale with q^4
        plotshift = 0.1;
        qdata0 = 4*pi/lambda*sind(alphai);  % q of the readback alphai
        errorbar(alphai,Idata.*qdata0.^4*plotshift^(iset-1),Ierr.*qdata0.^4*plotshift^(iset-1),'o','markersize',6,'CapSize',4,'LineWidth',0.5);  % data
        I_cal_by_q4 = I_cal.*qdata.^4;
        lb_I_cal = quantile(I_cal_by_q4,0.025,2);   % 2.5% quantile
        ub_I_cal = quantile(I_cal_by_q4,0.975,2);   % 97.5% quantile
        patch_x = [alphai; flip(alphai)];
        patch_y = [lb_I_cal;flip(ub_I_cal)];
        patch(patch_x,patch_y*plotshift^(iset-1),'m','EdgeColor','m','FaceAlpha',.4);
        
    elseif normFlag == 3  % --- no scale
        plotshift = 0.1;
        qdata0 = 4*pi/lambda*sind(alphai);  % q of the readback alphai
        errorbar(alphai,Idata*plotshift^(iset-1),Ierr*plotshift^(iset-1),'o','markersize',6,'CapSize',4,'LineWidth',0.5);  % data
        lb_I_cal = quantile(I_cal,0.025,2);   % 2.5% quantile
        ub_I_cal = quantile(I_cal,0.975,2);   % 97.5% quantile
        patch_x = [alphai; flip(alphai)];
        patch_y = [lb_I_cal;flip(ub_I_cal)];
        patch(patch_x,patch_y*plotshift^(iset-1),'m','EdgeColor','m','FaceAlpha',.4);
        
    end
end
if normFlag == 1
    set(gca,'yscale','linear');
    ylabel('R/R_F');
elseif normFlag == 2
    set(gca,'yscale','log');
    ylabel('R\timesq^4');
elseif normFlag == 3
    set(gca,'yscale','log');
    ylabel('R');
end
set(gca,'xminortick','on','YMinorTick','on');
xlabel('incident angle (deg)');
set(gca,'xlim',[0,1.8]);
set(gca,'TickLength',[0.02,0.02]);
legend({'Reflectivity','95% confidence'},'box','off');


%% Plot EDP (electron density profile in terms of dispersion)

% E_d = E_param(3:7);
% E_delta = [E_param(8:11),siox.dispersion];
E_param_transformed = transformx(E_param(1:end-1)',lbub(logical(fitFlag),:));
plotEDPFlag = 1;
[~,~] = fcn_xray_ref_hmcmc(E_param_transformed,alphai_set);   
close 
close
% --- plot EDP delta
figure('Position',[644,380,250,280]);

% --- plot EDP delta
% ha1 = subplot(2,1,1);
hold on; box on;

z_mean = linspace(min(cellfun(@(x,y)(x(1,1)-y(end-1,1))',edp_all,rawedp_all)),max(cellfun(@(x,y)(x(end,1)-y(end-1,1))',edp_all,rawedp_all)),size(edp_all{1},1));
z_mean = z_mean(:);

delta_interp = cell2mat(cellfun(@(x,y)interp1(x(:,1)-y(end-1,1),x(:,2),z_mean)',edp_all,rawedp_all,'UniformOutput',false))';
lb_delta = quantile(delta_interp,0.025,2);   % 2.5% quantile
ub_delta = quantile(delta_interp,0.975,2);   % 97.5% quantile
delta_median = nanmedian(delta_interp,2);
delta_mean = nanmean(delta_interp,2);
delta_mode = mode(delta_interp,2);
patch_x = [-z_mean; flip(-z_mean)];
patch_y = [lb_delta;flip(ub_delta)];
patch(patch_x,patch_y,'m','EdgeColor','none','FaceAlpha',1);
plot(-z_mean,delta_median,'-','LineWidth',1.5);
    plot(-edpstep(:,1)+d_total,edpstep(:,2),'-','linewidth',1);     % plot step-like profile
    plot(-edp_0(:,1)+d_total,slice_dispersion,'b--','linewidth',1);     % plot individual slice profile
 
%     plot(-z_mean,delta_mean,'-','LineWidth',1.5);
%     plot(-z_mean,delta_mode,'-','LineWidth',1.5);
xlabel(['Height (',char(197),')']);
% set(gca,'ylim',[2e-7,14e-7]);
set(gca,'xlim',[-10,70]);
ylabel('\delta');
set(gca,'TickLength',[0.02,0.02]);
% set(gca,'xlim',[-200,2500]);
set(gca,'xminortick','on','YMinorTick','on');
legend({'95% confidence','median','unsmeared profile','smeared profile'},'box','off')


  
return
