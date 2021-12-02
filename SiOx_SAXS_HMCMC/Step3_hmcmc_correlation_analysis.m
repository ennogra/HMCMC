load fresult_NUTS_NoDA_v1.mat

%% trace plot
% % --- clear NaN (in case of incomplete iteration loops);
nmcmc = nnz(all(~isnan(param_mcmc),2));
param_mcmc(nmcmc+1:end,:) = [];

% --- inverse transform to normal scale
x1_mcmc = inversetransformx(param_mcmc(:,1:end-1)',lbub(logical(fitFlag),:))';
Sigmavi_mcmc = exp(param_mcmc(:,end));
param_mcmc_inverted = [x1_mcmc,Sigmavi_mcmc];

% --- trace plot of parameters
% figure('position',[103,100,600,400]);
figure('position',[103,100,600,300]);

ha = tight_subplot(size(param_mcmc,2),1,[0.01 0.],[.1 .01],[.125 .025]);
axes(ha(1));
plot([0:nmcmc-1],param_mcmc_inverted(:,1));
set(gca,'ylim',[1.2,1.5]);
ylabel('I_0');

axes(ha(2));
plot([0:nmcmc-1],param_mcmc_inverted(:,2)*1e5);
%set(gca,'ylim',[1.1,1.5]);
ylabel({'I_b','(10^{-5})'});

axes(ha(3));
plot(0:nmcmc-1,param_mcmc_inverted(:,3));
set(gca,'ylim',[740,775]);
ylabel(['R (',char(197),')']);

axes(ha(4));
plot([0:nmcmc-1],param_mcmc_inverted(:,4));
set(gca,'ylim',[52,68]);
ylabel(['\sigma_R (',char(197),')']);

axes(ha(5));
plot([0:nmcmc-1],param_mcmc_inverted(:,5)*1e5);
%set(gca,'ylim',[44,59]);
ylabel({'\sigma_q', ['(10^{-5} ',char(197),'^{-1})']});

axes(ha(6));
plot([0:nmcmc-1],param_mcmc_inverted(:,6)*1e3);
set(gca,'ylim',[1,5]);
ylabel({'\sigma^2', '(10^{-3})'});
  
xlabel('iteration');
set(ha,'XMinorTick','on');
set(ha,'YMinorTick','on'); 
linkaxes(ha,'x');
%set(ha,'XTick',[0:200:1800]); 
%set(ha,'xlim',[0,2000]);
%set(ha,'ticklength',[0.015,0.015]);
set(ha(1:end-1),'XTickLabel','');
%return

%% statistics of parameters after burnin
% --- burn-in
nburnin=1000;
nstable = nmcmc-nburnin;
param_postburnin = param_mcmc(nburnin+1:end,:);
param_mcmc_inverted_postburnin = param_mcmc_inverted(nburnin+1:end,:);

% --- mean and variance of parameters
E_param = mean(param_mcmc_inverted_postburnin);
V_param = var(param_mcmc_inverted_postburnin);

array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'})


% --- MCMC standard error by effective sample size
[ESS,Sigma] = multiESS(param_mcmc_inverted_postburnin);    % book assume no correlation between thetas
% [ESS,Sigma] = multiESS(param_postburnin);    % book assume no correlation between thetas

ESS = round(ESS);
% MCMCSE = sqrt(diag(cov(param_mcmc))'/ESS);    % useing covariance of mcmc
MCMCSE = sqrt(diag(Sigma)'/ESS);    % useing covariance from mutliESS

array2table([E_param; sqrt(V_param);MCMCSE],'VariableNames',{'I0','bk','R','sigma_R','sigma_q','variance'},'RowNames',{'Mean','SD','MCMCSE'})
fprintf('Effective sample size is %d\n',ESS);

% --- matrix plot
param_mcmc_inverted_postburnin_scaled = param_mcmc_inverted_postburnin;
param_mcmc_inverted_postburnin_scaled(:,2) = param_mcmc_inverted_postburnin_scaled(:,2)*1e5;
param_mcmc_inverted_postburnin_scaled(:,5) = param_mcmc_inverted_postburnin_scaled(:,5)*1e5;
param_mcmc_inverted_postburnin_scaled(:,6) = param_mcmc_inverted_postburnin_scaled(:,6)*1e3;

figure('position',[573,125,400,350]);
%[S,AX,BigAx,H,HAx] = plotmatrix(param_postburnin);    % plot transformed parameters
[S,AX,BigAx,H,HAx] = plotmatrix(param_mcmc_inverted_postburnin_scaled); % plot parameters in normal scale
set(AX,'xminortick','on','yminortick','on','ticklength',[0.05,0.05]);
ylabel(AX(1,1),'I_0');
ylabel(AX(2,1),'I_b (10^{-5})');
ylabel(AX(3,1),['R (',char(197),')']);
ylabel(AX(4,1),['\sigma_R (',char(197),')']);
ylabel(AX(5,1),{'\sigma_q',['(10^{-5}',char(197),'^{-1})']});
ylabel(AX(6,1),'\sigma^2 (10^{-3})');
xlabel(AX(end,1),'I_0');
xlabel(AX(end,2),'I_b (10^{-5})');
xlabel(AX(end,3),['R (',char(197),')']);
xlabel(AX(end,4),['\sigma_R (',char(197),')']);
xlabel(AX(end,5),['\sigma_q (10^{-5}',char(197),'^{-1})']);
xlabel(AX(end,6),'\sigma^2 (10^{-3})');
set(H,'EdgeColor','none')
set(AX,'ticklength',[0.05,0.05]);
% for ii=1:6
%     delete(AX(ii,ii+1:end));
% end

%% autocorrelation
% --- calculate
acf_param_mcmc = nan(size(param_postburnin));
for ii=1:size(param_postburnin,2)
%     acf_param_mcmc(:,ii) = corrFFT(param_postburnin(:,ii));    
    acf_param_mcmc(:,ii) = corrFFT(param_mcmc_inverted_postburnin(:,ii));    

end

% --- plot
ylabel_list = {
    'I_0'
    'I_b'
    'R'
    '\sigma_R'
    '\sigma_q'
    '\sigma^2'
    };
% figure('position',[103,100,600,400]);
figure('position',[673,125,400,350]);

%ha = tight_subplot(size(param_mcmc,2),1,[0.01 0.],[.1 .01],[.125 .025]);
% xlist = 1:size(acf_param_mcmc,1)-1;
ha = nan(1,size(acf_param_mcmc,2));
for ii=1:size(acf_param_mcmc,2)
 %
    %axes(ha(ii));
    ha(ii) = subplot(size(param_postburnin,2),1,ii);
    hold on; box on;
    stem(acf_param_mcmc(2:end,ii),'marker','none');
%    legend(ylabel_list{ii},'Location','best','box','off','LineWidth',0.01);
    plot([0,size(acf_param_mcmc,1)],1.96/sqrt(size(acf_param_mcmc,1))*ones(1,2),'r--');
    plot([0,size(acf_param_mcmc,1)],-1.96/sqrt(size(acf_param_mcmc,1))*ones(1,2),'r--');

    hold off;
    if ii == round(size(acf_param_mcmc,2)/2)
        ylabel('Auto-correlation');
    end
%     ylabel(['ACF (',ylabel_list{ii},')']);
    %ylabel(ylabel_list{ii});
    text(90,0.21,ylabel_list{ii},'FontSize',8,'BackgroundColor','w');

end
xlabel('Lags');
set(ha,'XMinorTick','on');
set(ha,'YMinorTick','on'); 
linkaxes(ha,'xy');
set(ha,'xlim',[0,100]);
set(ha,'ylim',[-0.06,0.4]);
set(ha(1:end-1),'XTickLabel','');
% sgtitle('Autocorrelation');


%% effective sample size
ESS_by_acf = nan(1,size(acf_param_mcmc,2));
ESS_idx_by_acf = nan(1,size(acf_param_mcmc,2));
for ii=1:size(acf_param_mcmc,2)
    idx = find(abs(acf_param_mcmc(:,ii))<0.05,1,'first');
    ESS_idx_by_acf(ii) = idx;
    ESS_by_acf(ii) = round(size(acf_param_mcmc,1)/(1+2*sum( abs(acf_param_mcmc(2:idx,ii) ) )));    
end

ESS_by_multiESS = nan(1,size(acf_param_mcmc,2));
for ii=1:size(acf_param_mcmc,2)
%    ESS_by_multiESS(ii) = multiESS(param_postburnin(:,ii));
     ESS_by_multiESS(ii) = round(multiESS(param_mcmc_inverted_postburnin(:,ii)));
end

round([ESS_by_acf',ESS_by_multiESS'])


%% collect result
counter = nnz(diff(param_postburnin(:,end))~=0);
nparam=size(param_mcmc_inverted_postburnin,2);
[repmat(size(param_mcmc_inverted_postburnin,1),nparam,1),...
    ESS_by_acf(:),ESS_by_multiESS(:),...
    ESS_by_acf(:)/size(param_mcmc_inverted_postburnin,1), ...
    ESS_by_multiESS(:)/size(param_mcmc_inverted_postburnin,1), ...
    repmat(counter,nparam,1),...
    repmat(counter/size(param_postburnin,1),nparam,1)]

%% plot iterations