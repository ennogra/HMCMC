addpath(fullfile(pwd,'xrayrefraction'),fullfile(pwd,'xrayrefraction','AtomicScatteringFactor'));

%% ========================== post MCMC analysis ==========================
load fresult_Simple_NoDA_randoma_nontrivalmass_v1.mat

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

idx_plot = [1:8,30,size(param_mcmc,2)];
ylabel_str = {'I_0','z_{offset}','d_{ptba}','\sigma_{dtba}','d_{ps}','\sigma_{ps}','d_{au}','f_e','a_9','\sigma^2'};

nsubplots = length(idx_plot);
figure
for kk=1:length(idx_plot)
    subplot(nsubplots,1,kk);
    % --- plot transformed parameters
    % plot(param_mcmc(:,idx_plot(kk)));
    % --- plot inverse transformed (in true scale)
    plot(param_mcmc_inverted(:,idx_plot(kk)));
    ylabel(ylabel_str{kk});
end
% --- plot epsilon
% subplot(nsubplots,1,kk+1);
% plot(epsilon_mcmc)
set(gca,'yscale','log')
% ylabel('\epsilon');
xlabel('iteration');
linkaxes(findall(gcf,'type','axes'),'x');

%% --- plot d_au, d_ps, Sigma2 trace
figure('pos',[500,100,430,280]);
ha = tight_subplot(4,1,[0.03 0.],[.1 .05],[.125 .025]);
axes(ha(1));
plot(param_mcmc_inverted(:,3));
set(gca,'xticklabel',[]);
set(gca,'XMinorTick','on','YMinorTick','on');
ylabel(['d_{ptba} (',char(197),')']);

axes(ha(2));
plot(param_mcmc_inverted(:,5));
set(gca,'xticklabel',[]);
set(gca,'XMinorTick','on','YMinorTick','on');
ylabel(['d_{ps} (',char(197),')']);

axes(ha(3));
plot(param_mcmc_inverted(:,5)+param_mcmc_inverted(:,3));
set(gca,'xticklabel',[]);
set(gca,'XMinorTick','on','YMinorTick','on');
ylabel(['d_{ps+ptba} (',char(197),')']);


axes(ha(4));
plot(Sigmavi_mcmc);
set(gca,'yscale','log');
xlabel('Iteration');
set(gca,'XMinorTick','on','YMinorTick','on');
% set(gca,'ylim',[2e-4,1e-2]);
set(gca,'TickLength',[0.02,0.02]);
% set(gca,'xtick',[100:100:900]);
ylabel('log(\sigma^2)');
linkaxes(ha,'x');
set(ha,'xlim',[0,2000]);


%% statistics of parameters
% --- burn-in
nburnin= 601; 
nstable = nmcmc-nburnin;
param_postburnin = param_mcmc(nburnin+1:end,:);
param_mcmc_inverted_postburnin = param_mcmc_inverted(nburnin+1:end,:);

% --- mean and variance of parameters
E_param = mean(param_mcmc_inverted_postburnin);
V_param = var(param_mcmc_inverted_postburnin);

% --- MCMC standard error by effective sample size
[ESS,Sigma] = multiESS(param_mcmc_inverted_postburnin);    % book assume no correlation between thetas
%[ESS,Sigma] = multiESS(param_postburnin);    % book assume no correlation between thetas
ESS = round(ESS);
% MCMCSE = sqrt(diag(cov(param_mcmc))'/ESS);    % useing covariance of mcmc
MCMCSE = sqrt(diag(Sigma)'/ESS);    % useing covariance from mutliESS

array2table([E_param; sqrt(V_param);MCMCSE],'RowNames',{'Mean','SE','MCMCSE'})
fprintf('Effective sample size is %d\n',ESS);

%% % --- matrix plot
% --- for non-bspline coefficients
figure
[S,AX,BigAx,H,HAx] = plotmatrix(param_postburnin(:,[1:8,end]));    % plot transformed parameters
%[S,AX,BigAx,H,HAx] = plotmatrix(param_mcmc_inverted_postburnin(:,[1:8,end]));  
ylabel(AX(1,1),'I_0');
ylabel(AX(2,1),'z_{offset}');
ylabel(AX(3,1),'d_{ptba}');
ylabel(AX(4,1),'\sigma_{ptba}');
ylabel(AX(5,1),'d_{ps}');
ylabel(AX(6,1),'\sigma_{ps}');
ylabel(AX(7,1),'d_{au}');
ylabel(AX(8,1),'f_{e}');
ylabel(AX(9,1),'\sigma^2');
xlabel(AX(end,1),'I_0');
xlabel(AX(end,2),'z_{offset}');
xlabel(AX(end,3),'d_{ptba}');
xlabel(AX(end,4),'\sigma_{ptba}');
xlabel(AX(end,5),'d_{ps}');
xlabel(AX(end,6),'\sigma_{ps}');
xlabel(AX(end,7),'d_{au}');
xlabel(AX(end,8),'f_{e}');
xlabel(AX(end,9),'\sigma^2');
set(H,'EdgeColor','none')



%% autocorrelation
% --- calculate
acf_param_mcmc = nan(size(param_postburnin));
for ii=1:size(param_postburnin,2)
%     acf_param_mcmc(:,ii) = corrFFT(param_postburnin(:,ii));    
    acf_param_mcmc(:,ii) = corrFFT(param_mcmc_inverted_postburnin(:,ii));    

end


figure('position',[773,125,400,350]);
%ha = tight_subplot(size(param_mcmc,2),1,[0.01 0.],[.1 .01],[.125 .025]);
% xlist = 1:size(acf_param_mcmc,1)-1;
ha = nan(1,length(idx_plot));
% for ii=1:size(acf_param_mcmc,2)
for ii=1:length(idx_plot)
    %axes(ha(ii));
    ha(ii) = subplot(length(idx_plot),1,ii);
    hold on; box on;
    stem(acf_param_mcmc(2:end,idx_plot(ii)),'marker','none');
%    legend(ylabel_list{ii},'Location','best','box','off','LineWidth',0.01);
    plot([0,size(acf_param_mcmc,1)],1.96/sqrt(size(acf_param_mcmc,1))*ones(1,2),'r--');
    plot([0,size(acf_param_mcmc,1)],-1.96/sqrt(size(acf_param_mcmc,1))*ones(1,2),'r--');

    hold off;
%     ylabel(['ACF (',ylabel_list{ii},')']);
    ylabel(ylabel_str{ii});

end
xlabel('Lags');
set(ha,'XMinorTick','on');
set(ha,'YMinorTick','on'); 
linkaxes(ha,'xy');
set(ha,'xlim',[0,1500]);
set(ha(1:end-1),'XTickLabel','');
sgtitle('Autocorrelation');


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

round([ESS_by_acf',ESS_by_multiESS'])


%% collect result
counter = nnz(diff(param_postburnin(:,end))~=0);
nparam=size(param_mcmc_inverted_postburnin,2);
[repmat(size(param_mcmc_inverted_postburnin,1),nparam,1),...
    ESS_by_acf(:),ESS_by_multiESS(:),...
    ESS_by_acf(:)/size(param_mcmc_inverted_postburnin,1), ...
    ESS_by_multiESS(:)/size(param_mcmc_inverted_postburnin,1), ...
    repmat(counter,nparam,1),...
    repmat(counter/size(param_postburnin,1),nparam,1)];

%% collect x-ray results
Is_cal_all = nan(length(Is),nmcmc);
Y1_all = nan(length(Is),nmcmc);
Y2_all = nan(length(Is),nmcmc);
Y_elastic_all = nan(length(Is),nmcmc);
rho_all = cell(1,nmcmc);
nSlice = 800;  
for ii=1:nmcmc
    fprintf('Calculating %d of %d\n',ii,nmcmc);
    [Is_cal_all(:,ii),Y1_all(:,ii),Y2_all(:,ii),Y_elastic_all(:,ii)] = fcn_flu_spline_hmcmc(param_mcmc(ii,1:end-1),vzs); 
    rho_all{ii} = rho;
end
save msplot_hmcmc.mat


%% load resutl back and plot
load msplot_hmcmc.mat
 
%% plot initial iterations
vz0_mean = vz0;
Inorm_mean = Inorm;
alphaf_list = atand((vzs-vz0_mean)/sdd);

zlist = linspace(min(cellfun(@(x)x(end,1),rho_all)),max(cellfun(@(x)x(1,1),rho_all)),size(rho,1));
zlist = zlist(:);
edp_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,2),zlist),rho_all,'UniformOutput',false));
phi_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,3),zlist),rho_all,'UniformOutput',false));


idx_plot = [1,6,51,101,201,401,601,801];%,101,201];   % index to be plotted

nidx = length(idx_plot);

yplot_scale = (Is(end)-Ibk)/Inorm_mean;
figure('pos',[560,200,500,450]);
[ha, pos] = tight_subplot(nidx,2,[.01 .12],[.1 .01],[.1 .04]);
for ii=1:nidx
    axes(ha(2*ii-1));
    hold on; box on;
    plot(alphaf_list,(Is-Ibk)/Inorm_mean/yplot_scale,'o','markersize',2);
    plot(alphaf_list,(Is_cal_all(:,idx_plot(ii))-Ibk)/Inorm_mean/yplot_scale,'linewidth',1);    
    set(gca,'XTickLabelMode','auto')
    set(gca,'YTickLabelMode','auto')
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'xlim',[0,0.8]);
    set(gca,'ylim',[0,7]);    
    set(gca,'xminortick','on')
    set(gca,'yminortick','on')    
    %set(gca,'YTick',0:2:6);        

   
    %set(gca,'YTickLabel',0:1:4);    
    if ii~=nidx
        set(gca,'xticklabel',[]);        
    else
        xlabel('Exit Angle (deg)');
        set(gca,'XTick',0:0.2:0.8);
        set(gca,'XTickLabel',0:0.2:0.8);
    end
    if ii==round(nidx/2)
        %ylabel('Au L\alpha_{1,2} and L\beta_{2,15} Fluorescence (a.u.)');
        ylabel('XWFH (a.u.)');
    end    
    
    if ii==1
        legend({'data','calculation'},'Location','best','box','off')
    end    

    axes(ha(2*ii));
    plot(zlist,phi_interp(:,idx_plot(ii)),'linewidth',1);
    set(gca,'xminortick','on')
    set(gca,'yminortick','on')
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'xlim',[0,800]);
    set(gca,'ylim',[0,0.025]);
    
    if ii~=nidx
        set(gca,'xticklabel',[]);        
    else
        xlabel(['Height above Pd mirror (',char(197),')']);
        set(gca,'XTick',0:200:800);
        set(gca,'XTickLabel',0:200:800);
    end    
    if ii==round(nidx/2)
        ylabel(['Atomic number density (% ',char(197),'^{-1})']);
    end
    if ii==1
        text(400,0.012,'initial guess','FontSize',8);        
    else
        text(400,0.012,['iteration #',num2str(idx_plot(ii)-1)],'FontSize',8);        
    end
 
end

%% plot post-burnin
nburnin = 600;
Is_cal_post = Is_cal_all(:,nburnin+1:end);
Y1_post = Y1_all(:,nburnin+1:end);
Y2_post = Y2_all(:,nburnin+1:end);
Y_elastic_pos = Y_elastic_all(:,nburnin+1:end);

Y1_post_mode = mode(Y1_post,2);
Y2_post_mode = mode(Y2_post,2);
Y_elastic_pos_mode = mode(Y_elastic_pos,2);

Y1_post_median = nanmedian(Y1_post,2);
Y2_post_median = nanmedian(Y2_post,2);
Y_elastic_pos_median = nanmedian(Y_elastic_pos,2);

Inorm_post = param_mcmc_inverted_postburnin(:,1);
vz0_post = param_mcmc_inverted_postburnin(:,2);

Inorm_mean = mean(Inorm_post);
vz0_mean = mean(vz0_post);
Is_cal_mean = mean(Is_cal_post,2);

alphaf_list = atand((vzs-vz0_mean)/sdd);

yplot_scale = (Is(end)-Ibk)/Inorm_mean;
list_colororder = colororder;
figure('pos',[500,100,400,450])
subplot(3,1,2);
hold on; box on;
plot(alphaf_list,(Is-Ibk)/Inorm_mean/yplot_scale,'o','markersize',6);
lb_I_cal = quantile((Is_cal_post-Ibk)/Inorm_mean/yplot_scale,0.025,2);   % 2.5% quantile
ub_I_cal = quantile((Is_cal_post-Ibk)/Inorm_mean/yplot_scale,0.975,2);   % 97.5% quantile
patch_x = [alphaf_list; flip(alphaf_list)];
patch_y = [lb_I_cal;flip(ub_I_cal)];
% plot(alphaf_list,Y1_post_mode/yplot_scale,'color',list_colororder(3,:),'linewidth',1);
% plot(alphaf_list,Y2_post_mode/yplot_scale,'color',list_colororder(4,:),'linewidth',1);
% plot(alphaf_list,Y_elastic_pos_mode/yplot_scale,'color',list_colororder(5,:),'linewidth',1);
plot(alphaf_list,Y1_post_median/yplot_scale,'color',list_colororder(3,:),'linewidth',1);
plot(alphaf_list,Y2_post_median/yplot_scale,'color',list_colororder(4,:),'linewidth',1);
plot(alphaf_list,Y_elastic_pos_median/yplot_scale,'color',list_colororder(5,:),'linewidth',1);

patch(patch_x,patch_y,'m','EdgeColor','m','FaceAlpha',.4);
%plot(alphaf_list,(Is_cal_mean-Ibk)/Inorm_mean,'b-','linewidth',1);
ylabel('XWFH (a.u.)');
xlabel('Exit Angle (deg)');
set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'xlim',[0,0.8]);
set(gca,'ylim',[0,5]);
set(gca,'TickLength',[0.025,0.025]);
hplot = get(gca,'children');
legend(hplot([5,1,4,3,2]),{'data','Total (95% confidence)','Au L\alpha_{1,2} (median)','Au L\beta_{2,15} (median)','Elastic (median)'},'location','best','box','off')


% plot EDP curve
rho_post = rho_all(nburnin+1:end);

zlist = linspace(min(cellfun(@(x)x(end,1),rho_post)),max(cellfun(@(x)x(1,1),rho_post)),size(rho,1));
zlist = zlist(:);
edp_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,2),zlist),rho_post,'UniformOutput',false));
phi_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,3),zlist),rho_post,'UniformOutput',false));

lb_edp = quantile(edp_interp,0.025,2);   % 2.5% quantile
ub_edp = quantile(edp_interp,0.975,2);   % 97.5% quantile
edp_mean = nanmean(edp_interp,2);
edp_mode = mode(edp_interp,2);
edp_median = nanmedian(edp_interp,2);


lb_phi = quantile(phi_interp,0.025,2);   % 2.5% quantile
ub_phi = quantile(phi_interp,0.975,2);   % 97.5% quantile
phi_mean = nanmean(phi_interp,2);
phi_mode = mode(phi_interp,2);
phi_median = nanmedian(phi_interp,2);

%figure('pos',[500,500,400,200])
subplot(3,1,3);
hold on; box on;
patch_x = [zlist; flip(zlist)];
patch_y = [lb_phi;flip(ub_phi)];
patch(patch_x,patch_y,'m','EdgeColor','none','FaceAlpha',.5);
plot(zlist,phi_median,'b-','linewidth',1);
%plot(zlist,phi_mode,'b-','linewidth',1);
% plot(zlist,phi_mean,'k-','linewidth',1);
xlabel(['Heigth above mirror (',char(197),')']);
%ylabel(['Au atomic number density (%',char(197),'^{-1})']);
ylabel(['\phi_{au} (%',char(197),'^{-1})']);
set(gca,'xminorTick','on');
set(gca,'yminorTick','on');
set(gca,'xlim',[150,350]);
%set(gca,'xlim',[1,800]);
set(gca,'TickLength',[0.025,0.025]);
legend({'95% confidence','median'},'Location','best','box','off')

% --- Analyze statistics of parameters
acoeff = fliplr(param_mcmc_inverted_postburnin(:,end-30:end-1));
% figure
% boxplot(acoeff,'Notch','on','OutlierSize',.001);
% set(gca,'XMinorTick','off','YMinorTick','on');
% set(gca,'TickLength',[0.02,0.02]);
% %set(gca,'ylim',[20,110]);
% %set(gca,'xtick',[1:2:size(epsilon_list,2)],'XTickLabel',cellstr(num2str(epsilon_list(1,1:2:end)')));
% xlabel('index');
% ylabel('coefficient value');

E_param = mean(acoeff);
M_param = mode(acoeff);
V_param = var(acoeff);
tstat = E_param./sqrt(V_param); 
pvalue = (1-tcdf(E_param./sqrt(V_param),nmcmc))*2;
disp('CBS coefficients')
array2table([E_param; sqrt(V_param); tstat; pvalue ],'RowNames',{'Mean','SD','t-stat','p-value'})

subplot(3,1,1);
hold on; box on;
hstem = stem(pvalue,'o','markersize',4);
hstem.BaseLine.LineStyle = 'none';
set(gca,'ylim',[0,1]);
set(gca,'xlim',[0,length(pvalue)+1]);
plot([0,length(pvalue)+1],[0.05,0.05],'--','linewidth',1);
% plot([0,length(pvalue)+1],[0.1,0.1],'--','linewidth',1);
set(gca,'XMinorTick','off','YMinorTick','on');
set(gca,'TickLength',[0.02,0.02]);
set(gca,'xtick',0:2:30)
% set(gca,'xtick',0:2:30)
xlabel('Index of CBS coefficents');
ylabel('p-value');

%% statistics of other parameters
param_other = param_mcmc_inverted_postburnin(:,[1:8,end]);
param_other = [param_other, sum(param_other(:,[3,5]),2)];
E_param = mean(param_other);
M_param = mode(param_other);
V_param = var(param_other);
tstat = E_param./sqrt(V_param); 
pvalue = (1-tcdf(E_param./sqrt(V_param),nmcmc))*2;
table_var = {'I0','z_offset','d_ptba','sigma_dtba','d_ps','sigma_ps','d_au','f_e','sigma2','d_film'};
array2table([E_param; sqrt(V_param); tstat; pvalue ],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',table_var)


%% plt image
ms_vzs = (0:520)*d_ypixel;
alphaf = atand((ms_vzs-vz0_mean)/sdd);
alphaf = flip(alphaf);
img = img2(430:950,1:480) + img2(430:950,end-479-2:end-2);
[~,ind_alphafi] = min(abs(alphaf));
[~,ind_alphaff] = min(abs(alphaf-0.6));

img3 = img(ind_alphafi:-1:ind_alphaff,:);
img3 = img3/mean(img3(:),'omitnan')*10;
nx = size(img3,2);
pangle = atand(0.172/2/sdd)*2;
angle_list = pangle*((1:nx)-nx/2-1/2)+90;

figure('pos',[500,300,200,200]);     
imagesc(angle_list,(alphaf(ind_alphafi:-1:ind_alphaff)),img3,[1,30]);
colormap('jet');
set(gca,'ydir','norm');
set(gca,'TickDir','both')
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
set(gca,'TickLength',[0.025,0.025]);
ylabel('Exit angle (deg.)')
xlabel('In-plane angle (deg.)');


