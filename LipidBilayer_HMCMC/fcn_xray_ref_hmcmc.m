function varargout = fcn_xray_ref_hmcmc(x1,alphai_set)
global fitFlag x2 lbub
global lambda
global plotEDPFlag
global nsets nSlice nlayer
global edp_0 rawedp
global interface_profile
global conv_method
global binedp_delta_cutoff      % delta cutoff for edp binning

% --- Generate all the fitting parameters
x=zeros(length(fitFlag),1);
x(fitFlag==1) = x1;
x(fitFlag==0) = x2;
x = inversetransformx(x,lbub);

Inorm               = x(1:nsets);
Ibk                 = x(1*nsets+1:2*nsets);
resdawave           =  x(2*nsets+1);
resconstant         =  x(2*nsets+2);
foot_angle          =  x(2*nsets+3);
alphai_offset       =  x(2*nsets+4);
delta_Si            =  x(2*nsets+5);
beta_Si             =  x(2*nsets+6);
sigma_Si      =  x(2*nsets+7);
asymmetry_Si   =  x(2*nsets+8);
CDF_Si         =  x(2*nsets+9);
delta_buffer        =  x(2*nsets+10);
beta_buffer         =  x(2*nsets+11);
layeredp            = reshape(x(2*nsets+12:end),nlayer,[]);


%%  calcualte volume fraction profile for entire PS/D2O layer
rawedp = nan(2+nlayer,5);
rawedp(:,1:4) = [...
    0                           delta_buffer            beta_buffer         layeredp(1,4)
    cumsum(layeredp(:,1))       layeredp(:,2)           layeredp(:,3)       [layeredp(2:end,4);sigma_Si]
    NaN                         delta_Si                beta_Si             NaN];
rawedp(:,5) = interface_profile;               % add interface profile
if nnz(ismember(interface_profile,[2,3,4,5,7])) ~= 0
    rawedp(:,6) = [layeredp(:,5);asymmetry_Si;NaN];            % asymmetric factor
    rawedp(:,7) = [layeredp(:,6);CDF_Si;NaN];            % CDF roughness
end

% edp_0 = edprofile(rawedp,nSlice);
[edp_0,slice_dispersion,slice_absorption,edpstep] = edprofile(rawedp,nSlice);
edp = binedprofile(edp_0,binedp_delta_cutoff);  % bin the density profile

if plotEDPFlag == 1
    d_total = sum(layeredp(:,1));
    
    % --- plot edp
    figure
    haxes1 = subplot(2,1,1);
    hold on;  box on;
    plot(edp_0(:,1)-d_total,edp_0(:,2),'-','linewidth',2);
    plot(edp(:,1)-d_total,edp(:,2),'o');
    plot(edpstep(:,1)-d_total,edpstep(:,2),'g:','linewidth',1);     % plot step-like profile
    plot(edp_0(:,1)-d_total,slice_dispersion,'m:','linewidth',1);     % plot individual slice profile
    ylabel('\delta')    
    set(gca,'XMinorTick','on','YMinorTick','on');
    haxes2 = subplot(2,1,2);
    hold on; box on;
    plot(edp_0(:,1)-d_total,edp_0(:,3),'-','linewidth',2);
    plot(edp(:,1)-d_total,edp(:,3),'o');
    plot(edpstep(:,1)-d_total,edpstep(:,3),'b:','linewidth',1);          % plot step-like profile
    plot(edp_0(:,1)-d_total,slice_absorption,'m:','linewidth',1);         % plot individual slice profile
    ylabel('\beta')    
    xlabel(['depth (',char(197),')']);    
    set(gca,'XMinorTick','on','YMinorTick','on');    
    linkaxes([haxes1,haxes2],'x');
    % shift and output
    edp_0_shift = edp_0;
    edp_0_shift(:,1) = edp_0_shift(:,1) - d_total;
    assignin('base','edp_0_shift',edp_0_shift);
    assignin('base','edpstep',edpstep);
    assignin('base','d_total',d_total);
    assignin('base','slice_dispersion',slice_dispersion);
    assignin('base','slice_absorption',slice_absorption);
  
%     EDP = [edp_0(:,1)-d_total, edp_0(:,2), edp_0(:,3)];
%     save(strcat(fname(1:end-4),'_EDP.dat'),'EDP','-ascii','-tabs');
    
%     
%         z = edp(1,1);
%     edp_step = rawedp(1,2:3);
%     for ii=1:nlayer+2
%         z = [z; rawedp(ii,1); rawedp(ii,1)];
%         edp_step = [edp_step; rawedp(ii,2:3); rawedp(ii+1,2:3)];
%     end
%     z = [z; edp(end,1)];
%     edp_step = [edp_step; edp(end,2:3)];
%     
%     % plot step function    
%     figure
%     haxes1 = subplot(2,1,1);
%     hold on;  box on;
%     plot(edp_0(:,1)-d_total,edp_0(:,2),'-','linewidth',2);
%     plot(edp(:,1)-d_total,edp(:,2),'o');
%     plot(z-d_total,edp_step(:,1),'g:','linewidth',1);
%     ylabel('\delta')    
%     set(gca,'XMinorTick','on','YMinorTick','on');
%     haxes2 = subplot(2,1,2);
%     hold on; box on;
%     plot(edp_0(:,1)-d_total,edp_0(:,3),'-','linewidth',2);
%     plot(edp(:,1)-d_total,edp(:,3),'o');
%     plot(z-d_total,edp_step(:,2),'b:','linewidth',2); 
%     ylabel('\beta')    
%     xlabel(['depth (',char(197),')']);    
%     set(gca,'XMinorTick','on','YMinorTick','on');    
%     linkaxes([haxes1,haxes2],'x');
%     % shift and output
%     edp_0_shift = edp_0;
%     edp_0_shift(:,1) = edp_0_shift(:,1) - d_total;
%     assignin('base','edp_0_shift',edp_0_shift);
% %     EDP = [edp_0(:,1)-d_total, edp_0(:,2), edp_0(:,3)];
% %     save('EDP.dat','EDP','-ascii','-tabs')

    % --- plot step like EDP
    figure
    haxes1 = subplot(2,1,1);
    hold on;  box on;
%     plot(edp_0(:,1)-d_total,edp_0(:,2),'-','linewidth',2);
%     plot(edp(:,1)-d_total,edp(:,2),'o');
    plot(edpstep(:,1)-d_total,edpstep(:,2),'-','linewidth',1);     % plot step-like profile
    plot(edp_0(:,1)-d_total,slice_dispersion,':','linewidth',1);     % plot individual slice profile
    ylabel('\delta')    
    set(gca,'XMinorTick','on','YMinorTick','on');
    legend
    haxes2 = subplot(2,1,2);
    hold on; box on;
%     plot(edp_0(:,1)-d_total,edp_0(:,3),'-','linewidth',2);
%     plot(edp(:,1)-d_total,edp(:,3),'o');
    plot(edpstep(:,1)-d_total,edpstep(:,3),'-','linewidth',1);          % plot step-like profile
    plot(edp_0(:,1)-d_total,slice_absorption,':','linewidth',1);         % plot individual slice profile
    ylabel('\beta')    
    xlabel(['depth (',char(197),')']);    
    set(gca,'XMinorTick','on','YMinorTick','on');    
    linkaxes([haxes1,haxes2],'x');
end

%% calculate reflectivity for every set
I_cal_set = cell(1,nsets);
indcov_set = cell(1,nsets);
for iset = 1:nsets
    alphai = alphai_set{iset}+alphai_offset;    
    qz = 4*pi/lambda*sind(alphai);
    if resdawave == 0 && resconstant == 0 || conv_method == 0%--- no reoslution
        RR = parratt(qz,edp,lambda);
        I_cal = RR(:,2);
        I_cal = I_cal/Inorm + Ibk;
        index = alphai<foot_angle;
        I_cal(index) = I_cal(index).*sind(alphai(index))/sind(foot_angle);
        indconv = (1:length(alphai))';
    else
        % --- with resolution
        resfactor = 3;
        nconvpoint = 11;
        RR_conv =  refconv(qz,edp,lambda,resdawave,resconstant,resfactor,nconvpoint,conv_method);
        I_cal = RR_conv(:,2);
        indconv = RR_conv(:,3);       % indexs of the qz after convolution
        alphaiconv = alphai(indconv);
        I_cal = I_cal/Inorm + Ibk;
        index = alphaiconv<foot_angle;
        I_cal(index) = I_cal(index).*sind(alphaiconv(index))/sind(foot_angle);
    end
    I_cal_set{iset} = I_cal;
    indcov_set{iset} = indconv;
    qset{iset} = qz;
end

if nargout == 1
    varargout{1} = I_cal_set;
elseif nargout == 2
    varargout{1} = I_cal_set;
    varargout{2} = qset;
elseif nargout == 3
    varargout{1} = I_cal_set;
    varargout{2} = qset;
    varargout{3} = indcov_set;
elseif nargout == 4
    varargout{1} = I_cal_set;
    varargout{2} = qset;
    varargout{3} = indcov_set;
    fresnel_set = cell(1,nsets);
    for iset = 1:nsets
        tmp = fresnel(alphai_set{iset}+alphai_offset,delta_Si,beta_Si,lambda);
        fresnel_set{iset} = tmp(:,2);   % 2nd column is fresnel
    end
    varargout{4} = fresnel_set;
else
    error('Invalid number of output argument');
end
