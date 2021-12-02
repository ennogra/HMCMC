function varargout = edprofile(varargin)
% EDPROFILE Create electron density profile (EDP).
%   EDP = EDPROFILE(RAWPROFILE,NSLICE) returns continuous 
%   dispersion/electron density profile (EDP) using using effective-density
%   model.
%
%   [EDP,DISPERSION,ABSORPTION,EDPSTEP] = EDPROFILE(RAWPROFILE,NSLICE) returns EDP 
%   and depth distribution of the dispersion and absorption for each input
%   slice. It also returnsthe step-like edp.
%
%   Input argument:
%       RAWPROFILE: N+1 x (5 or 7), with each column representing as follows,
%               column 1: depth (z-coordinate, unit: A) of N physical 
%                   interfaces from top to bottom (last z can be any 
%                   arbitray value, only standing for substrate); 
%               column 2: dispersion (delta) above each interface. The last
%                   row means substrate;
%               column 3: absorption (beta) above each interface. The last
%                   row means substrat;
%               column 4: roughness of each interface (unit: A). The last 
%                   row can be arbitrary); 
%               column 5: PDF of each interface with 
%                   0 = norm (Gaussian)
%                   1 = sech2 (hyperbolic secant square)
%                   2 = skewed norm
%                   3 = skewed sech2
%                   4 = skewed norm-sech2
%                   5 = skewed sech2-norm
%                   6 = Cauchy profile
%                   7 = cumulative skewed Cauchy profile
%               column 6: asymmetry factor alpha. This column apply when 
%                   the interface profile is either 2-5 or 7.
%                   alpha=0 reduces the asymmetric profiles to symmetric
%                   ones, e.g., cumulative skewed normal  
%                   becomes error function profile, cumulative skewed 
%                   Cauchy becomes arctangential profile. Alpha<0 for
%                   negative skewness towards the interface above (medium
%                   above has less penetration to medium below the
%                   interface, a.k.a, medium below has more penetration to
%                   medium above); vice versa.
%               column 7: roughness of the CDF part.
%               Note: z increases from top to bottom. N is the total
%                  number of interfaces.
%       NSLICE: total number of slice used to discribe EDP. For example,
%           NSLICE could be 1000, if the totle sample thickness is 1000A. 
%           Each slice thickness had better be as smaller as the wavelength
%           of the probing beam, i.e 1A in case of x-ray.
%
%   Note: 
%   1. Final number of slices is determined automatically based on the
%       input number of slices NSLICE and the roughness parameters. It is
%       usually slightly larger than NSLICE.
%   2. The true electron density (Unit: A^-3) can be obtained by relation
%           rho = 2*pi*edp(:,2)/(re*lambda^2);
%       Here, re = 2.8179e-5 is the Thomson constant or electron radius 
%       (Unit: A), and lambda is the x-ray wavelength (Unit: A).
%
%   Output argument:
%       EDP: NSLICE+3 x 3 (with enforced top and bottom delta and beta)
%               column 1: the total number of interfaces plus last row for 
%                   the substrate; 
%               column 2: dispersion (delta);
%               column 3: absorption (beta).
%       DISPERSION/ABOSRPTION: NSLICE+1 x N+1
%               each column represnts the depth distribution for each slice
%               (from top to bottom).  
%       EDPSTEP: profile for sharp interfaces.
%
%   Example: A wafer of SiO2(10A thick, 4A rough) on Si (3A rough) in 
%       water; x-ray of 8kev; error function profile for roughness
%       
%       h2o  = refrac('H2O',8,1);       % require refrac package
%       sio2 = refrac('SiO2',8,2.65);
%       si   = refrac('Si',8,2.33);
%       rawedp = [0,   h2o.dispersion,  h2o.absorption,  4,  0;
%                 10,  sio2.dispersion, sio2.absorption, 3,   0;
%                 NaN, si.dispersion,   si.absorption,   NaN, NaN];
%       [edp,delta,beta,edpstep] = edprofile(rawedp,200);
%       figure; 
%       hold on
%       plot(edp(:,1),delta(:,1),'r-'); % plot dipsersion of water
%       plot(edp(:,1),delta(:,2),'g-'); % plot dipsersion of sio2
%       plot(edp(:,1),delta(:,3),'m-'); % plot dipsersion of si
%       plot(edp(:,1),edp(:,2),'b-');   % plot total dispersion
%       hold off; xlabel('depth (A)'); ylabel('dispersion');
%
% Reference: 
%   (1) Tolan, M. (1999) X-Ray Scattering from Soft-Matter Thin Films,
%       Springer, P26, Effective-density Model.
%   (2) S. Nadarajah and S. Kotz, Skewed distributions by the normal
%       kernel, Statistics and Probability Letters 65, 269 (2003)
%   (3) Z. Jiang and W. Chen, Z. Jiang and W. Chen, Generalized
%       Skew-symmetric Interfacial Probability Distribution in Reflectivity
%       and Small-angle Scattering Analysis
%
% By Zhang Jiang
% $Revision: 1.0 $  $Date: 2004/11/27 $
% $Revision: 1.1 $  $Date: 2011/11/20 $ 
%   (1) Add extra output argument for EDP of individual layers.
%   (2) Add an example to help info.
% $Revision: 2.0 $   
%   (1) Include cummulative skewness for asymmetric interface profiles to
%   account for near hard walls or incomplete layers. $Date: 2017/02/09 $
%   (2) Include cummulative hyperbolic tangent profile. $Date: 2017/03/18 $
%   (3) Include the 2nd roughness for the CDF of the skew-symmetric model.
%       $Data: 2017/03/30$
%
% See also BINEDPROFILE, REFRAC

rawProfile = varargin{1};               % N+1 x 5 or 6
nSlice = varargin{2};                   % number of total slices (M slices, M+1 slice interfaces)

zInterface = rawProfile(1:end-1,1)';    % 1 x N, physical interface positions
delta = rawProfile(:,2)';               % 1 x N+1, delta
beta = rawProfile(:,3)';                % 1 x N+1, beta
% sigma1 = abs(rawProfile(1:end-1,4))';         % 1 x N, roughness
sigma1 = (rawProfile(1:end-1,4))';         % 1 x N, roughness
profileFcn = rawProfile(1:end-1,5)';    % 1 x N, interface profile

skew_ind = [2,3,4,5,7];

if nnz(ismember(profileFcn,skew_ind)) && size(rawProfile,2) ~= 7    % need 7 column for skew-symmetric profiles
    error('Incorrect input raw profile');
end
if size(rawProfile,2) == 7              % for skew-symmetric profiles   
    alpha   = rawProfile(1:end-1,6)';      % 1xN, asymmetry factor
    sigma2  = rawProfile(1:end-1,7)';      % 1xN, sigma2 for CDF part of the profile
elseif size(rawProfile,2) == 5          % for symetric profiles
    alpha = zeros(size(profileFcn));    % 1xN zeros for symmetric profiles
    sigma2 = alpha;                     % sigma2 is zero
else
    error('Invalid input rawprofile...');
end
%zInterface = zInterface;               % reverse z-axis direction (pointing up)

% --- redefine roughness, preventing dividing by zero rounghess below
sigma1(sigma1==0) = eps;                
if size(rawProfile,2) == 7              % for skew-symmetric profiles   
    sigma2(sigma2==0) = eps;                
end

% --- physical number of interfaces and layers excluding infinite top and substrate
nInterface = size(zInterface,2);                % equal to N
%numLayer = nInterface-1;   % equal to N-1

% --- total thickness including roughness effect (10 times of roughness)
sigmaFactor = 10;   % roughness effect to the total thickness
%thick = abs(zInterface(end) - zInterface(1));   % physical thickness
%totalThick = thick + sigmaFactor*(sigma(1)+sigma(end));

% --- position of all the slice interfaces 1 x M+1
% zSliceInterface = linspace(...
%     min(zInterface-sigmaFactor*sigma1),...
%     max(zInterface+sigmaFactor*sigma1),nSlice+1);
zSliceInterface = linspace(...
    min(min(zInterface-sigmaFactor*sigma1),min(zInterface-sigmaFactor*sigma2)),...
    max(max(zInterface+sigmaFactor*sigma1),max(zInterface+sigmaFactor*sigma2)),nSlice+1);

% --- calculate the profile derivatives
% fcnh = {@f_norm,@f_tanh,@f_skewnorm,@f_skewtanh,@f_skewnormtanh,@f_skewtanhnorm,@f_cauchy,@f_skewcauchy};
pInterface = [zInterface(:),sigma1(:),alpha(:),sigma2(:)];
pd = nan(nInterface,nSlice+1);         % profile deriavatives
for ii=1:nInterface
    switch profileFcn(ii)+1
        case 1
            pd(ii,:) = f_norm(zSliceInterface,pInterface(ii,:));
        case 2
            pd(ii,:) = f_tanh(zSliceInterface,pInterface(ii,:));            
        case 3
            pd(ii,:) = f_skewnorm(zSliceInterface,pInterface(ii,:));            
        case 4
            pd(ii,:) = f_skewtanh(zSliceInterface,pInterface(ii,:));            
        case 5
            pd(ii,:) = f_skewnormtanh(zSliceInterface,pInterface(ii,:));            
        case 6
            pd(ii,:) = f_skewtanhnorm(zSliceInterface,pInterface(ii,:));            
        case 7
            pd(ii,:) = f_cauchy(zSliceInterface,pInterface(ii,:));            
        case 8
            pd(ii,:) = f_skewcauchy(zSliceInterface,pInterface(ii,:));            
    end
%    ftmp = fcnh{profileFcn(ii)+1};
%     pd(ii,:) = ftmp(zSliceInterface,pInterface(ii,:));    
%     pd(ii,:) = fcnh{profileFcn(ii)+1}(zSliceInterface,pInterface(ii,:));
end

% --- get cumulative profile for each interface
cpd = cumtrapz(zSliceInterface,pd,2);

% --- get weight proifile for each bounded layer (excluding the topmost and
% bottomost unbounded layers) at its two boundary interface; then for for
% each layer find the zeta position/ind where the two weight profiles merge 
if nInterface>=2     % if there is bounded layer
    cpd_bl_t = cpd(1:end-1,:);    % cpd for bounded layer at the interface above
    cpd_bl_b = 1-cpd(2:end,:);    % cpd for bounded layer at the interface below
    wp = cpd_bl_t;  % initialize wegith profile with  cpd_bl_t
    ind = (cpd_bl_t > cpd_bl_b);
    wp(ind) = cpd_bl_b(ind); 
else
    wp = [];    % no weight profile for single interface
end    

% --- get the weight profile for topmost and bottom-most phases
wp = [1-cpd(1,:);wp;cpd(end,:)];

% --- construct delta profile
sumwp = sum(wp,1);
allDelta = delta*wp./sumwp;
allBeta  = beta*wp./sumwp;


% % add top reference and substrate to both ends to force edp to converge
% sliceThick = abs(zSliceInterface(2)-zSliceInterface(1));
% zSliceInterface = [...
%     zSliceInterface(1)-sliceThick,...
%     zSliceInterface,...
%     zSliceInterface(end)+sliceThick];
% allDelta = [...
%     delta(1),...
%     allDelta,...
%     delta(end)];
% allBeta = [...
%     beta(1),...
%     allBeta,...
%     beta(end)];

% construct edp
edp = [zSliceInterface(:),allDelta(:),allBeta(:)];

% replace the 1st and last with with top and bottom delta and beta to
% enforce an converge
edp([1,end],2:3) = [delta([1,end])',beta([1,end])'];


% output
if nargout == 1
    varargout{1} = edp;
elseif nargout == 4
    deltaEachSlice = repmat(delta(:),1,nSlice+1).*wp./sumwp;
    betaEachSlice = repmat(beta(:),1,nSlice+1).*wp./sumwp;
    % get the box profile edp_step
    z_step = [edp(1,1);reshape(repmat(zInterface,[2,1]),[],1);edp(end,1)];
    delta_step = reshape(repmat(delta,[2,1]),[],1);
    beta_step = reshape(repmat(beta,[2,1]),[],1);
    edp_step = [z_step, delta_step, beta_step];
    
    varargout{1} = edp;
    varargout{2} = deltaEachSlice';
    varargout{3} = betaEachSlice';
    varargout{4} = edp_step;
else
    error('Invalid output argument');
end
% --- EOF

% % --- debugging plots
% figure
% hold on;
% for ii=1:nInterface
%     plot(zSliceInterface,pd(ii,:));
%     plot(zSliceInterface,cpd(ii,:));
% end
% if nInterface>=2
% for ii=1:size(cpd_bl_t,1)
%     plot(zSliceInterface,cpd_bl_t(ii,:),'k.-');
%     plot(zSliceInterface,cpd_bl_b(ii,:),'k.-');    
% end
% end
% for ii=1:nInterface
%     plot([zInterface(ii),zInterface(ii)],[0,1],'g--');
% end
% 
% for ii=1:nInterface+1
%     plot(zSliceInterface,wp(ii,:),'r--')
% end
% set(gca,'ylim',[-0.1,1.1]);
% figure
% plot(zSliceInterface,allDelta);
% hold on
% for ii=1:nInterface+1
%     plot(zSliceInterface,deltaEachSlice(ii,:),'--')
% end

function y = f_skewnorm(x,p)
x0 = p(1);
sigma1 = p(2);
alpha = p(3);
sigma2 = p(4);
g = 1/sqrt(2*pi)/sigma1*exp(-(x-x0).^2/2/sigma1^2);   
cg = 1/2*(1+erf(alpha*(x-x0)/sigma2/sqrt(2)));
y = 2*g.*cg;

function y = f_skewtanh(x,p)
x0 = p(1);
sigma1 = p(2);
alpha = p(3);
sigma2 = p(4);
t = pi/(4*sqrt(3)*sigma1)./ (cosh(pi/2/sqrt(3)*(x-x0)/sigma1)).^2;
ct = 1/2*(1+tanh(alpha*pi/2/sqrt(3)*(x-x0)/sigma2));
y = 2*t.*ct; 

function y = f_skewnormtanh(x,p)
x0 = p(1);
sigma1 = p(2);
alpha = p(3);
sigma2 = p(4);
g = 1/sqrt(2*pi)/sigma1*exp(-(x-x0).^2/2/sigma1^2);   
ct = 1/2*(1+tanh(alpha*pi/2/sqrt(3)*(x-x0)/sigma2));
y = 2*g.*ct; 

function y = f_skewtanhnorm(x,p)
x0 = p(1);
sigma1 = p(2);
alpha = p(3);
sigma2 = p(4);
t = pi/(4*sqrt(3)*sigma1)./ (cosh(pi/2/sqrt(3)*(x-x0)/sigma1)).^2;
cg = 1/2*(1+erf(alpha*(x-x0)/sigma2/sqrt(2)));
y = 2*t.*cg; 

function y = f_skewcauchy(x,p)
x0 = p(1);
gamma1 = p(2);
alpha = p(3);
gamma2 = p(4);
c = gamma1/pi./((x-x0).^2+gamma1^2);
cc = 1/2+atan(alpha*(x-x0)/gamma2)/pi;
y = 2*c.*cc;


function y = f_norm(x,p)
x0 = p(1);
sigma1 = p(2);
y = 1/sqrt(2*pi)/sigma1*exp(-(x-x0).^2/2/sigma1^2);   

function y = f_tanh(x,p)
x0 = p(1);
sigma1 = p(2);    
y = pi/(4*sqrt(3)*sigma1)./ (cosh(pi/2/sqrt(3)*(x-x0)/sigma1)).^2;

function y = f_cauchy(x,p)
x0 = p(1);
gamma1 = p(2);       % HWHM of Lorentzian
y = gamma1/pi./((x-x0).^2+gamma1^2);