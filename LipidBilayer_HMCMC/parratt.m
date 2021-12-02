function pardata = parratt(varargin)
% PARRATT Calculates reflectivity of a multilayer by Parratt's recursion
%   PARDATA = PARRATT(QZ,PROFILE,WAVELENGTH) calculates the refelctvity data
%       given qz span, sample profile and input wavelength (in vacuum).
%
%   Format of input:
%       QZ: 1 x M list of qz in the reference material in A-1.
%       PROFILE: N+1 x 3. 1st column z, 2nd delta (dispersion),3rd beta
%               (absorption). 1st row reference material. Last row substrate.
%               The last value in z-column can be any arbitrary number, just
%               stands for the semifinite substrate.
%       WAVELENGTH: wavelength of x-ray or neutron in the vacuum in A
%
%   Format of output:
%       PARDATA: M x 2 with 1st column qz, 2nd reflectivity.
%
% Reference:
% 1) Tolan, M. (1999) X-Ray Scattering from Soft-Matter Thin Films, Springer
% 2) Parratt, L.G. (1954). Phys. Rev. 95, 359-369
% 3) Als-Nielsen, J. (2001) Elements of Modern X-ray Physics, Wiley John & Sons
%
% Copyright 2004 Zhang Jiang
% $Revision: 1.0 $  $Date: 2004/11/18 17:49:30 $
% $Revision: 1.1 $  $Date: $ To be checked. Wrong code when top layer has larger EDP;

% See also REFCONV

qz = varargin{1};                   % qz list (M points)
profile = varargin{2};              % profile of sample; N+1 points with 
                                    % 1st the media above sample and last 
                                    % point the substrate
wavelength = varargin{3};           % wavelength of beam

nLayer = size(profile,1);           % N+1 layers including 'air' and substrate
z = profile(:,1);                  % slice position; the last z(N+1) dos not
                                    % count, simply stands for substrate
delta = profile(:,2);
beta = profile(:,3);
k = 2*pi/wavelength;                % wave vector in vacuum

% --- transfer qz (1 x M)
if size(qz,2) == 1
    qz = qz';
end

% --- construct refractive index (N+1 x 1)
n = 1-delta+1i*beta;

% --- get qz in vacuum (N+1 x 1)
qzInVac = sqrt(qz.^2+4*k^2-4*k^2*n(1)^2);

% --- Wavevector in each layer (N+1 x M)
QLayer = sqrt(4*(n*ones(1,size(qzInVac,2))).^2*k^2-4*k^2+(ones(size(n,1),1)*qzInVac).^2);

% --- Reflectance coefficients (no multiple scattering) (N x M)
r = (QLayer(1:end-1,:)-QLayer(2:end,:))./(QLayer(1:end-1,:) + QLayer(2:end,:));

% --- Reflectacne from more layers (with multiple scattering) (1 x M)
R = r(end,:);                           % reflectance from substrate
% if there are more than two layers, do recursion
if nLayer > 2
    thick = diff(z);    %  thickness of each slab (1 x N-1)
    thick = thick(1:end-1);
    for iInterface = nLayer-2:-1:1          % N-1 reflectances needed to be calculated
        R = (r(iInterface,:)+R.*exp(1i*QLayer(iInterface+1,:)*thick(iInterface)))./...
            (1+r(iInterface,:).*R.*exp(1i*QLayer(iInterface+1,:)*thick(iInterface)));
    end
end

% --- reflectivity
% RR = abs(R).^2;     % 1 x M
RR = (real(R).^2 + imag(R).^2);
pardata = [qz' RR'];     % M x 2