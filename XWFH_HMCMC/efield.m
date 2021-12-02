function varargout = efield(varargin)
% EFIELD Spatial distribution of x-ray excited complex electric field.
%
%   EF = EFIELD(ALPHA,PROFILE,WAVELENGTH) calculates spatial distribution 
%       of x-ray excited electric field (normalized to the incident beam)
%       in films by Parratt's recursion.
%
%   [EF,ER,ET,KZ] = EFIELD(ALPHA,PROFILE,WAVELENGTH) also returns R and T, 
%       which are the electric field amplitudes of reflected and transmitted
%       waves within  the film.
%
%   Input Argument:
%       ALPHA       : 1 x M list. Incident/exit angle in the reference
%               material. Unit: degree.
%       PROFILE     : N+1 x 3 matrix. Film depth, dispersion (delta) and
%               absorption (beta) profile. N is the total number of
%               interfaces. 1st column: depth (Unit: A) of each interface
%               with z-axis pointing down into the film. 2nd and 3rd 
%               columns: dispersion (delta) and absorption (beta) 
%               coefficients in the slab above the correponding interface 
%               defined in the 1st column. NOTE: PROFILE(1,2:3) are delta 
%               and beta of the reference material above interface 
%               PROFILE(1,1); PROFILE(N+1,2:3) are used as delta and beta 
%               of the substrate; PROFILE(N+1,1) can be any arbitrary 
%               number and is not used in the calculation.
%       WAVELENGTH  : Wavelength of x-ray in the vacuum. Unit: A.
%
%   Output Argument:
%       EF:     N+1 x M. Complex electric field.
%       ER:     N+1 x M. Complex electric field amplitude of reflected waves.
%       ET:     N+1 x M. Complex electric field amplitude of transmitted waves.
%       KZ:     N+1 x M. Complex kz within the film.
%   
%   Note:
%       |ER(1,:)|^2 is reflectivity
%
%   Reference:
%   1) Tolan, M. (1999) X-Ray Scattering from Soft-Matter Thin Films, Springer
%   2) Parratt, L.G. (1954). Phys. Rev. 95, 359-369
%   3) Als-Nielsen, J. (2001) Elements of Modern X-ray Physics, Wiley John & Sons
%   4) Jiang, Z. (2007) Ph.D thesis
%
%   Zhang Jiang
%   $Revision: 1.0 $  $Date: 2005/05/03$
%   $Revision: 1.1 $  $Date: 2005/06/25$ : Add options to output electric
%       field amplitudes of reflected and transmitted waves, and kz within the
%       film as a function of depth and angle.
%   $Revision: 2.0 $  $Date: 2007/09/21$ : Small deviations found in ER,
%       therefore ET and EF by Tolan's method from those obtained by
%       conventional method (parratt and Als-Nielsen). Rewrite the code
%       using formulas described in Als-Nielsen's book and Jiang's thesis.
%   $Revision: 3.0 $ $Date: 2008/09/29$: Use Tolan's recursive method.
%
% See also REFCONV, PARRATT, EDPROFILE, REFRAC

%% get input argument
if nargin ~=3 
    error('Invalid input argument.');
end
alpha   = varargin{1};              % alpha in the first medium (M points)
edp = varargin{2};                  % profile of sample; N+1 points with 
                                    % 1st the media above sample and last 
                                    % point the substrate
wavelength = varargin{3};           % wavelength of beam in vacuum
nLayer = size(edp,1);               % N+1 layers including 'vacuum' and substrate
z = -edp(:,1);                      % interface position 
delta = edp(:,2);
beta = edp(:,3);
k = 2*pi/wavelength;                % wave vector in vacuum

% --- construct refractive index (N+1 x 1)
n = 1-delta+1i*beta;

% --- kz in each layer (N+1 x M)
kz = k*sqrt( (n*ones(1,length(alpha))).^2 - (ones(length(n),1)*cos(alpha*pi/180)).^2);

%% --- Fresnel coefficients (no multiple scattering) (N x M)
r12 = ( kz(1:end-1,:) - kz(2:end,:) ) ./ ( kz(1:end-1,:) + kz(2:end,:) );
% ind_real = abs(real(r12)) < eps;
% r12(ind_real) =  1i*imag(r12(ind_real));
% ind_imag = abs(imag(r12)) < eps;
% r12(ind_imag) = real(r12(ind_imag)) ;
t12 = 1+r12;
r21 = -r12;
t21 = 1+r21;

%% --- calculate X in each layer
X = NaN*ones(nLayer,length(alpha));     % N+1 x M
X(end,:) = 0;
for jj = nLayer-1:-1:1
    X(jj,:) = r12(jj,:)+X(jj+1,:).*exp(2*1i*kz(jj+1,:)*z(jj));
    X(jj,:) = X(jj,:)./ (1+r12(jj,:).*X(jj+1,:).*exp(2*1i*kz(jj+1,:)*z(jj)));
    X(jj,:) = X(jj,:).*exp(-2*1i*kz(jj,:)*z(jj));
end

%% --- calculate R and T in each layer
R = NaN*ones(nLayer,length(alpha));        % N+1 x M
T = NaN*ones(nLayer,length(alpha));        % N+1 x M
R(1,:) = X(1,:);
T(1,:) = 1;
for jj=1:nLayer-1
    R(jj+1,:) = R(jj,:) .* exp( -1i * ( kz(jj+1,:) - kz(jj,:) ) * z(jj) );    
    R(jj+1,:) = R(jj+1,:) + T(jj,:) .* r21(jj,:) .* exp( -1i * ( kz(jj+1,:) + kz(jj,:) ) * z(jj) );
    R(jj+1,:) = R(jj+1,:) ./ t21(jj,:);
    T(jj+1,:) = T(jj,:) .* exp( 1i * (kz(jj+1,:) - kz(jj,:) ) * z(jj) );
    T(jj+1,:) = T(jj+1,:) + R(jj,:) .* r21(jj,:) .* exp( 1i * ( kz(jj+1,:) + kz(jj,:) ) * z(jj) );
    %T(jj+1,:) = T(jj+1,:) + R(jj,:) .* r21(jj,:) .* exp( 1i * real( kz(jj+1,:) + kz(jj,:) ) * z(jj) );
    T(jj+1,:) = T(jj+1,:) ./ t21(jj,:);
end

            
%% -- calcuate complex electric field
EF = T.*exp(-1i*kz.*(z*ones(1,length(alpha)))) + R.*exp(1i*kz.*(z*ones(1,length(alpha))));
%EF = T.*exp(-1i*kz.*(z*ones(1,length(alpha)))) + R.*exp(1i*real(kz).*(z*ones(1,length(alpha))));

%% --- output
varargout{1} = EF;
varargout{2} = R;
varargout{3} = T;
varargout{4} = kz;
return;