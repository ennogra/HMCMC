function convpardata = refconv(varargin)
% REFCONV Calculates resolution convoluted reflectivity of a multilayer by
%       Parratt's recursion.
%   CONVPARDATA =  REFCONV(QZ,PROFILE,WAVELENGTH,RESDWAVE,RESCONSTANT,RESFACTOR,NCONVPOINT,METHOD)
%
%   Format of input:
%       QZ: 1 x M list of qz in the reference material (usually air). Unit: A-1.
%       PROFILE: N+1 x 3. 1st column z, 2nd delta (dispersion),3rd beta
% %               (absorption). 1st row reference material. Last row substrate.
%               The last value in z-column can be any arbitrary number and
%               will not be used.
%       WAVELENGTH: wavelength. Unit: A
%       RESDWAVE: wavelength distribution width, d(wavelength)/wavelength.
%               A typical value is ~10^-3 for monochromatic beam.
%       RESCONSTANT: constant resolution. Unit: A^-1.
%       RESFACTOR (>=3): determines the qz range to do resolution
%               [<Q(i)>-RESFACOTR*resolution,<Q(i)>+RESFACTOR*resolution].
%               Typically at least three times of the resolution is needed
%               for a sufficiently accurate convolution.
%       NCONVPOINT (>=10): number of points in the above range for convolution.
%               Typically at least 10 points are needed for a sufficiently
%               accurate convolution. Also at least twice of RESFACTOR in
%               order to keep enough point density.
%       METHOD (1/2/3): 1/2 are function-based convolutions using matrix
%               vectorization for loop, respectivly. 1 is faster than 2.
%               3 is data-based convolution, which is the fastest of the
%               three. It uses Matlab built-in convolution function conv. So
%               it requires an evenly distributed data and after convoluiton
%               the beginning and end of the data will be removed because
%               the lack of the data out of the range. Default: 1
%
%   Format of output:
%       CONVPARDATA: M x 3, 1st column qz, 2nd the convoluted reflectivity,
%               3rd is the data index range (used for METHOD 3 to indicate
%               the data points kept after the convolution).
%
% Referece: Jan Sko Pedersen, J. Appl. Cryst. (1994). 27, 36-49
%
% Zhang Jiang
% $Revision: 1.0 $  $Date: 2004/11/18 $
% $Revision: 1.1 $  $Date: 2008/05/14 $
%
% See also PARRATT

qz = varargin{1};                   % qz list (M points)
profile = varargin{2};              % profile of sample; N+1 points with
% 1st the media above sample and last
% point the substrate
wavelength = varargin{3};           % wavelength of beam
resDWave = varargin{4};             % d(wavelength)/(wavelength) (unitless)
resConstant = abs(varargin{5});          % constant resolution in A^-1
resFactor = varargin{6};            % resolution factor
nConvPoint = varargin{7};           % total number of points

if nargin == 8
    conv_method = varargin{8};
else
    conv_method = 1;
end

% --- transfer qz (M x 1)
if size(qz,1) == 1
    qz = qz';
end

switch conv_method
    case 1
        % --- Method 1: vectorization
        resVar = sqrt((resDWave/(2*sqrt(2*log(2)))*qz).^2 + resConstant^2); %calculate the varaince of Guassian smearing function for each point (Mx1)
        resStep = 2*resVar*resFactor/(nConvPoint-1);                    % Mx1
        qzOffset = resVar*linspace(-resFactor,resFactor,nConvPoint);    % Mxn
        Qz = qz*ones(1,nConvPoint) + qzOffset;                          % Mxn
        rr = parratt(Qz(:),profile,wavelength);                         % Mnx2
        rr = reshape(rr(:,2),length(qz),nConvPoint);                    % Mxn
        g  = 1./(sqrt(2*pi)*(resVar*ones(1,nConvPoint))) .*...
            exp(-qzOffset.^2./(2*(resVar*ones(1,nConvPoint)).^2)); %Mxn
        convRR = sum(rr.*g,2).*resStep;     % Mx1
        indqz = 1:length(qz);
        convpardata = [qz convRR indqz'];
    case 2
        % ---Method 2: "for" loop
        resVar = sqrt((resDWave/(2*sqrt(2*log(2)))*qz).^2 + resConstant^2); %calculate the varaince of Guassian smearing function for each point (Mx1)
        convRR = ones(length(qz),1)*NaN;
        for iqz = 1:length(qz)
            resWindow = resVar(iqz)*resFactor;
            resStep   = 2*resWindow/(nConvPoint-1);
            Qz = linspace(qz(iqz)-resWindow,qz(iqz)+resWindow,nConvPoint)';
            rr = parratt(Qz,profile,wavelength);
            g  = 1./(sqrt(2*pi)*resVar(iqz))*exp(-(Qz-qz(iqz)).^2/(2*resVar(iqz)^2));
            convRR(iqz) = sum(rr(:,2).*g*resStep);
        end
        indqz = 1:length(qz);
        convpardata = [qz convRR indqz'];
    case 3
        % --- Method 3: (the fastest method) use conv only when resolution is constant (q-independece), i.e. resVar is a constant
        % Note:
        % 1. Points near the boundaries of the qz list will be removed.
        % 2. nconvpoint is not used
        rr = parratt(qz,profile,wavelength);
        resWindow = resConstant*resFactor;
        nconv = round(resWindow/(qz(2)-qz(1)));
        qzconv = (-nconv:nconv)*(qz(2)-qz(1));  % qzconv has to be odd number and step size identical to the data
        g = 1./(sqrt(2*pi)*resConstant)*exp(-qzconv.^2/(2*resConstant^2));
        r = [flipud(rr(:,2));rr(:,2)];      % mirror the reflectivity
        convr = conv(r,g);                  % convolute with resolution
        convRR = convr(length(convr)/2+(1:length(r)/2));
        if nconv == 0
            convRR = convRR / g;
        else
            convRR = convRR*(qz(2)-qz(1));      % rescale back using qz step
        end
        n_cut = round(3*resConstant/(qz(2)-qz(1)));     % calcuate the number points needed to be removed at the data boundaries
        indqz = n_cut+1:length(qz)-n_cut;
        qz_cut = rr(indqz,1);   % remove qz at boundaires
        convRR = convRR(indqz); % remove convRR at the boundaries
        convpardata = [qz_cut convRR indqz'];
end