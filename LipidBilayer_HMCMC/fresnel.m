function f = fresnel(alpha,delta,beta,lambda)
% FRESNEL Calculate fresnel reflectivity
%   F = FRESHNEL(ALPHA,DELTA,BETA,LAMBDA)
% 
%   Input argument:
%       ALPHA: list (Mx1) of incident angles
%       DELTA: dispersion of the subphase
%       BETA: absorption of the subphase
%       LAMBDA: beam wavelength (unit: A)
%
%   Output argument:
%       F: list (Mx2) with col 1 and 2 for qz and fresnel reflectivity. 
% 
%   See also parratt
%
%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2015/12/15 $

edp = [
    0           0           0
    NaN         delta       beta];
qz = 4*pi/lambda*sind(alpha);
f = parratt(qz,edp,lambda);
