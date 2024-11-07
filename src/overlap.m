function [gOver,g1,g2] = overlap(nu0, res, range, aL, dnuG1, dnuL1, dnuG2, dnuL2)
% overlap - This function computes two lineshapes and overlaps them.
%
%   [gOver,g1,g2] = overlap(nu0, res, range, aL, dnuG1, dnuL1, dnuG2, dnuL2)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function computes the spectral overlap or interaction between a 
%  laser profile and a %  transition profile. The profiles are modeled as 
%  combinations of Lorentzian and Gaussian functions, and this function 
%  quantifies their overlap. The function here uses the McLean model of 
%  calculating the Voigt line-shapes for laser and transition. Those Voigt 
%  line-shapes are then convoluted to obtain the overlap line-shape. The
%  returned vector is the higher half of the overlap function. 
%  
%  e.g.
%  gOver = overlap(nu0, res, range, aL, dnuG1, dnuL1, dnuG2, dnuL2);
%  pos = abs(laserPos-nu0-shift)/resolution;
%  gamma = interp1(0:length(gOver)-1, gOver, pos, 'linear', 'extrap');
%
% INPUTS:
%   nu0    -  The central frequency or reference frequency.
%   res    -  The resolution of the analysis.
%   range  - The spectral range. It can be an array of frequencies.
%   aL     -   The Lorentzian amplitude.
%   dnuG1  -  The Gaussian component of the transition profile.
%   dnuL1  -  The Lorentzian component of the transition profile.
%   dnuG2  - The Gaussian component of the laser profile.
%   dnuL2  - The Lorentzian component of the laser profile.
%
% OUTPUT:
%   gOver - The calculated spectral overlap between the signal 1 and signal 2.
%   g1    - The calculated lineshape 1.
%   g2    - The calculated lineshape 2.
%
% SEE ALSO:
%   voigtlineMcLean, dnuGFun, collisionalBroadening, excitationSpec, 
%   emissionSpec, absorptionSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments 
    nu0 (1,1) double
    res (1,1) double
    range (:,:) double
    aL (1,1) double
    dnuG1 (1,1) double
    dnuL1 (1,1) double
    dnuG2 (1,1) double = 0
    dnuL2 (1,1) double = 0
end

nu  = -fliplr(range-nu0);
g1 = voigtlineMcLean(nu,nu0,dnuG1,dnuL1,aL);
g2 = voigtlineMcLean(nu,nu0, dnuG2, dnuL2, aL);
gOver = conv(g2, g1)*res;

i = floor(length(gOver)/2);
if i == 0
    error('Error: overlap.m')
end
 
gOver = gOver(i+1:end);
gOver = gOver./max(gOver(:));
end


% G = @(v,FWHM) exp(-v.^2./(2*(FWHM./sqrt(8*log(2))).^2))./((FWHM./sqrt(8*log(2)))*sqrt(2*pi));
% L = @(v,FWHM) FWHM./(pi*(v.^2+FWHM.^2/4));
% gTrans = conv(G(range,dnuG),L(range,dnuL))*res;
% gLaser = conv(G(range,dnuGL),L(range,dnuLL))*res;