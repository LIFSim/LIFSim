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
%   gOver - The calculated spectral overlap between the lineshape 1 and
%           lineshape 2.
%   g1    - The calculated lineshape 1.
%   g2    - The calculated lineshape 2.
%
% SEE ALSO:
%   voigtlineMcLean, dnuGFun, collisionalBroadening, excitationSpec, 
%   emissionSpec, absorptionSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
%
% CHANGES:
% 2025-02-13: 
%   Change 1: Normalizing g1 and g2 according to Partridge and Normand, 
%             aka Y and L their paper: https://doi.org/10.1364/AO.34.002645
%   Change 2: Removed the normalization of gOver.
%             Changes by Abbas El Moussawi.
%
% 2025-02-28: 
%   Change 3: Added normalization for gOver through numerical integration
%             so that changing the linewidth does not affect the integrated
%             absorbance. 
%             By Weitian Wang.

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
%%% Start change 1 
g1 = g1./max(g1);
g2 = g2./max(g2);
%%% End change 1
gOver = conv(g2, g1).*res;

i = floor(length(gOver)/2);
if i == 0
    error('Error: overlap.m')
end


%%% Start change 3
gOver = gOver/(trapz(gOver)*res); 
% Norm the line shape by numerical intergration to make gamma independent of line strength. 
% A better solution can be to find the relationship between the numerical value
% of gamma and the parameters.  
%%% End change 3

gOver = gOver(i+1:end);
end


% G = @(v,FWHM) exp(-v.^2./(2*(FWHM./sqrt(8*log(2))).^2))./((FWHM./sqrt(8*log(2)))*sqrt(2*pi));
% L = @(v,FWHM) FWHM./(pi*(v.^2+FWHM.^2/4));
% gTrans = conv(G(range,dnuG),L(range,dnuL))*res;
% gLaser = conv(G(range,dnuGL),L(range,dnuLL))*res;
