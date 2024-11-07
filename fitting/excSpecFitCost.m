function cost = excSpecFitCost(...
    wnum,...
    measSpec,...
    linelist,...
    MM,...
    T,...
    dnuGL,...
    dnuLL,...
    specParams, ...
    fitParams,...
    collParam...
    )
% excSpecFitCost - Calculates cost between simulated and measured excitation spectrum.
%
%   cost = excSpecFitCost(...
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   The difference between a measured and simulated spectrum is calculated with 
%  a baseline correc-tion and a scaling factor. This function is used as an error function 
%  for fitting with lsqnonlin. A baseline is calculated from the user parameters and then 
%  subtracted from the measured spectrum. The excitation parameters are used to simulate 
%  the spectrum, which is scaled by the given factor. 
%
% INPUTS:
%   measSpec     - The measurement data points, which should have the same length as wnum. 
%   The inputs (wnum, linelist, MM, T, dnuGL, dnuLL, resFactor, Z, limit) are forwarded to
%   the function excitationSpec. Also, if colls is not available, then dnuSh and dnuL of 
%   specParams are forwarded to excitationSpec, if passed to the current function.
%
% OPTIONAL INPUTS:
%   a          - The slope of the baseline (default = 1).
%   offset     - The base line offset (default = 0). 
%   scale       - The scaling factor for the simulated spectrum to match the measured data (default = 1).
%   collParam   
%   -	colls and gas - The collision model and the gas composition. If colls is empty the 
%       quenching and collisional broadening are not calculated, then the given or default 
%       dnuSh and dnuL value would not be overridden. If colls and gas are available and dnuSh 
%       is given by user, then this given shift will be added to the collisional shift. 
%   -	P             - Pressure in bar. This value is only considered if colls are available.
%   
%
% OUTPUT:
%   cost     - The difference between the measured and simulated spectra with respect to the given parameters.
%   
%
% SEE ALSO:
%   fitExcitationSpec, lsqnonlin, fitImageDataset, calBaseline, fitExcSpec-trumImage.mlx, 
%   fitExcSpectrumSingle.mlx, lineShiftCost
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    wnum (:,:) double
    measSpec (:,1) double
    linelist (:,:) cell
    MM (1,1) double
    T (1,1) double
    dnuGL (1,1) double
    dnuLL (1,1) double
    
    specParams.dnuL (1,1) double
    specParams.dnuSh (1,1) double
    specParams.resFactor (1,1) double 
    specParams.Z (:,:) double
    specParams.limit (1,1) double

    fitParams.a (1,1) double = 0
    fitParams.offset (1,1) double = 0
    fitParams.scale (1,1) double = 1

    collParam.colls = {}
    collParam.gas = {}
    collParam.P (1,1) double = 1
end

if ~isempty(collParam.colls)

    [specParams.dnuL, dnuSh] = collisionalBroadening(...
        collParam.gas, collParam.colls, collParam.P, T);
    [~,linelist] = quenchRate(...
        collParam.gas, collParam.colls, T, collParam.P, MM, linelist);
    specParams.dnuSh = dnuSh + specParams.dnuSh;
end

p =  namedargs2cell(specParams);

baseline = calBaseline(fitParams.a,fitParams.offset,length(measSpec));

spec = excitationSpec(wnum, linelist, MM, T, dnuGL, dnuLL, p{:}, normalize=true);

cost =  (measSpec-baseline-fitParams.scale.*spec);

end

