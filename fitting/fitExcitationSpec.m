function [solutionVals, residue]= fitExcitationSpec(...
    wnum, measSpec, linelist, MM, bestGuess, options, fit, specParams, collParam)
% fitExcitationSpec - Fit optimization function for excitation spectrum.
%
%   fitExcitationSpec(wnum, measSpec, linelist, MM, bestGuess, options, fit, specParams, collParam)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This is the fit optimization function for a measured excitation spectrum. This function uses 
% the optimization function lsqnonlin with a user-given options, see getFitOptions. Here, the 
% expected fit parameters are seven with lower and upper boundaries vectors
%
% INPUTS:
%   wnum           - The wavenumber region in 1/cm
%   measSpec       - The measured spectrum data point, with the same length as wnum.
%   linelist       - he linelist obtained from function selectLines.
%   MM             - Molar mass of the species, obtained from function molMass or selectLines.
%   bestGuess      - A matrix vector with seven places, which are the initial fit parameters to start the optimization. 
%   options        - The options objects needed by the optimization function lsqnonlin. Some presets can be obtained from the function getFitOptions.
%
%   fit.lb         - Lower boundaries for the fit
%   fit.ub         - Upper boundaries for the fit
%   specParams (resFactor, Z, limit) - These parameters are forwarded to the function excitationSpec.
%   collParam (colls, gas and P)     - These are the parameters to calculate the collisional broadening and quench-ing rate. These are forwarded to the function excSpecFitCost. 
%
%
% OUTPUT:
%   solutionVals   - The solution result for the optimization problem. This is vector with same varia-ble positions as bestGuess.
%   residue        - The normalized residual for the fit function returned by lsqnonlin
%
% SEE ALSO:
%  excSpecFitCost, lsqnonlin, fitImageDataset, collisions, quenchRate, col-lisionalBroadening, 
%   getFitOptions, fitExcSpectrumImage.mlx, fitExcSpectrumSingle.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    wnum (:,1) double
    measSpec (:,1) double
    linelist (:,:) cell
    MM (1,1) double
    bestGuess (:,:) double
    options (1,1)
    fit.lb (:,:) double = [250  0.001 0.001 0.01 -2 -02.01 -.01]
    fit.ub (:,:) double = [3200 3     3     100     2  02.01  .01]

    specParams.resFactor (1,1) double
    specParams.Z (:,:) double
    specParams.limit (1,1) double

    collParam.colls = {}
    collParam.gas = {}
    collParam.P (1,1) double = 1
end

if fit.lb(1)<=0

   fit.lb(1) = 0.01;
end
if fit.lb(2)<=0

   fit.lb(2) = 0.0000000000001;
end
if fit.lb(3)<=0

   fit.lb(3) = 0.0000000000001;
end
if fit.lb(4)<=0

   fit.lb(4) = 0.0000000000001;
end


p =  namedargs2cell(specParams);


c = namedargs2cell(collParam);
% measSpec = measSpec./max(measSpec);



fitHandle = @(params) excSpecFitCost(...
    wnum,...
    measSpec,...
    linelist,...
    MM,...
    params(1),...%T
    params(2),... %dnuGL
    params(3),... %dnuLL 
    p{:}, c{:},...
    scale = params(4), ...
    dnuSh= params(5), ...
    offset= params(6), ...
    a = params(7) ...... % baseline slope
    );

[solutionVals,residue] =  lsqnonlin(fitHandle, bestGuess , ...
    fit.lb,fit.ub, options);





end