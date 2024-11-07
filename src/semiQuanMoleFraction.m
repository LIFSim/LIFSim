function x = semiQuanMoleFraction(wnum, params, line, Z, MM,  varargin)
% semiQuanMoleFraction - Calculate the mole fraction semi-quantitatively.
%
%   x = semiQuanMoleFraction(wnum, params, line, Z, MM,  varargin)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   The semi quantitative mole fraction of a probed species can be inferred 
%  with this function, for an isolated transition, i to k. The number density 
%  Nk can be substituted with, px/(kT), according to the ideal gas law. the 
%  overlap function is integrated along the excitation wavenumber. A qualitative 
%  mole fraction detection is valuable for a relative analysis, where concentrations
%   at multiple conditions can be investigated and compared.
%
% INPUTS:
%   wnum            - Wavenumber region in cm–1.  
%   solutionValues  - A matrix vector with seven places, which are the resulting 
%                     fit parameters. See output of fitExcitationSpec.
%   line            - A struct, lineStruct, for the line to be used for this calculation. 
%   Z               - Partition function from selectLines.
%   MM              - Molar mass in g/mol.
%
% OUTPUT:
%   x     - The semi quantitative mole fraction. 
%   
%
% SEE ALSO:
%   fitExcitationSpec, lineStruct, selectLines, fluorTransm, excitation-SpecSolution,
%   semiQuantMoleFract.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
%   
%   Based on: Chrystie et al. Applied Physics B (2019) 125:29
%           https://doi.org/10.1007/s00340-019-7137-8
%           Section 2.1

arguments
    wnum (1,:) double
    params (1,:) double
    line (1,1) struct
    Z (:,:) double 
    MM (1,1) double
end
arguments (Repeating)
    varargin
end

linelist = {line};

T= params(1);

% p =  namedargs2cell(specParams);
% 
% 
% c = namedargs2cell(collParam);
% 



% removes the baseline correction.
params(6) = 0;% offset
params(7) = 0;% slope

fitSpec = excitationSpecSolution(...
    wnum, linelist, MM, params,varargin{:}  );

% Scaling is taken care of in excitationSpecSolution


res = (max(wnum)-min(wnum))/length(wnum);
specToGamma = sum(fitSpec)./res;


fb = calcFb(line.jLo, T, Z, line.EGr); % Boltzmann fraction


yield = line.emSumTransm/(line.emSum+line.P+line.Q+line.W); % tao-eff / tao-nat

x = specToGamma*T/(fb*yield);



end