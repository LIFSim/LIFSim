function [spec, collEffect] = excitationSpecSolution(...
    wnum, linelist, MM, solutionParams, varargin, collParam)
% excitationSpecSolution - Computes the excitation spectrum based on a
% given fit solution
%
%   [spec, collEffect] = excitationSpecSolution(wnum, linelist, MM, solutionParams, varargin, collParam)
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
%   This function computes the excitation spectrum based on a given
%   solution parameters, typically the result of a fit. This is to save
%   time writing reduntant code to do so using the simulation function
%   excitationSpec. The collision parameters can be passed here
%
% INPUTS:
%   wnum           - The wavenumber region in 1/cm
%   linelist       - he linelist obtained from function selectLines.
%   MM             - Molar mass of the species, obtained from function molMass or selectLines.
%   solutionParams - A matrix vector with seven places, which are the solution
%                   parameters to start the optimization. 
%
% OPTIONAL INPUTS:
%   varargin            - This is forwarded to excitationSpec
%   collParam (colls, gas and P) - These are the parameters to calculate the 
% collisional broadening and quench-ing rate. These are forwarded to the 
% function excSpecFitCost. 
%
% OUTPUT:
%   spec       - The calculated spectrum based on the solution values
%   collEffect: Struct
%       - dnuShC: Calculated colisional shift.
%       - dnuL: Calculed colisional Lorentian linewidth.
%       - linelist: The line list with the collision rate in the Q variable
%       of the struct in each cell.
%
% SEE ALSO:
%   semiQuanMoleFraction, fitExcSpectrumImage.mlx, fitExcSpectrumSingle.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments (Input)
    wnum (:,:) double
    linelist (:,:) cell
    MM (1,1) double
    solutionParams (:,:) double
end

arguments(Input, Repeating)
    varargin
end

arguments (Input)
    collParam.colls = {}
    collParam.gas = {}
    collParam.P (1,1) double = 1
end

collEffect = {};

if ~isempty(collParam.colls)

    [dnuL, dnuShC] = collisionalBroadening(...
        collParam.gas, collParam.colls, collParam.P, solutionParams(1));
    [~,linelist] = quenchRate(...
        collParam.gas, collParam.colls, solutionParams(1), collParam.P, MM, linelist);

    collEffect.dnuShC = dnuShC;
    collEffect.dnuL = dnuL;
    collEffect.linelist = linelist;

    solutionParams(5) = solutionParams(5) + dnuShC;
    varargin = [varargin(:)', {'dnuL'}, {dnuL}];
end


spec = excitationSpec(wnum,...
    linelist,...
    MM,...
    solutionParams(1),...%T
    solutionParams(2),... %dnuGL
    solutionParams(3),... %dnuLL
    varargin{:},... % to do test colls
    dnuSh=solutionParams(5) ...
    );


baseline = calBaseline(solutionParams(7),solutionParams(6),length(spec));

spec = baseline+spec.*solutionParams(4) ;

end