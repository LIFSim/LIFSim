function [specs, noisy] = generateTRange(TRange, wnum, linelist, MM, specParams,collParam, noise)
% lifsimTRange - Calculates and plots lif spectrum for a given wavenumbers range with a
% given temperature range.
%
%   spec = lifsimTRange(wnum, linelist, MM, T, dnuGL, dnuLL, params)
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   Calculates the Laser-Induced Fluorescence (LIF) excitation spectra for a
%   specified wavenumber range (in 1/cm), as a function of a given
%   temperature range. This function uses 'excitationSpec'.
%
% INPUTS:
%   TRange     - Temperature range, e.g. 600:100:1400
%   wnum       - Spectrum in wavenumbers, e.g. 44400:0.025:44425 1/cm
%   MM         - Molar mass for species.
%   linelist   - Line list. See fn: selectLines
%
% OPTIONAL PARAMETERS:
%   options    - doPlot: Plot the spectra (true, false)
%                plotRate: e.g. 2, will plot every other spectrum,
%                however returns all spectra.
%                dnuGL: Laser Guassian linewidth
%                dnuLL: Laser Lorentzian linewidth
%   specParams - Forwarded to excitationSpec (documentation there).
%
% OUTPUT:
%   specs      - The calculated spectra
%   fig        - Figure handle, if doPlot is false returns 0.
%
% SEE ALSO:
%  excitationSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    TRange (:,1) double
    wnum (:,1) double
    linelist (:,:) cell
    MM (1,1) double
    specParams.dnuGL (1,1) double = 0.20
    specParams.dnuLL (1,1) double = 0.01
    specParams.dnuL (1,1) double = 0.01
    specParams.dnuSh (1,1) double = 0 
    specParams.resFactor (1,1) double
    specParams.normalize logical = true 
    specParams.Z (:,:) double
    specParams.limit (1,1) double
    collParam.colls = {}
    collParam.gas = {}
    collParam.P (1,1) double = 1
    noise.snr (1,1) double
    noise.repetitions (1,1) double
end


specs = double(zeros(length(TRange), length(wnum)));

params = rmfield(specParams, 'dnuGL');
params = rmfield(params, 'dnuLL');
params = rmfield(params, 'dnuL');
params = rmfield(params, 'dnuSh');

params =  namedargs2cell(params);




hndl = @(T, LL, dnuL, dnuSh)  excitationSpec(...
    wnum,...
    LL,...
    MM, ...
    T,...%T
    specParams.dnuGL,specParams.dnuLL, params{:}, dnuL=dnuL, dnuSh=dnuSh ...
    );
 

lTRange = length(TRange);
dnuL = specParams.dnuL*ones(lTRange,1);
dnuSh = specParams.dnuSh*ones(lTRange,1);
LL = cell(lTRange,1);

for t = 1:length(TRange)
        LL{t} = linelist;

end



if ~isempty(collParam.colls)

    for t = 1:length(TRange)
        [dnuL(t),  dnuSh(t)] = collisionalBroadening(...
            collParam.gas, collParam.colls, collParam.P, TRange(t));

        [~,LL{t}] = quenchRate(...
            collParam.gas, collParam.colls, TRange(t), collParam.P, MM, linelist);
    end 
end


parfor t = 1:lTRange



    st = hndl(TRange(t),LL{t}, dnuL(t), dnuSh(t));

    specs(t, :) = st(:);

end

if ~isfield(noise,"snr") || ~isfield(noise,"repetitions")
    noisy = {};
    return;
end

noisy = cell(length(TRange), 1);





for t = 1:length(TRange)

    sp.spec = repmat(specs(t, :),noise.repetitions,1);
    sp.spec = awgn(sp.spec, noise.snr, "measured");
    % sp.mean = mean(sp.spec,1);
    sp.T = TRange(t);
    noisy{t} = sp;
end

end
