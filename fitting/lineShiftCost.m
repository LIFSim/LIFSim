function cost = lineShiftCost(...
    wnum, measSpec, linelist, lineIndex, pos, fitParams, MM, varargin)
% lineShiftCost - This is necessary for the optimization function linesShiftCorrection.
%
%   cost = lineShiftCost(...
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function uses the cost function excSpecFitCost for calculating the
%   difference between a measured and a simulated spectrum with given parameters. 
%  However, it uses two extra arguments to adjust a single line position. These 
%  arguments are the index of the line in the array linelist and the new position 
%  in wavenumbers. This is necessary for the optimization function linesShiftCorrection.
%
% INPUTS:
%   wnum     - Wavenumber region in cm–1.
%   measSpec     - The measured spectrum data point, with the same length as wnum.
%   linelist     - The linelist obtained from function selectLines.
%   bestGuess     - A matrix vector with seven places, which are the initial fit parameters to start the optimization. 
%   lineIndex     - A matrix vector with the index of each line to be shifted in the linelist cell vec-tor. The positions can be adjusted by specifying the new position through pos. The length of pos should be identical to lineIndex.
%   pos     - The adjusted line position wavenumber in cm–1.
%   params     - A matrix vector with seven places, which are the initial fit parameters to start the optimi-zation. 
%   MM     - Molar mass of the species, obtained from function molMass or selectLines.
%   
%
% OPTIONAL INPUTS:
%   varargin     - Further variable arguments to forward to function excSpecFitCost.
%
% OUTPUT:
%   cost     - The difference between the measured and simulated spectra, with respect to the given parameters.
%
% SEE ALSO:
%   linesShiftCorrection, excSpecFitCost, lineStruct, selectLines
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    wnum
    measSpec
    linelist
    lineIndex
    pos
    fitParams
    MM

end
arguments (Repeating)
    varargin
end


for s = 1:length(lineIndex)
    linelist{lineIndex(s)}.nu0 = pos(s);
end

cost = excSpecFitCost(...
    wnum,...
    measSpec,...
    linelist,...
    MM,...
    fitParams(1),...%T
    fitParams(2),... %dnuGL
    fitParams(3),... %dnuLL
    varargin{:},...
    scale= fitParams(4), ...
    dnuSh= fitParams(5)+pos(end), ...
    offset= fitParams(6), ...
    a=fitParams(7) ... % baseline slope
    );
end

