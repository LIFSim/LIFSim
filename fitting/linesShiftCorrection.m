function [linelist] = linesShiftCorrection(...
    linelist, wnum, measSpec, fitParams, MM, varargin, opts)
% linesShiftCorrection - This function optimizes the individual line positions. 
%
%   [linelist] = linesShiftCorrection(linelist, wnum, measSpec, fitParams, MM, varargin, opts)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   The lines in the array linelist will be looped over and each corresponding position 
%   will be ad-justed based on a fit with the line position as free parameter. If no 
%   indices are specified, all lines in the array linelist will be shifted if needed. 
%   Furthermore, a boundary can be set as well to limit the shifting of lines. The 
%   variable of nu0 will be adjusted for each line in the array linelist and this array 
%   will be returned.
%
% INPUTS:
%   linelist    - The linelist obtained from function selectLines.
%   wnum        - Wavenumber region in cm–1.  
%   measSpec    - The measurement data points, which should have the same length as wnum.
%   fitParams   - A matrix vector with seven places, which are the initial fit parameters
%                to start the optimization. 
%   MM          - Molar mass in g/mol.
%
% OPTIONAL INPUTS:
%   fitRegions   - A matrix of size (n,2) defining the regions to consider for position fitting.
%   Each row of this matrix defines a region to search for lines in it.
%   E.g. A matrix, fitRegions = [44415.4 44416.8; 44410.3 44412.3], will allow the function to 
%   adjust the positions of lines within 44415.4–44416.8 cm–1 and 44410.3–44412.3 cm–1. If not 
%   specified, all lines will be con-sidered.
%   boundary     - The boundary for shifting the lines. This will be considered as upper and 
%   lower bound-ary around the line position nu0 (default = 1 1/cm).
%   fitOpts      - The options objects needed by the optimization function lsqnonlin. Some presets
%   can be obtained from the function getFitOptions.
%   verbose      - A logical flag, if true, the progress and a final report for the fit details 
%   will be displayed. (default = false).
%   varargin     - Further variable arguments to forward to function excSpecFitCost. Use named 
%   argument cell.
%
% OUTPUT:
%   linelist     - The given array linelist with the adjusted variable nu0 for each line.
%   
%
% SEE ALSO:
%   lineShiftCost, excSpecFitCost, getFitOptions, lineStruct, selectLines
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    linelist
    wnum
    measSpec
    fitParams
    MM
    
end
arguments(Repeating)
    varargin
end
arguments
    opts.verbose (1,1) logical = false
    opts.fitRegions (:,2) double = [] %all
    opts.fitOpts = {}
    opts.boundary (1,1) double = 1
end

if isempty(opts.fitOpts)
    opts.fitOpts = getFitOptions(id="lineShift");
end
if isempty(opts.fitRegions)
    opts.fitRegions = [-Inf Inf] ;
end

% lines2Shift = 1:length(linelist);
lines2Shift = [];
for i = 1:length(linelist)
    nu0 = linelist{i}.nu0;
    for j = 1:height(opts.fitRegions)
        if nu0 < opts.fitRegions(j,2) && nu0 > opts.fitRegions(j,1)

            lines2Shift(end+1) = i;
            break;
        end
    end
end



nu0s = zeros(length(lines2Shift),1);
lbs = nu0s;
ubs = nu0s;

for s = 1:length(lines2Shift)
    nu0s(s) = linelist{lines2Shift(s)}.nu0;
    lbs(s) = nu0s(s) - opts.boundary;
    ubs(s) = nu0s(s) + opts.boundary;
end
nu0s(end+1) = 0;
lbs(end+1) = nu0s(end) - opts.boundary;
ubs(end+1) = nu0s(end) + opts.boundary;

objLinePos = @(pos) lineShiftCost(...
    wnum,...
    measSpec,...
    linelist,...
    lines2Shift,...
    pos,...% fitting this
    fitParams,...
    MM,...
    varargin{:} ...
    );

[FitValuesLinFit, res] = lsqnonlin(objLinePos, nu0s, lbs,ubs, opts.fitOpts);

for s = 1:length(lines2Shift)
    i = lines2Shift(s);
    if opts.verbose
        disp([num2str(i) '. '])
        disp(['Shifted ' linelist{i}.branch ...
            ', delta: ' num2str(linelist{i}.nu0) ...
            ' - ' num2str(FitValuesLinFit(s)) ...
            ' = ' num2str(linelist{i}.nu0-FitValuesLinFit(s))]);
    end
    linelist{i}.nu0=FitValuesLinFit(s);

end
if opts.verbose
    disp(['Residue: ' num2str(res)]);
end
end

