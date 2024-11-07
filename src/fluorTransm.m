function [linelist,emList] = fluorTransm(linelist, emList, filter)
% fluorTransm - This function interpolates the transmission at the emission
% wavenumbers, by a given  filter transmission curve.
%
%   linelist = fluorTransm(linelist, emList, filter)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function interpolates the transmission at the emission wavenumbers,
% by a given filter transmission curve. Furthermore, it calculates the sum
% of Einstein A coefficients with the considered transmission.
%
% INPUTS:
%   linelist    - The list of lines as return by function selectLines
%   emList      - The emission lists of the given transitions, as returned
%                 by function selectLines
%   filter      - The transmission curve of the used filter, as a table with
%                 wnum and intens as headers. See function loadSpectrum.
%                 The intens column is the transmission of the filter
%
% OUTPUT:
%   linelist    - The updated input linelist, with the field emSumTrapped from
%                 the lineStruct updated.
% SEE ALSO:
%   lineStruct, loadSpectrum, excitationSpec, emissionSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    linelist (:,:) cell
    emList (:,:) cell
    filter (:,2) table

end



for i = 1:length(linelist)
    line = linelist{i};

    emSumFiltered = 0;
    emLineList = emList{i};
    for j = 1:length(emLineList)

        % filter
        fTransm = 1;
        if length(filter.wnum)>1
            fTransm = ...
                interp1(filter.wnum, filter.intens, emLineList{j}.nu0,...
                'linear', 'extrap');
        end

        emLineList{j}.transm =fTransm; % for function emission Spectra, if needed

        emSumFiltered = emSumFiltered + emLineList{j}.A*fTransm;
    end
    line.emSumTransm = emSumFiltered;
    emList{i} = emLineList;
    linelist{i} = line;

end
end
