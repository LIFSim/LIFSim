function availableSpecies = getAvailableSpecies
% getAvailableSpecies - Returns the available species
%
%   availableSpecies = getAvailableSpecies;
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
%   Returns the available species
%
% INPUTS:
%   None.
%
% OUTPUT:
%   availableSpecies: Available species under input-data\lines\*_lines*.mat
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen


sp = dir("input-data\lines\*_lines*.mat");
availableSpecies = cell(length(sp),1);
for i = 1: length(sp)
    s = split(sp(i).name,'_lines');
    availableSpecies{i} = s{1};
end
availableSpecies = unique(availableSpecies);
availableSpecies = string(availableSpecies);


end