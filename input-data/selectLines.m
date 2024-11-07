function [linelist, Z, n, emRange, MM, emList, listCell] = selectLines(species, wnrange, o)
% selectLines - Loads the lines database and the partition function for a
% given species.
%
% [linelist, Z, n, emRange, MM, emList, listCell] = selectLines(species, wnrange, o)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   Loads the lines database and the partition function for a given species
%   (NO, SiO, OH, O2). If further species are needed, the lines database
%   should be added in the given format. Then this script should be
%   updated. e.g. Add the following under the switch statement after adding
%   the .mat files.
%
% INPUTS:
%   species        - A char verctor, e.g species = 'NO'. 
%
% OPTIONAL PARAMETERS:
%   wnrange        - Wavenumber range in 1/cm, leave empty to get all.
%   outside        - Get more lines in an extended region around wnrange.
%   
%
% OUTPUT:
%   linelist       - Cell array of contain lines data. Each cell has the
%                    line info as struct, see lineStruct.
%
%   Z              - Partition function as matrix (T, Z).
%   n              - Number of lines found in the region
%   emRange        - The region of emission for the given wavenumber range.
%   MM             - Species molar mass
%   emList         - Cell list containing at each index a list of emissions 
%                    for a given line in linelist. The lists emList and linelist 
%                    are of the same order, so emList{i} is the list emission 
%                    lines for the line in linelist{i}, (i is any given index).
%   listCell       - This returns the lines in a cell matrix as load from the
%                    selected MATLAB workspace. The indices are appended to 
%                    the last column, with a consistent order as in linelist.
%
% SEE ALSO:
%  lineStruct, fluorTransm, importLines.mlx, importPartitionFunc.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    species
    wnrange (:,1) double = 0
    o.outside (1,1) uint32 = 0
    o.listVersion (1,:) char = ''
    o.partVersion (1,:) char = ''
end


outside = o.outside;
Z = 0;
linelist = cell(0,1);
emList= cell(0,1);


if isstring(species)
    species = species.char;
end


if ~isempty(o.listVersion)
    o.listVersion = ['_' o.listVersion];
end
if ~isempty(o.partVersion)
    o.partVersion = ['_' o.partVersion];
end

load([species '_lines' o.listVersion '.mat']);
load([species '_partfunc' o.partVersion '.mat']);  % exomol
listCell = eval([species '_lines']);
Z =  eval([species '_partfunc']);

MM = molMass(species);

minEm = Inf;
maxEm = 0;
emRange = [minEm maxEm];


t = cell2table(listCell);
t = sortrows(t,4);
listCell = table2cell(t);

linesReadInAllMat = cell2mat(listCell(:,2:10));
linesReadInAllCell = listCell;
n = size(listCell,1);
m = (1:n)';
m = mat2cell(m, ones(n, 1));

listCell = [listCell m];


if wnrange ~= 0
    outside = double(outside);
    if length(wnrange)==1 && outside ==0
        outside = 1;
    end
    
    wnmin=min(wnrange)-outside;
    wnmax=max(wnrange)+outside;

    linesSelect=cell2mat(listCell(:,4));
    linesSelect(linesSelect<wnmin)=0;
    linesSelect(linesSelect>wnmax)=0;
    linesSelect=find(linesSelect);

    listCell=listCell(linesSelect,:);

    n = size(listCell,1);
    if n == 0
        return;
    end
end



linelist = cell(n,1);
emList = cell(n,1);



for i = 1:n

    line = lineStruct(listCell(i,:));

    % Find rows with with same upper bands (v', j')
    hits = find(linesReadInAllMat(:, 2) == line.jUp & linesReadInAllMat(:, 9) == line.vUp);
    AList = cell(length(hits),1);
    emSum = 0;

    for l = 1:length(hits)
        h = hits(l);
        emissionLine = lineStruct(linesReadInAllCell(h,:));
        emissionLine.transm = 1;
        AList{l} =emissionLine;
        emSum = emSum + emissionLine.A;

        minEm = min(minEm, emissionLine.nu0);
        maxEm = max(maxEm, emissionLine.nu0);
    end
    line.emSum = emSum;
    line.emSumTransm = emSum;

    emList(i) = {AList};
    linelist(i) = {line};
end
emRange = [minEm maxEm];

end

