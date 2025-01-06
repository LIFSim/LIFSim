function gas = loadGasComposition(path)
% loadGasComposition - Loads a gas composition from a csv file.
%
%   gas = loadGasComposition(path)
%
% Author:    Abbas El Moussawi
% 
% DESCRIPTION:
% This function considers mole fractions of a given gas composition. The given
% composition should be complemented by a model in collisions.m. The output
% will be a table with columns 'molecule' and 'fraction'. The elements with
% 0 will not be imported. The fractions will be normalized to sum 1.
%
% INPUTS:
%   path     -  System path to your spreadsheet, e.g. "D:/path/to/file.csv", if
%               this parameter is not passed, a prompt window will open to
%               browse to file location and select it. E.g. file: 
%                   molecule,fraction
%                   CH4,0.545
%                   O2,0.455
%                   Ar,1.090
%                   NO,0.010
%               
% OUTPUT:
%   gas      -  A table with two columns, molecule and fraction
%               corresponding to your spreadsheet.
%
% SEE ALSO:
%  collisions, quenchRate, collisionalBroadening
% 
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    path (1,:) char = ''
end

if isempty(path)
    [f,p] = uigetfile('./input-data/*.csv',"Gas composition file");
    path = [p f];
end
gas = readtable(path);
gas.fraction=gas.fraction./sum(gas.fraction);
rowsToRemove = gas.fraction == 0;
gas(rowsToRemove, :) = [];


end