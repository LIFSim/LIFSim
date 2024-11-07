function [dnuL, dnuSh, missing] = collisionalBroadening(gas, colls, P, T)
% collisionalBroadening - 
% Returns the collisional broadening, 
% dnuL, and shifting parameters, dnuSh
%
%   [dnuL, dnuSh, missing] = collisionalBroadening(gas, colls, P, T)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%  The function calculates the collisional broadening, dnuL, and shifting parameters, 
%  dnuSh, for a given gas mixture at specified pressure and temperature conditions. It 
%  calculates the total collisional broadening and shifting resulting from the contributing 
%  collision partners. These partners and their corresponding fraction are provided to 
%  the function as input, along with temperature and pressure. If a partner is added in 
%  the gas composition, where this molecule has no model provided, this will be notified 
%  within a third output.
%
% INPUTS:
%   gas    - A table containing two columns (molecule, fraction) listing the 
%            composition. See function loadGasComposition 
%   colls  - The collision model based on the function collisions.
%   P      - Pressure in bar.
%   T      - Temperature in Kelvin.
%
% OUTPUT:
%   C         - The collisional model as a structure
%   dnuL      - The total collisional broadening from the calculated model, 
%               see section ‎4.5.1. 
%   dnuSh     - The total collisional shift from the calculated model, see 
%               section ‎4.5.1.
%   missing   - A cell list of missing collisional model for a given molecule/atom 
%               in gas was not present. Otherwise, this will be empty.
%
% SEE ALSO:
%   loadGasComposition, collisions, quenchRate, excSpecFitCost, fitExcitationSpec,
%   excitationSpec, emissionSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
 
dnuSh = 0;
dnuL = 0;
missing = {};

for i=1:length(gas.molecule)
    fract = gas.fraction(i);
    
    try
        cc = colls.(gas.molecule{i});
        dnuL  = dnuL+fract*cc.br(T);
        dnuSh = dnuSh+fract*cc.sh(T);
    catch
        missing{end+1} = gas.molecule{i};
        continue; 
    end

    
 
end

dnuL  = dnuL.*P;
dnuSh = dnuSh.*P;  
end