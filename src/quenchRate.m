function [quen,linelist] = quenchRate(gas, colls, T, P, MM, linelist)
% quenchRate - Computes the total quenching rate, Q, for a given molecule.
%
%   [quen,linelist] = quenchRate(gas, colls, T, P, MM, linelist)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   Computes the total quenching rate, Q, for a given molecule with respect to 
%  a gas composition of collision partners. The quenching rate is calculated based 
%  on the pressure in bar, Temperature in Kelvin, the gas composition of different 
%  elementsas fractions, the molar masses of the collision partners to calculate the 
%  reduced mass in the relative velocity. The Avogadro’s constant, bar to Pascal (1e5), 
%  and Angstroms^2 to m^2 (1e–20), are used to have Q in 1/s. 
%
% INPUTS:
%   gas       - A table containing two columns (molecule, fraction) listing the 
%               composition. See loadGasComposition
%   colls     - The collision model based on the function collisions.
%   T         - Temperature in Kelvin.
%   P         - Pressure in bar.
%   MM        - The molar mass of the diatomic species in g/mol.
%
% OPTIONAL INPUTS:
%   linelist  - A loaded linelist can be optionally passed to the function. 
%               This will update the field Q in each lineStruct in the list.
%
% OUTPUT:
%   quen      - The quenching rate in 1/s.
%   linelist  - The updated input linelist, with the field Q from the 
%               lineStruct updated.
%
% SEE ALSO:
%   loadGasComposition, collisions, collisionalBroadening, 
%   excSpecFitError, fitExcitationSpec, excSpecFitCost, emissionSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    gas (:,2) table
    colls (1,1) struct
    T (1,1) double
    P (1,1) double
    MM (1,1) double
    linelist (:,1) = {}

end

% qrate +=
% GetMolefraction(qparm->qname[j],bparm)
% * sqrt(1 + GetMolemass(bparm->lifmolename)
% / GetMolemass(qparm->qname[j])) * cross;
quen = 0;

for i=1:length(gas.molecule)
    fract = gas.fraction(i);
    gasMolec = gas.molecule{i};

    try
        cross = colls.(gasMolec).cross(T);
        %%% reduced mass from Molar mass
        % 1/u_ab =  Na * (1/MMa + 1/MMb), Na is Avagadros constant
        % 
        quen = quen + cross * fract * sqrt(1/MM + 1/ molMass(gasMolec));
    catch
        continue;
    end



end

k = CONSTANTS("k");
conv = 0.0245400504481959; %
% conv
% changes bar*angstrom^2*sqrt(mol/g)...
% to
% 1e5 Pa * 1e-20 m^2 * sqrt(6.02214076e23 1/g) = 0.0245400504481959
% 
quen = quen * (P/(k*T)) * sqrt(8*k*T/pi) * conv;

if isempty(linelist)
    return;
end

for i = 1:height(linelist)

    linelist{i}.Q = quen;

end

end