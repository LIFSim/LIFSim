function fB = calcFb(jLo,T, Z, EGr)
% calcFb - Calculates the Boltzmann’s factor 
%
%   fB = calcFb(jLo,T, Z, EGr)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function calculates the Boltzmann’s equation for the fraction of 
%  molecules in energy level i for a given ground state energy, temperature,
%  partition function, lower rotation quantum number.
%
% INPUTS:
%   jLo  - Lower rotation quantum number
%   T    - Temperature in Kelvin
%   Z    - Partition function for the temperature T
%   EGr  - Ground state energy in cm–1.
%
% OUTPUT:
%   fB   - The calculated Boltzmann factor 
%
% SEE ALSO:
%   excitationSpec, emissionSpec, absorptionSpec, fluorTransm
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
arguments
    jLo (1,1) double
    T (1,1) double
    Z (:,:) double
    EGr (1,1) double
end

Z=interp1(Z(:,1),Z(:,2),T, 'linear', 'extrap');

fB = (2*jLo+1)*exp(-1.4387863*EGr/(T))/Z;

end

