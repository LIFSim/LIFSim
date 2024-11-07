function dnuG = dnuGFun(T, MM, nu0) 
% dnuGFun - This is a short hand function to calculate the doppler Gaussian broadening.
%
%   dnuG = dnuGFun(T, MM, nu0) 
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This is a short hand function to calculate the doppler Gaussian broadening
%   of a transition based on temperature, mole mass, and center wavenumber. 
%  
%
% INPUTS:
%   T     - The temperature in Kelvin
%   MM     - The molar mass of the diatomic species in g/mol.
%   nu0     - The center wavenumber.
%
% OUTPUT:
%   dnuG     - The the doppler broadening line-width for the given inputs.
%
% SEE ALSO:
%   excitationSpec, emissionSpec, absorptionSpec
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    T (1,1) double
    MM (1,1) double
    nu0 (1,1) double
end
    dnuG = 7.1623e-7*sqrt(T/MM)*nu0;
end