function const = CONSTANTS(inp)
% CONSTANTS - Store and load constants in this function
%
%   const = CONSTANTS(inp)
%
% Author:   Abbas El Moussawi 
%
% DESCRIPTION:
%   Store and load constants in this function, such as c, h, k.
%
% INPUTS:
%   inp: Constant name a char, supported: speed of light 'c', Planck's
%   constant 'h',  Boltzmann's constant 'k'
%
% OUTPUT:
%   const: The value of the constant
%   
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
     inp {mustBeMember(inp,["c", "h", "k"])}
end

switch inp
    case 'c'
        %% Speed of light m/s
        const = 299792458;
    case 'h'
        %% Planck's constant J/Hz
        const = 6.62607015e-34; 
    case 'k'
        %% Boltzmann constant J/K
        const = 1.380662e-23;
end
end