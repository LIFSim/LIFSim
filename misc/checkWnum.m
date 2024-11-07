function checkWnum(unit,value)
% checkWnum - Warns if a value and unit do not match.
%
%   checkWnum(unit,value)
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
%  The function checkWnum warns if a given value in a specified unit (nm 
% or 1/cm) is matching the unit. If the value exceeds 1000 and is in nm
% then the function warns and suggests a converted value in 1/cm. And vice
% versa is true for 1/cm.
%
% INPUTS:
%   unit:'nm' or'1/cm'
%   value: The value to be checked.
%   
% OUTPUT:
%   None
%
% SEE ALSO:
%   fitExcSpectrumImage.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen


suggestion = num2str(convWnumWlen(value));
if strcmpi(unit,'nm') && value > 1000
    warning(['Unit selected in nm and input value ' num2str(value) ', suggestion: ' suggestion]);
elseif strcmpi(unit,'1/cm') && value < 1000
    warning(['Unit selected in 1/cm and input value ' num2str(value) ', suggestion: ' suggestion]);
end

end

