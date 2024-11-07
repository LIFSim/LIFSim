function output = convWnumWlen(value, lambda)
% convWnumWlen - Converts wavenumbers to wavelength and vice versa.
%
%   output = convWnumWlen(value, lambda)
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
% This function converts the input value based on the parameter lambda. 
% If lambda is 0, it assumes that the value is a single position and inverts 
% it by 1e7/value. For non-zero lambda, it assumes that value is a delta then
% computes the target delta based on the given lambda and lambda+value, if
% the target delta is greater than 0.009, it is then rounded to two decimal places.
%
% INPUTS:
%   value: A wavelength or wavenumber position or delta.
%   lambda: If value is delta, then lambda must be specified as a reference
%   position.
%   
% OUTPUT:
%   output: The converted value.
%   
%
% SEE ALSO:
%   loadSpectrum, fitExcSpectrumSingle.mlx, ...
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen


arguments 
  value (:,1) double
  lambda (1,1) double = 0 
end

if lambda == 0
    
    value = 1e7./value;
    output = flipud(value);
else

    output = 1e7/lambda - 1e7/(lambda + value);
    output = abs(output);
    if output>0.009
    output = round(output,2);
    end
end



end

