function baseline = calBaseline (a,b,L)
% calBaseline - This function calculates the linear baseline with the equation y = ax +b for a given length.
%
%   baseline = calBaseline (a,b,L)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function calculates the linear baseline with the equation y = ax +b for a given length.
%
% INPUTS:
%   a     - The slope of the line.
%   b     - The y-intercept or the offset. 
%   L     - The length of the baseline vector needed.
%
% OUTPUT:
%   baseline     - The baseline vector. 
%
% SEE ALSO:
%   excSpecFitCost, fitExcitationSpec, fitExcSpectrumImage.mlx, fitExcSpectrumSingle.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

baseline = (1:L)';
baseline = a*(baseline)+b;
end