function [wndata, nmData] = loadSpectrum(filepath, unit)
% loadSpectrum - Loads a spectrum curve.
%
%   [wndata, nmData] = loadSpectrum(filepath, unit)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function loads a spectrum curve from a spread sheet file
%   and outputs a table with wavenumbers and intensity, data.wnum and
%   data.intens.
%
% INPUTS:
%   filepath  - System path to your file
%
% OPTIONAL INPUTS:
%   unit      - Set the unit of the file, if not defined then 1/cm is assumed.
%
% OUTPUT:
%   wndata    - The wavenumbers (if input is in nm, it will be converted to 1/cm),
%               and transmission. As wndata.wnum and wndata.intens
%   nmData    - The wavelength (if input is in 1/cm, it will be converted to nm),
%               and transmission. As nmData.wnum and nmData.intens
%
% SEE ALSO:
%   fluorTransm
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    filepath  (1,:) char = ""
    unit  {mustBeMember(unit,["nm", "1/cm"])} = "1/cm"
end

if ~isfile(filepath)
    error("File does not exist");
end



data = readmatrix(filepath);

data = array2table(data, 'VariableNames', {'wnum', 'intens'});


nmData = data;
wndata = data;

if unit == "nm"
    wndata.wnum = convWnumWlen(wndata.wnum);
    wndata.intens = flipud(data.intens);
else
    nmData.wnum = convWnumWlen(data.wnum);
    nmData.intens = flipud(data.intens);
end

end