function spec = excitationSpec(wnum, linelist, MM, T, dnuGL, dnuLL, params) 
% excitationSpec - Calculates the Laser-Induced Fluorescence (LIF) excitation spectrum
%
%   spec = excitationSpec(wnum, linelist, MM, T, dnuGL, dnuLL, params) 
%
% Author:    Abbas El Moussawi and Torsten Endres
%
% DESCRIPTION:
%   Calculates the Laser-Induced Fluorescence (LIF) excitation spectrum for a specified
%  wavenumber range (in 1/cm). This function computes two Voigt functions for the laser 
%  and line, based on the Voigt McLean model (see function "voigtlineMcLean"). It then 
%  calculates the overlap function, which is weighted with the given physical parameters. 
%
% INPUTS:
%   wnum      - The wavenumber range to be simulated, in cm-1, e.g. wnum = 44407:0.1:44417;.
%   linelist  - A cell containing the list of lines, each line data is stored 
%               in a line struct, see selectLines and linesStruct.
%   MM        - Molar mass g/mol
%   T         - Temperature in Kelvin
%   dnuGL     - The Gaussian FWHM of the laser profile.
%   dnuLL     - The Lorentzian FWHM of the laser profile. 
%
% OPTIONAL INPUTS:
%   params    - This is a set of optional parameters used as a name-value 
%               cell or struct
%     -	Z        - The partition function acquired from selectLines. This may 
%                  be omitted if not available, however it is not recommended.
%     -	dnuL     - The Lorentzian linewidth of the transition, resulting from 
%                  collisional broadening. This parameter is set by using 
%                  collisions and collisionalBroadening. If the latter is not 
%                  available and no value for dnuL is set, then a minimal default
%                  value is assumed, 0.0001 cm-1, and thereby the overall 
%                  Lorentzian feature of the overlap function is controlled
%                   mainly through the laser parameter dnuLL.
%     -	dnuSh     - The collisional shifting parameter, to be set using collisions
%                   and collisionalBroadening. In case no collision data available,
%                   this value is zero and maybe used as a free parameter to 
%                   adjust the shift.
%     -	resFactor - A divider factor to the step value of wnum. This is used
%                   for defining the resolution for the overlap function generation.
%                   See overlap. If not given, the default value is 5. 
%     -	normalize - A Boolean flag to normalize the spectra. The default value is true. 
%     -	limit     - This is a limiting factor used as multiplier to the sum 
%                   of the linewidth (dnuG, dnuL, dnuGL, dnuLL). The result 
%                   multiplication is considered as a range limit around a 
%                   given center frequency. The higher the limit is, the 
%                   further lines from the center frequency would be considered 
%                   for interpolation. This is set to 12 by default, if the 
%                   user did not pass the parameter. 
%
% OUTPUT:
%   spec     - The calculated excitation spectrum with the same length as wnum.
%
% SEE ALSO:
%   selectLines, collisions, quenchRate, collisionalBroadening, overlap, 
%   voigtlineMcLean, simulateExcSpectra.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    wnum (:,1) double
    linelist (:,:) cell
    MM (1,1) double
    T (1,1) double
    dnuGL (1,1) double = 0.105
    dnuLL (1,1) double = 0.005 
    params.dnuL (1,1) double = 0.0001
    params.dnuSh (1,1) double = 0  
    params.resFactor (1,1) double = 3
    params.normalize logical = true
    params.Z (:,:) double = [T 1;T+1 1]
    params.limit (1,1) double = 12
end
 
res = round((wnum(2)-wnum(1))/params.resFactor, 2, "significant");

nu0 = linelist{1}.nu0;
dnuG = dnuGFun(T, MM, nu0); % doppler broadening for transition

aL = 1;
rangeLimit = round((params.dnuL+dnuG+dnuGL+dnuLL)*params.limit, 2, "significant");
lwnum = length(wnum);
llines = size(linelist,1);

range  = -rangeLimit:res:rangeLimit;

% use any transition for nu0 >> faster computation,
% the difference is negligible in dnuG if calculated for all transition
% in a limited range

gOver = overlap(nu0, res, range,  aL, dnuG, params.dnuL, dnuGL, dnuLL);

lgOver = length(gOver)-1;
spec = double(zeros(lwnum,1));

for ln = 1:llines
    line = linelist{ln};

    DD = abs(wnum-line.nu0-params.dnuSh); % distances of each wavenumber to line(ss)
    xq = DD./res; % convert to number within range, for interpolation
    withinDistance = rangeLimit>DD; % calc all conditions, if in range limit
    f_B = calcFb(line.jLo, T, params.Z, line.EGr);
    const = f_B*line.B*line.emSumTransm/(line.emSum+line.P+line.Q+line.W);

    for wn=1:lwnum
        if withinDistance(wn) 
            gamma = interp1(0:lgOver,gOver,xq(wn),'linear', 'extrap');
            spec(wn) = spec(wn)+gamma*const;
        end
    end


end


if params.normalize
    spec= spec./max(spec(:));
end




end