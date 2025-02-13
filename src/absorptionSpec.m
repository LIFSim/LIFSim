function spec = absorptionSpec(wnum, linelist, MM, T, dnuGL, dnuLL, params)
% absorptionSpec - Calculates the absorption spectrum using the absorbance
%
%   spec = absorptionSpec(wnum, linelist, MM, T, dnuGL, dnuLL, params)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%  Calculates the absorption spectrum using the absorbance for a given 
%  wavenumber region of the incident light (laser), molecule, temperature, 
%  and instrument characteristics. An overlap function is calculated, based 
%  on two Voigt functions of the laser and the transition, see function 
%  overlap, with a given laser linewidth parameter (Gaussian and Lorentzian
%  linewidth contributions).
%
% INPUTS:
%   wnum     - The wavenumber range to be simulated, in cm-1, e.g. wnum = 
%              44407:0.1:44417;.
%   linelist - A cell containing the list of lines, each line data is stored 
%              in a line struct, see selectLines and linesStruct.
%   MM       - Molar mass g/mol
%   T        - Temperature in Kelvin
%   dnuGL    - The Gaussian FWHM of the laser profile.
%   dnuLL    - The Lorentzian FWHM of the laser profile.
%
% OPTIONAL INPUTS:
%   params - This is a set of optional parameters used as a name-value cell 
%            or struct
%     -	Z         - The partition function acquired from selectLines. This may
%                   be omitted if not available, however it is not recommended.
%     -	dnuL      - The Lorentzian linewidth of the transition, resulting from 
%                   collisional broadening. This parameter is set by using collisions
%                   and collisionalBroadening. If the latter is not available 
%                   and no value for dnuL is set, then a minimal default value
%                   is assumed, 0.0001 cm-1, and thereby the overall Lorentzian
%                   feature of the overlap function is controlled mainly through
%                   the laser parameter dnuLL.
%     -	dnuSh     - The collisional shifting parameter, to be set using collisions
%                   and collisionalBroadening. In case no collision data available, 
%                   this value is zero and maybe used as a free parameter to 
%                   adjust the shift.
%     -	resFactor - A divider factor to the step value of wnum. This is used 
%                   for defining the resolution for the overlap function generation. 
%                   See overlap in section ‎4.6.1. If not given, the default value is 5. 
%     -	normalize - A Boolean flag to normalize the spectra. The default value 
%                   is true. 
%     -	limit     - This is a limiting factor used as multiplier to the sum of 
%                   the linewidth (dnuG, dnuL, dnuGL, dnuLL). The result multiplication
%                   is considered as a range limit around a given center frequency.
%                   The higher the limit is, the further lines from the center 
%                   frequency would be considered for interpolation. This is set
%                   to 12 by default, if the user did not pass the parameter. 
%     -	N         - The number density of the absorbing molecules, m^-3.
%     -	L         - The length of the absorption path, m.
%
% OUTPUT:
%   spec  - The calculated absorption spectrum with the same length as wnum.
%   
%
% SEE ALSO:
%   selectLines, collisions, quenchRate, collisionalBroadening, overlap,
%   voigtlineMcLean, simulateAbsorptionSpectra.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
% 
% CHANGES:
% 2025-02-13: 
%   Change 1: Dividing gOver by total FWHM according to Partridge and
%             Normand, gOver is capital Gamma in https://doi.org/10.1364/AO.34.002645
%   Change 2: Dividing by speed of light, c, to change gOver from cm to s.

arguments
    wnum (1,:) double
    linelist (:,:) cell
    MM (1,1) double
    T (1,1) double
    dnuGL (1,1) double
    dnuLL (1,1) double
    params.dnuL (1,1) double = 0.005
    params.dnuSh (1,1) double = 0
    params.resFactor (1,1) double = 5
    params.normalize logical = false
    params.limit (1,1) double = 12
    params.Z (:,:) double = [T 1;T+1 1]
    params.N (1,1) double = 1
    params.L (1,1) double = 1
end


aL = 1;
 
res = round((wnum(2)-wnum(1))/params.resFactor, 2, "significant");

nu0 = linelist{1}.nu0;
dnuG = dnuGFun(T, MM, nu0); % doppler broadening for transition


rangeLimit = round((params.dnuL+dnuG+dnuGL+dnuLL)*params.limit, 2, "significant");
lwnum = length(wnum);
llines = size(linelist,1);

range  = -rangeLimit:res:rangeLimit;

% use any transition for nu0 >> faster computation,
% the difference is negligible in dnuG if calculated for all transition
% in a limited range
gOver = overlap(nu0, res, range, aL, dnuG, params.dnuL, dnuGL, dnuLL);
lgOver = length(gOver)-1;

%%% Start change 1
halfMax = max(gOver) / 2;
dnu = interp1(gOver, 0:lgOver,halfMax,'linear', 'extrap')*res*2; 
% *2 because gOver represents the overlap function from the peak to the end.
% the value is based on the index of the matrix gOver, thus *res

% Change the dimensionless gamma to cm. Gamma = dnu*g. ref Partridge and
% Normand, https://doi.org/10.1364/AO.34.002645
gOver = gOver./dnu; 
%%% End change 1


spec = double(zeros(1,lwnum));


for ln = 1:llines
    line = linelist{ln};

    DD = abs(wnum-line.nu0-params.dnuSh); % distances of each wavenumber to line(s)
    xq = DD./res; % convert to number within range, for interpolation
    withinDistance = rangeLimit>DD; % calc all conditions
    
    for wn=1:lwnum
        if withinDistance(wn) % if distance larger  7 time linewidth then include
            gamma = interp1(0:lgOver,gOver,xq(wn),'linear', 'extrap');
            f_B = calcFb(line.jLo, T, params.Z, line.EGr);
            overlapWeighed = gamma*line.B*f_B*wnum(wn); 
            spec(wn) = spec(wn) + overlapWeighed;
        end
    end
    
end


%%% Start change 2
% J*s*(1/cm)*(m^3)*(1/J)*(s^-2)*cm*(m*^-3)*m --> m/s
% so devide by c to have unitless absorbance
h = CONSTANTS("h"); 
c = CONSTANTS("c");
hNLg = (h*params.N*params.L./c);

spec= spec.*hNLg;
%%% End change 2

if params.normalize
    spec= spec./max(spec);
end
end