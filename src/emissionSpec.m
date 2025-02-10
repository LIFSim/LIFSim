function [EmSpectra, excitedLines, sumFluor] = ...
    emissionSpec(laserPos, emWnum, linelist, MM, T, emList, em, inst, ex)
% emissionSpec - This function calculates the spectrally resolved fluorescence emissions spectra.
%
% [EmSpectra, excitedLines, sumFluor] = emissionSpec(laserPos, emWnum, linelist, MM, T, emList, em, inst, ex)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function calculates the spectrally resolved fluorescence emissions 
%  based on a laser excitation at a given wavenumber position. Here an overlap
%  function will be computed based on the laser and transition line-shape 
%  functions. The laser function is computed from the FWHMs of the Gaussian 
%  and Lorentzian contributions and then convoluted with the transitions. 
%
% INPUTS:
%   laserPos     - The laser position in wavenumbers. 
%   emWnum       - The emission wavenumber range to be simulated, in 1/cm, 
%                  e.g. emWnum = 33000:1:47000. 
%   linelist     - The list of lines, each line data is stored in a line 
%                  struct, see selectLines and linesStruct.
%   MM           - Molar mass, g/mol
%   T            - Temperature, Kelvin
%   emList       - Cell list containing at each index a list of emissions 
%                  for a given line in linelist. The lists emList and linelist
%                  are of the same order, so emList{i} is the list emission 
%                  lines for the line in linelist{i}, (i is any given index).
%                  See function selectLines.
%
% OPTIONAL INPUTS:
%   dnuL       - The Lorentzian linewidth of the emission transition, resulting 
%                from collisional broadening. This parameter is set by using 
%                collisions and collisionalBroadening. If the latter is not 
%                available and no value for dnuL is set, then a minimal default 
%                value is assumed, 0.0001 1/cm, and thereby the overall Lorentzian
%                feature of the overlap function is controlled mainly through 
%                the laser parameter dnuLL.
%   dnuSh      - The collisional shifting parameter, to be set using collisions
%                and collisionalBroadening. In case no collision data is available,
%                this value is zero and maybe used as a free parameter to 
%                adjust the shift.
%   normalize  - A Boolean flag to normalize the spectra. The default value 
%                is true. This normalization considers the maximum of the 
%                output matrix EmSpectra.
%   Z          - The partition function acquired from selectLines. This may 
%                be omitted if not available, however it is not recommended.
%   limit      - This is a limiting factor used as multiplier to the sum of 
%                the linewidth (em.dnuL, emdnuG, inst.dnuGsm,  inst.dnuLsm). 
%                The result multiplication is considered as a range limit around 
%                a given center frequency. The higher the limit is, the further 
%                lines from the center frequency would be considered for interpolation. 
%                This is set to 100 by default, if the user did not pass the parameter.
%   dnuGsm     - The Gaussian linewidth of the spectrometer line-shape, 1/cm.
%   dnuLsm     - The Lorentzian linewidth of the spectrometer line-shape, 1/cm.
%   resEx      - The resolution for calculating the laser lineshape.
%   dnuGL      - The Gaussian linewidth of the laser line-shape, 1/cm, default = 0.1.
%   dnuLL      - The Lorentzian linewidth of the laser line-shape, 1/cm, default = 0.01.
%   rangeWidth - The width of wavenumber range for generating the excitation overlap. 
%   dnuLex     - Same as dnuL, unless otherwise specified.
%   dnuShex    - Same as dnuSh, unless otherwise specified.
%   exLimit    - The threshold limit of the excitation overlap function when 
%                to consider lines to excite. A higher limit will consider 
%                a narrower region around the peak, hence less exited lines. 
%                This is important to omit regions of the overlap when the 
%                intensity is negligible, which makes the computation faster.
%
% OUTPUT:
%   EmSpectra     - The emission spectra for every excited transition within 
%                   the line-shape function of the laser. 
%   excitedLines  - The indices of the lines found in the region of excitation. 
%                   The indices can be used in the cell list linelist.
%   sumFluor      - The sum of the emission spectra EmSpectra 
%
% SEE ALSO:
%   selectLines, collisions, quenchRate, collisionalBroadening, overlap, 
%   voigtlineMcLean, simulateEmSpectra.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    laserPos (:,:) double
    emWnum (:,:) double
    linelist (:,:) cell
    MM (1,1) double
    T (1,1) double
    emList (:,1) cell
    em.dnuL (1,1) double = 0.0001
    em.dnuSh (1,1) double = 0
    em.normalize logical = true
    em.Z (:,:) double = [T 1;T+1 1]
    em.limit (1,1) double = 100
    inst.dnuGsm (1,1) double = 0.1
    inst.dnuLsm (1,1) double = 0.005
    ex.resEx (1,1) double = 0.01
    ex.dnuGL (1,1) double = 0.1
    ex.dnuLL (1,1) double = 0.01
    ex.exLimit (1,1) double = 0.00001
    ex.dnuLex (1,1) double = 0
    ex.dnuShex (1,1) double = NaN
end
em_aL = 1;
ex_aL = 1;

if ex.dnuLex == 0
    ex.dnuLex = em.dnuL;
end
if isnan(ex.dnuShex)
    ex.dnuShex= em.dnuSh;
end

nu0 = linelist{1}.nu0;

dnuG = dnuGFun(T, MM, nu0);

tripleLinewidth = round((dnuG+ex.dnuGL+ex.dnuLL+ex.dnuLex)*10, 2, "significant");

rangeWidth = tripleLinewidth*2;

exRange = -rangeWidth:ex.resEx:rangeWidth;
gOver = overlap(nu0,ex.resEx, exRange,ex_aL, dnuG, ex.dnuLex, ex.dnuGL, ex.dnuLL);

gOver = gOver./max(gOver(:));
llines = size(linelist,1);

% excited = zeros(llines,1),
% laserRange = laserRange+exWnum;

emRes = emWnum(2)-emWnum(1);
lwnum = length(emWnum);

EmSpectra = double(zeros(llines, lwnum));


excitedLines = cell(size(linelist));

for ln = 1:llines

    line = linelist{ln};

    pos = abs(laserPos-line.nu0-ex.dnuShex)/ex.resEx;

    gamma = interp1(0:length(gOver)-1, gOver, pos, 'linear', 'extrap');

    linelist{ln,2} = 0;
    linelist{ln,3} = gamma;
    if gamma<ex.exLimit
        continue;
    end

    linelist{ln,2} = line.nu0;

    fB = calcFb(line.jLo, T, em.Z, line.EGr);

    const =   fB*gamma*line.B...
        /...
        (line.emSum+line.Q+line.P + line.W);


    lineEmissions = emList{ln};
    for e = 1:length(lineEmissions)
        emLine = lineEmissions{e};

        if emLine.nu0 < emWnum(1) || emLine.nu0 > emWnum(end)
            continue;
        end

        emdnuG = dnuGFun(T, MM, emLine.nu0);
        rangeLimit = round((em.dnuL+emdnuG+inst.dnuGsm +inst.dnuLsm)*em.limit, 2, 'significant');
        emRange  = -rangeLimit:emRes:rangeLimit;

        %spectrometer emission interaction
        % To speed up:
        % This can be assumed the same and move outside the for loop, if
        % nu0 is in the same range, then the change in emdnuG is
        % insignificant. 
        gOverEm = overlap(emLine.nu0,ex.resEx, emRange, em_aL,...
            emdnuG, em.dnuL, inst.dnuGsm, inst.dnuLsm);

        range = 0:length(gOverEm)-1;

        % see fluorTransm, there you can preset a filter for initializing
        % a transmission value per emission line.
        const2 = const*emLine.A*emLine.transm;

        DD = abs(emWnum-emLine.nu0-em.dnuSh);
        xq = DD./emRes;

        within = rangeLimit>DD;
        within = find(within);

        for wn = within
            gammaEm = interp1(range, gOverEm, xq(wn), "linear", "extrap");
            EmSpectra(ln, wn) = EmSpectra(ln, wn) + gammaEm'.*const2;
        end
    end


    excitedLines{ln} = linelist{ln,1};
end

sumFluor = sum(EmSpectra,1);
if em.normalize
    EmSpectra =  EmSpectra./max(EmSpectra(:));
    sumFluor = sumFluor./max(sumFluor(:));
end
excitedLines = find(~cellfun('isempty',excitedLines));
end