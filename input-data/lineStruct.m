function line = lineStruct(lineRaw)
% lineStruct - Converts raw cell date of line to the a structure.
%
%   line = lineStruct(lineRaw)
%
% Author:    Abbas El Moussawi  
% 
% DESCRIPTION:
%   Converts raw cell date of line to the a structure, from a cell vector:
%   Branch, j_low, j_upp, E(trans), E(ground), Ein A, Ein B, Pred, v_low, v_up 
%
% INPUTS:
%   lineRaw     - cell vector:
% Branch, j_low, j_upp, E(trans), E(ground), Ein A, Ein B, Pred, v_low, v_up 
%
% OUTPUT:
%   line        - The line info as struct:
%        line.A    - Einstein A emission coefficient (1/s)
%        line.B    - Einstein B absorption coefficient (m^3/(Js^2))
%        line.EGr  - Ground state energy (1/cm)
%        line.nu0  - Line position / Transition energy (1/cm)
%        line.jLo  - Lower rotational quantum number 
%        line.jUp  - Upper rotational quantum number 
%        line.vLo  - Lower vibrational quantum number 
%        line.vUp  - Upper vibrational quantum number 
%        
%        To include these edit create a new version of this function or a
%        postprocessing function to handle this after selectLines return
%        line.P    - Predissociation rate (0)
%        line.W = 0;
%        line.Q = 0;
%        line.emSum = 1;
%        line.emSumTransm = 1;
%        line.transm = 1;
%
% SEE ALSO:
%  selectLines
% 
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen  

    line.A = lineRaw{6};
    line.B = lineRaw{7};
    line.EGr = lineRaw{5};
    line.nu0 = lineRaw{4};
    line.jLo = lineRaw{2};
    line.jUp = lineRaw{3};
    line.vLo = lineRaw{9};
    line.vUp = lineRaw{10};

    %vvvv modify for inclusion vvvv
    line.P = 0; 
    line.W = 0;
    %^^^^ modify for inclusion ^^^^

    line.Q = 0;

    line.emSum = 1;

    line.emSumTransm = 1;
    line.transm = 1;
    % Here assumed 1, however can be post this function estimated using
    % functions: fluorTransm.

    if ~ischar(lineRaw{1})
        line.branch = char(lineRaw{1});
    else
        line.branch = lineRaw{1};
    end

end