function options = getFitOptions(opts)
% getFitOptions - This function stores preset for fitting parameters for the function lsqnonlin.
%
%   options = getFitOptions(opts)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   This function stores preset for fitting parameters for the function lsqnonlin 
%  that is used for fitting the spectra in LIFSim. The user may add here more presets as 
%  needed. 
%
% INPUTS:
%   opts     - The name of the options preset to be returned.
%   
%
% OUTPUT:
%   options     - The options object for fitting functions. 
%   
%
% SEE ALSO:
%   fitExcitationSpec, lsqnonlin, fitImageDataset, fitExcSpectrumImage.mlx, fitExcSpectrumSingle.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

arguments
    opts.id {mustBeMember(opts.id,[ ...
        "default", ...
        "coarse", ...
        "shift", "lsqnonlin", "lineShift"])} = 'default'
end

 
switch opts.id
    case 'default'
        options = optimoptions('lsqnonlin');
        options.Algorithm="trust-region-reflective";
        options.Display='final-detailed';
        options.TolFun = 1e-5;
        options.TolX = 1e-4;
        options.MaxIter = 2200;
        options.StepTolerance = 1e-5;
    case 'lsqnonlin'
        options = optimoptions('lsqnonlin');
        options.Algorithm="trust-region-reflective";
        options.Display='final-detailed';
        options.TolFun = 1e-5;
        options.TolX = 1e-4;
        options.MaxIter = 2200;
        % options.FiniteDifferenceStepSize= 1e-9;
        % options.FiniteDifferenceType="forward";
        options.StepTolerance = 1e-6;
        % options.UseParallel = true;
        return;
    case 'shift'
        options = optimoptions('lsqnonlin');
        options.Algorithm="trust-region-reflective";
        options.Display='final-detailed';
        options.TolFun = 1e-1;
        options.TolX = 1e-1;
        options.MaxIter = 2000;
        options.FiniteDifferenceStepSize= 1e-9;
        options.StepTolerance = 1e-5;
        return;
    case 'lineShift'
        options = optimoptions('lsqnonlin');
        options.Algorithm="trust-region-reflective";
        options.Display='final-detailed';
        options.TolFun = 1e-5;
        options.TolX = 1e-4;
        options.MaxIter = 2200;
        options.StepTolerance = 1e-5;
        return;
    case 'coarse'
        options = optimoptions('lsqnonlin');
        options.Algorithm="trust-region-reflective";
        options.Display='final-detailed';
        options.TolFun = 1e-1;
        options.TolX = 1e-2;
        options.MaxIter = 1000;
        options.StepTolerance = 1e-2;

end



end