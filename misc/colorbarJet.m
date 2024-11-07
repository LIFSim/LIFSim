function c = colorbarJet(figureHandle, opts)
% colorbarJet - Converts a figure wiht scale image to jet colors with bar option on.
%
%   c = colorbarJet(figureHandle, opts)
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
%   Converts a figure wiht scale image to jet colors with bar option on.
%
% INPUTS:
%   figureHandle:
%   
%
% OPTIONAL INPUTS:
%   opts.limits: Limits for bar limits
%   opts.label: Label for the side bar
%
% OUTPUT:
%   c: colorbar handle.
%   
%
% SEE ALSO:
%   fitExcImageDataset, fitExcSpectrumImage.mlx, semiQuantMoleFractImage.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
arguments
    figureHandle (1,1) 
    opts.limits (1,2) = [0 0]
    opts.label (1,:) char = ''
end

if opts.limits(1) ~= opts.limits(2)
clim(opts.limits);
end

ax = gca(figureHandle);
c = colorbar(ax);
c.Label.String = opts.label;

jt=jet(255);
jt(1,:) = [0 0 0];
jt(end,:)=[1 1 1];

set(ax,'Colormap',jt)


end
