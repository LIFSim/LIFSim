function l = LABELS(inp, legendShow)
% LABELS - Runs a preset for an active figure
%
%   l = LABELS(inp, legendShow)
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
%   Runs preset for an active figure, to adjust the labels. More presets
%   can be added. This can be called multiple times for a figure, to label
%   the x and y axis.
%
% INPUTS:
%   inp: "1/cm", "nm", "fluorIntArb", "flourIntNorm", "em_nm", "em_cm", "abs"
%
% OPTIONAL INPUTS:
%   legendShow: Boolean for enabling or disable the legend in a figure.
%   
%
% OUTPUT:
%   l: Label handle.
%   
%
% SEE ALSO:
%   sensitivityAnalysis.mlx, ... 
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen


arguments
    inp {mustBeMember(inp,[  "1/cm", "nm", "fluorIntArb", "flourIntNorm", "em_nm", "em_cm", "abs"])}
    legendShow (1,1) logical = false
end

fN = 'Arial';
fS = 12.5; 
ax = gca;

switch inp

    case '1/cm'
        l = xlabel('Excitation wavenumber / cm^{−1}');
        
        ax.XAxis.Exponent = 0;
        xtickangle(45);
        
    case 'nm'
        l = xlabel('Excitation wavelength / nm');
    case 'flourIntNorm'
        l = ylabel('Fluorescence intensity / normalized');
    case 'fluorIntArb'
        l = ylabel('Fluorescence intensity / Arb.u.');
    case 'abs'
        l = ylabel( 'Absorbance –ln(I/I0)');
    case 'em_nm'
        l = ylabel('Emmission wavelength / nm');
    case 'em_cm'
        l = ylabel('Emmission wavenumber / cm^{−1}  ');
        % l.Interpreter="latex";
        ax = gca;
        ax.YAxis.Exponent = 0;
    otherwise
        return;

end
% l.FontSize = fS;
% l.FontName = fN; 
% set(findall(gca,'-property','FontSize'),'FontSize',fS)
% set(findall(gca,'-property','FontSize'),'FontSize',fN)

fontname(fN)
fontsize(fS,"points")
ax.FontSize=fS-1; 

if legendShow
legend1 = legend(ax, "show" );
set(legend1,'FontSize',fS-2.5);
end

end