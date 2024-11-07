function MM = molMass(formula)
% molMass - Calculates an approximate molar mass of a given formula in simple form.
%
%   MM = molMass(formula)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
%   Calculates an approximate molar mass of a given formula in simple form.
%   This function supports no brackets support. e.g. inputs   'NO', 'SiO', 'OH','O2'.
%
% INPUTS:
%   formula  - Simple molecular formula or structural; case sensitive and no brackets supported.
%
% OUTPUT:
%   MM        - The calculated molar mass in g/mol.
%
% SEE ALSO:
%   selectLines, excitationSpec, emissionSpec, absorptionSpec, dnuGFun, quenchRate
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen

elements = {
    'H';'He';'Li';'Be';'B';'C';'N';'O';'F';'Ne';'Na';'Mg';'Al';'Si';'P';
    'S';'Cl';'Ar';'K';'Ca';'Sc';'Ti';'V';'Cr';'Mn';'Fe';'Co';'Ni';'Cu';
    'Zn';'Ga';'Ge';'As';'Se';'Br';'Kr';'Rb';'Sr';'Y';'Zr';'Nb';'Mo';'Tc';
    'Ru';'Rh';'Pd';'Ag';'Cd';'In';'Sn';'Sb';'Te';'I';'Xe';'Cs';'Ba';'La';
    'Ce';'Pr';'Nd';'Pm';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';
    'Hf';'Ta';'W';'Re';'Os';'Ir';'Pt';'Au';'Hg';'Tl';'Pb';'Bi';'Po';'At';
    'Rn';'Fr';'Ra';'Ac';'Th';'Pa';'U';'Np';'Pu';'Am';'Cm';'Bk';'Cf';'Es';
    'Fm';'Md';'No';'Lr'
    };
masses = [
    1.0079,4.0026,6.941,9.0121,10.811,12.0107,14.0067,15.9994,...
    18.9984,20.1797,22.9897,24.305,26.9815, 28.0855,30.9737,...
    32.065,35.453,39.948,39.0983,40.078,44.9559,47.867,50.9415,51.9961,...
    54.9380,55.845,58.9332,58.6934,63.546,65.39,69.723,72.64,74.9216,...
    78.96,79.904,83.8,85.4678,87.62,88.9058,91.224,92.9063,95.94,98,...
    101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,121.76,127.6,...
    126.9044,131.293,132.9054,137.327,138.9055,140.116,140.9076,...
    144.24,145,150.36,151.964,157.25,158.9253,162.5,164.9303,167.259,...
    168.9342,173.04,174.967,178.49,180.9479,183.84,186.207,190.23,...
    192.217,195.078,196.9665,200.59,204.3833,207.2,208.9803,209,210,222,...
    223,226,227,232.0381,231.0358,238.0289,237,244,243,247,247,251,252,...
    257,258,259,262 ...
    ];

if isstring(formula)
    formula = formula.char;
end

splitted = {formula(1)};

eln = 1;
for i = 2:length(formula)
    if isstrprop(formula(i),'upper')
        el = formula(i);
        eln = eln +1;
    else
        el = strcat(splitted{eln}, formula(i));
    end

    splitted{eln} = el;

end
MM = 0;
for i = 1:length(splitted)
    el = splitted{i};
    nm = regexp(el,'\d*','Match');
    
    if isempty(nm) 
        nm = 1;
    else
        nm = nm{1};
        el = split(el, nm);
        el = el{1};
        nm = str2double(nm);
    end

    index = find(strcmp(el, elements));
    if isempty(index)
        error(['Invalid element in formula: ' el]);
    end

    MM = MM + nm*masses(index);
end

end