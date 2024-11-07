function C = collisions(species)
% collisions - Returns the collisional model as a structure.
%
%   C = collisions(species)
%
% Author:    Jan Menser and Abbas El Moussawi
%
% DESCRIPTION:
%   In this function the models for collisional broadening, shift, and quenching are 
%  implemented as structure. This version supports the NO species based on the Harpoon
%  model, by Paul et al. {Paul, 1995, DOI: 10.2172/10111713}. 
%  Broadening:
%  The supported partners for collisional broadening and shift are: 
%      N2, O2, H2O, Ar, CO2, CO, CH4
%  Quenching:
%  The supported collision partners are: 
%      N2, O2, CO2, CO, H2O, CH4, C2H6, C3H8, C2H4, C2H2, NO,
%      NO2, N2O, NH3, NH, H2, O, H, OH, CH, He, Ne, Ar, Kr, Xe
%
% INPUTS:
%   species     - The species to which the model should be returned, e.g. ‘NO’. 
%
% OUTPUT:
%   C   - The collisional model as a structure
%          C.(partner).br(T)
%          C.(partner).sh(T)
%          C.(partner).cross(T)
%   e.g.
%   
%   The broadening, shift and collisions, related to a corresponding partner, 
%   at a temperature T can be acquired as.
%          T = 1500;
%          C.O2.br(T)
%          C.O2.sh(T)
%          C.O2.cross(T)
%          C.Ar.br(T)
%          C.Ar.sh(T)
%          C.Ar.cross(T)
%
% SEE ALSO:
%   overlap, dnuGFun, collisionalBroadening, excitationSpec, emissionSpec, absorptionSpec
%
% SOURCE:
%   P.H. Paul et al., Office of Scientific and Technical Information (OSTI),
%   1995-01-01, 1995.

switch species
    case 'NO'
        %% Broadening of NO with collision partner
        C.N2.br     = @(T) 0.585*(295/T).^0.75;
        C.N2.sh     = @(T) -0.18*(295/T).^0.56;
        
        C.O2.br     = @(T) 0.527*(295/T).^0.66;
        C.O2.sh     = @(T) -0.16*(295/T).^0.52;
        
        C.H2O.br    = @(T) 0.792*(295/T).^0.79;
        C.H2O.sh    = @(T) -0.21*(295/T).^0.71;
        
        C.Ar.br     = @(T) 0.505*(295/T).^0.65;
        C.Ar.sh     = @(T) -0.16*(295/T).^0.58;
        
        C.CO2.br    = @(T) 0.61*(295/T).^0.70;
        C.CO2.sh    = @(T) -0.19*(295/T).^0.55;
        
        C.CO.br     = @(T) 0.58*(295/T).^0.70;
        C.CO.sh     = @(T) -0.18*(295/T).^0.55;
        
        C.CH4.br    = @(T) 0.71*(295/T).^0.70;
        C.CH4.sh    = @(T) -0.22*(295/T).^0.55;
        
        %% Collisional quenching by other species
        % Harpoon model -> functional representation of quenching rates
        % 6 different approaches
        % (see Paul, Correlations for the NO A�E+ Electronic Quenching Cross-section, Sandia 1995)
        C.TRef       = 300;
        
        C.eta      = @(c2,c3,T) c2*(C.TRef/T)+c3*(C.TRef/T)^2;
        
        Harpoon1 = @(c0,c1,c2,c3,c4,T) 0;
        Harpoon2 = @(c0,c1,c2,c3,c4,T) c0+c1*exp(-c2*C.TRef/T)+c3*exp(-c4*C.TRef/T);
        Harpoon3 = @(c0,c1,c2,c3,c4,T) c0*((1+C.eta(c2,c3,T))*exp(-C.eta(c2,c3,T))+c1*C.eta(c2,c3,T).^(1/3)*gammainc(C.eta(c2,c3,T),5/3));
        Harpoon4 = @(c0,c1,c2,c3,c4,T) c0;
        Harpoon5 = @(c0,c1,c2,c3,c4,T) c0+c1(C.TRef/T).^c2;
        Harpoon6 = @(c0,c1,c2,c3,c4,T) c0*((1+C.eta(c2,c3,T))*exp(-C.eta(c2,c3,T))+c1*C.eta(c2,c3,T).^(1/3)*gammainc(C.eta(c2,c3,T),5/3));
        
        
        C.N2.cross  = @(T) Harpoon2(0,0.88,4.9,48,32,T);
        C.O2.cross  = @(T) Harpoon4(25.1,0,0,0,0,T);
        C.CO2.cross = @(T) Harpoon3(54.2,0.95,3.24,0.18,0,T);
        C.CO.cross  = @(T) Harpoon2(5.9,5.3,7,22.1,14,T);
        C.H2O.cross = @(T) Harpoon3(28.2,3.39,0.15,2.95,0,T);
        
        C.CH4.cross = @(T) Harpoon1(0,0,0,0,0,T);
        C.C2H6.cross= @(T) Harpoon1(0,0,0,0,0,T);
        C.C3H8.cross= @(T) Harpoon1(0,0,0,0,0,T);
        C.C2H4.cross= @(T) Harpoon6(5,7.1,1.07,0.72,0,T);
        C.C2H2.cross= @(T) Harpoon6(28,0.81,13.1,14.85,0,T);
        
        C.NO.cross  = @(T) Harpoon4(43,0,0,0,0,T);
        C.NO2.cross = @(T) Harpoon5(82,9,0.54,0,0,T);
        C.N2O.cross = @(T) Harpoon3(59,0.99,3.98,0.16,0,T);
        C.NH3.cross = @(T) Harpoon6(50.7,0.91,10.74,6.5,0,T);
        C.NH.cross  = @(T) Harpoon4(58,0,0,0,0,T);
        
        C.H2.cross  = @(T) Harpoon1(0,0,0,0,0,T);
        C.O.cross   = @(T) Harpoon4(32,0,0,0,0,T);
        C.H.cross   = @(T) Harpoon4(12,0,0,0,0,T);
        C.OH.cross  = @(T) Harpoon4(82,0,0,0,0,T);
        C.CH.cross  = @(T) Harpoon4(101,0,0,0,0,T);
        
        C.He.cross  = @(T) Harpoon1(0,0,0,0,0,T);
        C.Ne.cross  = @(T) Harpoon1(0,0,0,0,0,T);
        C.Ar.cross  = @(T) Harpoon1(0,0,0,0,0,T);
        C.Kr.cross  = @(T) Harpoon1(0,0,0,0,0,T);
        C.Xe.cross  = @(T) Harpoon1(0,0,0,0,0,T);
    case 'add new Species'
        C =  struct([]);
    otherwise
        C = struct([]);
        
end
end