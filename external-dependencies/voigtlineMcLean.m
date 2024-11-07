function V = voigtlineMcLean(nu,nu0,dnuG,dnuL,aL)
% voigtlineMcLean - Calculates a Voigt profile for modeling a single spectral transition.
%
%   V = voigtlineMcLean(nu,nu0,dnuG,dnuL,aL)
%
% Authors:    Björn Daniel Brans and Abbas El Moussawi
%
% DESCRIPTION:
%   Calculates a Voigt profile for modeling a single spectral transition based on the 
%  implementation of McLean et al. (1994). The line shapes are calculated
%  in non-dimensional form, for usage in modelling absorption, emission and excitation
%  spectra, where subsequently each line-shape is weighed accordeningly.
%
% INPUTS:
%   nu     - Frequency range 
%   nu0    - The central frequency or reference frequency.
%   dnuG   - The Gaussian component of the transition profile.
%   dnuL   - The Lorentzian component of the transition profile.
%   aL     - The Lorentzian amplitude.
%
% OUTPUT:
%   V     - The Voigt profile as a vector.
%
% SEE ALSO:
%   overlap, dnuGFun, collisionalBroadening, excitationSpec, emissionSpec, absorptionSpec
%
% SOURCE:
% Adapted from C implementation: McLean et al., J. Electron Spectrosc. Rel.
% Phen. 69, 125 (1994) 
% 
 
n_spec_pos = length(nu);     
% number of spectral positions where to calculate Voigt function

A = repmat([-1.215; -1.3509; -1.215; -1.3509],1,n_spec_pos);
B = repmat([1.2359; .3786; -1.2359; -.3786],1,n_spec_pos);
C = repmat([-.3085; .5906; -.3085; .5906],1,n_spec_pos);
D = repmat([.0210; -1.1858; -.021; 1.1858],1,n_spec_pos);

X=2*sqrt(log(2))*(nu-nu0)./dnuG;     % == X (Eq. 10) in McLean
X=repmat(X,size(A,1),1);
Y=sqrt(log(2))*dnuL./dnuG;        % == Y (Eq. 11) in McLean
constant = aL*dnuL*sqrt(pi)*sqrt(log(2))/dnuG;
alphaV = C.*(Y-A)+D.*(X-B);
betaV = (Y-A).^2+(X-B).^2;
%
V = sum(alphaV./betaV);

V = constant*V;


% dVdX = sum(D./betaV-2*(X-B).*alphaV./betaV.^2);
% dVdY = sum(C./betaV-2*(Y-A).*alphaV./betaV.^2);

% Use this line if you want a real voigt profile

% % %g=voigt(x,a);g=g./trapz(v,g);%g=g./mytrapzz(v,g); 			
% % % Use this line if you have to renormalize the area anyway
% % % Form derivatives of line intensities --> for Jacobian in fitting routine:
% dgdaL = constant*V/aL;
% dgdv0 = -constant*dVdX*2*sqrt(log(2))/delvD;
% dgddelvL = constant*(V/delvL+dVdY*sqrt(log(2))/delvD);
% dgddelvD = constant*(V+(sqrt(log(2))/delvD)*(2*(X(1,:)-v0).*dVdX+delvL*dVdY))/delvD;    
% % % x is a matrix!