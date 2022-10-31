% FUNCTION [deemedWhite] = checkIfWhite( data, [K], [alpha] )
%
% The function computes the Monti test (using K terms, with default K = 2,
% and confidence alpha, with default = 0.05) to determine if the data is
% white nor not, returning the decision. 
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function deemedWhite = checkIfWhite( data, K, alpha )

if nargin<2
    K=20;
end
if nargin<3
    alpha = 0.05;
end

[deemedWhite, Q, chiV] = montiTest( data, K, alpha );
if deemedWhite
     fprintf('The residual is deemed to be WHITE according to the Monti-test (as %5.2f < %5.2f).\n', Q, chiV );
else
    fprintf('The residual is NOT deemed to be white according to the Monti-test (as %5.2f > %5.2f).\n', Q, chiV );
end
