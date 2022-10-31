% FUNCTION [ signPerFt, T, g ] = fisherTest( data, [alpha], [indF] )
%
% The function computes the Fisher test to determine if data contains a
% significant periodicty with significance alpha (default 0.05). If given,
% the test is computed for frequency index indF (given as an index from 1 
% to floor(N/2-0.5).
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ signPerFt, T, g ] = fisherTest( data, alpha, indF )

if nargin<2
    alpha = 0.05;
end
N = length( data );
X = abs( fft( data ).^2 )/N;

N2 = floor( N/2 - 0.5 );
X2 = X(1:N2);

if nargin<3
    T = max( X2 ) / sum( X2 );
else
    if indF<1 || indF>N2
        error('fisherTest: index out of bound.');
    end
    T = X2( indF ) / sum( X2 );
end
g = 1 - ( alpha/N2 ).^( 1/(N2-1) ); 

signPerFt = T > g;
