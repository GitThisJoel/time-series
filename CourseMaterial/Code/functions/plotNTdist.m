% FUNCTION plotNTdist( data, [titleStr] )
%
% The function shows how well the data fits a Normal and a t-distribution.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function plotNTdist( data, titleStr )

if nargin<2
    titleStr = 'Probability plot';
end
probplot( data )
p = mle( data, 'dist', 'tlo' );
t = @(data,mu,sig,df) cdf('tlocationscale', data, mu, sig, df );
h = probplot( gca, t, p );
%set( h, 'color', 'r', 'linestyle', '-' );
set( h, 'linestyle', ':' );
title( titleStr )
legend( 'Normal dist', 'Data', 't dist', 'Location', 'NW' );
