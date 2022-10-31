% FUNCTION [ deemedWhite, Q, chiV ] = mlTest( x, [K], [alpha] )
%
% The function computes the McLeod-Li statistic using K considered 
% correlations. With significance alpha, one may reject the hypothesis 
% that the residual is white if Q > \chi^2_{1-alpha}(K). Thus, for a 
% 95% confidence, alpha = 0.05. Unless specified, the function uses K=20 
% and alpha=0.05 as default.
%
% The function returns deemedWhite = 1 if the sequence is deemed white, 
% together with the Q value and the chi2 significance level.
%

% Reference: 
%   "An Introduction to Time Series Modeling" by Andreas Jakobsson
%   Studentlitteratur, 2013
%
function [ deemedWhite, Q, chiV ] = mlTest( x, K, alpha )

N = length(x);
if nargin<2
    K=20;
end
if nargin<3
    alpha = 0.05;
end

r = xcorr( x.^2 - mean( x.^2 ), K, 'biased' );
r = r/max(r);
r = r( K+2:2*K+1 );

Q = N*(N+2)*sum( r.^2 ./ (N-(1:K)') );

chiV = chi2inv( 1-alpha, K );
deemedWhite = Q < chiV;

