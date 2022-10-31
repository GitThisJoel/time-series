% FUNCTION [ deemedWhite, Q, chiV ] = montiTest( x, [K], [alpha] )
%
% The function computes the Monti statistic using K considered 
% correlations. With significance alpha, one may reject the hypothesis 
% that the residual is white if Q > \chi^2_{1-alpha}(K). Thus, for a 
% 95% confidence, alpha = 0.05. Unless specified, the function uses K=20 
% and alpha=0.05 as default.
%
% The function returns deemedWhite = 1 if the sequence is deemed white, 
% together with the Q value and the chi2 significance level.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ deemedWhite, Q, chiV ] = montiTest( x, K, alpha )

if nargin<2
    K=20;
end
if nargin<3
    alpha = 0.05;
end
N = length(x);

r = pacf( x, K );
Q = N*(N+2)*sum( r(2:end).^2 ./ (N-(1:K)') );

chiV = chi2inv( 1-alpha, K );
deemedWhite = Q < chiV;
