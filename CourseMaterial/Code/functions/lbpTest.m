% FUNCTION [ deemedWhite, Q, chiV ] = lbpTest( x, [K], [alpha] )
%
% The function computes the modified Ljung-Box-Pierce statistic using K
% considered correlations. The test forms the Q statistics and compares
% this to the corresponding chi^2 value. 
%
% For a univariate process, one may reject the hypothesis that the residual 
% is white if Q > chi^2_{1-alpha}(K). For a m-variate process, the test is
% Q > chi^2_{1-alpha}(Km^2).
%
% Here, alpha denotes the significance of the test, such that for a 95% 
% confidence, alpha = 0.05. Unless specified, the function uses K=20 and 
% alpha=0.05 as default.
%
% The function returns deemedWhite = 1 if the sequence is deemed white, 
% together with the Q value and the chi2 significance level.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ deemedWhite, Q, chiV ] = lbpTest( x, K, alpha )

if nargin<2
    K=20;
end
if nargin<3
    alpha = 0.05;
end

% Arrange x to have each source data as a column vector.
[m,N] = size(x);
if m>N
    x = x.';
    [m,N] = size(x);
end

if m>1
    % Form the multivariate Portmanteau statistic
    r = corrM( x, K+1 );
    Q = 0;
    for k=2:K+1
        tmp = r{1}\( r{k}.' );
        Q = Q + trace( tmp*tmp )/(N-k-1);
    end
    Q = N^2*Q;
    chiV = chi2inv( 1-alpha, K*m^2 );

else
    % Form the univariate Portmanteau statistic
    r = xcorr( x, K, 'biased' );
    r = r/max(r);
    r = r( K+2:2*K+1 );

    Q = N*(N+2)*sum( r.^2 ./ (N-(1:K) ) );
    chiV = chi2inv( 1-alpha, K );
end
deemedWhite = Q < chiV;
