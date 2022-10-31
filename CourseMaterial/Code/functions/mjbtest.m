% FUNCTION [ isNormal ] = mjbtest( data, signLvl )
%
% The function computes the multivariant version of the Jarque-Bera test of 
% normality with significance signLvl (default 0.05). 
%
% The function returns isNormal == 1 if the data sequence is deemed to be
% normal distributed, otherwise 0.

%
% Reference: 
%   "An Introduction to Time Series Modeling" by Andreas Jakobsson
%   Studentlitteratur, 2013
%
function [ isNormal ] = mjbtest( data, signLvl )

if nargin<2
    signLvl = 0.05;
end
if (signLvl <= 0 || signLvl >= 1)
   error('dagostinoK2test: significance level must be between 0 and 1\n');
end

% Arrange x to have each source data as a column vector.
[m,N] = size(data);
if m>N
    data = data.';
    [m,N] = size(data);
end

% Form Sx.
mx = mean( data.' ).';
Sx = zeros(m,m);
for t=1:N
    Sx = Sx + ( data(:,t)-mx )*( data(:,t)-mx )';
end
Sx = Sx./(N-1);

% Form Vx, b1, and b2.
b1 = zeros(m,1); 
b2 = b1;
Ps = inv( chol( Sx )' );
for t=1:N
    V(:,t) = Ps*( data(:,t)-mx );
end
for t=1:m
    b1(t) = sum( V(t,:).^3 ) / N;
    b2(t) = sum( V(t,:).^4 ) / N;
end

% Form test statistics, lambda_s and lambda_k, and form the test.
lambda_s = N*b1'*b1/6;
lambda_k = N*( b2-3*ones(m,1) )'*( b2-3*ones(m,1) )/24;

P = 1-chi2cdf(lambda_s+lambda_k,2*m); % Joint test statistics.
isNormal = P>signLvl;

