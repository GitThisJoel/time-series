% FUNCTION [ esacfM, esacfX, condInt ] = esacf( x, pMax, qMax )
%
% The function computes the ESACF estimate. See Tsay & Tiao (1984) for 
% further details on the algorithm.
%
% Inputs:  data vector to be estimated, maximum AR(p) and MA(q) orders.
% Outputs: the ESACF matrix as well as the corresponding indicator matrix, 
% and the zero truncation level.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ esacfM, esacfX, condInt ] = esacf( x, pMax, qMax )

% Initialize.
x = x(:) - mean(x);
N = length( x );
esacfM = zeros( pMax, qMax );
residM = zeros( N-1, pMax );

% Form the ESACF matrix. For p=0, simply use the ACF of the process.
r = xcorr( x, qMax, 'biased' );
r = r/max(r);
r = r( qMax+2:end );
esacfM(1,:) = r';

% Then, for each order, compute the 
for k = 1:pMax
    for j = 0:qMax-1 
        % Compute the iterated AR coefficients.
        A = arIterEsacf( x, j+1, k );
           
        % Using the estimated AR-coeffients, estimate the MA-part. Remove
        % the samples corrupted due to the lack of initialization of the
        % filter.
        z = filter( A, 1, x );
        z = z( length(A):end );
        
        % Form the ESACF as the ACF of the resulting process.
        r = xcorr( z, qMax, 'biased' );
        r = r/max(r);
        r = r( qMax+2:end );

        esacfM(k+1,j+1) = r(j+1);
    end
end

condInt = 2/sqrt(length(x)-pMax-qMax);
esacfX = abs(esacfM) > condInt;



