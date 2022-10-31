% FUNCTION Rx = covM( x, L, [fb] )
%
% The function estimates the L x L covariance matrix
% of the data x. If fb is set, forward-backward averaging
% is used.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function Rx = covM( x, L, fb )

N = length(x);
if nargin<3
    fb = 0;
end
x = x - mean(x);

yv = hankel( x(1:L), x(L:N) );
if fb 
    yb = conj( hankel( x(N:-1:N-L+1), x(N-L+1:-1:1) ) );    
	Rx = ( yv*yv' + yb*yb' )/( 2*(N-L+1) );
else
	Rx = ( yv*yv' )/( N-L+1 );    
end
