% FUNCTION A = arIterEsacf( x, j, k )
%
% The function computes the iterative AR estimates used in the ESACF
% estimate. See Tsay & Tiao (1984) for further details.
%
% Inputs: data vector to be estimated, iteration j, model order k.
% Output: the iterated AR coefficients.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function A = arIterEsacf( x, j, k )

if j<0
    error('arIterEsacf: Illegal function call.');
    
elseif j==0
    A = arcov( x, k );
    
else
    A_k  = arIterEsacf( x, j-1, k );
    A_k1 = arIterEsacf( x, j-1, k+1 );
    
    A(1) = 1;
    for ell=1:k
        A(ell+1) = A_k1(ell+1) - A_k(ell) * A_k1(k+2) / A_k(k+1);
    end
    
end
