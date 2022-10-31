% FUNCTION [ signPerInd, signPerVal ] = bolvikenTest( data, Na, [alpha] )
%
% The function computes the Bolviken test to determine if data contains
% significant periodicties with significance alpha (default 0.05). The
% function returns all the frequency indices found to be significant in the 
% variable signPerInd, as well as the decision values for these indices, in 
% signPerVal. 
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ signPerInd, signPerVal ] = bolvikenTest( data, Na, alpha )

if nargin<3
    alpha = 0.05;
end
N  = length( data );
N2 = floor( N/2 - 0.5 );

rho2 = abs( fft( data ).^2 )/N;
rho2 = rho2(1:N2);
ordR = sort( rho2 );
sigmaEst = sum( ordR(1:N2-Na) );

signPerInd = zeros( N2, 1 );
signPerValVec = signPerInd;
for k=1:length( rho2 )
    ratioT = rho2(k) / sigmaEst;
    
    prodQ = 1;
    for k2=1:N2-Na
        prodQ = prodQ / ( 1 + k2*ratioT / (Na+k2-1) );
    end
    signPerValVec(k) = prodQ*N2;
end

signPerInd = find( signPerValVec<alpha )-1;
signPerVal = signPerValVec( signPerInd+1 );

