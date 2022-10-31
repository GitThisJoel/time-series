% FUNCTION [ nRatio, confInt ] = countSignChanges( data, [confLev] )
%
% The function counts the number of sign changes in the data and returns
% the ratio of sign changes. 
%
% For larger datasets, this ratio should then approximately follow a Normal 
% distribution. If the confidence level, confLev, is given, the function 
% also returns the resulting confidence interval for the resulting whiteness 
% test.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ nRatio, confInt ] = countSignChanges( data, confLev )

N = length( data );
data( data > 0 )  = 1;
data( data <= 0 ) = 0;
nRatio = length( find( diff( data ) ~= 0 ) ) / (N-1);

if nargin>1
    if ( confLev > 1 ) || ( confLev < 0 )
        error('countSignChanges: not a legal probability.');
    end
    normConfLev = norminv( (1-(1-confLev)/2), 0, 1 );
    confInt = ( (N-1)/2 + normConfLev*sqrt( (N-1)/4 )*[-1 1] ) / (N-1);    
end