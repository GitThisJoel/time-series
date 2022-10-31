% FUNCTION [deemedNormal] = checkIfNormal( data, titleStr, [whichTest], [alpha] )
%
% The function computes a normality test to determine if the data is normal
% distributed or not. If whichTest is set to 'D' (default), the
% D'Agostino-Pearson's K2 test is computed, otherwise the Jarque-Bera test.
% 
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function deemedNormal = checkIfNormal( data, titleStr, whichTest, alpha )

if nargin<3
    whichTest = 'D';
end
if nargin<4
    alpha = 0.05;
end

if strcmp(whichTest,'D') 
    deemedNormal = dagostinoK2test(data, alpha);
    if ~deemedNormal
        fprintf( 'The D''Agostino-Pearson''s K2 test indicates that the %s is NOT normal distributed.\n', titleStr )
    else
        fprintf( 'The D''Agostino-Pearson''s K2 test indicates that the %s is NORMAL distributed.\n', titleStr) 
    end
else
    deemedNormal = jbtest(data);
    if deemedNormal
        fprintf( 'The Jarque-Bera test indicates that the %s is NOT normal distributed.\n', titleStr)
    else
        fprintf( 'The Jarque-Bera test indicates that the %s is NORMAL distributed.\n', titleStr)
    end
end

figure
normplot( data )
title( sprintf('Normal probability plot for %s', titleStr ))
