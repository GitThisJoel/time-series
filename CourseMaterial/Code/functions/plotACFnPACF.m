% FUNCTION plotACFnPACF( data, noLags, titleStr, [signLvl] )
%
% The function plots the ACF and the PACF for data using noLags lags. If
% not set, the significance level (signLvl) is set to 0.05.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [acfEst, pacfEst] = plotACFnPACF( data, noLags, titleStr, signLvl )

% If the fourth argument is not given, set the value to its default.
if nargin<4
    signLvl = 0.05;
end

figure
subplot(211)
acfEst = acf( data, noLags, signLvl, 1 );
title( sprintf('ACF (%s)',titleStr))
subplot(212)
pacfEst = pacf( data, noLags, signLvl, 1 );
title( sprintf('PACF (%s)',titleStr))
