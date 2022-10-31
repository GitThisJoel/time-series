%
% Time series analysis
% Lund University
%
% Example code 10: Residual analysis (see also code 4)
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Load the data.
dataTobacco;
data = data(:);
noLags = 25;                                    % How many lags should we trust? Recall: not more than N/4.

% Recall the model from code 4. 
Cm = 1;
Am = conv([1 1], [1 -1]);

% Lets estimate the model again and check if it white.
[~, ey] = estimateARMA( data, Am, Cm, 'Modelling the tobacco data', noLags );               % Note that this function is now updated. Download it again!
figure
whitenessTest( ey )


%% What happened - why is it not white now? 
% Can you see the problem? Lets fix it.
foundModel = pem( data, idpoly(Am, [], Cm) );
present( foundModel )
ey = filter( foundModel.A, foundModel.C, data );  ey = ey(length(foundModel.A):end );
[acfEst, pacfEst] = plotACFnPACF( ey, noLags, 'Modelling the tobacco data, version 2' );     % Note that this function is now updated. Download it again!

figure
whitenessTest( ey )


%% Ok, so now it is white.
% The issue is with the use of the free command - using it, you get a
% different model in this case. Are the ACF and PACF normal distributed?

% What does the D'Agostino-Pearson's K2 test indicate?
checkIfNormal( acfEst(2:end), 'ACF', 'D' );
checkIfNormal( pacfEst(2:end), 'PACF', 'D' );


%% What does the Jarque-Bera test indicate?
checkIfNormal( acfEst(2:end), 'ACF', 'J' );
checkIfNormal( pacfEst(2:end), 'PACF', 'J' );
