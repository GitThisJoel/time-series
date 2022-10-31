%
% Time series analysis
% Lund University
%
% Example code 11: transfer function example (see also Ex. 4.22)
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Load the data.
dataLydiaPinkham;
y1 = data_adver(:);                             % Advertising
y2 = data_sales(:);                             % Sales
time   = 1907:1960;
noLags = 20;

figure
plot( time, [y2 y1] )
xlabel('Year')
ylabel('Dollars (thousands)')
axis([1905 1962 0 3500 ])
legend('Sales','Advertising')
title('The Lydia-Pinkham data set')


%% Should we transform the data?
figure
hold on
lambda_max = bcNormPlot( y1 );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
lambda_max = bcNormPlot( y2 );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
hold off
ylim([-380 -310])
xlabel('Lambda')
ylabel('Log-likelihood')

% Both seems to indicate that a log-transform could be helpful.
y1=log(y1);
y2=log(y2);


%% Lets model the sales as input.
plotACFnPACF( y1, noLags, 'Modelling the sales data' );


%% Looks like it might need to be differentiated - lets try that.
AS = [1 -1]; 
y1_diff = filter(AS, 1, y1);    y1_diff = y1_diff(length(AS):end );
[rhoEst, phiEst] = plotACFnPACF( y1_diff, noLags, 'Differentiated sales data' );


%% Seems that this was enough, the differentiated input is white.
checkIfWhite( y1_diff );    
checkIfNormal( rhoEst(2:end), 'ACF' );
checkIfNormal( phiEst(2:end), 'PACF' );


%% Lets look at the cross correlation.
% Note that A3 is the differentiation, i.e., AS, and C3 = 1.
ey = filter( AS, 1, y2 );   ey = ey(length(AS):end );

figure;
[Cxy,lags] = xcorr( ey, y1_diff, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between filtered in- and output')


%% Form the BJ model.
% It seems we get away with just a scaling, b0.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
estimateBJ( y2, y1, [1], [1], [1], [1], 'BJ model 1', noLags );         % Note that this function is now updated. Download it again!


%% There is a dependence at lag 1 in the PACF, lets add that.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
[foundModel, ey] = estimateBJ( y2, y1, [1], [1 1], [1], [1], 'BJ model 2', noLags );


%% That seems fairly white. Lets check.
[acfEst, pacfEst] = plotACFnPACF( y1_diff, noLags, 'BJ residual' );
checkIfNormal( acfEst(2:end), 'ACF' );
checkIfNormal( pacfEst(2:end), 'PACF' );


%% Did we manage to extract all the information from the input?
tilde_et = y2 - filter( foundModel.B, foundModel.F, y1 );      

% Note that we now have to remove samples from x as well.
tilde_et  = tilde_et(length(foundModel.B):end );
filter_xt = y1(length(foundModel.B):end );

figure
[Cxy,lags] = xcorr( filter_xt, tilde_et, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(y1) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between input and residual without the influence of the input')
