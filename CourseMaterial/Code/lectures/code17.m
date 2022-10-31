%
% Time series analysis
% Lund University
%
% Example code 17: Revisting the Lydia-Pinkham data set (see also code 11
% as well as examples 4.22 and 7.11).
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Load the data.
dataLydiaPinkham;
y1     = data_adver;
y2     = data_sales;
time   = 1907:1960;
noLags = 10;
N      = length(time);

figure
plot( time, [y2(:) y1(:)] )
xlabel('Year')
ylabel('Dollars (thousands)')
axis([1905 1962 0 3500 ])
legend('Sales','Advertising')
title('The Lydia-Pinkham data set')


%% Estimate the multivariate ACF and PACF. 
% Using code 11, we first differentiate the log of the data.
y1 = filter([1 -1], 1, log(y1));    y1 = y1(2:end);
y2 = filter([1 -1], 1, log(y2));    y2 = y2(2:end);

Y = [y1' y2'];                                  % Construct data matrix.
[Ry, rhoY] = corrM( Y, noLags );                % Estimate the ACF and PACF.
S = pacfM( Y, noLags );

Ra = [ rhoY{2} rhoY{3} rhoY{4} rhoY{5} rhoY{6} rhoY{7} rhoY{8} rhoY{9} ]
Sa = [ S{1} S{2} S{3} S{4} S{5} S{6} S{7} S{8}]


%% Which coefficients are significant?
fprintf('The estimated variance of the ACF and PACF are %4.2f.\n', 1/N )
abs(Ra)> 2 / sqrt(N)
abs(Sa)> 2 / sqrt(N)

fprintf('The PACF suggests a VAR(2) model might be appropriate.\n')


%% Lets try to estimate the parameters for the VAR(2) model.
% The lsVAR function returns the estimated parameters, the estimated
% variance of the residual, and the residual.
[ thEst, sigEst, resEst] = lsVAR( Y.', 2 );

fprintf('The estimated VAR(2) model is: \n')
thEst

fprintf('The estimated covariance of the noise residual is: \n')
sigEst
