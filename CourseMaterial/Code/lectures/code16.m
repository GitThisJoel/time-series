%
% Time series analysis
% Lund University
%
% Example code 16: Multivariate identification (see also example 7.10).
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

rng(1)                                          % Set the seed (just done for the lecture!)
N      = 1000;                                  % Try using a shorter realization.
m      = 2;                                     % Dimension of the VAR-process.
noLags = 10;

% Simulate a VAR(2) process.
A1 = [.5 .4 ; .1 .8];                           % Try changing the AR-polynomials (remember to check for stability!)
A2 = [-.2 -.1 ; .3 .6 ];
e = randn(m, N);
Y = zeros(m, N);
for k=3:N
    Y(:,k) = -A1*Y(:,k-1) -A2*Y(:,k-2) + e(:,k);
end

% Prot the resulting processes.
figure; 
subplot(211)
plot( Y(1,:) )
xlabel('Time')
ylabel('Amplitude')
title('Simulated VAR(2) process')
subplot(212)
plot( Y(2,:) )
xlabel('Time')
ylabel('Amplitude')


%% Estimate the ACF and PACF.
[Ry, rhoY] = corrM( Y, noLags );
S = pacfM( Y, noLags );

fprintf('The estimated covariance sequence is: \n')
Rya = [ Ry{1} Ry{2} Ry{3} Ry{4} Ry{5} ] 

fprintf('The estimated ACF is: \n')
Ra  = [ rhoY{1} rhoY{2} rhoY{3} rhoY{4} rhoY{5} ] 

fprintf('The estimated PACF is: \n')
Sa  = [ S{1} S{2} S{3} S{4} S{5} ]

signBound = 2 / sqrt(N);
fprintf('The (upper) confidence interval is %6.4f.\n', signBound );

% Test is the ACF and PACF are normal distributed.
if mjbtest( [ rhoY{2} rhoY{3} rhoY{4} rhoY{5} rhoY{6} ]  )
    fprintf('The ACF is deemed NORMAL distributed.\n')
else
    fprintf('The ACF is NOT deemed normal distributed.\n')
end
if mjbtest( [ S{2} S{3} S{4} S{5} S{6} ] )
    fprintf('The PACF is deemed NORMAL distributed.\n')
else
    fprintf('The PACF is NOT deemed normal distributed.\n')
end


%% Which values are outside the confidence interval?
fprintf('The ACF values outside the confidence interval are:')
abs(Ra)>signBound
fprintf('The PACF values outside the confidence interval are:')
abs(Sa)>signBound


%% Lets estimate the AR polynomials for varying orders.
% To deal with that Matlab does not allow vectors with index 0, store the
% vectors shifted one step.
signLvl = 0.05;
for k=1:6
    [ thEst{k}, sigEst{k}, resEst{k}] = lsVAR( Y, k-1 );
end

for p=2:6
    M(p-1) = -(N-p-p*m-0.5)*log( det(sigEst{p})/det(sigEst{p-1}) );
end
chiV = chi2inv( 1-signLvl, m^2 );
fprintf('The threshold for the likelihood ratio test is %6.4f.\n', chiV);
fprintf('The likelihood ratio test, for k=1, 2, ... (the model order estimate is the last non-zero ratio).\n')

M
M>chiV


%% Lets estimate the AIC, BIC, and FPE as well.
for k=1:6
    AIC(k) = N*log( det( sigEst{k} ) ) + 2*k*m^2;
    BIC(k) = N*log( det( sigEst{k} ) ) + k*m^2 * log(N);
    FPE(k) = ( (N+m*k+1)/(N-m*k-1) )^m *det( sigEst{k} ) ;
end

figure
plot( [AIC(:)/max(AIC) BIC(:)/max(BIC) FPE(:)/max(FPE)])
line( [2 2], [0 1 ], 'Color','red','LineStyle',':' )
title('Normalized AIC, BIC, and FPE')
legend('AIC', 'BIC', 'FPE', 'True order')
xlabel('k')
ylabel('Normalized cost function')


%% Form the LBP test for the different residuals.
for k=1:6
    [ deemedWhite(k), Q(k), chiV ] = lbpTest( resEst{k} ); 
end
fprintf('The threshold for the LBP test is %6.4f.\n', chiV);

figure
plot( Q )
line( [1 6], [chiV chiV], 'Color','red','LineStyle',':' )
legend('Test statistics', 'Threshold')
title('The LBP whiteness test')
xlabel('k')
ylabel('Test statistics')
deemedWhite
