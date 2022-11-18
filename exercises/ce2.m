clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

%%

n = 500; % Number of samples
A3 = [ 1 .5 ] ;
C3 = [ 1 -.3 .2 ] ;
w = sqrt ( 2 ) * randn( n + 100 , 1 ) ;
x = filter(C3 , A3,w ) ; % Create the input
A1 = [ 1 -.65];
A2 = [ 1 .90 .78 ] ;
C = 1 ;
B = [ 0 0 0 0 .4 ] ;
e = sqrt ( 1.5 ) * randn( n + 100 , 1 ) ;
y = filter(C, A1, e ) + filter(B, A2, x ) ; % Create the output
x = x ( 101 : end); 
y = y ( 101 : end); % Omit initial samples
clear A1 A2 C B e w A3 C3 

%% Examine the data.

noLags = 32;

figure; 
subplot(211); 
plot( x );
ylabel('Input signal')
title('Measured signals')
subplot(212); 
plot( y ); 
ylabel('Output signal')
xlabel('Time')

figure
[Cxy,lags] = xcorr( x, y, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(y) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between in- and output')

%% Construct model for input signal
ACFnPACFnNormplot(x,32);

% Test an AR1 first
ar_model_1 = estimateARMA(x,[1 1], [1] , "AR1",32)

% Works ok, maybee remove something at lag 4
ar_model_4 = estimateARMA(x,[1 1 0 0 1], [1], "AR4",32)

% the model choosen for x is an AR(4), which is not too close to the 
% ARMA(1,2) that generated the data.
