%
% Time series analysis
% Lund University
%
% Example code 2: examining roots, spectra, and ACF.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate some processes. Ignore the initial samples to avoid any
% initialization effects.
extraN = 100;
N  = 100;
C  = [ 1 0.6 0 0 0.4 ];
A  = [ 1 -0.4 0.6 ];                            % Note: the roots must be inside the unit circle! Try putting them outside...
e  = randn( N+extraN, 1 );
y1 = filter( C, 1, e );     y1 = y1(extraN:end);
y2 = filter( 1, A, e );     y2 = y2(extraN:end);
y3 = filter( C, A, e );     y3 = y3(extraN:end);

% Plot realisations and the roots of the generating polynomials.
figure
subplot(311)
plot(y1)
title('Time-domain')
ylabel('MA process')
subplot(312)
plot(y2)
ylabel('AR process')
subplot(313)
plot(y3)
ylabel('ARMA process')
xlabel('Time')

% Estimate the ACF. 
noLags = 20;
figure
subplot(311)
acf( y1, noLags, 0.05, 1 );
title('ACF')
subplot(312)
acf( y2, noLags, 0.05, 1 );      
subplot(313)
acf( y3, noLags, 0.05, 1 );      


%% Examine the roots of the characteristic polynomials. 
% Given these, what should the spectrum look like?
figure
subplot(211)
zplane( C )
title('C polynomial')
subplot(212)
zplane( A )
title('A polynomial')

% Which frequencies should we look a bit more at?
rootsA = angle(roots(A))/pi/2;
rootsC = angle(roots(C))/pi/2;
fprintf('The (positive) angle of the root of the A polynomial is %3.2f degrees.\n', rootsA(1));
fprintf('The (positive) angles of the roots of the C polynomial are %3.2f and %3.2f degrees.\n', rootsC(3), rootsC(1));


%% Estimate the power spectral density using the periodogram.
Padd = 1024;
Y1 = fftshift( abs( fft(y1, Padd) ).^2 / N );
Y2 = fftshift( abs( fft(y2, Padd) ).^2 / N );
Y3 = fftshift( abs( fft(y3, Padd) ).^2 / N );

% Compute the true spectra. 
w  = exp( 1i* linspace(-.5,.5,Padd)*2*pi )';
X1 = abs( polyval(C,w) ).^2;
X2 = 1 ./ abs( polyval(A,w) ).^2;
X3 = X1 .* X2;

% Plot the resulting spectra
ff = (0:Padd-1)'/Padd-0.5;
figure
subplot(311)
semilogy(ff, [Y1 X1])
title('Frequency-domain')
ylabel('MA process')
subplot(312)
semilogy(ff, [Y2 X2])
ylabel('AR process')
subplot(313)
semilogy(ff, [Y3 X3])
ylabel('ARMA process')
xlabel('Frequency')
legend('Periodogram', 'True spectrum')

% % Estimate the PACF. Note how the PACF can be used to identify AR-processes.
% figure
% subplot(311)
% pacf( y1, noLags, 0.05, 1 );
% title('PACF')
% subplot(312)
% pacf( y2, noLags, 0.05, 1 );      
% subplot(313)
% pacf( y3, noLags, 0.05, 1 );      
