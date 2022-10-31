%
% Time series analysis
% Lund University
%
% Example code 1: examining realisations, spectra, and ACF.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear;
close all;

% Simulate some processes.
N  = 100;                                       % How does the ACF estimate change if you increase N to 1000?
w0 = 2*pi*0.1;                                  % Periodicity of 1/0.1 = 10 samples.
x  = 2 * sin( w0*(1:N)' + rand*2*pi );          % Try reducing the amplitude to 1. Can you see the periodicity in the realisation?
w  = randn(N,1);
y  = x + w;

% Plot realisations
figure
subplot(311)
plot(x)
title('Time-domain')
ylabel('x')
subplot(312)
plot(w)
ylabel('w')
subplot(313)
plot(y)
ylabel('y')
xlabel('Time')


%% Estimate ACF. Notice the periodicity as compared to the realisation.
noLags = 20;
figure
subplot(311)
acf( x, noLags, 0.05, 1 );                      % How to call acf? Use help acf. Also, try type acf
title('Correlation-domain')
subplot(312)
acf( w, noLags, 0.05, 1 );      
subplot(313)
acf( y, noLags, 0.05, 1 );      


%% Lets estimate the power spectral density.
Padd = 100;                                     % Try increasing the zero-padding to 1024.
X = fftshift( abs( fft(x, Padd) ).^2 / N );
W = fftshift( abs( fft(w, Padd) ).^2 / N );
Y = fftshift( abs( fft(y, Padd) ).^2 / N );

ff = (0:Padd-1)'/Padd-0.5;
figure
subplot(311)
semilogy(ff, X)
title('Frequency-domain')
ylabel('x')
subplot(312)
semilogy(ff, W)                                 % Theoretically, this spectrum should be flat... 
ylabel('w')
subplot(313)
semilogy(ff, Y)
ylabel('y')
xlabel('Frequency')


%% What about windowing? 
% Notice the loss of power, the increasing with the mainlobe, as well as
% the weaker sidelobes when using the windowing. 
Y2 = fftshift( abs( fft(y.*hamming(N), Padd) ).^2 / N );
Y3 = fftshift( abs( fft(y.*blackman(N), Padd) ).^2 / N );

figure
subplot(211)
plot( ff, [Y Y2 Y3])                
title('Frequency-domain')
legend('Rectangular','Hamming','Blackman')
subplot(212)
semilogy( ff, [Y Y2 Y3])
ylabel('Log scale')
xlabel('Frequency')
%axis([-.5 .5 .1 max(Y)])                       % What does this command do?


%% Lets try do this over and over...
% Notice how the variance increases proportionally to the true spectrum.
% This is obviously not good. Windowing reduce the sidelobes, but at the
% price of resolution and loss of power - but it does not reduce the
% variance.
figure; hold on;
for k=1:100
    w  = randn(N,1);
    y  = x + w;
    %Y = fftshift( abs( fft(y, Padd) ).^2 / N );
    Y = fftshift( abs( fft(y.*hamming(N), Padd) ).^2 / N );     % Try using windowing. 
    Y = pyulear( y, 20, Padd, 'centered' );                     % Try using the Yule-Walker estimate.
    plot(ff, Y)
end
title('Frequency-domain')
ylabel('Magnitude')
xlabel('Frequency')
xlim([0.07 0.18])                               % What does this command do?
hold off



