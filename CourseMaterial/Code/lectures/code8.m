%
% Time series analysis
% Lund University
%
% Example code 8: Estimation using least-squares.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear;
close all;
rng(1)                                          % Set the seed (just done for the lecture!)

% Simulate some processes. To avoid mess with the LS implementation below,
% we use complex sinusoids.
N  = 100; 
w0 = 2*pi*[0.1 0.11];                           % Try with the second frequency at 0.11.
A0 = [2 0.4];
x1 = A0(1)*exp( 1i*w0(1)*(1:N)' + 1i*rand*2*pi );
x2 = A0(2)*exp( 1i*w0(2)*(1:N)' + 1i*rand*2*pi );
e  = ( randn(N,1) + 1i*randn(N,1) )/sqrt(2);    % For complex noise, you have to scale with sqrt(2).
y  = x1 + x2 + e;

% Plot realisations
figure
subplot(211)
plot([real(x1) real(x2)])
title('Time-domain')
ylabel('Signal')
legend('x1', 'x2')
subplot(212)
plot( real(y) )
ylabel('Data')
xlabel('Time')


%% Lets examine the periodogram estimates.
Padd = 512;
X1 = fftshift( abs( fft(x1, Padd) ).^2 / N );
X2 = fftshift( abs( fft(x2, Padd) ).^2 / N );
Y  = fftshift( abs( fft(y, Padd) ).^2 / N );

ff = (0:Padd-1)'/Padd-0.5;
figure
subplot(211)
semilogy(ff, [ X1 X2])
title('Frequency-domain')
ylabel('Signal')
legend('x1', 'x2')
ylim([1e-4 2e2])
subplot(212)
semilogy(ff, Y)
ylabel('Data')
xlabel('Frequency')
ylim([1e-2 2e2])


%% Try estimating the amplitudes using LS.
X = exp( 1i*(1:N)'*w0 );
theta = inv( X'*X )*X'*y;
fprintf('Estimation using LS, with known frequencies.\n')
fprintf('  True amplitudes:      %5.3f and %5.3f.\n', sort(A0) )
fprintf('  Estimated amplitudes: %5.3f and %5.3f.\n', sort(abs(theta)) )


%% What if we do not know the frequencies? 
% Lets compute the cost function for all candidate frequencies. The below
% is a very inefficient implementation - why is that?
ff = ( (1:Padd)/Padd - 0.5 );
C  = zeros(Padd, Padd);
tic;
for k1=1:Padd
    for k2 = 1:Padd
        if k1~=k2       % If it is the same frequencies, the matrix Xw'*Xw is not invertible. 
            Xw = exp( 1i*(1:N)'*[ 2*pi*ff(k1) 2*pi*ff(k2) ] );
            C(k1,k2) = y' * Xw * inv( Xw'*Xw ) * Xw' * y;
        else
            Xw = exp( 1i*(1:N)'*[ 2*pi*ff(k1) ] );
            C(k1,k1) = y' * Xw * inv( Xw'*Xw ) * Xw' * y;
        end
    end
end
fprintf('Computing the grid search took %4.2f s.\n', toc )

figure
mesh( ff, ff, real(C) )
xlabel('Frequency, f1')
ylabel('Frequency, f2')
xlim([0.08 w0(1)/2/pi+0.02])
ylim([0.08 w0(2)/2/pi+0.02])

figure
mesh( ff, ff, real(C) )
xlabel('Frequency, f1')
ylabel('Frequency, f2')
