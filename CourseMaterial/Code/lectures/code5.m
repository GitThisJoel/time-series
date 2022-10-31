%
% Time series analysis
% Lund University
%
% Example code 5: why should one remove the initial samples when filtering?
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear;
close all;

% Simulate a sinusoidal processes.
N  = 100;
T  = 10;                                % The periodicity in the data.
w0 = 2*pi*1/T;
x  = 2 * sin( w0*(1:N)' + rand*2*pi );
w  = randn(N,1);
y  = x + w;

% Remove a periodicity of T.
AS = [1 zeros(1, T-1) -1];
z1 = filter(AS, 1, y);
z2 = filter(AS, 1, x);

% Plot the resulting filtered signals. Notice that if you would just plot
% the differentiated data with noise, it is not obvious that the initial T
% samples are corrupted.
subplot(211)
plot( [x z2] )
title('Without noise')
subplot(212)
plot( [y z1] )
title('With noise')
xlabel('Time')


%% Compute the residual using resid instead 
% Note that the initial samples are now just set to zero! 
eh = resid( iddata(y), idpoly( AS,[],1 ) );
figure
plot( [z1 eh.y] )
title('Differentiated data')
legend('Filter','Resid')
xlabel('Time')

