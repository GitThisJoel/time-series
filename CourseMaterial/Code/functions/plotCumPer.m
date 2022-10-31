% FUNCTION [ C, x1, x2, Ka ] = plotCumPer( x, alpha, plotIt )
%
% The function plots the cumulative periodogram for some typical settings
% of alpha, and returns the cumulative periodogram together with the 
% corresponding confidence interval. The function currently supports alpha 
% values of 0.1, 0.05, 0.1, and 0.25.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ C, x1, x2, Ka ] = plotCumPer( x, alpha, plotIt )

if nargin<3
    plotIt = 0;
end

N = length(x);
X = abs(fft(x) ).^2;
S = sum(X)/2;

n = floor(N/2);
C = zeros(n,1);
for k=1:n
    C(k) = sum( X(1:k) )/S;
end

% Find confidence interval.
switch alpha
    case 0.01, Ka = 1.63;
    case 0.05, Ka = 1.36;
    case 0.10, Ka = 1.22;
    case 0.25, Ka = 1.02;
    otherwise
        % This is an unknown alpha value
        Ka = 0;
        C  = -1;
end

x1 = (0:n-1)/n - Ka / sqrt( floor( (n-1)/2 ) );
x2 = (0:n-1)/n + Ka / sqrt( floor( (n-1)/2 ) );

if plotIt && all(C~=-1)
    ff = 0.5*(0:n-1)/n;
    plot( ff, C )
    hold on
    plot( ff, x1, '--')
    plot( ff, x2, '--')
    hold off
    %legend( sprintf('%d%% significance test', alpha*100) )
    axis([0 0.5 0 1])
    xlabel('Frequency')
    ylabel('C(w)')
end
