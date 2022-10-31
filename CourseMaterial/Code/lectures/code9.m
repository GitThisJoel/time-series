%
% Time series analysis
% Lund University
%
% Example code 9: Model order estimation.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate some processes. Ignore the initial samples to avoid any
% initialization effects.
rng(0)                                          % Set the seed (just done for the lecture!)
extraN = 100;
N  = 100;                                      % Try using only N=100 samples.
C  = [ 1 -0.6 0.5 ];
A  = [ 1 0.2 0.3 0.8 ];
p  = length(A)-1;
q  = length(C)-1;
e  = randn( N+extraN, 1 );
y1 = filter( 1, A, e );     y1 = y1(extraN:end);
y2 = filter( C, A, e );     y2 = y2(extraN:end);

% Examine the data.
figure(1); 
subplot(211); 
plot( y1 );
ylabel('AR process')
title('Simulated data')
xlim([1 N])
subplot(212); 
plot( y2 ); 
ylabel('ARMA process')
xlabel('Sample')
xlim([1 N])


%% Form the least squares estimate using an AR model of varying order.
tic; pMax = 8;
S_th = zeros(pMax,1);
for p=1:pMax
    modelAR = arx(y1,p);
    ey = filter( modelAR.A, modelAR.C, y1 );  ey = ey(length(modelAR.A):end );
    S_th(p) = N*log( var(ey) );
end
    
figure(2)
bic_k = (1:pMax).'*log(N);
plot( [ S_th bic_k S_th+bic_k ])
xlabel('Order, p')
ylabel('Sum of squared residual')
legend('S_\theta', 'BIC penalty', 'BIC(k)')
title( sprintf('Estimating the model order of an AR(%i)', length(A)-1) )
fprintf('Average computational cost per estimate: %4.2f ms.\n', 1000*toc/pMax)


%% Lets try the same for the ARMA process. 
% Note the difference in computational cost!
tic; qMax = 6;
S_th = zeros(pMax,qMax);
for p=1:pMax
    for q=1:qMax
        modelARMA = armax(y2,[p q]);
        ey = filter( modelARMA.A, modelARMA.C, y2 );  ey = ey(length(modelARMA.A):end );
        S_th(p,q) = N*log( var(ey) );
    end
end

figure(3)
bic_pq = ( (1:pMax).'+(1:qMax) )*log(N);
mesh( log(S_th+bic_pq) )
xlabel('q')
ylabel('p')
zlabel('Cost')
title( sprintf('BIC(p,q) of an ARMA(%i,%i)', length(A)-1, length(C)-1) )
fprintf('Average computational cost per estimate: %4.2f ms.\n', 1000*toc/pMax/qMax)
