%
% Time series analysis
% Lund University
%
% Example code 3: examining the ACF and PACF.
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
C  = [ 1 0.6 0.5 -0.8 ];
A  = [ 1 -0.6 0.8 ];
e  = randn( N+extraN, 1 );
y1 = filter( C, 1, e );     y1 = y1(extraN:end);
y2 = filter( 1, A, e );     y2 = y2(extraN:end);
y3 = filter( C, A, e );     y3 = y3(extraN:end);

% Plot realisations as well as roots of the AR polynomial.
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

% Check that the polynomial is stable. The roots should be inside the unit
% cirle.
figure
zplane(A)


%% Examine ACF and PACF. 
% Can you determine which model that could be suitable? Which model order
% to try? Note how similar the results are - does this mean either models
% is just as good? 
noLags = 20;    p = length(A)-1;    q = length(C)-1;
figure
subplot(311)
acf( y1, noLags, 0.05, 1 );
title( sprintf('ACF (p=%i, q=%i)', p, q ))
subplot(312)
acf( y2, noLags, 0.05, 1 );      
subplot(313)
acf( y3, noLags, 0.05, 1 );      

figure
subplot(311)
pacf( y1, noLags, 0.05, 1 );
title( sprintf('PACF (p=%i, q=%i)', p, q ))
subplot(312)
pacf( y2, noLags, 0.05, 1 );      
subplot(313)
pacf( y3, noLags, 0.05, 1 );      


%% Examine the ESACF. This is not always all that helpful.
[ esacfM, esacfX, condInt ] = esacf( y3, 6, 8 );
fprintf('This is an ARMA(%i,%i). The noise threshold is %2.4f.\n', p, q, condInt );
    esacfM
    esacfX

    
%% What happens if there is a season in the data?
s = 5;
A = conv( A, [ 1 zeros(1,s-1) -1 ] );
y1 = filter( C, 1, e );     y1 = y1(extraN:end);
y2 = filter( 1, A, e );     y2 = y2(extraN:end);
y3 = filter( C, A, e );     y3 = y3(extraN:end);

% Plot realisations, ACF, and PACF.
p = length(A)-1;    q = length(C)-1;
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

figure
subplot(311)
acf( y1, noLags, 0.05, 1 );
title( sprintf('ACF (p=%i, q=%i)', p, q ))
subplot(312)
acf( y2, noLags, 0.05, 1 );      
subplot(313)
acf( y3, noLags, 0.05, 1 );      

figure
subplot(311)
pacf( y1, noLags, 0.05, 1 );
title( sprintf('PACF (p=%i, q=%i)', p, q ))
subplot(312)
pacf( y2, noLags, 0.05, 1 );      
subplot(313)
pacf( y3, noLags, 0.05, 1 );      
