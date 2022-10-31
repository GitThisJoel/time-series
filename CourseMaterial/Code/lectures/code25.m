%
% Time series analysis
% Lund University
%
% Example code 25: Example of predicting transformed data.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;
rng(0)                                          % Set the seed (just done for the lecture!)
extraN   = 400;
N        = 1200;
valStart = 900;                                 % Determine where the validation data starts.
noLags   = 100;
time     = ( 1:N )/N*45*12*2;                   % The lectures have lasted 12 x 45 x 2 min.

% Simulate and plot a transformed process,
C  = [ 1 ];
A = conv([1 -1],[ 1 0.37 ]);
e  = 0.1*randn( N+extraN, 1 );
data = exp( filter( C, A, e )+8 );     data = data(extraN+1:end);

figure
plot(time, data)
line( [time(valStart) time(valStart)], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
title('Number of Dante references during the course (or?)')
xlabel('Time (minutes)')
axis([time(1) time(N) 0 max(data)*1.2])
legend('Data', 'Prediction starts','Location','NW')
fprintf('The mean of the data is %4.2f.\n', mean(data))


%% Should we transform the data?
% We examine all the data when making this decision. See also table 4.6.
figure; 
lambda_max = bcNormPlot(data,1);
fprintf('The Box-Cox curve is maximized at %4.2f. This suggests that a log-transform might be helpful.\n', lambda_max)

figure
normplot( data )


%% Lets transform the data. 
logData = log( data );

figure
normplot( logData )

figure
plot(time, logData)
title('log-transformed data')
xlabel('Time (minutes)')
xlim([time(1) time(N)])


%% Extract the modelling data and examine the ACF and PACF.
logModellingData = logData(1:valStart);
plotACFnPACF( logModellingData, noLags, 'log-transformed data' );


%% The ACF is decaying very slowly, suggesting we need to differentiate the data.
y = filter([1 -1], 1, logModellingData);   
y = y(2:end);
plotACFnPACF( y, noLags, 'Differentiated log-transformed data' );


%% Looks like we need an a1 term. 
% Great, this gives us a white residual.
dataModel = estimateARMA( y, [ 1 1 ], [1], 'Differentiated log-transformed data', noLags );

% Add the differentation to the model.
dataModel.A = conv([1 -1], dataModel.A);   % Add the differentiation to the model.


%% Lets predict the signal.
k = 4;
[F, G] = polydiv( dataModel.C, dataModel.A, k );
yhatk = filter(G, dataModel.C, logData );
 
% Compute the average group delay.
shiftK = round( mean( grpdelay(G, 1) ) );
fprintf('The average group delay is %i.\n', shiftK)

% This is in the wrong domain!
figure
plot(time(1:end-shiftK), [logData(1:end-shiftK) yhatk(shiftK+1:end)] )
line( [time(valStart) time(valStart)], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
legend('log-transformed data', 'Predicted log-data', 'Prediction starts','Location','SW')
title( sprintf('Shifted predicted log-transformed data, x_{t+%i|t}', k) )
axis([time(1) time(N) 4 10])

% Lets examine the residuals.
ehat = logData - yhatk;
ehat = ehat(valStart:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step input prediction log-transformed residual', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF (log-data)' );


%% Plot the data in the correct domain.
yhatkE = exp( yhatk );
figure
plot(time(1:end-shiftK), [data(1:end-shiftK) yhatkE(shiftK+1:end)] )
line( [time(valStart) time(valStart)], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
legend('Data', 'Predicted data', 'Prediction starts','Location','SW')
title( sprintf('Shifted predicted data, x_{t+%i|t}', k) )
axis([time(1) time(N) 0 max(data)*1.2])
legend('Data', 'Prediction starts','Location','NW')


%% Examine the prediction residual.
ehat = data - yhatkE;
ehat = ehat(valStart:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step input prediction residual', k) )
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF (data)' ); 
