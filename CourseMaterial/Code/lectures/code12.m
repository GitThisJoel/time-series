%
% Time series analysis
% Lund University
%
% Example code 12: Predicting the tobacco production (see also code 4)
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Plot the examined data. 
dataTobacco;

% Divide data into model and validation data set.
data = data(:);
modelLim  = 90;
modelData = data(1:modelLim);

figure(1)
time = 1871:1984;
plot( time, data, 'r' )
hold on
plot(time(1:modelLim), modelData,'b')
line( [time(modelLim) time(modelLim)], [-40 3000 ], 'Color','red','LineStyle',':' )
hold off
axis([time(1) time(end) 0 2500 ])
xlabel('Year')
ylabel('Production')
legend('Validation data', 'Model data', 'Location','SE')
title('US tobacco production')
noLags = 20;

% For debugging, try generating data that has similar characteristics.
% N = 10000; extraN = 100; time = 1:N;
% C  = [ 1 ];
% A = conv([1 -1],[ 1 0.37 ]);
% e  = randn( N+extraN, 1 );
% data = filter( C, A, e );     data = data(extraN+1:end) + 1000;
% figure
% plot(data)
% noLags = 100;


%% Estimate unknown parameters with the model structure found in code 4.
% Differentiate to remove trend.
y = filter([1 -1], 1, data);   
y = y(2:end);

% Note that the found model is in the differentiated domain!
foundModel = estimateARMA( y, [1 1], [1], 'Differentiated data, \nabla y_t', noLags );


%% Lets form the k-step prediction.
% Predict future values. Some important things to note:
%   1) We start the prediction long before the validation data - here from
%      the beginning of the modeling data. This to save data; recall that
%      we need to omit ord(G) samples. 
%   2) The G polynomial will start with k zeros. It is (in general) not monic and is of order max( p-1, q-k ).
%   3) The F polynomial is monic of order k-1.
k = 1;                                              % Try a larger k, e.g., k=5.
foundModel.A = conv([1 -1], foundModel.A);          % Form the A polynomial taking the differentiation into account.
[F, G] = polydiv( foundModel.C, foundModel.A, k )   % Compute the G and F polynomials.
yhatk  = filter( G, foundModel.C, data );           % Form the predicted data.


%% Lets have a look at the prediction, comparing it to the true data.
% Looking at the figure, note that:
%   1) The prediction appears to be shifted k step to the data (it should be, so this is correct).
%   2) The inital k predicted values are zeros due to the zeros in the G polynomial.
figure
plot(time, [data yhatk] )
hold on
line( [time(modelLim) time(modelLim)], [-40 max(data)*2 ], 'Color','red','LineStyle',':' )
hold off
axis([time(1) time(end) 0.8*min(data(10:end)) 1.2*max(data(10:end)) ])
legend('True data', sprintf('%i-step predition', k), 'Validation starts', 'Location','NW')
title('US tobacco production')


%% To simplify the comparison, one can shift the predicted data k steps.
% This is only done for visual representation - the prediction is correct as it is! 
figure
grpdelay(G, 1)
title('Estimated group delay')
shiftK = round( mean( grpdelay(G, 1) ) );    % Compute the average group delay of the filter.
fprintf('The average group delay is %i.\n', shiftK)


%% Plot the prediction shifted by the average group delay.
figure
plot(time(1:end-shiftK), [data(1:end-shiftK) yhatk(shiftK+1:end)] )
hold on
line( [time(modelLim) time(modelLim)], [-40 max(data)*2 ], 'Color','red','LineStyle',':' )
hold off
axis([time(10) time(end) 0.8*min(data(10:end)) 1.2*max(data(10:end)) ])
legend('True data', sprintf('%i-step predition', k), 'Validation starts', 'Location','NW')
title('US tobacco production with shifted predictions')


%% How does the model work?
% We focus only on the validation data. Note that we now also remove the
% inital corrupted samples due to the filtering.
ehat = data - yhatk;
ehat = ehat(modelLim+1:end);

% Note that the prediction residual should only be white if k=1.
figure
acfEst = acf( ehat, noLags, 0.05, 1 );
title( sprintf('Prediction residual, %i-step prediction', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
checkIfWhite( ehat );


%% What does the D'Agostino-Pearson's K2 test indicate?
% As the PACF should ring for an MA(k-1) process, we only check the ACF.
checkIfNormal( acfEst(k+1:end), 'ACF' );
