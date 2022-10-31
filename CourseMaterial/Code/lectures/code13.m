%
% Time series analysis
% Lund University
%
% Example code 13: A bit more on linear prediction.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate some processes.
noLags = 40;
extraN = 100;
pstart  = 800;                                  % Start the prediction at this sample.
rng(0)                                          % Set the seed (just done for the lecture!)
N  = 1000;                                      % Try reducing the number of samples.
C  = [ 1 -0.6 0.5 ];                            % Try other polynomials.
A  = conv([1 -.96],[ 1 -0.4 .3 ]);
data = filter( C, A, randn(N+extraN,1) );     data = data(extraN:end);

f1 = figure;

plot( data )
figProp = get(f1);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
legend('Measured data', 'Prediction starts','Location','NW')
title('Simulated data to predict')


%% Part 1: Form the prediction using the true polynomials. 
% These are of course not known, but when examining an issue, you want to
% have as few complicating factors as possible. 
k = 4;                                          % Try varying k, say k=3.
[F, G] = polydiv( C, A, k )                     % Compute the G and F polynomials.
yhatk  = filter( G, C, data );                  % Form the predicted data.
 
% Compute the average group delay.
shiftK = round( mean( grpdelay(G, 1) ) );
fprintf('The average group delay is %i.\n', shiftK)


%% Plot the prediction
% Recall that the shifting is just for visualization. The predictions are
% correct as they are.
f1 = figure;
plot([data(1:end-shiftK) yhatk(shiftK+1:end)] )
figProp = get(f1);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
xlim([600 900])
title( sprintf('Shifted %i-step predictions',k))
legend('Measured data', 'Predicted data', 'Prediction starts','Location','NW')


%% Form the residual. Is it behaving as expected? Recall, no shift here!
ehat = data - yhatk;
ehat = ehat(100+1:end);

figure
acfEst = acf( ehat, noLags, 0.05, 1 );
title( sprintf('Prediction residual, %i-step prediction', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
checkIfWhite( ehat );



%% Part 2: Lets try to predict the future, as if it was not known (this is actually often the case :-)
% Predict the future values, i.e., form the prediction for k=1, 2, 3, ...

% Begin with estimating the model - to simplify, we use the true model.
[foundModel, ey] = estimateARMA( data, A, C, 'Model residual', noLags );
q = length(foundModel.C);

% We proceed to form the predictions, one for each k.
kMax    = 20;                                   % How far ahead do we wish to predict. 
pred_yk = data(1:pstart-1);                     % This is the predicted vector; it is known up to the start of the prediction.
var_yk  = zeros(pstart,1);                      % This is the variance of the prediction; the variance is zero up to the start of the prediction.
var_et  = var(ey);                              % Estimate the noise variance using the model residual.
for k=1:kMax
    % Compute the prediction polynomials.
    [F, G] = polydiv( foundModel.C, foundModel.A, k ); 

    % Form the prediction - note that you have to use the earlier predicted
    % values, not the true values, for the C \hat{y}_{t+1|t} part. It is
    % easy to get errors on this part. Use simulate data to check, at first
    % using C=1. 
    pred_yk(pstart-1+k) = G*data(pstart-1+k:-1:pstart+k-length(G)) - ...    % This is G y_t. Note that you have to reverse the time order.
        foundModel.C(2:end)*pred_yk(pstart+k-2:-1:pstart+k-q);              % This is all but the first term of C \hat{y}_{t+k|t}.

    % The theorerical variance is estimated from the F polynomial.
    var_yk(pstart-1+k) = sum( F.^2 )*var_et;
end


% Plot the prediction and compare to the validation data.
f1 = figure;
indV = pstart:pstart+kMax-1;
plot( [data(1:pstart+kMax-1)])
xlim([pstart-50 pstart+kMax])
hold on
indV = pstart-1:pstart+kMax-1;
plot( indV, pred_yk(indV) )
plot( indV, pred_yk(indV) + 2*sqrt(var_yk(indV)),'b:' )
plot( indV, pred_yk(indV) - 2*sqrt(var_yk(indV)),'b:' )
figProp = get(f1);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
hold off
legend('Measured data', 'Predicted data', 'Prediction starts','Location','NW')
title( 'Predictng the signal for k=1, 2, ...')
xlabel('k')


%% Plot the theoretical variance.
f2 = figure;
plot(indV-pstart+1, var_yk(indV) )
figProp = get(f2);
line( [1 1], figProp.CurrentAxes.YLim, 'Color','red','LineStyle',':' )
xlim([-4 kMax])
title('Theoretical k-step prediction variance')
xlabel('k')
ylabel('Variance')
legend('Variance','Prediction starts','Location','SE')
