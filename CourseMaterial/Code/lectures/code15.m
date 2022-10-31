%
% Time series analysis
% Lund University
%
% Example code 15: A bit more on BJ prediction (see also code 14).
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;
rng(0)                                          % Set the seed (just done for the lecture!)

extraN = 1000;
N      = 1300;
noLags = 30;
pstart = 1000;                                  % Determine where the validation data starts.

% Simulate some process (same as in code 14).
sX = 24;
A1 = [ 1 -1.8 0.82 ];
C1 = [ 1 0 -0.8 ];
A3 = conv( [ 1 zeros(1,sX-1) -1 ], [ 1 -0.4 ] );
C3 = [ 1 0.8 2.1 ];   
B  = [ 1.2 ];
A2 = 1;

% Generate the noise and the input signals.
z = filter( C1, A1, randn(N+extraN, 1) );
x = filter( C3, A3, randn(N+extraN, 1) );

% Form the output using the filtered input. Remove the initial samples.
y = filter( B, A2, x ) + z;
y = y(extraN:end);
x = x(extraN:end);
xM = x(1:pstart);
yM = y(1:pstart);

% Examine the data.
figure; 
subplot(211); 
plot( x );
line( [pstart pstart], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(x)*1.8 max(x)*1.8])
ylabel('Input signal')
title('Measured signals')
legend('Input signal', 'Prediction starts','Location','SW')
subplot(212); 
plot( y ); 
line( [pstart pstart], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(y)*1.5 max(y)*1.5])
ylabel('Output signal')
xlabel('Time')
legend('Output signal', 'Prediction starts','Location','SW')


%% Estimate the BJ model (this is the model we found in code 14).
diff_xM = filter([ 1 zeros(1,sX-1) -1 ], 1, xM);   diff_xM = diff_xM(sX+1:end);
inputModel   = estimateARMA( diff_xM, [ 1 1 1 ], [1 0 1 1], 'Differentiated input model', noLags );
inputModel.A = conv([1 zeros(1, sX-1) -1], inputModel.A);
[ foundModel, ey ] = estimateBJ( yM, xM, [1 0 1], [1 1 1], [1], [1], 'Final BJ model', noLags );
ex = filter( inputModel.A, inputModel.C, xM );
var_ex = var( ex(100:pstart));
var_ey = var( ey(100:pstart));


%% Lets try to predict the future, as if it was not known.
% Predict the future values, i.e., form the prediction for k=1, 2, 3, ...
kMax    = 20;                                 % How far ahead do we wish to predict. 
pred_yk = x(1:pstart-1);                      % This is the predicted yt vector; it is known up to the start of the prediction.
pred_xk = y(1:pstart-1);                      % This is the predicted xt vector; it is known up to the start of the prediction.
var_yk  = zeros(pstart,1);                    % This is the variance of the y prediction; the variance is zero up to the start of the prediction.
var_xk  = zeros(pstart,1);                    % This is the variance of the x prediction; the variance is zero up to the start of the prediction.

KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );    

lC3 = length(inputModel.C);
lKC = length(KC);

for k=1:kMax
    % Compute the prediction polynomials for the input signal.
    [Fx, Gx] = polydiv( inputModel.C, inputModel.A, k ); 

    % Compute the prediction polynomials for the output signal.
    [Fy, Gy]   = polydiv( foundModel.C, foundModel.D, k );
    [Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

    % Predict the input signal.
    pred_xk(pstart-1+k) = Gx*x(pstart-1+k:-1:pstart+k-length(Gx)) - ... % G x_t
        inputModel.C(2:end)*pred_yk(pstart+k-2:-1:pstart+k-lC3);        % All but the first term of C3 \hat{x}_{t+k|t}.

    % Predict the output signal.
    FhhKC = conv( Fhh, KC );
    pred_yk(pstart-1+k) = ...
        FhhKC*pred_xk(pstart-1+k:-1:pstart+k-length(FhhKC)) + ...       % KC * \hat\hat{F} \hat{x}_{t+k | t }
        Ghh*x(pstart-1+k:-1:pstart+k-length(Ghh) ) + ...                % \hat\hat{G} x_t 
        Gy*y(pstart-1+k:-1:pstart+k-length(Gy)) - ...                   % Gy y_t. 
        KC(2:end)*pred_yk(pstart+k-2:-1:pstart+k-lKC);                  % All but the first term of KC \hat{y}_{t+k|t}.

    % The theorerical variances.
    var_xk(pstart-1+k) = sum(Fx.^2)*var_ex;
    var_yk(pstart-1+k) = sum(Fy.^2)*var_ey + sum(Fhh.^2)*var_xk(pstart-1+k);
end


% Plot the output prediction and compare to the validation data.
f1 = figure;
indV = pstart:pstart+kMax-1;
plot( [y(1:pstart+kMax-1)])
xlim([pstart-30 pstart+kMax])
hold on
plot( indV, pred_yk(indV) )
indV = pstart-1:pstart+kMax-1;
plot( indV, pred_yk(indV) + 2*sqrt(var_yk(indV)),'b:' )
plot( indV, pred_yk(indV) - 2*sqrt(var_yk(indV)),'b:' )
figProp = get(f1);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
hold off
legend('Output data', 'Predicted output', 'Prediction starts','Location','NW')
title( 'Predictng the output for k=1, 2, ...')
xlabel('Time')


%% Plot the theoretical output variance.
f2 = figure;
plot(indV-pstart+1, var_yk(indV) )
hold on
plot(indV-pstart+1, var_xk(indV) )
hold off
figProp = get(f2);
line( [1 1], figProp.CurrentAxes.YLim, 'Color','red','LineStyle',':' )
xlim([-4 kMax])
title('Theoretical k-step prediction variance')
xlabel('k')
ylabel('Variance')
legend('Variance output prediction', 'Variance input prediction', 'Prediction starts', 'Location', 'SE')
