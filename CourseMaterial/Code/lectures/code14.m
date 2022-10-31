%
% Time series analysis
% Lund University
%
% Example code 14: Predicting a BJ model.
%
% Note that Matlab uses a slightly different notation for the BJ model as
% compared to that used in the course. In Matlab's notation:
%
%   A(z) y(t) = [B(z)/F(z)] u(t) + [C(z)/D(z)] e(t)
%
% This means that:
%   A(z) = 1,       B(z) = B(z),    F(z) = A2(z)
%   C(z) = C1(z),   D(z) = A1(z)
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;
rng(1)                                          % Set the seed (just done for the lecture!)

extraN   = 1000;                                % Lets check again why we do this; try removing this.
N        = 1300;
noLags   = 30;
modelLim = 1000;                                % Determine where the validation data starts.

% Simulate some process.
sX = 24;
A1 = [ 1 -1.8 0.82 ];
C1 = [ 1 0 -0.8 ];
A3 = conv( [ 1 zeros(1,sX-1) -1 ], [ 1 -0.4 ] );
C3 = [ 1 0.8 2.1 ];   
B  = [ 1.2 ];
A2 = 1;

% Generate the noise and the input signals.
z = filter( C1, A1, randn(N+extraN, 1) );       % This is the noise model.
x = filter( C3, A3, randn(N+extraN, 1) );       % This is the input signal.

% Form the output using the filtered input. Remove the initial samples.
y = filter( B, A2, x ) + z;
y = y(extraN:end);
x = x(extraN:end);
xM = x(1:modelLim);                             % Extract model data.
yM = y(1:modelLim);

% Examine the data.
figure; 
subplot(211); 
plot( x );
line( [modelLim modelLim], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(x)*1.8 max(x)*1.8])
ylabel('Input signal')
title('Measured signals')
legend('Input signal', 'Prediction starts','Location','SW')
subplot(212); 
plot( y ); 
line( [modelLim modelLim], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(y)*1.5 max(y)*1.5])
ylabel('Output signal')
xlabel('Time')
legend('Output signal', 'Prediction starts','Location','SW')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a model for the input.
plotACFnPACF( xM, noLags, 'Input, x_t' );


%% There seems to be a strong periodicity at 24, suggesting that a differentiation might help.
diff_xM = filter([ 1 zeros(1,sX-1) -1 ], 1, xM);   diff_xM = diff_xM(sX+1:end);
plotACFnPACF( diff_xM, noLags, sprintf('Differentiated input, \\nabla_{%i} x_t', sX) );


%% Lets try add a1.
estimateARMA( diff_xM, [ 1 1 ], [ 1 ], 'Differentiated input, version 2', noLags );


%% We need c3, and perhaps c2?
estimateARMA( diff_xM, [ 1 1 ], [1 0 0 1], 'Differentiated input, version 3', noLags );


%% Ok, add c2 too... Maybe an a2 term as well...
% Good, now it is white and all coefficients are significant.
inputModel = estimateARMA( diff_xM, [ 1 1 1 ], [1 0 1 1], 'Differentiated input, version 4', noLags );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the BJ model.
inputModel.A = conv([1 zeros(1, sX-1) -1], inputModel.A);   % Add the differentiation to the model.
ex = filter( inputModel.A, inputModel.C, xM );              % Lets look at the residaual

figure
plot( ex )                                                  % As there are lots of roots on the unit circle, the ringing lasts longer than expected!
ylabel('Pre-filtered input signal')
xlabel('Time')
title('Prolonged ringing due to root on the unit circle')


%% Lets remove some additional samples just to be safe...
ey = filter( inputModel.A, inputModel.C, yM );   
ex = ex(length(inputModel.A)+50:end );                      % Remove some more samples given the ringing (this is much more than needed).
ey = ey(length(inputModel.A)+50:end );
var_ex = var(ex);

figure;
[Cxy,lags] = xcorr( ey, ex, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between filtered in- and output')
% Seems like we only need b0.


%% Lets form an initial model.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
estimateBJ( yM, xM, [1], [1], [1], [1], 'BJ model 1', noLags );


%% Lets try to add a1 and a2
estimateBJ( yM, xM, [1], [1 1 1], [1], [1], 'BJ model 2', noLags );


%% Better... Maybe add a c2 term?
% Yes, now it is white.
[ foundModel, ey, ~, pacfEst ] = estimateBJ( yM, xM, [1 0 1], [1 1 1], [1], [1], 'BJ model 3', noLags );
var_ey = var(ey);


%% We now have a white residual; can we trust the Monti test?
checkIfNormal( pacfEst(2:end), 'PACF' );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lets predict the input first.
k = 4;
[Fx, Gx] = polydiv( inputModel.C, inputModel.A, k )
xhatk = filter(Gx, inputModel.C, x);
 
% Compute the average group delay.
shiftK = round( mean( grpdelay(Gx, 1) ) );
fprintf('The average group delay is %i.\n', shiftK)

figure
plot([x(1:end-shiftK) xhatk(shiftK+1:end)] )
line( [modelLim modelLim], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
legend('Input signal', 'Predicted input', 'Prediction starts')
title( sprintf('Shifted predicted input signal, x_{t+%i|t}', k) )
axis([1 N min(x)*1.5 max(x)*1.5])

std_xk = sqrt( sum( Fx.^2 )*var_ex );
fprintf( 'The theoretical std of the %i-step prediction error is %4.2f.\n', k, std_xk)


%% Form the residual. Is it behaving as expected? Recall, no shift here!
ehat = x - xhatk;
ehat = ehat(modelLim:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step input prediction residual', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proceed to predict the data using the predicted input.
% Form the BJ prediction polynomials. In our notation, these are
%   A1 = foundModel.D
%   C1 = foundModel.C
%   A2 = foundModel.F
% 
% The KA, KB, and KC polynomials are formed as:
%   KA = conv( A1, A2 );
%   KB = conv( A1, B );
%   KC = conv( A2, C1 );
%
KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

% Form the ARMA prediction for y_t (note that this is not the same G
% polynomial as we computed above (that was for x_t, this is for y_t).
[Fy, Gy] = polydiv( foundModel.C, foundModel.D, k );

% Compute the \hat\hat{F} and \hat\hat{G} polynomials.
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

% Form the predicted output signal using the predicted input signal.
yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, x) + filter(Gy, KC, y);

% A very common error is to forget to add the predicted inputs. Lets try
% that to see what happens.
% yhatk  = filter(Ghh, KC, x) + filter(Gy, KC, y);

% Plot the prediction
shiftK = round( mean( grpdelay(Fhh, 1) ) );     % Beware: As we have three filters, this might be quite a poor estimate.

figure
plot([y(1:end-shiftK) yhatk(shiftK+1:end)] )
line( [modelLim modelLim], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
legend('Output signal', 'Predicted output', 'Prediction starts')
title( sprintf('Predicted output signal, y_{t+%i|t}', k) )
axis([1 N min(y)*1.5 max(y)*1.5])


%% What does the ACF look like - what should it look like?
ehat = y - yhatk;
ehat = ehat(modelLim:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step output prediction residual', k) )
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );
