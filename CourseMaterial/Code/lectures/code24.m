%
% Time series analysis
% Lund University
%
% Example code 24: Example of dual-input transfer function.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;
rng(0)                                          % Set the seed (just done for the lecture!)
extraN = 400;
N      = 1300;
pstart = 1000;
noLags = 40;

% Noise process.
C1 = [ 1 0 0 0.6 ];
A1 = [ 1 zeros(1,7) -.8 ];

% First input process, x1.
C3_1 = [ 1 -0.4 0.9 ]; 
A3_1 = [ 1 zeros(1,23) -.98 ];
B_1  = [ 4 ];
A2_1 = 1;

% Second input process, x2.
C3_2 = [ 1 0 0 0.8 ]; 
A3_2 = conv([1 -1],[ 1 -0.6 ]);
B_2  = [0 0 0.8 0 -1];
A2_2 = 1;

% Generate the noise and the input signals.
z = filter( C1, A1, randn(N+extraN, 1) );       % This is the noise model.
x1 = filter( C3_1, A3_1, randn(N+extraN, 1) );  % This is the first input signal.
x2 = filter( C3_2, A3_2, randn(N+extraN, 1) );  % This is the second input signal.

% Form the output using the filtered input. Remove the initial samples.
y  = filter( B_1, A2_1, x1 ) + filter( B_2, A2_2, x2 ) + z;
y  = y(extraN:end);
x1 = x1(extraN:end);
x2 = x2(extraN:end);

% Plot the processes.
figure(1)
subplot(311)
plot( x1 )
title('Input signal, x_1(t)')
line( [pstart pstart], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(x1)*1.8 max(x1)*1.8])
legend('x_1(t)', 'Location','SW')

subplot(312)
plot( x2 )
line( [pstart pstart], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(x2)*1.8 max(x2)*1.8])
legend('x_2(t)', 'Location','SW')
title('Input signal, x_2(t)')
xlim([1 N])

subplot(313)
plot( y )
line( [pstart pstart], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
axis([1 N min(y)*1.8 max(y)*1.8])
legend('y(t)', 'Prediction starts', 'Location','SW')
title('Output signal, y(t)')
xlabel('Time')



%% Create the dual-input idpoly and iddata structures for the output data. 
% See the documentation for idpoly and iddata for further details. Recall
% that the initial values might be modified (as here).
B{1,1}  = [1];                                  % B_1;
B{1,2}  = [0 0 1 0 0.3];                        % B_2;
A2{1,1} = [1];                                  % A2_1;
A2{1,2} = [1];                                  % A2_2;
A1_0    = [1 0 0 0 0 0 0 0 0.3];                % A1
C1_0    = [1 0 0 0.3];                          % C1
polyContainer = idpoly( 1, B, C1_0, A1_0, A2 );
dataContainer = iddata( y(1:pstart), [x1(1:pstart) x2(1:pstart)] );

% Determine the parameters to estimate.
polyContainer.Structure.B(1,1).Free = B_1;
polyContainer.Structure.B(1,2).Free = B_2;
polyContainer.Structure.C.Free = C1;
polyContainer.Structure.D.Free = A1;

% Estimate the polynomials from data.
foundModel = pem( dataContainer, polyContainer );
present( foundModel );


%% Estimate the models of the input signals.
inputModel1 = estimateARMA( x1(1:pstart), [ 1 zeros(1,23) 1 ], [1 1 1], 'Input model, x_1(t)', noLags );

% To handle the slowly decaying ACF, the data is first differentiated,
% after which the differentiation is again added to the model.
x2_diff = filter([1 -1], 1, x2(1:pstart));      x2_diff = x2_diff(2:end); 
inputModel2 = estimateARMA( x2_diff, [ 1 1 ], [1 0 0 1], 'Differentiated input model, x_2(t)', noLags );
inputModel2.A = conv([1 -1], inputModel2.A);


%% How do you predict a dual-input BJ?
k = 3;

% We need to form the polynomials:
%
% KA = A_1 A_2^1 A_2^2
% KB = B^1 A_1 A_2^2
% KC = C_1 A_2^1 A_2^2
% KD = B^2 A_1 A_2^1
%
% In Matlab's notation:
%   A(z) y(t) = [B(z)/F(z)] u(t) + [C(z)/D(z)] e(t)
%
% Here, u(t) is two-dimensional. This means that:
%   A(z) = 1,       B(z) = B(z),    F(z) = A2(z)
%   C(z) = C1(z),   D(z) = A1(z)
%
KA = conv( conv( foundModel.D, foundModel.F{1}), foundModel.F{2} );
KB = conv( conv( foundModel.D, foundModel.B{1}), foundModel.F{2} );
KC = conv( conv( foundModel.F{1}, foundModel.F{2}), foundModel.C );
KD = conv( conv( foundModel.D, foundModel.B{2}), foundModel.F{1} );

% Form the polynomial divisions. Note that KC/KA is the same as C1/A1. The
% two latter polynomials will also contain shared components that will
% cancel. Can you see what the resulting divisions are?
[Fy, Gy]   = polydiv( foundModel.C, foundModel.D, k );  % KC/KA = C1/A1
[Fh1, Gh1] = polydiv( conv(Fy, KB), KC, k );            % [F KB]/KC = ?
[Fh2, Gh2] = polydiv( conv(Fy, KD), KC, k );            % [F KD]/KC = ?

% Predict the input signals.
[Fx1, Gx1] = polydiv( inputModel1.C, inputModel1.A, k );
xhatk1 = filter(Gx1, inputModel1.C, x1);
[Fx2, Gx2] = polydiv( inputModel2.C, inputModel2.A, k );
xhatk2 = filter(Gx2, inputModel2.C, x2);

% Form the predicted output signal using the predicted input signals.
yhatk  = filter(Fh1, 1, xhatk1) + filter(Gh1, KC, x1) + ...
         filter(Fh2, 1, xhatk2) + filter(Gh2, KC, x2) + ...
         filter(Gy, KC, y);

% Plot the resulting predictions
f1 = figure;
subplot(211)
shiftK_x1 = round( mean( grpdelay(Gx1, 1) ) );
plot([x1(1:end-shiftK_x1) xhatk1(shiftK_x1+1:end)] )
figProp = get(f1);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
xlim([800 N-100])
title( sprintf('Shifted %i-step predictions of x_1(t)',k))
legend('x_1(t)', 'Predicted data', 'Prediction starts','Location','NW')

subplot(212)
shiftK_x2 = round( mean( grpdelay(Gx2, 1) ) );
plot([x2(1:end-shiftK_x2) xhatk2(shiftK_x2+1:end)] )
figProp = get(f1);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
xlim([800 N-100])
title( sprintf('Shifted %i-step predictions of x_2(t)',k))
legend('x_2(t)', 'Predicted data', 'Prediction starts','Location','NW')

f2 = figure;
shiftK_y = round( mean( grpdelay(Fh1, 1) ) );
plot([y(1:end-shiftK_y) yhatk(shiftK_y+1:end)] )
figProp = get(f2);
line( [pstart pstart], figProp.CurrentAxes.YLim, 'Color', 'red', 'LineStyle', ':' )
xlim([800 N-100])
title( sprintf('Shifted %i-step predictions of y(t)',k))
legend('y(t)', 'Predicted data', 'Prediction starts','Location','NW')


%% Check the resulting residuals.
ehat = y - yhatk;
ehat = ehat(pstart:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step output prediction residual', k) )
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );
