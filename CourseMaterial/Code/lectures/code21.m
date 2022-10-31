%
% Time series analysis
% Lund University
%
% Example code 22: Estimating the unknown parameters of an ARMA process
% using the Kalman filter (see also example 8.12).
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate a process.
rng(1)                                          % Set the seed (just done for the lecture!)
extraN = 100;
N  = 1000;
A0 = [1 -0.95];
C0 = [1 0.5 0 -0.2];                            
e  = randn( N+extraN, 1 );
y  = filter( C0, A0, e );   y = y(extraN+1:end);    e = e(extraN+1:end);

% Plot realisation.
figure
plot(y)
title( 'Realisation of an ARMA-process' )
ylabel('Amplitude')
xlabel('Time')
xlim([1 N])


%% Estimate the unknown parameters using a Kalman filter.
p0 = 1;                                         % Number of unknowns in the A polynomial.
q0 = 2;                                         % Number of unknowns in the C polynomial.

A     = eye(p0+q0);
Rw    = 1;                                      % Measurement noise covariance matrix, R_w
Re    = 1e-6;                                   % System noise covariance matrix, R_e
Rx_t1 = eye(p0+q0);                             % Initial covariance matrix, R_{1|0}^{x,x}
h_et  = zeros(N,1);                             % Estimated one-step prediction error.
xt    = zeros(p0+q0,N);                         % Estimated states. Intial state, x_{1|0} = 0.
yhat  = zeros(N,1);                             % Estimated output.
for t=4:N                                       % We use t-3, so start at t=4.
    % Update the predicted state and the time-varying state vector.
    x_t1 = A*xt(:,t-1);                         % x_{t|t-1} = A x_{t-1|t-1}
    C    = [ -y(t-1) h_et(t-1) h_et(t-3) ];     % Use earlier prediction errors as estimate of e_t.
    
    % Update the parameter estimates.
    Ry = C*Rx_t1*C' + Rw;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    Kt = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;                           % One-step prediction, \hat{y}_{t|t-1}.
    h_et(t) = y(t)-yhat(t);                     % One-step prediction error, \hat{e}_t = y_t - \hat{y}_{t|t-1}
    xt(:,t) = x_t1 + Kt*( h_et(t) );            % x_{t|t}= x_{t|t-1} + K_t ( y_t - Cx_{t|t-1} ) 

    % Update the covariance matrix estimates.
    Rx_t  = Rx_t1 - Kt*Ry*Kt';                  % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} - K_t R_{t|t-1}^{y,y} K_t^T
    Rx_t1 = A*Rx_t*A' + Re;                     % R^{x,x}_{t+1|t} = A R^{x,x}_{t|t} A^T + Re
end

% Examine the estimate of the driving noise.
figure
plot( [e h_et])
xlabel('Time')
title('Estimating the driving noise process')
legend('Noise process, e_t', 'Prediction error, \epsilon_{t|t-1}', 'Location','SW')


%% Examine the estimated parameters.
Q0 = [A0(2) C0(2) C0(4)];                       % These are the true parameters we seek.
figure
plot( xt' )
for k0=1:length(Q0)
    line([1 N], [Q0(k0) Q0(k0)], 'Color','red','LineStyle',':')
end
title(sprintf('Estimated parameters, with Re = %7.6f and Rw = %4.3f', Re, Rw))
xlabel('Time')
ylim([-1.5 1.5])


%% Show the one-step prediction. 
figure
plot( [y yhat] )
title('One-step prediction using the Kalman filter')
xlabel('Time')
legend('Realisation', 'Kalman estimate', 'Location','SW')


%% Examine one-step prediction residual.
ey = y-yhat;
ey = ey(N-100:N);                               % Ignore the initial values to let the filter converge first.
checkIfWhite( ey );

[~, pacfEst] = plotACFnPACF( ey, 40, 'One-step Kalman prediction' );
checkIfNormal( pacfEst(2:end), 'PACF' );
