%
% Time series analysis
% Lund University
%
% Example code 20: State space model of an AR(1) process with unknown a1
% parameter (see also examples 8.3 and 8.9). 
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate an AR(1) process.
rng(1)                                          % Set the seed (just done for the lecture!)
N  = 500;
A0_start = -0.95;
A0_stop  = -0.45;
A1 = linspace( A0_start,A0_start, N);
%A1(N/2:end) = A0_stop;                         % Abruptly change a_1.
y = zeros(N,1);
e = randn( N, 1 );
for k=2:N                                       % Implement filter by hand to allow a1 to change.
    y(k) = e(k) - A1(k)*y(k-1);
end

% Show the a_1 value over time.
figure 
plot( A1 )
title('True value of a_1')
xlabel('Time')
ylim([-1 0])

% Set the state as the "unknown" parameter a1 (for now, treated as known).
A  = 1; 
B  = 1;
xt = A0_start;
ys = zeros(N,1);
for t=2:N
    C     = [ -y(t-1) ];                        % The C matrix is now time-varying.    
    xt    = A*xt;
    ys(t) = C*xt + e(t);
end

figure
plot( [y ys] )
title('Realisation of an AR(1) process')
xlabel('Time')
legend('Transfer function form', 'State space form', 'Location','SW')


%% What if we do not know a1?
% Construct a Kalman filter to estimate a_1 (see theorem 8.7). Try
% modifying Re and Rw to see how these affect the estimate (make large 
% changes). Also try different intial covariance matrices - what does it
% mean if you select this "small" or "large"?
Rw    = 1e-1;                                   % Measurement noise covariance matrix, R_w
Re    = 1e-6;                                   % System noise covariance matrix, R_e
Rx_t1 = 1;                                      % Initial covariance matrix, R_{1|0}^{x,x}
xt    = zeros(1,N);                             % Estimated states. Intial state, x_{1|0} = 0.
yhat  = zeros(N,1);                             % Estimated output.
for t=2:N                                       % We use t-1, so start at t=2.
    % Update the predicted state and the time-varying state vector.
    x_t1 = A*xt(:,t-1);                         % x_{t|t-1} = A x_{t-1|t-1}
    C    = [ -y(t-1) ];     
    
    % Update the parameter estimates.
    Ry = C*Rx_t1*C' + Rw;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    Kt = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;                           % \hat{y}_{t|t-1} = Cx_{t|t-1}
    xt(:,t) = x_t1 + Kt*( y(t)-yhat(t) );       % x_{t|t}= x_{t|t-1} + K_t ( y_t - \hat{y}_{t|t-1} ) 

    % Update the covariance matrix estimates.
    Rx_t  = Rx_t1 - Kt*Ry*Kt';                  % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} - K_t R_{t|t-1}^{y,y} K_t^T
    Rx_t1 = A*Rx_t*A' + Re;                     % R^{x,x}_{t+1|t} = A R^{x,x}_{t|t} A^T + Re
end

% Show the estimated a1 parameter.
figure
plot( xt )
line([1 N/2], [A0_start A0_start], 'Color','red','LineStyle',':')
line([N/2 N/2+1], [A0_start A1(N/2+1)], 'Color','red','LineStyle',':')
line([N/2+1 N], [A1(N/2+1) A1(N/2+1)], 'Color','red','LineStyle',':')
title(sprintf('Estimated parameter, a_1, with Re = %7.6f and Rw = %4.3f', Re, Rw))
legend('Estimated a_1', 'True a_1', 'Location','NW')
xlabel('Time')


%% Show the one-step prediction. 
figure
plot( [y yhat] )
title('One-step prediction using the Kalman filter')
xlabel('Time')
legend('Realisation', 'Kalman estimate', 'Location','SW')


%% How good is this prediction?
ey = y-yhat;
ey = ey(100:end);                               % Ignore the initial values to let the filter converge first.
checkIfWhite( ey );

[~, pacfEst] = plotACFnPACF( ey, 40, 'Kalman prediction' );
checkIfNormal( pacfEst(2:end), 'PACF' );

