%
% Time series analysis
% Lund University
%
% Example code 21: Estimating the unknown parameters of an ARMA process
% using the Kalman filter with missing samples.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate a process.
rng(0)                                          % Set the seed (just done for the lecture!)
extraN = 100;
N  = 1000;
A0 = [1 -0.95];
C0 = [1 0.5 0 -0.2];                            
e  = randn( N+extraN, 1 );
y  = filter( C0, A0, e );   y = y(extraN+1:end);
tt = 1:N;

% Insert missing samples, marked with NaN.
noNan = 30;                                     % Number of missing samples.
noVal = randi( [100 N-100], noNan, 1);          % Lets stress the system a bit.
%noVal = randi( [700 850], noNan, 1);
y0 = y;                                         % Retain the original.
y(noVal) = NaN;                                % Replace missing samples.
y1 = y;                                         % Replica to allow for modifications.

% Plot realisation.
figure
plot(y)
ylabel('Amplitude')
xlabel('Time')
if sum(isnan(y))
    hold on
    plot( tt(noVal), y0(noVal), 'r*')
    hold off
    legend('Realisation', 'Missing sample', 'Location','SW')
    title('Realisation of an ARMA-process with missing samples')
else 
    legend('Realisation', 'Location','SW')
    title('Realisation of an ARMA-process')
end
xlim([1 N])


%% Lets form the one-step prediction using the Kalman filter.
% Construct a Kalman filter that assumes the model parameters to be known.
% Note how this differs from the setup in code20, where the parameters
% where instead treated as unknown.
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
    C    = [ -y1(t-1) h_et(t-1) h_et(t-3) ];    % Use earlier prediction errors as estimate of e_t.
    
    % Update the parameter estimates.
    Ry = C*Rx_t1*C' + Rw;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    Kt = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;                           % One-step prediction, \hat{y}_{t|t-1}.

    % If a sample is missing, just retain the earlier state.
    if isnan( y(t) )
        xt(:,t) = x_t1;                         % x_{t|t} = x_{t|t-1} 
        Rx_t    = Rx_t1;                        % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} 
        y1(t)   = yhat(t);                      % Replace the missing sample with the estimated value. 
    else
        h_et(t) = y(t)-yhat(t);                 % One-step prediction error, \hat{e}_t = y_t - \hat{y}_{t|t-1}
        xt(:,t) = x_t1 + Kt*( h_et(t) );        % x_{t|t} = x_{t|t-1} + K_t ( y_t -  \hat{y}_{t|t-1} ) 
        Rx_t    = Rx_t1 - Kt*Ry*Kt';            % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} - K_t R_{t|t-1}^{y,y} K_t^T
    end
    Rx_t1 = A*Rx_t*A' + Re;                     % R^{x,x}_{t+1|t} = A R^{x,x}_{t|t} A^T + Re
end


% Show the one-step prediction. We compare with y, not y1, as the former
% retains the missing samples (this was the reason to use the y1 vector);
% this is just to illustrate the missing samples in the figure.
figure
plot(tt(1:N), [y(1:N) yhat(1:N)] )
xlabel('Days')
if sum(isnan(y))
    hold on
    plot( tt(noVal), y0(noVal), 'b*')
    hold off
    legend('Realisation', 'Kalman estimate', 'Missing sample', 'Location','SW')
    title('One-step prediction using the Kalman filter with missing samples')
else 
    legend('Realisation', 'Kalman estimate', 'Location','SW')
    title('One-step prediction using the Kalman filter')
end


%% Examine one-step prediction residual.
ey = y0-yhat;
ey = ey(600:N);                                 % Ignore the initial values to let the filter converge first.
checkIfWhite( ey );

[~, pacfEst] = plotACFnPACF( ey, 40, 'One-step Kalman prediction' );
checkIfNormal( pacfEst(2:end), 'PACF' );
