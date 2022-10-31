%
% Time series analysis
% Lund University
%
% Example code 23: Form a k-step prediction of an ARMA process with unknown
% parameters using the Kalman filter (see also example 8.12).
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


%% Estimate the unknown parameters using a Kalman filter and form the k-step prediction.
k  = 5;                                         % k-step prediction.
p0 = 1;                                         % Number of unknowns in the A polynomial.
q0 = 2;                                         % Number of unknowns in the C polynomial.

A     = eye(p0+q0);
Rw    = 1;                                      % Measurement noise covariance matrix, R_w
Re    = 1e-6;                                   % System noise covariance matrix, R_e
Rx_t1 = eye(p0+q0);                             % Initial covariance matrix, R_{1|0}^{x,x}
h_et  = zeros(N,1);                             % Estimated one-step prediction error.
xt    = zeros(p0+q0,N);                         % Estimated states. Intial state, x_{1|0} = 0.
yhat  = zeros(N,1);                             % Estimated output.
yhatk = zeros(N,1);                             % Estimated k-step prediction.
for t=4:N-k                                     % We use t-3, so start at t=4. As we form a k-step prediction, end the loop at N-k.
    % Update the predicted state and the time-varying state vector.
    x_t1 = A*xt(:,t-1);                         % x_{t|t-1} = A x_{t-1|t-1}
    C    = [ -y(t-1) h_et(t-1) h_et(t-3) ];     % C_{t|t-1}
    
    % Update the parameter estimates.
    Ry = C*Rx_t1*C' + Rw;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    Kt = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;                           % One-step prediction, \hat{y}_{t|t-1}.
    h_et(t) = y(t)-yhat(t);                     % One-step prediction error, \hat{e}_t = y_t - \hat{y}_{t|t-1}
    xt(:,t) = x_t1 + Kt*( h_et(t) );            % x_{t|t}= x_{t|t-1} + K_t ( y_t - Cx_{t|t-1} ) 

    % Update the covariance matrix estimates.
    Rx_t  = Rx_t1 - Kt*Ry*Kt';                  % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} - K_t R_{t|t-1}^{y,y} K_t^T
    Rx_t1 = A*Rx_t*A' + Re;                     % R^{x,x}_{t+1|t} = A R^{x,x}_{t|t} A^T + Re

    % Form the k-step prediction by first constructing the future C vector
    % and the one-step prediction. Note that this is not yhat(t) above, as
    % this is \hat{y}_{t|t-1}.
    Ck = [ -y(t) h_et(t) h_et(t-2) ];           % C_{t+1|t}
    yk = Ck*xt(:,t);                            % \hat{y}_{t+1|t} = C_{t+1|t} A x_{t|t}

    % Note that the k-step predictions is formed using the k-1, k-2, ...
    % predictions, with the predicted future noises being set to zero. If
    % the ARMA has a higher order AR part, one needs to keep track of each
    % of the earlier predicted values.
    for k0=2:k
        Ck = [ -yk h_et(t+k0-1) h_et(t+k0-3) ]; % C_{t+k|t}
        yk = Ck*A^k*xt(:,t);                    % \hat{y}_{t+k|t} = C_{t+k|t} A^k x_{t|t}
    end
    yhatk(t+k) = yk;                            % Note that this should be stored at t+k.
end


%% Examine the estimated parameters.
figure
Q0 = [A0(2) C0(2) C0(4)];                       % These are the true parameters we seek.
plot( xt' )
for k0=1:length(Q0)
    line([1 N], [Q0(k0) Q0(k0)], 'Color','red','LineStyle',':')
end
title(sprintf('Estimated parameters, with Re = %7.6f and Rw = %4.3f', Re, Rw))
xlabel('Time')
ylim([-1.5 1.5])
xlim([1 N-k])


%% Show the one-step prediction. 
figure
plot( [y yhat] )
title('One-step prediction using the Kalman filter')
xlabel('Time')
legend('Realisation', 'Kalman estimate', 'Location','SW')
xlim([1 N-k])


%% Show the k-step prediction. 
figure
plot( [y(1:N-k) yhatk(k+1:N)] )
title( sprintf('%i-step prediction using the Kalman filter (shifted %i steps)', k, k) )
xlabel('Time')
legend('Realisation', 'Kalman estimate', 'Location','SW')
xlim([1 N-k])


%% Examine k-step prediction residual.
ey = y-yhatk;
ey = ey(N-200:N-k);                             % Ignore the initial values to let the filter converge first.
plotACFnPACF( ey, 40, sprintf('%i-step prediction using the Kalman filter', k)  );
