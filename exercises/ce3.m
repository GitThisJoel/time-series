clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

%% 2.1 load data

load('tar2.dat')
load('thx.dat')

%% 2.1 plot tar2 and thx

close all
figure
subplot(211)
plot(tar2)
title("tar2")
subplot(212)
plot(thx)
title("thx")

%% 2.1 recusiveAR
close all
lambdas = [1 0.95 0.9];

plot_results = 1;

if plot_results
    figure
end

for i = 1:length(lambdas)
    X = recursiveAR(2);
    X.ForgettingFactor = lambdas(i);
    X.InitialA = [1 0 0];

    N = length(tar2);
    Aest = zeros(N, 3);
    yhat = zeros(1, N);

    for kk = 1:N
        [Aest(kk, :), yhat(kk)] = step(X, tar2(kk));
    end

    if plot_results
        subplot(1, length(lambdas), i)
        plot(Aest);
        title(sprintf("lambda = %f", lambdas(i)))
        % subplot(2, length(lambdas), i + length(lambdas))
        hold on
        plot(thx);
        legend("est1", "est2", "est3", "thx1", "thx2", "Location", "southeast")
    end

end

%% 2.1 choose lambda
close all

n = 100;
lambda_line = linspace(0.85, 1, n);
ls2 = zeros(n, 1);
yhat = zeros(n, 1);

for i = 1:length(lambda_line)
    reset(X);
    X.ForgettingFactor = lambda_line(i);
    X.InitialA = [1 0 0];

    for kk = 1:length(tar2)
        [~, yhat(kk)] = step(X, tar2(kk));
    end

    ls2(i) = sum((tar2 - yhat) .^ 2);
end

min_ind = find(ls2 <= min(ls2));
opt_lambda = lambda_line(min_ind);

plot(lambda_line, ls2);
title("lambda estimate (min = " + opt_lambda + ")")

%% 2.2 Kalman filtering
% Example of Kalman filter
% Simulate N samples of a process to test your code.

y = ? % Simulated data

% Define the state space equations.
A = [1 0; 0 0]; % ?? ar2
Re = [.004 0; 0 0]; % State covariance matrix
Rw = 1.25; %          Observation variance

% Set some initial values
Rxx1 = 10 * eye(2); % Initial state variance
xtt1 = [0 0]'; %      Initial statevalues

% Vectors to store values in
Xsave = zeros(2, N); % Stored states
ehat = zeros(1, N); %  Prediction residual
yt1 = zeros(1, N); %   One step prediction
yt2 = zeros(1, N); %   Two step prediction

% The filter use data up to time t −1 to predict value at t,
% then update using the prediction error. Why do we start
% from t = 3? Why stop at N−2?
for t = 3:N - 2
    Ct = [? ?]; % C { t | t-1}
    yhat(t) = ? % y { t | t-1}
    ehat(t) = y(t) - yhat(t); % et = yt - y{ t | t-1 }

    % Update
    Ryy = Rxx1 + eye(2); % Rˆ{ yy } { t | t-1 }
    Kt = ? %  Kt
    xtt = ? % x{ t | t }
    Rxx = ? % R{ xx } { t | t }

    % Predict the next state
    xtt1 = ? % x{ t+1 | t }
    Rxx1 = A * Rxx * A' + Re; % Rˆ{ xx }{ t+1 | t }

    % Form 2− step prediction. Ignore this part at first.
    Ct1 = [? ?]; %   C{ t+1 | t }
    yt1(t + 1) = ? % y{ t+1 | t } = C{ t+1 | t } x{ t | t }
    Ct2 = [? ?]; %   C{ t+2 | t }
    yt2(t + 2) = ? % y{ t+2 | t } = C{ t+2 | t } x{ t | t }

    % Store the state vector
    Xsave(:, t) = xtt;
end
