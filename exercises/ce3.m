clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

% rng(0);

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

%% 2.1 opt lambda
close all

n = 100;
lambda_line = linspace(0.85, 1, n);
ls2 = zeros(n, 1);
yhat = zeros(n, 1);

X = recursiveAR(2);

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

X = recursiveAR(2);
X.ForgettingFactor = opt_lambda;
X.InitialA = [1 0 0];

N = length(tar2);
Aest = zeros(N, 3);
yhat = zeros(1, N);

for kk = 1:N
    [Aest(kk, :), yhat(kk)] = step(X, tar2(kk));
end

figure
plot(Aest(:, 2:3));
title(sprintf("lambda = %f", opt_lambda))
hold on
plot(thx);
legend("est1", "est2", "thx1", "thx2", "Location", "southeast")

disp("Estimated coeffs = " + Aest(end, 2) + "   " + Aest(end, 3));
% Aest(end) =  1.0000    1.4283    0.6737

%% 2.2 Kalman filtering
close all

N = length(tar2);
y = tar2;

% Define the state space equations.
A = [1 0; 0 1]; % ?? ar2
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

% The filter use data up to time t−1 to predict value at t,
% then update using the prediction error. Why do we start
for t = 3:N
    Ct = [-y(t - 1) -y(t - 2)]; % C{ t | t-1}
    yhat(t) = Ct * xtt1; % y{ t | t-1}
    ehat(t) = y(t) - yhat(t); % et = yt - y{ t | t-1 }

    % Update
    Ryy = Ct * Rxx1 * Ct' + Rw; % Rˆ{ yy } { t | t-1 }
    Kt = Rxx1 * Ct' / Ryy; %  Kt
    xtt = xtt1 + Kt * (ehat(t)); % x{ t | t }
    Rxx = Rxx1 - Kt * Ryy * Kt'; % R{ xx } { t | t }

    % Predict the next state
    xtt1 = A * xtt; % x{ t+1 | t }
    Rxx1 = A * Rxx * A' + Re; % Rˆ{ xx }{ t+1 | t }

    % Store the state vector
    Xsave(:, t) = xtt;
end

figure
plot(Xsave(:, 3:end)')
hold on
plot(thx)

disp("sum of square residuals = " + norm(ehat) .^ 2);

%% 2.3 Kalman filtering
% Example of Kalman filter
% Simulate N samples of a process to test your code.

close all
clear

rng(0);

% simulate data
extraN = 100;
N = 10000;
A0 = [1 -0.8 0.2];
ee = 0.1 * randn(N + extraN, 1);

y = filter(1, A0, ee);
y = y(extraN + 1:end);
ee = ee(extraN + 1:end);

% Define the state space equations.
A = [1 0; 0 1]; % ?? ar2
Re = [1e-6 0; 0 1e-6]; % State covariance matrix
Rw = 0.1; %          Observation variance

% Set some initial values
Rxx1 = 10 * eye(2); % Initial state variance
xtt1 = [0 0]'; %      Initial statevalues

% Vectors to store values in
Xsave = zeros(2, N); % Stored states
ehat = zeros(1, N); %  Prediction residual
yt1 = zeros(1, N); %   One step prediction
yt2 = zeros(1, N); %   Two step prediction

% The filter use data up to time t−1 to predict value at t,
% then update using the prediction error. Why do we start
% from t = 3? Why stop at N−2?
for t = 3:N - 2
    Ct = [-y(t - 1) -y(t - 2)]; % C{ t | t-1}
    yhat(t) = Ct * xtt1; % y{ t | t-1}
    ehat(t) = y(t) - yhat(t); % et = yt - y{ t | t-1 }

    % Update
    Ryy = Ct * Rxx1 * Ct' + Rw; % Rˆ{ yy }{ t | t-1 }
    Kt = Rxx1 * Ct' / Ryy; %  Kt
    xtt = xtt1 + Kt * (ehat(t)); % x{ t | t }
    Rxx = Rxx1 - Kt * Ryy * Kt'; % R{ xx }{ t | t }

    % Predict the next state
    xtt1 = A * xtt; % x{ t+1 | t }
    Rxx1 = A * Rxx * A' + Re; % Rˆ{ xx }{ t+1 | t }

    % Form 2−step prediction. Ignore this part at first.
    Ct1 = [-y(t) -y(t - 1)]; %   C{ t+1 | t }
    yt1(t + 1) = Ct1 * xtt; % y{ t+1 | t } = C{ t+1 | t } x{ t | t }
    Ct2 = [-yt1(t + 1) -y(t)]; %   C{ t+2 | t }
    yt2(t + 2) = Ct2 * xtt; % y{ t+2 | t } = C{ t+2 | t } x{ t | t }

    % Store the state vector
    Xsave(:, t) = xtt;
end

figure
plot(y(end - 100 - 2:end - 2))
hold on
plot(yt1(end - 100 - 1:end - 1), 'g')
plot(yt2(end - 100:end), 'r')
hold off
legend('y', 'k = 1', 'k = 2')

err_resid = norm(ehat(end - 200:end)) .^ 2
disp("sum pred residuals = " + err_resid)

%% 2.4 quality control
clear
close all

N = 500;

b = 20;
vare = 1;
varv = 4;

P = [7/8 1/8; 1/8 7/8];
ut = simulate_ut(P, N);
et = sqrt(vare) * randn(1, N);
vt = sqrt(varv) * randn(1, N);

xt = zeros(1, N);
yt = zeros(1, N);

y(1) = xt(1) + b * ut(1) + vt(1);

for t = 2:N
    xt(t) = xt(t - 1) + et(t);
    yt(t) = xt(t) + b * ut(t) + vt(t);
end

% Define the state space equations.
A = eye(2); % ?? ar2
Re = [vare 0; 0 0]; % State covariance matrix
Rw = varv; %          Observation variance

% Set some initial values
Rxx1 = 20 * eye(2); % Initial state variance
xtt1 = [0 0]'; %      Initial statevalues

[xtt, ehat, Xsave] = kalmanfilter(yt, A, eye(2), N, Rw, Re, xtt1, Rxx1, ut);

figure
plot(Xsave(:, 3:end - 2)')
hold on
plot(xt(3:end - 2))
yline(b)
legend("est x", "est b", "real x", "real b")

%% 2.5
clear
close all

load svedala94.mat

T = linspace(datenum(1994, 1, 1), datenum(1994, 12, 31), length(svedala94));

S = 6;
AS = [1 zeros(1, S - 1) -1];
ydiff = myFilter(AS, 1, svedala94, 2);

figure
plot(T, svedala94)
datetick('x')
hold on
plot(T(2:end), ydiff)

th = armax(ydiff, [2 2]);
th_winter = armax(ydiff(1:540), [2 2]);
th_summer = armax(ydiff(907:1458), [2 2]);
%
% figure
% plot(th.A, '*')
% hold on
% plot(th_summer.A, '*')
% plot(th_winter.A, '*')
% legend('th', 'th\_summer', 'th\_winter')
%
% figure
% plot(th.C, '*')
% hold on
% plot(th_summer.C, '.')
% plot(th_winter.C, '.')
% legend('th', 'th\_summer', 'th\_winter')

figure
subplot(1, 3, 1)
pzplot(th)
subplot(1, 3, 2)
pzplot(th_summer)
subplot(1, 3, 3)
pzplot(th_winter)

% part 2.5.3
X = recursiveARMA([2 2]);
X.InitialA = [1 th_winter.A(2:end)];
X.InitialC = [1 th_winter.C(2:end)];
X.ForgettingFactor = 0.99;

for k = 1:length (ydiff)
    [Aest(k, :), Cest(k, :), yhat(k)] = step(X, ydiff(k));
end

figure
subplot(3, 1, 1)
plot (T, svedala94);
datetick('x')
subplot(3, 1, 2)
plot (Aest (:, 2:3))
hold on
plot(repmat(th_winter.A(2:end), [length(ydiff) 1]), 'g:');
plot(repmat(th_summer.A(2:end), [length(ydiff) 1]), 'r:');
axis tight
hold off
subplot(3, 1, 3)
plot (Cest (:, 2:3))
hold on
plot(repmat(th_winter.C(2:end), [length(ydiff) 1]), 'g:');
plot(repmat(th_summer.C(2:end), [length(ydiff) 1]), 'r:');
axis tight
hold off

% The A- coeficients of the recursive process correspond well to the
% non-recursive process during winter or during summer, while spring and
% autum is transfer periods between the non-recursive processes. For the
% C-coefficients, and especially the c1-coefficient, the difference to the
% summerprocess is bigger for the recursive process conmpared to the
% non-rcursive.

%% 2.6
clear
close all

load svedala94.mat

%y = svedala94(850:1100);
%y = y - mean(y);
y = svedala94;
y = y -y(1);

t = (1:length(y))';
U = [sin(2 * pi * t / 6) cos(2 * pi * t / 6)];
Z = iddata(y, U);
model = [3 [1 1] 4 [0 0]]; %[ na [ nb_1 nb_2 ] nc [ nk_1 nk_2 ] ] ;

thx = armax(Z, model);

figure
plot(y)
hold on
plot(U * cell2mat(thx.b)')
legend('y', 'Seasonal function')

% Looks like the period matches well, might be more flexible if it is done
% recursivly

U = [sin(2 * pi * t / 6) cos(2 * pi * t / 6) ones(size(t))];
Z = iddata(y, U);
m0 = [thx.A(2:end) cell2mat(thx.B) 0 thx.C(2:end)];
Re = diag ([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re, m0);

m = thr(:, 6);
a = thr(end, 4);
b = thr(end, 5);
y_mean = m + a * U(:, 1) + b * U(:, 2);
y_mean = [0; y_mean(1:end - 1)];

figure
plot(y)
hold on
plot(y_mean)
legend('y', 'y\_mean')

mse = rms(y - y_mean)

% För hela datasettet diffar det en del i början men följer rätt fint sen
