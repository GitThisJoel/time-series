%% setting up data etc.
close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");
addpath("functions");

load('raw_data.mat')

halt_konc = raw_data(:, 1);
[modeling_set, validation_set, test_set, index_validation, index_test] = load_data(halt_konc);

halt_ing_rep = create_input_signal(raw_data);
[inp_modeling_set, inp_validation_set, inp_test_set, inp_index_validation, inp_index_test] = load_data(halt_ing_rep);

noLags = 50;

run_test = 1;

if run_test
    x = halt_ing_rep;
    y = halt_konc;
else
    x = [inp_modeling_set; inp_validation_set];
    y = [modeling_set; validation_set];
end

% models from part B
load('model_part_B_final.mat');
load('input_model.mat');

% %% Recursive model

%% input model
Ax_inp = input_model.A;
Cx_inp = input_model.C;

degA_inp = poly_degree(Ax_inp);
degC_inp = poly_degree(Cx_inp);

% non zeros elements
nzA_inp = zeros(size(Ax_inp)); nzA_inp(Ax_inp ~= 0) = 1;
nzC_inp = zeros(size(Cx_inp)); nzC_inp(Cx_inp ~= 0) = 1;

polys = [Ax_inp Cx_inp];

% removes 1 since we do not want the const. term in A
n_inp = length(find(polys)) - 1;

N = length(y);
A_inp = eye(n_inp);

ReA = 0.006;
ReC = 0.008;
Re_inp = diag([ReA * ones(1, sum(nzA_inp) - 1), 0 ReC * ones(1, sum(nzC_inp) - 1)]);

Rw_inp = 3;

% Initial input values
Rxx1_inp = 10 * eye(size(A_inp));
xtt1_inp = [Ax_inp(2:end) Cx_inp(find(Cx_inp))]'; % start with initial model parameters

% Vectors to store values in
ehat_inp = zeros(1, N);
xhat = zeros(1, N);
xt1 = zeros(1, N);
xtk = zeros(1, N);
Xsave_inp = zeros(n_inp, N);

%% output model
% Form the BJ prediction polynomials. In our notation, these are
%   A1 = foundModel.D
%   C1 = foundModel.C
%   A2 = foundModel.F
%
% The KA, KB, and KC polynomials are formed as:
%   KA = conv( A1, A2 );
%   KB = conv( A1, B );
%   KC = conv( A2, C1 );

Ax = conv(modelB.D, modelB.F);
Bx = conv(modelB.D, modelB.B);
Cx = conv(modelB.F, modelB.C);

% some poly properties
% degrees of the polynomials
degA = poly_degree(Ax);
degB = poly_degree(Bx);
degC = poly_degree(Cx);

degmax = max([degA degB degC degA_inp degC_inp]);

% non zeros elements
nzA = zeros(size(Ax)); nzA(Ax ~= 0) = 1;
nzB = zeros(size(Bx)); nzB(Bx ~= 0) = 1;
nzC = zeros(size(Cx)); nzC(Cx ~= 0) = 1;

indC = sum(nzA) + 1;
indB = indC + sum(nzC) + 1;

polys = [Ax Bx Cx];

% removes 1 since we do not want the const. term in A
n = length(find(polys)) - 1;

N = length(y);
A = eye(n);

ReA = 0.00006;
ReB = 6e-5; % where does B begin, should there be a value?
ReC = 0.002; % where does C begin, should there be a value?
Re = diag([ReA * ones(1, sum(nzA) - 1), 0 ReC * ones(1, sum(nzC) - 1), ReB * ones(1, sum(nzB))]);

Rw = 3;

% Initial values
Rxx1 = 10 * eye(size(A));
xtt1 = [Ax(2:end) Cx Bx(find(Bx))]';

% Vectors to store values in
Xsave = zeros(n, N);
ehat = zeros(1, N);
yhat = zeros(1, N);
yt1 = zeros(1, N);
ytk = zeros(1, N);

k = 9; % k-step predictor

%% kalman
for t = degmax + 1:N - k
    % constructing the input prediction
    Ct_inp = create_ct_inp(t, x, ehat_inp, Ax_inp, Cx_inp, 0);

    xhat(t) = Ct_inp * xtt1_inp;
    ehat_inp(t) = x(t) - xhat(t);

    Ryy = Ct_inp * Rxx1_inp * Ct_inp' + Rw_inp;
    Kt = Rxx1_inp * Ct_inp' / Ryy;
    xtt_inp = xtt1_inp + Kt * (ehat_inp(t));
    Rxx = Rxx1_inp - Kt * Ryy * Kt';

    % Predict the next state
    xtt1_inp = A_inp * xtt_inp;
    Rxx1_inp = A_inp * Rxx * A_inp' + Re_inp;

    % Form 1-step and k-step prediction.

    xt_inp = x(t - degA_inp + 1:t)';
    C1_inp = create_ct_inp(t, xt_inp, ehat_inp, Ax_inp, Cx_inp, 1);

    xk = C1_inp * xtt_inp;
    xt1(t + 1) = xk;

    xt = x(t - degB + 1:t + 1)';
    xt = [xt(2:end) xt1(t)];

    Xsave_inp(:, t) = xtt_inp;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % output prediction
    % construct the Ct vector.
    Ct = create_ct(t, y, x, ehat, Ax, Bx, Cx, 0);

    yhat(t) = Ct * xtt1;
    ehat(t) = y(t) - yhat(t);

    % Update
    Ryy = Ct * Rxx1 * Ct' + Rw;
    Kt = Rxx1 * Ct' / Ryy;
    xtt = xtt1 + Kt * (ehat(t));
    Rxx = Rxx1 - Kt * Ryy * Kt';

    % Predict the next state
    xtt1 = A * xtt;
    Rxx1 = A * Rxx * A' + Re;

    % Form 1-step and k-step prediction.
    C1 = create_ct(t, y, xt, ehat, Ax, Bx, Cx, 1);

    yk = C1 * xtt;
    yt1(t + 1) = yk;

    yt = y(t - degA + 1:t + 1)';
    yt = [yt(2:end) yk];

    for k0 = 2:k
        % predict input
        Ck_inp = create_ct_inp(t, xt_inp, ehat_inp, Ax_inp, Cx_inp, k0);

        xk = Ck_inp * A_inp ^ k0 * xtt_inp;
        xt = [xt(2:end) xk];

        % predict output
        Cky = flip(yt(find(flip(nzA(2:end))')));
        Ck = create_ct(t, yt, xt, ehat, Ax, Bx, Cx, k0, Cky);

        yk = Ck * A ^ k0 * xtt;
        yt = [yt(2:end) yk];
    end

    xtk(t + k) = xk;
    ytk(t + k) = yk;

    % Store the state vector
    Xsave(:, t) = xtt;
end

%% evalute the prediction
% we should probably update our evaluation here, something more is needed.
close all

% input prediction
figure
plot(x)
hold on
plot(xt1)
plot(xtk)
legend('x', 'k = 1', sprintf("k = %d", k))
title('input prediction')

figure
plot(Xsave_inp')
legend('1', '2', '3', '4', '5', '6')
title('parameters for input prediction')

% output prediction
additional_plot_start = 500; % max 2000, min 0

figure
plot(y(end - additional_plot_start - 100 - k:end - k))
hold on
plot(yt1(end - additional_plot_start - 100 - k + 1:end - k + 1), 'g')
plot(ytk(end - additional_plot_start - 100:end), 'r')
hold off
legend('y', 'k = 1', sprintf("k = %d", k))
title('Output prediction')

plot_hist_param = 0;

if plot_hist_param
    figure
    plot(y)
    hold on
    plot(yt1)
    plot(ytk)
    legend('x', 'k = 1', sprintf("k = %d", k))
    title('input prediction')
end

% err_resid = norm(ehat(end - 200:end)) .^ 2;
% disp("sum pred residuals = " + err_resid)

figure
plot(Xsave')
legend('1', '2', '3', '4', '5', '6', '7', '8', '9')
title('parameters for output prediction')

% one step prediction
ehat_inp_1 = xt1 - x;
ehat_1 = yt1 - y;

[varx_val, varx_test] = evaluate_performance(ehat_inp_1, index_validation, index_test);
[vary_val, vary_test] = evaluate_performance(ehat_1, index_validation, index_test);

fprintf('one step prediction ')
fprintf('Validation set: Variance of y is %s, variance of x is %d \n', vary_val, varx_val)
fprintf('Test set: Variance of y is %s, variance of x is %d \n', vary_test, varx_test)

%% nine step prediction

ehat_inp_9 = xtk - x;
ehat_9 = ytk - y;
[varx_val, varx_test] = evaluate_performance(ehat_inp_9, index_validation, index_test);
[vary_val, vary_test] = evaluate_performance(ehat_9, index_validation, index_test);

fprintf('Validation set: Variance of y is %s, variance of x is %d \n', vary_val, varx_val)
fprintf('Test set: Variance of y is %s, variance of x is %d \n', vary_test, varx_test)

%% Figures to report

figure
plot(y(index_validation(1) - 1:index_validation(2) - 1))
hold on
plot(yt1(index_validation(1):index_validation(2)))
title("1-step output prediction")
legend('Real value', 'Prediction')

figure
plot(y(index_validation(1) - 9:index_validation(2) - 9))
hold on
plot(ytk(index_validation(1):index_validation(2)))
title("9-step output prediction")
legend('Real value', 'Prediction')

figure
plot(x(index_validation(1) - 1:index_validation(2) - 1))
hold on
plot(xt1(index_validation(1):index_validation(2)))
title("1-step input prediction")
legend('Real value', 'Prediction')

figure
plot(x(index_validation(1) - 9:index_validation(2) - 9))
hold on
plot(xtk(index_validation(1):index_validation(2)))
title("9-step input prediction")
legend('Real value', 'Prediction')

%% parametrar for validation
% byt notation i legenderna timm typ a_1, B_1 osv

figure
plot(Xsave(:, index_validation(1):index_validation(2))')
legend('1', '2', '3', '4', '5', '6', '7', '8', '9')
title('parameters for output prediction')
ylim([-1.2 1.2])

figure
plot(Xsave_inp(:, index_validation(1):index_validation(2))')
legend('1', '2', '3', '4', '5', '6')
title('parameters for input prediction')
ylim([-1.2 1.2])

%% ACF av residualerna, är det ens relevant för kalman?!
figure
acf(ehat_1(index_validation(1):index_validation(2)), 50, 0.05, 1);
figure
acf(ehat_9(index_validation(1):index_validation(2)), 50, 0.05, 1);

figure
acf(ehat_inp_1(index_validation(1):index_validation(2)), 50, 0.05, 1);
figure
acf(ehat_inp_9(index_validation(1):index_validation(2)), 50, 0.05, 1);
