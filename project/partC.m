close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");
addpath("functions");

load('raw_data.mat')

halt_konc = raw_data(:, 1);
[modeling_set, validation_set, test_set] = load_data(halt_konc);

halt_ing_rep = create_input_signal(raw_data);
[inp_modeling_set, inp_validation_set, inp_test_set] = load_data(halt_ing_rep);

noLags = 50;

y = modeling_set;
x = inp_modeling_set;

% model from part B
load('model_part_B.mat')

%% Recursive model
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
degmax = max([degA degB degC]);

% non zeros elements
nzA = zeros(size(Ax)); nzA(Ax ~= 0) = 1;
nzB = zeros(size(Bx)); nzB(Bx ~= 0) = 1;
nzC = zeros(size(Cx)); nzC(Cx ~= 0) = 1;

polys = [Ax Bx Cx];
n = sum(nonzeros(polys) ~= 1);

N = length(y);
A = eye(n);
Re = zeros(size(A)); % think about this one...
Re(1) = 0.04;
% Re(...) = ... % where does B begin, should there be a value?
% Re(...) = ... % where does C begin, should there be a value?
Rw = 1.25;

% Initial values
Rxx1 = 10 * eye(size(A));
xtt1 = zeros(size(A, 1), 1);

% Vectors to store values in
Xsave = zeros(n, N);
ehat = zeros(1, N);
yhat = zeros(1, N);
yt1 = zeros(1, N);
yt2 = zeros(1, N);

k = 1; % k-step predictor

%% kalman
for t = degmax + 1:N - k
    Cty = (y(t - degA:t - 1, :) .* flip(nzA(2:end))');
    Ctx = (y(t - degB:t, :) .* flip(nzB)'); % b0 != 1
    Cte = (y(t - degC:t - 1, :) .* flip(nzC(2:end))');
    Ct = nonzeros([Cty' Cte' Ctx'])';

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

    % Store the state vector
    Xsave(:, t) = xtt;
end
