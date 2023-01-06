%% setting up data etc.
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

y = [modeling_set; validation_set];
x = [inp_modeling_set; inp_validation_set];

% model from part B
load('model_part_B_final.mat')

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

indB = sum(nzA) ;
indC = indB + sum(nzB);

polys = [Ax Bx Cx];
n = sum(nonzeros(polys) ~= 0)-1;

N = length(y);
A = eye(n);
ReA = 0.00001;          % delem som beror av tidigare värden på y
ReB = 0.0001;           % delen som beror av insignalen
ReC = 0.00001;          % delen som beror av tidigare fel
Re = diag([ReA, ReB*ones(1,4), 0 ReC*ones(1,3)]); 
Rw = 0.1;               % "mätbrus"

% Initial values
Rxx1 = 10 * eye(size(A));
xtt1 = [Ax(2) Bx(4:end) Cx]'; % start with initial model parameters

% Vectors to store values in
Xsave = zeros(n, N);
ehat = zeros(1, N);
yhat = zeros(1, N);
yt1 = zeros(1, N);
ytk = zeros(1, N);

k = 1; % k-step predictor

%% kalman
% Frida, kolla tidsförskjutningen, är C skiftad ett steg för långt åt vänster?
for t = degmax + 1:N - k
    % construct the Ct vector.
    % start by cheeting and using the real input instead of the prediction
    %Ct = [-yhat(t+8) x(t+6) x(t+5) x(t+4) x(t+3) 0 0 0 0];
%    Ct = [-yhat(t+8) xhat(t+6) xhat(t+5) xhat(t+4) xhat(t+3) 0 0 0 0];
    %Skapa och skatta xhat? hur? eget kalmanfilter?- tror det....
    Ct = [-y(t-1) x(t-3) x(t-4) x(t-5) x(t-6) 0 ehat(t-1) ehat(t-2) ehat(t-3)];
    yhat(t) = Ct * xtt1; % y{ t | t-1}
    ehat(t) = y(t) - yhat(t); % et = yt - y{ t | t-1 }
    
if t == 12
    Ct1=Ct
end

        % Update
        Ryy = Ct * Rxx1 * Ct' + Rw; % Rˆ{ yy }{ t | t-1 }
        Kt = Rxx1 * Ct' / Ryy; %  Kt
        xtt = xtt1 + Kt * (ehat(t)); % x{ t | t }
        Rxx = Rxx1 - Kt * Ryy * Kt'; % R{ xx }{ t | t }
    
        % Predict the next state
        xtt1 = A * xtt; % x{ t+1 | t }
        Rxx1 = A * Rxx * A' + Re; % Rˆ{ xx }{ t+1 | t }
    
        % Store the state vector
        Xsave(:, t) = xtt;

end

%% evalute the prediction
% we should probably update our evaluation here, something more is needed.
close all

figure

if k > 0
    plot(y) %(end - 100 - k:end - k))
    hold on
    plot(yhat) %(end - 100 - k + 1:end - k + 1), 'g')
    hold off
    legend('y', 'k = 1')
end

err_resid = norm(ehat(end - 200:end)) .^ 2;
disp("sum pred residuals = " + err_resid)

figure
plot(Xsave(:, 1:end)')