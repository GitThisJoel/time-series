clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

%%
rng(0);
n = 500; % Number of samples
A3 = [1 .5];
C3 = [1 -.3 .2];
w = sqrt (2) * randn(n + 100, 1);
x = filter(C3, A3, w); % Create the input
A1 = [1 -.65];
A2 = [1 .90 .78];
C = 1;
B = [0 0 0 0 .4];
e = sqrt (1.5) * randn(n + 100, 1);
y = filter(C, A1, e) + filter(B, A2, x); % Create the output
x = x (101:end);
y = y (101:end); % Omit initial samples
clear A1 A2 C B e w A3 C3

%% Examine the data.

noLags = 32;

figure;
subplot(211);
plot(x);
ylabel('Input signal')
title('Measured signals')
subplot(212);
plot(y);
ylabel('Output signal')
xlabel('Time')

figure
[Cxy, lags] = xcorr(x, y, noLags, 'coeff');
stem(lags, Cxy)
hold on
condInt = 2 * ones(1, length(lags)) ./ sqrt(length(y));
plot(lags, condInt, 'r--')
plot(lags, -condInt, 'r--')
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between in- and output')

%% Construct model for input signal
ACFnPACFnNormplot(x, 32);

%% Test an AR1 first
ar_model_1 = estimateARMA(x, [1 1], [1], "AR1", 32)

%% Works ok, maybee remove something at lag 4
ar_model_4 = estimateARMA(x, [1 1 0 0 1], [1], "AR4", 32)

% the model choosen for x is an AR(1), since adding the 4:th component does
% not add significantly to them odel acuracy in our opinion.
% Data genererat med ARMA(1,2)

%%
arma_model = estimateARMA(x, [1 1], [1 1 1], "ARMA(1,2)", 32);

%% Pre-whiten y_t to create eps_t
w_t = filter(arma_model.A, arma_model.C, x);
w_t = w_t(length(arma_model.A):end);
eps_t = filter(arma_model.A, arma_model.C, y);
eps_t = eps_t(length(arma_model.A):end);

%%
M = 40; 
figure
stem(-M:M, crosscorr(w_t, eps_t, M));
title(' Cross\_correlation\_function'), xlabel (' Lag ')
hold on
plot(-M:M, 2 / sqrt (n) * ones(1, 2 * M + 1), '--')
plot(-M:M, -2 / sqrt (n) * ones(1, 2 * M + 1), '--')
hold off

%% Approximated model order (4,2,0), create model and find the residual of x
% försiktig med initialvärden här, sätter vi bara 1:or så riskerar vi att
% kanske inte få nån lösning
d = 4; r = 2; s = 0;
A2 = [1 zeros(1,r)];
B = [zeros(1, d) ones(1, s + 1)];
Mi = idpoly([1], [B], [], [], [A2]);
z = iddata(y, x);
Mba2 = pem(z, Mi); present(Mba2)
etilde = resid(Mba2, z);

%% Plot CCF of x and e_tilde 
CCF(x, etilde.OutputData, "(x vs. etilde)");

%% analyse e_tilde
ACFnPACFnNormplot(etilde.y,32)

% e_tilde is not white and it should not be since the delay is 4, so we are
% expecting an MA(3), which is roughly what we have.

%% Now use the model to estimate e_tilde
% This is equvalent to what is done when etilde is formed with the
% resid-comand

e_tilde_2 = y - filter(Mba2.B, Mba2.F, x)

ACFnPACFnNormplot(e_tilde_2 , 32)

%% Lets try to model e_tilde
etilde_ar_1 = estimateARMA(etilde.y,[1 1], [1], 'ar1\_etilde', 32)
present(etilde_ar_1)


% This model looks good :)

%% Create the complete model
A1 = [1 0];
A2 = A2;
B = B;
C = [1];
Mi = idpoly( 1 ,B,C,A1 ,A2 ) ;
z = iddata (y , x ) ;
MboxJ = pem( z ,Mi ) ;
present (MboxJ)
ehat = resid (MboxJ , z ) ;

%% Analyse the model residual

% all parameters in the model are significant

% looks white
ACFnPACFnNormplot(ehat.y,32);

% test if white, ok
checkIfNormal(ehat.y,'');
checkIfWhite(ehat.y);

% crosscorrelation between x and ehat, insignificant for low lags. some
% small significant componenets for relatively big lags, can probably be
% ignored.
CCF(x, ehat.y, 'x vs ehat')





%% 2.2 Hairdryer data
close all
clear all


ts = 0.08;              % Sampling time [s]
load ( 'tork.dat' )
tork = tork - repmat (mean(tork), length(tork), 1) ;
y = tork (:,1);
x = tork (:,2);
z = iddata (y,x) ;
plot (z(1:300))

%% Analyse input ata and construct model for it

ACFnPACFnNormplot(x, 32)

%% Test an AR(1)

model_input_AR1 = estimateARMA(x,[1 1], [1], "AR(1) model of input", 32)
present(model_input_AR1)

% This model seams to work well, both ACF and PACF looks white

%% Pre-whiten y and x with the abowe model
w_t = filter(model_input_AR1.A, model_input_AR1.C, x);
w_t = w_t(length(model_input_AR1.A):end);
eps_t = filter(model_input_AR1.A, model_input_AR1.C, y);
eps_t = eps_t(length(model_input_AR1.A):end);

%% Look at cross correlation of pre-whitened x and y
CCF(w_t,eps_t, "cross correlation of w\_t and eps\_t")

% From the table in the book, first guess:
d = 3;
r = 1;
s = 2;

%% Create filter for input
A2 = [1 zeros(1,r)];
B = [zeros(1, d) ones(1, s + 1)];
Mi = idpoly([1], [B], [], [], [A2]);
z = iddata(y, x);
Mba2 = pem(z, Mi); present(Mba2)
etilde = resid(Mba2, z);