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

%% plot the (modelling) data
close all

figure
subplot(211)
plot(x);
title("Input signal")
subplot(212)
plot(y);
title("Output signal")

ACFnPACFnNormplot(x, noLags);

%% cross correlation between input and output data

close all
% plot_crosscorr(x, y, noLags);
CCF(x, y, "", noLags);


%% pre white

close all

%[input_model, input_model_res ] = estimateARMA(x, [1 1 1 1], [1 0 1 1 zeros(1, 7) 1], "ARMA(3,11) of input", noLags);
[input_model, input_model_res ] = estimateARMA(x, [1 1 1 1], [1 0 1 1], "ARMA(3,3) of input", noLags);
present(input_model)
checkIfWhite(input_model_res);

% tried diff, did not help.
[w_t, eps_t] = pre_white(input_model, x, y);
cc = CCF(w_t, eps_t, "", noLags);

%% choose model order
close all

r = 1;
[~, d] = max(cc); d = d - noLags - 1;
s = 1;

[Mba2, etilde] = create_input_model(d, r, s, x, y);


ACFnPACFnNormplot(etilde.OutputData,noLags)
na = 1;
nc = 3;
% etilde_model = estimateARMA(etilde.y, [1 ones(1, na)], [1 ones(1, nc)], "ar1\_etilde", noLags);
etilde_model = estimateARMA(etilde.y, [1 1], [1 0 1 1], "ar1\_etilde", noLags);
present(etilde_model)

 %% entire model
A1 = [1 zeros(1, etilde_model.na)];
A2 = [1 zeros(1, Mba2.nf)];
B = [zeros(1, Mba2.nk) ones(1, Mba2.nb)];
C = [1 zeros(1, etilde_model.nc)]; % ones or zeros here?
Mi = idpoly(1, B, C, A1, A2);
z = iddata (y, x);
MboxJ = pem(z, Mi);
present(MboxJ)
ehat = resid(MboxJ, z);

ACFnPACFnNormplot(ehat.y, 32);

checkIfNormal(ehat.y, ''); close
checkIfWhite(ehat.y);

CCF(x, ehat.y, 'x vs ehat');

present(MboxJ)
%% Do predictions with model:
k = 1;  % prediction horizon
% The data we do the predictions on:
% x = [x ; inp_validation_set];
% y = [y ; validation_set];

% first predict the input:
[Fx, Gx] = polydiv( input_model.C, input_model.A, k );
xhatk = filter(Gx, input_model.C, x);

figure
plot(x)
hold on
plot(xhatk)


ehat = x - xhatk;
ehat = ehat(10:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step input prediction residual', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );

%% predict output

% Form the BJ prediction polynomials. In our notation, these are
%   A1 = foundModel.D
%   C1 = foundModel.C
%   A2 = foundModel.F
% 
% The KA, KB, and KC polynomials are formed as:
%   KA = conv( A1, A2 );
%   KB = conv( A1, B );
%   KC = conv( A2, C1 );
%

foundModel = MboxJ;

KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

% Form the ARMA prediction for y_t (note that this is not the same G
% polynomial as we computed above (that was for x_t, this is for y_t).
[Fy, Gy] = polydiv( foundModel.C, foundModel.D, k );

% Compute the \hat\hat{F} and \hat\hat{G} polynomials.
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

% Form the predicted output signal using the predicted input signal.
yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, x) + filter(Gy, KC, y);

%%
figure
plot(y)
hold on
plot(yhatk)


ehat = y - yhatk;
ehat = ehat(10:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step output prediction residual', k) )
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );
