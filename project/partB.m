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

% input_model = estimateARMA(x, [1 1 1 1], [1 0 1 1 zeros(1, 7) 1], "ARMA(3,11) of input", noLags);
input_model = estimateARMA(x, [1 1 1 1], [1 0 1 1], "ARMA(3,3) of input", noLags);
present(input_model)

% tried diff, did not help.
[w_t, eps_t] = pre_white(input_model, x, y);
cc = CCF(w_t, eps_t, "", noLags);

%% choose model order
close all

r = 2;
[~, d] = max(cc); d = d - noLags - 1;
s = 1;

[Mba2, etilde] = create_input_model(d, r, s, x, y);

na = 1;
nc = 3;
% etilde_model = estimateARMA(etilde.y, [1 ones(1, na)], [1 ones(1, nc)], "ar1\_etilde", noLags);
etilde_model = estimateARMA(etilde.y, [1 1], [1 0 0 1], "ar1\_etilde", noLags);
present(etilde_model)

% %% entire model
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
