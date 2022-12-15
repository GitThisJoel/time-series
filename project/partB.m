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
plot_crosscorr(x, y, noLags);

%% pre white
close all

input_model = estimateARMA(x, [1 1 1 1], [1 0 1 1 zeros(1, 7) 1], "ARMA(3,11) of input", noLags);
% input_model = estimateARMA(x, [1 1 1 1], [1 0 1 1], "ARMA(3,3) of input", noLags);
present(input_model)
% ACFnPACFnNormplot(ex, noLags);

S = 13; % from looking at estimated arma of x
AS = [1 zeros(1, S - 1) -1];
xdiff = myFilter(AS, 1, x, 2);
ACFnPACFnNormplot(xdiff, noLags);
input_model = estimateARMA(x, [1 1 1 1], [1 0 1 1], "ARMA(3,3) of input", noLags);
present(input_model)
