close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");
addpath("functions");

load('raw_data.mat')

halt_konc = raw_data(:, 1);
% halt_konc = halt_konc - mean(halt_konc);
[modeling_set, validation_set, test_set] = load_data(halt_konc);

halt_ing_rep = create_input_signal(raw_data);
[inp_modeling_set, inp_validation_set, inp_test_set] = load_data(halt_ing_rep);

noLags = 32;

%% plot the (modelling) data
close all

y = modeling_set;
x = inp_modeling_set;

figure
subplot(211)
plot(x);
title("Input signal")
subplot(212)
plot(y);
title("Output signal")
