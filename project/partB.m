close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");

load('raw_data.mat')

halt_konc = raw_data(:, 1);
[modeling_set, validation_set, test_set] = load_data(raw_data)

halt_ing_rep = (halt_raw_konc .* flow_raw_remill + halt_scav_konc .* flow_scav_remill) ./ flow_ing_rep;
%%
