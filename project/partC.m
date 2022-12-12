close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");

load('raw_data.mat')
[modeling_set, validation_set, test_set] = load_data(raw_data)

%%
