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

t = [Ax Bx Cx];

N = length(y);
A = eye(sum(nonzeros(t) ~= 1));
Re = zeros(size(A)); % think about this one...
Re(1) = 0.04;
% Re(...) = ... % where does B begin, should there be a value?
% Re(...) = ... % where does C begin, should there be a value?
Rw = 1.25;

% Initial valuess
Rxx1 = 10 * eye(size(A));
xtt1 = zeros(size(A, 1));
