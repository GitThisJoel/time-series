addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");

close all
clear

%% t1

close all

N = 1000;
s = 12;

A = [1 -0.5 0.7]; B = []; C = [1 zeros(1, 11) -0.7];
A12 = [1 zeros(1, s - 1) -1];
A_star = conv(A, A12);

e = randn(N, 1);
y = filter(C, A_star, e); y = y(50:end);

% plot(y)
% ACFnPACFnNormplot(y, 32);
% ACF shows strong seasoning.

y_s = filter(A12, 1, y); y_s = y_s(length(A12):end);
data = iddata(y_s);

A = [1 0 0]; B = []; C = [];
model_init = idpoly(A, B, C);
model_arma = pem(data, model_init);

% figure
% resid(model_arma, data);
% rarma = resid(model_arma, data);
% ACFnPACFnNormplot(rarma.OutputData, 32);
% strong season at 12

% add c_12 factor.

A = [1 0 0]; B = []; C = [1 zeros(1, 12)];
model_init = idpoly(A, B, C);
model_init.Structure.c.Free = [zeros(1, 12) 1];
model_arma = pem(data, model_init);

figure
resid(model_arma, data);
rarma = resid(model_arma, data);
ACFnPACFnNormplot(rarma.OutputData, 32);

% no more seasoning :)

%% t2
