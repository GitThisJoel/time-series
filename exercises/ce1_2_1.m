clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);


%% 2.1

A1 = [1 -1.79 0.84];
C1 = [1 -0.18 -0.11];

% ARMA process A(z)y_t = C(z)e_t => y_t = C(z)/A(z)e_t
arma1 = idpoly(A1, [], C1);

rng(0)
sigma2 = 1.5;
N = 300;
e = sqrt(sigma2) * rand(N, 1);

y1 = filter(arma1.c, arma1.a, e);

y1 = y1(101:end);

y = iddata(y1);

na = 2;
nc = 2;
ar_model = arx(y, [na]);
arma_model = armax(y1, [na nc]);

present(ar_model)
present(arma_model)
% the ACF/PACF plot points to a AR(2) process.
% normplot = ?
ACFnPACFnNormplot(y1, 32);

e_hat = myFilter(ar_model.a, ar_model.c, y1,3);
e_hat = myFilter(arma_model.a, arma_model.c, y1,5);

figure
plot(e_hat)

[eacf, epacf]  = ACFnPACFnNormplot(e_hat, 32);

stem(eacf)

checkIfWhite(e_hat)

figure
resid(arma_model, e_hat)


% FPE??
% model est. ...

