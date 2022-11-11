clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

%% 2.1

A1 = [1 -1.79 0.84];
C1 = [1 -0.18 -0.11];

A2 = [1 -1.79];
C2 = [1 -0.18 -0.11];

% ARMA process A(z)y_t = C(z)e_t => y_t = C(z)/A(z)e_t
arma1 = idpoly(A1, [], C1);
arma2 = idpoly(A2, [], C2);

rng(0)
sigma2 = 1.5;
N = 300;
e = sqrt(sigma2) * rand(N, 1);

y1 = filter(arma1.c, arma1.a, e);
y2 = filter(arma2.c, arma2.a, e);

y1 = y1(101:end);
y2 = y2(101:end);

close all
subplot(221)
pzmap(arma1);
title("pole zero - arma1")
subplot(222)
pzmap(arma2);
title("pole zero - arma2")

subplot(223)
plot(y1)
title("filtered signal y1")
subplot(224)
plot(y2)
title("filtered signal y2")

% arma2 is unsable since one pole is outside of the unit circle.
% why...

m = 20;
r_theo = kovarians(arma1.c, arma1.a, m);
figure
stem(0:m, r_theo * sigma2)
hold on
r_est = covf(y1, m + 1);
stem(0:m, r_est, 'r')

%% 2.1 cont.
close all

y = iddata(y1);

na = 2;
nc = 2;
ar_model = arx(y, [na]);
arma_model = armax(y, [na nc]);

present(ar_model)
present(arma_model)
% the ACF/PACF plot points to a AR(2) process.
% normplot = ?
ACFnPACFnNormplot(y1, 32);

% e_hat = filter(ar_model.a, ar_model.c, y1);
e_hat = filter(arma_model.a, arma_model.c, y1);

figure
subplot(121)
plot(e_hat(1:20))
title("first samples corrupted")

subplot(122)
e_hat = e_hat(100:end);
plot(e_hat(1:20))
title("no corruption")
ACFnPACFnNormplot(e_hat, 32);

% FPE??
% model est. ...

%% 2.1 loop
close all

y = iddata(y1);

lim = 4;
minfpe = inf;
disp("starting loop")

for na = 1:lim

    for nc = 1:lim
        arma_model = armax(y, [na nc]);
        % e_hat = filter(arma_model.a, arma_model.c, y1);
        % e_hat = e_hat(100:end);
        f = fpe(arma_model);
        % if f < minfpe
        % minfpe = f;
        % minna = na;
        % minnc = nc;
        fprintf("na = %d, nc = %d, fpe = %d\n", na, nc, f)
        % e_hat = myFilter(arma_model.a, arma_model.c, y1);
        % ACFnPACFnNormplot(e_hat, 32);
        % hold on
        % pause;
        % end
    end

end

disp("done")

%% 2.2

clear
close all

load data.dat % arma(1,1)
load noise.dat

data = iddata(data);

% ACFnPACFnNormplot(data, 32);
% Some ringing in the ACF, 4 non zero values in PACF => AR(4)?

%% 2.2 AR(p) est.
% This shows that a AR(3) is enough
for i = 1:5
    arx_model = arx(data, i);
    resid(arx_model, data);
    title(sprintf("AR(%d) process", i))
    pause;
end

%% 2.2 arp-model
ar3_model = arx(data, 3);
present(ar3_model)

resid(ar3_model, data);
rar3 = resid(ar3_model, data);

na = 3;
figure
plot(noise(na:end))
hold on
plot(rar3.y(na:end))
legend("noise", "residuals")
title("noise vs. residuals")

%% 2.2 arma-model

arma12 = armax(data, [1 2]);
resid(arma12, data);
rarma12 = resid(arma12, data);

na = 3;
figure
plot(noise(na:end))
hold on
plot(rarma12.y(na:end))
legend("noise", "residuals")
title("noise vs. residuals")

% arma(1,2) had a lower MSE than ar(3).

% show with plot?

%% 2.3

close all
clear

N = 10000; % 600

rng(0)
A = [1 -1.5 0.7];
C = [1 zeros(1, 11) -0.5];
A12 = [1 zeros(1, 11) -1];
Astar = conv(A, A12);
e = randn(N, 1);
y = filter(C, Astar, e);
y = y(101:end);
plot(y)

ACFnPACFnNormplot(y, 32);

y_s = filter(A12, 1, y);
y_s = y_s(length(A12):end);
data = iddata(y_s);

A = [1 0 0]; B = []; C = [];
model_init = idpoly(A, B, C);
model_armax = pem(data, model_init);

resid(model_armax, data);
rarma = resid(model_armax, data);

ACFnPACFnNormplot(rarma.OutputData, 32);

A = [1 0 0]; B = []; C = [1 zeros(1, 12)];
model_init = idpoly(A, B, C);
model_init.Structure.c.Free = [zeros(1, 12) 1];
model_armax = pem(data, model_init);

resid(model_armax, data);
rarma = resid(model_armax, data);
ACFnPACFnNormplot(rarma.OutputData, 32);

checkIfNormal(rarma.OutputData, "arma");
checkIfWhite(rarma.OutputData);

%% 2.3 no season

arma12 = armax(y, [1 12]);
resid(arma12, data);
rarma12 = resid(arma12, data);

ACFnPACFnNormplot(arma12.y, 32);

%% 2.4

load svedala
data = svedala;

N = length(data);
X = [ones(N, 1) (1:N)'];
mk = (X' * X) \ X' * data;

plot_trend = 0;

if plot_trend
    figure
    subplot(211)
    plot(data)
    hold on
    plot(mk(1) + mk(2) * (1:N)', 'r')
    legend("data", "linear trend")
    title("data and linear trend")

    subplot(212)
    plot(data - mk(1) - mk(2) * (1 * N)')
    title("data without trend")
end

% differenting data
data = filter([1 -1], 1, data);
data = data(2:end);

% ACFnPACFnNormplot(data, 32);
% AR(2)?

A = [1 0 0]; B = []; C = [];
model_init = idpoly(A, B, C);
model_armax = pem(data, model_init);
figure
resid(model_armax, data);
rarma = resid(model_armax, data);
ACFnPACFnNormplot(rarma.OutputData, 32);

A = [1 0 0]; B = []; C = [1 zeros(1, 24)];
model_init = idpoly(A, B, C);
model_armax = pem(data, model_init);

figure
resid(model_armax, data);
rarma = resid(model_armax, data);
