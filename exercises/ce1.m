clear
close all
% hej
addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

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

% the ACF/PACF plot points to a AR(2) process.
% normplot = ?
% ACFnPACFnNormplot(y1, 32);

e_hat = filter(ar_model.a, ar_model.c, y1);

figure 
subplot(121)
plot(e_hat(1:20))
title("first samples corrupted")

subplot(122)
e_hat = e_hat(100:end);
plot(e_hat(1:20))
title("no corruption")

% FPE??
% model est. ...

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
for i=1:5
    arx_model = arx(data, i);
    resid(arx_model, data);
    title(sprintf("AR(%d) process", i))
    pause;
end

%% 2.2 cont.
ar3_model = arx(data, 3);
resid(ar3_model, data);
rar3 = resid(ar3_model, data);

na = 3; 
figure
plot(noise(na:end))
hold on 
plot(rar3.y(na:end)) 
legend("noise", "residuals")
title("noise vs. residuals")

%% 2.3

%% 2.4
