clear
close all

addpath("../CourseMaterial/Code/data");

%% 2.1

A1 = [1 -1.79 0.84];
C1 = [1 -0.18 -0.11];

A2 = [1 -1.79];
C2 = [1 -0.18 -0.11];

% ARMA process A(z)y_t = C(z)e_t => y_t = C(z)/A(z)e_t
arma1 = idpoly(A1, [], C1);
arma2 = idpoly(A2, [], C2);

rng(0)
sigma2 = 4;
N = 300;
e = sqrt(sigma2) * rand(N, 1);

y1 = filter(arma1.c, arma1.a, e);
y2 = filter(arma2.c, arma2.a, e);

y1 = y1(101:end);
y2 = y2(101:end);

% arma2 is unsable since one pole is outside of the unit cirle.
% why...

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
