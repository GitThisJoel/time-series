clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

%% 2.1 load data

load('tar2.dat')
load('thx.dat')

%% 2.1 plot tar2 and thx

close all
figure
subplot(211)
plot(tar2)
title("tar2")
subplot(212)
plot(thx)
title("thx")

%% 2.1 recusiveAR
close all
lambdas = [1 0.95 0.9];

plot = 1;

if plot
    figure
end

for i = 1:length(lambdas)
    X = recursiveAR(2);
    X.ForgettingFactor = lambdas(i);
    X.InitialA = [1 0 0];

    N = length(tar2);
    Aest = zeros(N, 3);
    yhat = zeros(1, N);

    for kk = 1:N
        [Aest(kk, :), yhat(kk)] = step(X, tar2(kk));
    end

    if plot
        subplot(1, length(lambdas), i)
        plot(Aest)
        title(sprintf("lambda = %f", lambdas(i)))
        % subplot(2, length(lambdas), i + length(lambdas))
        hold on
        plot(thx)
        legend("est1", "est2", "est3", "thx1", "thx2", "Location", "southeast")
    end

end

%% 2.1 choose lambda
close all

n = 100;
lambda_line = linspace(0.85, 1, n);
ls2 = zeros(n, 1);
yhat = zeros(n, 1);

for i = 1:length(lambda_line)
    reset(X);
    X.ForgettingFactor = lambda_line(i);
    X.InitialA = [1 0 0];

    for kk = 1:length(tar2)
        [~, yhat(kk)] = step(X, tar2(kk));
    end

    ls2(i) = sum((tar2 - yhat).^2);
end

min_ind = find(ls2 <= min(ls2));
opt_lambda = lambda_line(min_ind);

plot(lambda_line, ls2);
title("lambda estimate (min = " + opt_lambda + ")")
