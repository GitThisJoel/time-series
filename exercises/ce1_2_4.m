clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

load svedala
data = svedala;

residOpts = residOptions("MaxLag", 32);

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
% AR(1)?

A = [1 0 0]; B = []; C = [];
model_init = idpoly(A, B, C);
model_armax = pem(data, model_init);

figure
resid(model_armax, data, residOpts);
rarma = resid(model_armax, data);
ACFnPACFnNormplot(rarma.OutputData, 32);

A = [1 0 0]; B = []; C = [1 zeros(1, 24)];
model_init = idpoly(A, B, C);
model_armax = pem(data, model_init);

figure
resid(model_armax, data, residOpts);
rarma = resid(model_armax, data);
