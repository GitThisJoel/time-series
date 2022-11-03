addpath("../../CourseMaterial/Code/data");

%% t1
[y, Fs] = audioread("fa.wav");

% sound(y, Fs)

time = (0:size(y, 1) - 1) * 1 / Fs;
plot(time, y')
xlabel("time [s]")
ylabel("Amplitude")
title("Sound wave")
% each peak represent a vowel?

%% t2

ind = 0.7 * Fs;
samples = y(0.7 * Fs:0.7 * Fs + 200); % 201 samples
sample_time = (ind:ind + 200) * 1 / Fs;

subplot(411)
plot(sample_time, samples)
title("sound wave")

nolags = 31
lags = acf(samples, nolags);
covar = cat(1, flipud(lags), lags);

subplot(412)
plot(covar, '*');
title("covariance function")

N = size(samples, 1);
padd = 2^9 - N;

sd = fftshift(abs(fft(samples .* hamming(N))).^2 / N);
sd_zp = fftshift(abs(fft(samples .* hamming(N), padd)).^2 / N);
ff = (0:padd - 1)' / padd - 0.5;
f = Fs * (0:(N / 2)) / N;

subplot(413)
% semilogy(f, sd(1:N / 2 + 1));
semilogy(sd);
title("estimation of spectral density with out zero padding")

subplot(414)
semilogy(ff * Fs, sd_zp);
title("estimation of spectral density with zero padding")

% for this time span the fundamental freq seems to be ~250 Hz

%% t3

% done above :)
