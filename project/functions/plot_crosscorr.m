function [Cxy, lags] = plot_crosscorr(x, y, noLags)
    % FUNCTION plots the cross correlation between x, y.
    % returns Cxy and lags

    if nargin < 3
        noLags = 32;
    end

    figure
    [Cxy, lags] = xcorr(x, y, noLags, "coeff");
    stem(lags, Cxy);
    hold on
    condInt = 2 * ones(1, length(lags)) ./ sqrt(length(y));
    plot(lags, condInt, "r--")
    plot(lags, -condInt, "r--")
    hold on
    xlabel("Lag")
    ylabel("Amplitude")
    title("Crosscorrelation")
end
