function [acfEst, pacfEst] = ACFnPACFnNormplot(data, noLags, signLvl, titleStr, includeZeroLag)
    % FUNCTION ACFnPACFnNormplot (data, noLags, signLvl, title)
    %
    % Plot the estimated ACF and PACF as well as norm plot in a single figure.

    if nargin < 5
        includeZeroLag = 1;
    end
    if nargin < 4
        titleStr = "";
    end

    if nargin < 3
        signLvl = 0.05;
    end

    figure
    subplot(311)
    acfEst = acf(data, noLags, signLvl, 1, 0,includeZeroLag);
    title(sprintf("ACF %s", titleStr))

    subplot(312)
    pacfEst = pacf(data, noLags, signLvl, 1, includeZeroLag);
    title(sprintf("PACF %s", titleStr))

    subplot(313)
    normplot(data);
end
