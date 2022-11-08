function [acfEst, pacfEst] = ACFnPACFnNormplot(data, noLags, signLvl, titleStr)
    % FUNCTION ACFnPACFnNormplot (data, noLags, signLvl, title)
    %
    % Plot the estimated ACF and PACF as well as norm plot in a single figure.

    if nargin < 4
        titleStr = "";
    end

    if nargin < 3
        signLvl = 0.05;
    end

    figure
    subplot(311)
    acfEst = acf(data, noLags, signLvl, 1);
    title(sprintf("ACF %s", titleStr))

    subplot(312)
    pacfEst = pacf(data, noLags, signLvl, 1);
    title(sprintf("PACF %s", titleStr))

    subplot(313)
    normplot(data);
end
