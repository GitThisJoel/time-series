function [ci1, ci2, delta] = ci(data, perc)

    if nargin < 2
        perc = 0.95;
    end

    if perc == 0.95
        alpha = 1.96;
    elseif perc == 0.99
        alpha = 2.5758;
    end

    m = mean(data);
    sigma = sqrt(var(data));
    ci1 = m - alpha * sigma;
    ci2 = m + alpha * sigma;
    delta = alpha * sigma;
end
