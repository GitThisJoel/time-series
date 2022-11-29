function CCF(w_t, eps_t, titleStr, M, n)

    if nargin < 5
        n = length(w_t);
    end

    if nargin < 4
        M = 40;
    end

    if nargin < 3
        titleStr = "";
    end

    stem(-M:M, crosscorr(w_t, eps_t, M), "b");
    title(sprintf('Cross correlation function %s', titleStr)), xlabel('Lag')
    hold on
    plot(-M:M, 2 / sqrt(n) * ones(1, 2 * M + 1), '--')
    plot(-M:M, -2 / sqrt (n) * ones(1, 2 * M + 1), '-- ')
    hold off
end
