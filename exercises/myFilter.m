function [filtered] = myFilter(a, c, data, na)
    %myFilter (a, c, data, na)
    %   Fitlers data and removes na samples, if na is omitted 100 samples are removed.

    if nargin < 4
        na = 100;
    end

    filtered = filter(a, c, data);
    filtered = filtered(na:end);
end
