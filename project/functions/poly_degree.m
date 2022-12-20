function deg = poly_degree(coeffs)
    % FUNCTION deg = poly_degree(coeffs)
    %
    % The function finds the degree of a polynomial that is represented as
    % a vector, e.g. poly_deg([1 2 3]) = 2.
    % This function will not take trailing zeros in to account.

    deg = length(coeffs(1:find(coeffs, 1, 'last'))) -1;
end
