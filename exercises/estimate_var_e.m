function vare = estimate_var_e(C, A, y)
    k = 1;
    [Fk, Gk] = polydiv(C, A, k);

    filter_skip = length(A);
    yhat_k = myFilter(Gk, C, y, filter_skip);
    ehat = y(filter_skip:end) - yhat_k;

    vare = var(ehat);
end
