function k_step_prediction(k, A, C, y, var_e)
    [Fk, Gk] = polydiv(C, A, k);

    % theo_e = myFilter(A, C, y, length(C));

    filter_skip = length(C); % max(length(Gk), length(C));
    yhat_k = myFilter(Gk, C, y, filter_skip);
    ehat = y(filter_skip:end) - yhat_k;

    mean_e = mean(ehat);
    % var_e_k = var(ehat);
    theo_var_e = norm(Fk, 2).^2 * var_e;

    disp("k: " + k)
    disp("mean ehat: " + mean_e)
    disp("var_e: " + var_e)
    disp("theo var e: " + theo_var_e)

    conf_int = 1.96 * sqrt(var_e) * norm(Fk, 2);
    ci = [(mean_e - conf_int) (mean_e + conf_int)];
    disp("95% conf interval: (" + ci(1) + ", " + ci(2) + ")")

    num_outside = sum(abs(ehat) > ci(2));
    disp("[%] outside conf int: " + num_outside / length(ehat))

    figure
    plot(y(filter_skip:end));
    hold on
    plot(yhat_k);
    legend("process", "prediction")
    title(sprintf("Process vs. prediction for k = %d", k))

    ACFnPACFnNormplot(ehat, 32, 0.05, "residual for k = " + k);

    M = length(ehat);
    figure
    plot(ehat)
    hold on
    plot(0:M, ci(1) * ones(1, M + 1), '--')
    plot(0:M, ci(2) * ones(1, M + 1), '--')
    hold off
    title(sprintf("conf interval for k = %d", k))
end
