function Ct_inp = create_ct_inp(t, x, ehat, Ax, Cx, k0)
    % create_ct_inp - Calculates the Ct matrix for a kalman filter.
    %
    % Syntax: Ct_inp = create_ct_inp(t,x,ehat,Ax,Cx,k0)

    degA = poly_degree(Ax);
    degC = poly_degree(Cx);

    nzA = zeros(size(Ax)); nzA(Ax ~= 0) = 1;
    nzC = zeros(size(Cx)); nzC(Cx ~= 0) = 1;

    if length(x) ~= degA
        xt = x(t - degA + k0:t + k0)';
    else
        xt = x;
    end

    Ctx = flip(xt(find(flip(nzA(2:end))))); % b0 != 1

    et = ehat(t - degC + k0:t + k0);
    Cte = flip(et(find(flip(nzC))));

    Ct_inp = [-Ctx Cte];

end
