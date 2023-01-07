function Ct = create_ct(t, y, x, ehat, Ax, Bx, Cx, k0, Cty)
    %create_ct - Description
    %
    % Syntax: Ct = create_ct(y,x,ehat,degA,degB,degC,k0)
    %
    % Calculates the Ct matrix for a kalman filter.
    degA = poly_degree(Ax);
    degB = poly_degree(Bx);
    degC = poly_degree(Cx);

    nzA = zeros(size(Ax)); nzA(Ax ~= 0) = 1;
    nzB = zeros(size(Bx)); nzB(Bx ~= 0) = 1;
    nzC = zeros(size(Cx)); nzC(Cx ~= 0) = 1;

    if nargin < 9

        Cty = [];

        if length(y) ~= 0
            yt = y(t - degA + k0:t - 1 + k0)';
            Cty = flip(yt(find(flip(nzA(2:end)))));
        end

    end

    if length(x) ~= (degB + 1)
        xt = x(t - degB + k0:t + k0)';
    else
        xt = x;
    end

    Ctx = flip(xt(find(flip(nzB)))); % b0 != 1

    et = ehat(t - degC + k0:t + k0);
    Cte = flip(et(find(flip(nzC))));

    Ct = [-Cty Cte Ctx];
end
