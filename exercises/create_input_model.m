function [Mba2, etilde] = create_input_model(d, r, s, x, y)

    if nargin < 5
        % TODO:
        % split 'x' in new x and y
    end

    A2 = [1 zeros(1, r)];
    B = [zeros(1, d) ones(1, s + 1)];
    Mi = idpoly([1], [B], [], [], [A2]);
    z = iddata(y, x);
    Mba2 = pem(z, Mi); present(Mba2)
    etilde = resid(Mba2, z);

    checkIfNormal(etilde.y, ""); % close
    checkIfWhite(etilde.y);

    ACFnPACFnNormplot(etilde.y, 32);
end
