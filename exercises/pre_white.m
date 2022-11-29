function [w_t, eps_t] = pre_white(model, x, y)
    w_t = myFilter(model.A, model.C, x, length(model.A));
    eps_t = myFilter(model.A, model.C, y, length(model.A));
end
