function [modeling_set, validation_set, test_set] = load_data(raw_data, start, l_modelling, l_test)

    halt_konc = halt_konc - mean(halt_konc);

    if nargin < 2
        start = 6800;
    end

    if nargin < 3
        l_modelling = 2200;
    end

    if nargin < 4
        l_test = 700;
    end

    modeling_set = halt_konc(start:start + l_modelling);
    validation_set = halt_konc(start + l_modelling + 1:start + l_modelling + 1 + l_test);
    test_set = halt_konc(round(start + length(halt_konc) * 0.47):round(start + length(halt_konc) * 0.47) + l_test);
end
