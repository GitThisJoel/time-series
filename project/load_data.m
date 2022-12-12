function [modeling_set, validation_set, test_set] = load_data(raw_data, start, l_modelling, l_test)

    raw_data = raw_data - mean(raw_data);

    if nargin < 2
        start = 6800;
    end

    if nargin < 3
        l_modelling = 2200;
    end

    if nargin < 4
        l_test = 700;
    end

    modeling_set = raw_data(start:start + l_modelling);
    validation_set = raw_data(start + l_modelling + 1:start + l_modelling + 1 + l_test);
    test_set = raw_data(round(start + length(raw_data) * 0.47):round(start + length(raw_data) * 0.47) + l_test);
end
