function [modeling_set, validation_set, test_set, index_validation, index_test] = load_data(data, start, l_modelling, l_test)

    data = data - mean(data);

    if nargin < 2
        start = 6800;
    end

    if nargin < 3
        l_modelling = 2200;
    end

    if nargin < 4
        l_test = 700;
    end

    modeling_set = data(start:start + l_modelling);
    validation_set = data(start + l_modelling + 1:start + l_modelling + 1 + l_test);
    test_set = data(round(start + length(data) * 0.47):round(start + length(data) * 0.47) + l_test);
    index_validation = [start + l_modelling + 1 start + l_modelling + 1 + l_test];
    index_test = [round(start + length(data) * 0.47) round(start + length(data) * 0.47) + l_test];
end
