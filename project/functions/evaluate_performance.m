function [var_validation, var_test] = evaluate_performance(ehat, index_validation, index_test)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    validation_ehat = ehat(index_validation(1):index_validation(2));
    ACFnPACFnNormplot(validation_ehat, 50, 0.05, 'Validation set', 0);

    checkIfWhite(validation_ehat);
    var_validation = var(validation_ehat);

    var_test = -1;

    if  (nargin > 2)
        test_ehat = ehat(index_test(1):index_test(2));
        ACFnPACFnNormplot(test_ehat, 50, 0.05, 'Test set');

        checkIfWhite(test_ehat);
        var_test = var(test_ehat);
    end

end
