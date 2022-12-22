function [var_validation, var_test] = evaluate_performance(ehat, index_validation, index_test)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
validation_ehat = ehat(index_validation(1):index_validation(2));
test_ehat = ehat(index_test(1):index_test(2));


ACFnPACFnNormplot(validation_ehat, 50);
ACFnPACFnNormplot(test_ehat, 50);

checkIfWhite(validation_ehat);
var(validation_ehat)

checkIfWhite(test_ehat);
var(test_ehat)

end