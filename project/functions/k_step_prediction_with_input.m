function [yhatk, ehaty, ehatx] = k_step_prediction_with_input(input_model, model, k, x, y, noLags, index_validation, index_test)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% first predict the input:
[Fx, Gx] = polydiv( input_model.C, input_model.A, k );
xhatk = filter(Gx, input_model.C, x);

shiftK = round( mean( grpdelay(Gx, 1) ) );    % Compute the average group delay of the filter.
figure
plot(x(index_validation(1):index_validation(2)-shiftK))
hold on
plot(xhatk(index_validation(1)+shiftK+1:index_validation(2)))
title(sprintf('inputput %i-step prediction',k))
legend('Real value','Prediction')




ehatx = x - xhatk;

figure
acf( ehatx(index_validation(1):index_validation(2)), noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step input prediction residual', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
ehatx = ehatx(10:end);
checkIfWhite( ehatx );
pacfEst = pacf( ehatx, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );

%% predict output

% Form the BJ prediction polynomials. In our notation, these are
%   A1 = foundModel.D
%   C1 = foundModel.C
%   A2 = foundModel.F
% 
% The KA, KB, and KC polynomials are formed as:
%   KA = conv( A1, A2 );
%   KB = conv( A1, B );
%   KC = conv( A2, C1 );

KA = conv( model.D, model.F );
KB = conv( model.D, model.B );
KC = conv( model.F, model.C );

% Form the ARMA prediction for y_t (note that this is not the same G
% polynomial as we computed above (that was for x_t, this is for y_t).
[Fy, Gy] = polydiv( model.C, model.D, k );

% Compute the \hat\hat{F} and \hat\hat{G} polynomials.
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

% Form the predicted output signal using the predicted input signal.
yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, x) + filter(Gy, KC, y);

%%

shiftK = round( mean( grpdelay(Gy, 1) ) ); 
figure
plot(y(index_validation(1):index_validation(2)-shiftK))
hold on
plot(yhatk(index_validation(1)+1+shiftK:index_validation(2)))
title(sprintf('output %i-step prediction',k))
legend('Real value','Prediction')


ehaty = y - yhatk;

figure
acf( ehaty(index_validation(1):index_validation(2)), noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step output prediction residual', k) )
ehaty = ehaty(10:end);
checkIfWhite( ehaty );
pacfEst = pacf( ehaty, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );


%% add zeros to get index-position right for future analysis
ehatx = [zeros(length(yhatk)-length(ehatx),1);ehatx];
ehaty = [zeros(length(yhatk)-length(ehaty),1);ehaty];
end