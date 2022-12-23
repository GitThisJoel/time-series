% To do:

% testa fler modellordningar
% skriv funktion för att utvärdera ressultaten

close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");
addpath("functions");

load('raw_data.mat')

halt_konc = raw_data(:, 1);
[modeling_set, validation_set, test_set, index_validation, index_test] = load_data(halt_konc);

halt_ing_rep = create_input_signal(raw_data);
[inp_modeling_set, inp_validation_set, inp_test_set] = load_data(halt_ing_rep);

noLags = 50;

y = modeling_set;
x = inp_modeling_set;

%% plot the (modelling) data
close all

figure
subplot(211)
plot(x);
title("Input signal")
subplot(212)
plot(y);
title("Output signal")

ACFnPACFnNormplot(x, noLags);

% Ser vad som händer med inputten av att differentiera den
AS = [1 -1];
xdiff = filter(AS, 1, x);
xdiff = xdiff(2:end);
ACFnPACFnNormplot(xdiff, noLags);

%% cross correlation between input and output data

close all
% plot_crosscorr(x, y, noLags);
CCF(x, y, "y vs x", noLags);

%% pre white

% Test with differentiated data:
%[input_model, input_model_res] = estimateARMA(xdiff, [1 1 1 ], [1 0 1 ], "ARMA(3,11) of input", noLags);
% har svårt att generera modeller som ens är i närheten av att vara vita
% här, så återgår till det odifferentierade datasettet.
%% close all

[input_model, input_model_res] = estimateARMA(x, [1 1 1 1], [1 0 1 1 zeros(1, 7) 1], "ARMA(3,11) of input", noLags);
%[input_model, input_model_res ] = estimateARMA(x, [1 1 1 1], [1 0 1 1], "ARMA(3,3) of input", noLags);
present(input_model)
checkIfWhite(input_model_res);

% tried diff, did not help.
[w_t, eps_t] = pre_white(input_model, x, y);
cc = CCF(w_t, eps_t, "", noLags);

% modellen blir såpass mycket bättre att jag behåller 11-komponenten i
% inputmodellen
% men prediktionene blir inte direkt bättre och får trevligare ACF utan
% 11:an så kanske kör utan ändå.
% x-prediktionene får lägre varians med 11-an och y-prediktionen påverkas
% inte så 11:an får vara med! 
%% choose model order
close all
% r=1, d=3, s=2 funkar men ger insignifikanta parametrar
r = 1; 
[~, d] = max(cc); d = d - noLags - 1; % denna behåller vi!
s = 2;



[Mba2, etilde] = create_input_model(d, r, s, x, y);

ACFnPACFnNormplot(etilde.OutputData, noLags);
na = 1;
nc = 3;
% etilde_model = estimateARMA(etilde.y, [1 ones(1, na)], [1 ones(1, nc)], "ar1\_etilde", noLags);
etilde_model = estimateARMA(etilde.y, [1 1], [1 1 1 1], "ar1\_etilde", noLags);
%present(etilde_model)


%% entire model

% A1 = [1 zeros(1, etilde_model.na)];                 % D
% A2 = [1 zeros(1, Mba2.nf)];                         % F
% B = [zeros(1, Mba2.nk) ones(1, Mba2.nb)];
% C = [1 zeros(1, etilde_model.nc)]; % ones or zeros here?
% Mi = idpoly(1, B, C, A1, A2);
% z = iddata(y, x);
% MboxJ = pem(z, Mi);
% present(MboxJ)
% ehat = resid(MboxJ, z);


D = [1 1];                          % A1
F = [1 ];                          % A2
B = [0 0 0 1 1 1 ];
C = [1 0 1 1];
z = iddata(y, x);
MboxJ = estimateBJ(y,x,C,D,B,F,'',50)
%MboxJ = pem(z, Mi);
present(MboxJ)
ehat = resid(MboxJ, z);

ACFnPACFnNormplot(ehat.y, 32);

checkIfNormal(ehat.y, ''); close
checkIfWhite(ehat.y);

CCF(x, ehat.y, 'x vs ehat');

present(MboxJ)
modelB = MboxJ;
%save('model_part_B.mat','modelB')
%% Do predictions with model:
k = 9; % prediction horizon
% The data we do the predictions on:
x = halt_ing_rep - mean(halt_ing_rep);
y = halt_konc - mean(halt_konc);

[yhat, ehaty, ehatx] = k_step_prediction_with_input(input_model, MboxJ, k, x, y, noLags);

% Analyse results:
varx = evaluate_performance(ehatx, index_validation, index_test);
vary = evaluate_performance(ehaty, index_validation, index_test);

fprintf('Variance of y is %s, variance of x is %d \n', vary, varx)

