close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");
addpath("functions");

load('raw_data.mat')
halt_konc = raw_data(:, 1);
[modeling_set, validation_set, test_set, index_validation, index_test] = load_data(halt_konc);

%% Någr observationer om datan
% spolvatten verkar vara rätt konstant så går nog bra att bara kvotera in
% om vi behöver det, men är sannorlikt försumbart då mätningarna är på % av
% fast material

% De signaler vi kommer jobba med är halten i koncentratet "halt_konc" och
% som insignal kommer vi sedan att använda "halt_ing_rep" (och eventuellt 
% "flow_ing_rep", men sanorlikt inte) som insignaler.

%% A. modeling without an external input

%% Utreda om outliers behöver behandlas
figure
subplot(1,2,1)
acf(modeling_set, 32, 0.05, 1);
subplot(1,2,2)
tacf(modeling_set, 32, 0.04, 0.05, 1);

% Relativt lika, testar utan out-lier-behandlign till att börja med

%% Check if wee need a transform
figure; 
lambda_max = bcNormPlot(modeling_set,3);
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)

figure
normplot( modeling_set )

% Har högt lambda redan vid 0, men ännu högre vid 2, väljer att inte göra
% någon transform till att börja med, men kanske eventuellt testar att 
% logga sen.

%% Inspect data
ACFnPACFnNormplot(modeling_set,50);

% Tar Lång tid innan ACF:en klingar av, men differentierar man datat så
% blir det ingen struktur kvar, så undvikaer att göra det och testar att
% modelera istället då AR(1)-komponenten är mycket tydlig

%% Log(data) and inspect again

modeling_set_mod = modeling_set + 15;
min(modeling_set_mod)
modeling_set_mod = log(modeling_set_mod);

ACFnPACFnNormplot(modeling_set_mod, 50);
figure
normplot(modeling_set_mod)

% Logga datat gör det något bättre i hur fort ACF:en klingar av, men inte
% jättemycket så testar utan transform till att börja med
% Normalfördelningen som approximation är dock tveksamm, det händer en hel
% del märkliga grejjer i svansarna och delvis lite nära mitten med på övre
% halvan. Värt att ta med sig vidare i analysen.

%% Långt beroende i ACF: bör differentiera datan, bara z^(-1) först
%AS = [1 -1];
%modeling_set = filter(AS, 1, modeling_set);
%modeling_set = modeling_set(2:end);

% att göra detta slår ihjäl typ all struktur i datan, så en direkt
% differentiering verkar inte vara lämpligt

%% Start with testing an AR(1)
[model_AR1, res_model_AR1] = estimateARMA(modeling_set, [1 1], [1], 'AR(1) model',50);
present(model_AR1)
checkIfWhite(res_model_AR1);
%% Test adding an MA(2)-component
[model_ARMA12, res_model_ARMA12] = estimateARMA(modeling_set, [1 1], [1 0 1], "ARMA(1,2) model", 50);
checkIfWhite(res_model_ARMA12);
%% Test adding an AR(3)-component
% removed MA(2) again since it became insignificant
[ model_AR3, res_model_AR3] = estimateARMA(modeling_set, [1 1 0 1], [1], "AR(3) model", 50);
checkIfWhite(res_model_AR3);
checkIfNormal(res_model_AR3, '');
% Den här är väldigt nära att vara vit
%% Test adding an MA(3)-component
[model_ARMA33, res_model_ARMA33] = estimateARMA(modeling_set, [1 1 0 1], [1 0 0 1], "ARMA(3,3) model", 50);
present(model_ARMA33)
checkIfWhite(res_model_ARMA33);
checkIfNormal(res_model_ARMA33, '');
% den här modellen tar bort alla låga signifikanta parametrar, utvärdera
% den vidare, jämför med vanlig AR(1), hur mycket vinner man på den högre
% ordningen?.. (alla parametrar är signifikanta), residualen är vit men
% inte normalfördelad, ser mer t-fördelad ut.


%% Gör predictions med modellen!
% gör det här på hela datasettet så att vi inte behöver filtrera i början
% av valideringsdata/testdata

close all

halt_konc = halt_konc - mean(halt_konc);

y = halt_konc;
model = model_AR1
C = model.C;
A = model.A;
k = 9; %Prediction horizon

[Fk , Gk] = polydiv ( C, A, k ) ;

filter_skip = length(A)
yhat_k = myFilter(Gk, C, y, filter_skip);
ehat = y(filter_skip:end) - yhat_k;
ehat = [zeros(filter_skip-1,1) ; ehat];
figure
plot(ehat)

figure
plot(y)
hold on
plot(yhat_k)
legend('y', 'y\_hat')

validation_ehat = ehat(index_validation(1):index_validation(2));
test_ehat = ehat(index_test(1):index_test(2));


ACFnPACFnNormplot(validation_ehat, 50);
ACFnPACFnNormplot(test_ehat, 50);

checkIfWhite(validation_ehat);
var(validation_ehat)

checkIfWhite(test_ehat);
var(test_ehat)
