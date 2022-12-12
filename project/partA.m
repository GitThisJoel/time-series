close all
clear

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");

load('raw_data.mat')
halt_konc = raw_data(:, 1);
[modeling_set, validation_set, test_set] = load_data(halt_konc);

%%
%% Vi behöver hantera outliers, finns gott om dom...
% spolvatten verkar vara rätt konstant så går nog bra att bara kvotera in
% om vi behöver det, annars kanske försumbart

% De signaler vi kommer jobba med är halten i koncentratet "halt_konc" och
% som insignal kommer vi sedan att använda "halt_ing_rep" och eventuellt 
% "flow_ing_rep" som insignaler.

%% A. modeling without an external input

%% Utreda om outliers behöver behandlas
figure
acf(modeling_set, 32, 0.05, 1);
figure
tacf(modeling_set, 32, 0.04, 0.05, 1);

% Relativt lika, testar utan out-lier-behandlign till att börja med

%% Check if wee need a transform
figure; 
lambda_max = bcNormPlot(modeling_set,3);
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)

figure
normplot( modeling_set )

% Har högt lambda redan vid 0, men ännu högre vid 2, väljer att inte göra
% någon transform till att börja med.

%% Inspect data
ACFnPACFnNormplot(modeling_set,50);

%% Log(data) and inspect again

modeling_set_mod = modeling_set + 15;
min(modeling_set_mod)
modeling_set_mod = log(modeling_set_mod);

ACFnPACFnNormplot(modeling_set_mod, 50)

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
model_AR1 = estimateARMA(modeling_set, [1 1], [1], 'AR(1) model',50)
present(model_AR1)
%% Test adding an MA(2)-component
model_ARMA12 = estimateARMA(modeling_set, [1 1], [1 0 1], "ARMA(1,2) model", 50)

%% Test adding an AR(3)-component
% removed MA(2) again since it became insignificant
model_AR3 = estimateARMA(modeling_set, [1 1 0 1], [1], "AR(3) model", 50)

%% Test adding an MA(3)-component
model_ARMA33 = estimateARMA(modeling_set, [1 1 0 1], [1 0 0 1], "ARMA(3,3) model", 50)
present(model_ARMA33)
% den här modellen tar bort alla låga signifikanta parametrar, utvärdera
% den vidare, jämför med vanlig AR(1), hur mycket vinner man på den högre
% ordningen?.. (alla parametrar är signifikanta)
