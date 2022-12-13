close all
clear all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");

%% Load data and split it into sets

load('raw_data.mat')


time = (1:length(raw_data(:,1)))* 6*60;   % Time in [s], concider minutes instead 
time = time';

% extract signals in raw data
halt_konc = raw_data(:,1);
halt_ing_flot_1 = raw_data(:,2);
halt_ing_flot_2 = raw_data(:,3);
halt_raw_konc_1 = raw_data(:,4);
halt_raw_konc_2 = raw_data(:,5);
halt_scav_konc_1 = raw_data(:,6);
halt_scav_konc_2 = raw_data(:,7);
flow_raw_remill_1 = raw_data(:,8);
flow_raw_remill_2 = raw_data(:,9);
flow_scav_remill_1 = raw_data(:,10);
flow_scav_remill_2 = raw_data(:,11);
flow_ing_1 = raw_data(:,13);
flow_ing_2 = raw_data(:,12);

% construct signals of interest (higher abstraction level)
flow_ing = flow_ing_1 + flow_ing_2 ; % [ton/h], kanske konvertera till m3/s?
halt_ing_flot = (halt_ing_flot_1.*flow_ing_1+halt_ing_flot_2.*flow_ing_2 )./flow_ing;
halt_raw_konc = (halt_raw_konc_1.*flow_ing_1+halt_raw_konc_2.*flow_ing_2 )./flow_ing;
flow_raw_remill = flow_raw_remill_1 + flow_raw_remill_2;
halt_scav_konc = (halt_scav_konc_1.*flow_ing_1+halt_scav_konc_2.*flow_ing_2 )./flow_ing;
flow_scav_remill = flow_scav_remill_1 + flow_scav_remill_2; 
flow_ing_rep = (flow_raw_remill + flow_scav_remill);
halt_ing_rep = (halt_raw_konc.*flow_raw_remill + halt_scav_konc.*flow_scav_remill) ./flow_ing_rep;

figure 
plot(halt_raw_konc)
title('halt rå konc')

figure
plot(flow_raw_remill)
hold on
plot(flow_scav_remill)

time = time/60/60/24;
figure
subplot(2,1,1)
plot(time,halt_konc)
title('Grade in final concentarte (y)')
ylabel('Grade [%]')
xlim([0 122])
xlabel('Time [days]')
%subplot(3,1,2)
%plot(time, halt_ing_flot)
%title('Concentration of incoming ore (x)')
subplot(2,1,2)
plot(time, halt_ing_rep)
title('Grade of middle product (x)')
xlabel('Time [days]')
ylabel('Grade [%]')
xlim([0 122])

%%

%% Vi behöver hantera outliers, finns gott om dom...
% spolvatten verkar vara rätt konstant så går nog bra att bara kvotera in
% om vi behöver det, annars kanske försumbart

% De signaler vi kommer jobba med är halten i koncentratet "halt_konc" och
% som insignal kommer vi sedan att använda "halt_ing_rep" och eventuellt 
% "flow_ing_rep" som insignaler.

%% A. modeling without an external input
close all

halt_konc = halt_konc - mean(halt_konc);
% för att hålla ungefär samma datastorlek som i orginalprojektet bör
% modeleringsdatat innehålla ca 1700 samples, this corresponds to roughly 7
% days with our sampling intervall. - inte viktigt, gjorde dom lite större
start = 6800;
l_modelling = 2200;
l_test = 700;
modeling_set = halt_konc(start:start+l_modelling);
validation_set = halt_konc(start+l_modelling+1:start+l_modelling+1+l_test );
test_set = halt_konc(round(start+length(halt_konc)*0.47):round(start+length(halt_konc)*0.47)+l_test);

figure
plot(modeling_set)

figure
plot(validation_set)

figure
plot(test_set)

figure
plot(halt_konc)
hold on
xline(start)
xline(start+l_modelling)
xline(start+l_modelling+l_test)
xline(start+length(halt_konc)*0.47)
xline(start+length(halt_konc)*0.47+l_test)

% I think this split looks ok, lets see if the data-size is sufficent

%% Naive predictors:

% A suitable naive predictor for the one step rediction is the current
% value of the grade.

% For the 9-step prediction, which corresponds to the prediction almost an
% hour from now (how about doing 10 step since that is an hour?), the mean
% value for the last 4 hours is a suitable naive predictor.

halt_mean = movmean(halt_konc, [40 0]);

figure
plot(halt_konc)
hold on
plot(halt_mean)
legend('halt', 'mean')


%% Calculate the variance of the naive predictors prediction errors

%The indecies of the testset in the real dataset
ind = [round(start+length(halt_konc)*0.47) round(start+length(halt_konc)*0.47)+l_test];

% one step:
res_one_step = halt_konc(2:end)-halt_konc(1:end-1);

ACFnPACFnNormplot(res_one_step,50);
% Looks very white

% extract part corresponding to the test set: (one step lag in the res-vector)

