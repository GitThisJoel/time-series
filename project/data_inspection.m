close all
clear all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

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
subplot(3,1,1)
plot(time,halt_konc)
title('Concentration in final concentarte (y)')
subplot(3,1,2)
plot(time, halt_ing_flot)
title('Concentration of incoming ore (x)')
subplot(3,1,3)
plot(time, halt_ing_rep)
title('Concentration of middle product, (alternative x)')
xlabel('Time [days]')

%% Vi behöver hantera outliers, finns gott om dom...
% spolvatten verkar vara rätt konstant så går nog bra att bara kvotera in
% om vi behöver det, annars kanske försumbart

% De signaler vi kommer jobba med är halten i koncentratet "halt_konc" och
% som insignal kommer vi sedan att använda "halt_ing_rep" och eventuellt 
% "flow_ing_rep" som insignaler.

%% A. modeling without an external input
close all
figure 
plot(halt_konc)

% för att hålla ungefär samma datastorlek som i orginalprojektet bör
% modeleringsdatat innehålla ca 1700 samples, this corresponds to roughly 7
% days with our sampling intervall. 
start = 6800;
modeling_set = halt_konc(start:start+1700);
validation_set = halt_konc(start+1701:start+1710+350 );
test_set = halt_konc(start+length(halt_konc)*0.5:start+length(halt_konc)*0.5+350);

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
xline(start+1700)
xline(start+1700+350)
xline(start+length(halt_konc)*0.5)
xline(start+length(halt_konc)*0.5+350)

% I think this split looks ok, lets see if the data-size is sufficent

%% Utreda om outliers behöver behandlas
figure
acf(modeling_set, 32, 0.05, 1);
figure
tacf(modeling_set, 32, 0.04, 0.05, 1);

% Relativt lika, testar utan out-lier-behandlign till att börja med

