close all
clear all

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
halt_ing_rep = (halt_raw_konc.*flow_raw_remill + halt_scav_konc.*flow_scav_remill) ./(flow_raw_remill + flow_scav_remill);

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