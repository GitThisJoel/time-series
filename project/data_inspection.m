close all

load('raw_data.mat')

time = (1:length(raw_data(:,1)))* 6*60;   % Time in [s], concider minutes instead 

konc = raw_data(:,1);
ing_malm = 0.5* (raw_data(:,2) + raw_data(:,3));
ra_konc = 0.5* (raw_data(:,4) + raw_data(:,5));
scav_konc = 0.5* (raw_data(:,6) + raw_data(:,7));


time = time/60/60/24;
figure
subplot(3,1,1)
plot(time,konc)
title('Concentration in final concentarte (y)')
subplot(3,1,2)
plot(time, ing_malm)
title('Concentration of incoming ore (x)')
subplot(3,1,3)
plot(time, ra_konc)
title('Concentration of middle product, (alternative x)')
xlabel('Time [days]')