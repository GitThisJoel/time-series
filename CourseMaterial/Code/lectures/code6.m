%
% Time series analysis
% Lund University
%
% Example code 6: examining the powerload data from Helsinki (see also Ex 5.17).
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Load data and examine it.
load dataHelsinki;

figure
plot(powerload),          
title('Electricity consumption in Helsinki');


%% Should we transform the data?
% Note that we examine all the data when making this decision. 
figure; 
lambda_max = bcNormPlot(powerload,1);   % See also table 4.6
fprintf('The Box-Cox curve is maximized at %4.2f. This suggests that a log-transform might be helpful.\n', lambda_max)

figure
normplot( powerload)

% Lets transform the data. 
% The difference is small... perhaps we do not need this... Try without!
original  = powerload;
powerload = log( powerload );


%% Extract a modelling set.
first = 800;
last  = first + 168*4 - 1;              % We select 4 weeks. Does the data seem stable enough?
powerload = powerload(first:last);
figure
plot(powerload),		
title('log-transformed powerload');


%% Examine the ACF and PACF.
% Note the clear periodicities at 24 and 168 hours. Why is this?
noLags = 200;
plotACFnPACF( powerload, noLags, 'Data' );


%% Differentiate to remove the periodicities. Remember to remove the initial samples!
sweek = 168; 
sday  = 24;
dayPoly  = [1 zeros(1,sday-1) -1];
weekPoly = [1 zeros(1,sweek-1) -1];
powerload = filter(dayPoly,1,powerload);      powerload = powerload(sday+1:end);
powerload = filter(weekPoly,1,powerload);     powerload = powerload(sweek+1:end);

plotACFnPACF( powerload, noLags, 'Differentiated data' );


%% Make a first model
% We begin modelling the AR-part. Note that the first PACF lag is
% significant. Lets try using just that.
dataContainer = iddata( powerload );
Am = [ 1 1 ];                           % Specify assumed model orders.
Cm = [ 1 ];
foundModel = pem( dataContainer, idpoly( Am,[],Cm ) );
present( foundModel );                  % Note the confidence interval on the estimated parameter!

% Compute the residual. Remember to remove the initial samples.
ey = filter( foundModel.A, foundModel.C, powerload );  ey = ey(length(foundModel.A):end );

% Is the residual white?
figure
plot( ey )
title('Residual for model 1. Is it white?')


%% Lets look at the ACF and PACF.
plotACFnPACF( ey, noLags, 'Residual, model 1' )
checkIfWhite( ey );


%% Make a second model
% The residual is clearly not white. After includeing an AR-part, include
% also an MA-part. Note the strong 24-season in the ACF. Lets add that too.
Cm = conv( dayPoly, weekPoly );
Cm(2:end) = 0.3*Cm(2:end);              % Force the model to be stable. Try this if you get odd results.
polyContainer = idpoly( Am,[],Cm );
polyContainer.Structure.c.Free = Cm;    % Ensure that only the relevant coefficients are estimated. Try commenting this line. What happens?  
foundModel = pem( dataContainer, polyContainer );
present( foundModel );

ey = filter( foundModel.A, foundModel.C, powerload );  ey = ey(length(foundModel.A):end);
plotACFnPACF( ey, noLags, 'Residual, model 2' )


%% Is the resulting model residual white?
% It is very close... Perhaps this is a better model than the one below? 
checkIfWhite( ey );


%% Make a third model, adding a further periodicity also for the AR-part.
% Almost... there seems to be some periodicity at lag 11, lets add that...
Am = conv([1 1], [1 zeros(1,11-1) -1]);
polyContainer = idpoly( Am,[],Cm );
polyContainer.Structure.c.Free = Cm;
polyContainer.Structure.a.Free = Am;
foundModel = pem( dataContainer, polyContainer );
present( foundModel );

ey = filter( foundModel.A, foundModel.C, powerload );  ey = ey(length(foundModel.A):end);
plotACFnPACF( ey, noLags, 'Residual, model 3' )


%% Is the resulting model residual white now?
% Notice that although the residual is now white, the last coefficient in
% the AR polynomial is not significant. Try removing this - is the residual
% still white? Do we need to keep it?
checkIfWhite( ey );
