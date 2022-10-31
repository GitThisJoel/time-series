%
% Time series analysis
% Lund University
%
% Example code 4: examining the US Tobacco data (see also Ex 4.17).
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Plot the examined data. 
signLvl = 0.05;     noLags = 20;
dataTobacco;
data = data(:);

figure(1)
time = 1871:1984;
plot( time, data )
axis([time(1) time(end) 0 2500 ])
xlabel('Year')
ylabel('Production')
title('US tobacco production')

plotACFnPACF( data, noLags, 'Data' );    % We do this so often, lets add a function for it.


%% There seems to be a linear trend in the data. Estimate this and remove its influence.
N = length(data);
X = [ ones(N,1) (1:N)' ];
theV = inv( X'*X )*X'*data;             % This is the least-squares estimate of the trend. More on this later.
z = data - theV(1) - theV(2)*(1:N)';    % Subtract the estimated trend.

% Add the estimated trend to the plot.
figure(1)
hold on
plot( time, theV(1) + theV(2)*(1:N)', 'r')
legend('Data', 'Linear trend', 'Location', 'SE')
hold off

% Plot the resulting de-trended data and its ACF and PACF.
figure
plot( time, z )
xlim([time(1) time(end)])
xlabel('Year')
ylabel('Relative production')
title('US tobacco production, without trend')
plotACFnPACF( z, noLags, 'without trend' );


%% Try differentiating the data instead. Which seems to remove the trend best?
y = filter([1 -1], 1, data);
y = y(2:end);                           % Why should one remove the p first samples when filtering with a p:th order filter? Look at code5!

% Plot the resulting de-trended data and its ACF and PACF.
figure
plot( y )
title('US tobacco production, differentiated')
plotACFnPACF( y, noLags, 'differentiated' );


%% Try estimating a model for the differentiated data.
dataContainer = iddata( y );            % Create an iddata struct for the data. Try using z instead; does this give a better model?
Am = [ 1 1 ];                           % Specify assumed model orders.
Cm = [ 1 ];

% You can also use the non-differentiated data. How does the resulting polynomials compare?
%dataContainer = iddata( data );
%Am = conv([ 1 1 ], [1 -1] );

% Estimate unknown parameters using PEM (more on this later).
foundModel = pem( dataContainer, idpoly( Am,[],Cm ) );                   
present( foundModel );                  % Present the found model.

% Compute the residual. Remember to remove the initial samples.
ey = filter( foundModel.A, 1, y );  ey = ey(length(foundModel.A):end );

% Plot the ACF and PACF.
figure
plot( ey )
title('Model residual. Is it white?')
plotACFnPACF( ey, noLags, 'Residual' );

% With the differentiation, the found A-polynomial is:
conv(foundModel.A, [1 -1])

% The variance of the residual is a measure for how well the model fits the
% data. The variance of the data is computed as a sanity check.
fprintf('The variance of the model residual is %2.3f.\nThe variance of the data is %2.3f.\n\n', var(ey), var(data) )


%% Does the data seem to be white? More on this later.
% Is the residual white? Lets examine the Monti test.
[deemedWhite, Q, chiV] = montiTest( ey );
if deemedWhite
     fprintf('The residual is deemed to be white according to the Monti-test (as %5.2f < %5.2f).\n', Q, chiV );
else
    fprintf('The residual is not deemed to be white according to the Monti-test (as %5.2f > %5.2f).\n', Q, chiV );
end

% Does the residual have a mean that is different from zero?
[ rejectH0, tRatio, tLimit ] = testMean( ey );
if rejectH0
    fprintf('The data is not deemed to be zero-mean (as %6.4f > %6.4f). \n', tRatio, tLimit )
else
    fprintf('The data is deemed to be zero-mean (as %6.4f < %6.4f). \n', tRatio, tLimit )
end
