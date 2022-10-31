%
% Time series analysis
% Lund University
%
% Example code 7: transfer function example.
%
% Note that Matlab uses a slightly different notation for the BJ model as
% compared to that used in the course. In Matlab's notation:
%
%   A(z) y(t) = [B(z)/F(z)] u(t) + [C(z)/D(z)] e(t)
%
% This means that:
%
%   A(z) = 1,       B(z) = B(z),    F(z) = A2(z)
%   C(z) = C1(z),   D(z) = A1(z)
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;
rng(0)                      % Set the seed (just done for the lecture!)

extraN = 100;
N      = 1000;              % Try reducing the number of available samples; how many samples do you need to get a white residual even with the true model?
noLags = 30;

% Simulate some processes.
C1 = [ 1 0.4 0 0.6 ];
A1 = [ 1 -0.4 0.8 ];
C3 = [ 1 0 0 0.8 ]; 
A3 = [ 1 .9 ];
B  = [ 0 0 2 0 -0.8 ];
A2 = 1;

% Generate the noise and the input signals.
z = filter( C1, A1, randn(N+extraN, 1) );       % This is the noise model.
x = filter( C3, A3, randn(N+extraN, 1) );       % This is the input signal.

% Form the output using the filtered input. Remove the initial samples.
y = filter( B, A2, x ) + z;
y = y(extraN:end);
x = x(extraN:end);

% Examine the data.
figure; 
subplot(211); 
plot( x );
ylabel('Input signal')
title('Measured signals')
subplot(212); 
plot( y ); 
ylabel('Output signal')
xlabel('Time')

figure
[Cxy,lags] = xcorr( x, y, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(y) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between in- and output')


%% Begin with constructing a model for the input signal.
plotACFnPACF( x, noLags, 'Input' );


%% Estimate a first model for the input
% The AR part has a strong dependency at lag 1 and 3, lets begin with that. 
estimateARMA( x, [ 1 1 0 1 ], [ 1 ], 'Input model 1', noLags );


%% Estimate a second model for the input
% The MA part seems to have dependencies at lag 3 and 4. Lets add those.
estimateARMA( x, [ 1 1 0 1 1], [ 1 0 0 1 1], 'Input model 2', noLags );


%% Estimate a third model for the input
% It is now white, but notice that two AR coefficienta - and the last MA
% coefficient - do not seem to be significant, lets try to remove c4. Is
% the residual still white? 
estimateARMA( x, [ 1 1 0 1 1], [ 1 0 0 1  ], 'Input model 3', noLags );


%% Estimate a fourth model for the input
% Can we cut the two AR components too? Note how the resulting model
% compares with the model we used to simulate the input.
foundModel = estimateARMA( x, [ 1 1 ], [ 1 0 0 1  ], 'Input model 4', noLags );


%% Form the filtered signals and compute their cross-correlation.
% Word of warning; if you get odd result, plot ex and ey and check to there
% are no transient effects, such as ringdown from the filtering, still in
% the residual. If so, just cut a few more samples.
ex = filter( foundModel.A, foundModel.C, x );   ex = ex(length(foundModel.A):end );
ey = filter( foundModel.A, foundModel.C, y );   ey = ey(length(foundModel.A):end );

figure;
[Cxy,lags] = xcorr( ey, ex, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between filtered in- and output')
% Suggests a delay of d=2 and a filter of length 3 (but missing the second
% lag). This implies s=2 (recall d+s). There does not seem to be much of a
% decay, so lets try r=0.


%% Form a first model using the transfer model.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
estimateBJ( y, x, [1], [1], [0 0 1 0 1], [1], 'BJ model 1', noLags );
% There seems to be strong dependencies for the first 3 AR lag, lets add them.


%% Form a second model using the transfer model.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
estimateBJ( y, x, [1], [1 1 1 1], [0 0 1 0 1], [1], 'BJ model 2', noLags );
% Both B coefficients are statistically significant. There is some MA
% dependency at lag 3 and possibly 4. Lets try with just the third lag.


%% Form a third model using the transfer model.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
estimateBJ( y, x, [1 0 0 1], [1 1 1 1], [0 0 1 0 1], [1], 'BJ model 3', noLags );
% Seems that this was not enough - try adding the fourth as well.


%% Form a fourth model using the transfer model.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
foundModel = estimateBJ( y, x, [1 0 0 1 1], [1 1 1 1], [0 0 1 0 1], [1], 'BJ model 4', noLags );
% Yes, now it seems to be white!!! :-)


%% Check correlation of the resulting model.
% Ideally, the residual formed as tilde_et = yt - [ B / A2 ] xt should be
% uncorrelated with xt. Lets check! All seems fine! Compare the resulting
% model with how the signal was generated.
tilde_et = y - filter( foundModel.B, foundModel.F, x );      

% Note that we now have to remove samples from x as well.
tilde_et  = tilde_et(length(foundModel.B):end );
filter_xt = x(length(foundModel.B):end );

figure
[Cxy,lags] = xcorr( filter_xt, tilde_et, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(y) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between input and residual without the influence of the input')


%% What if we had used the correct model?
% We did not reach the model we used to simulate the data - how does our
% model compare to the true model? Reflections?
foundModel = estimateBJ( y, x, C1, A1, B, A2, 'True BJ model', noLags );

