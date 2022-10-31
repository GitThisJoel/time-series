% FUNCTION rho = acf( y, maxOrd, [signLvl], [plotIt], [maOrder], [includeZeroLag] )
%
% The function estimates the ACF of y for lags up to maxOrd. If the
% optional parameter plotIt is given the ACF is plotted together with a
% (Gaussian) confidence interval with significance signLvl (default 0.05), 
% assuming that y can be modeled as an MA(maOrder) process. 
% 
% If maOrder is unspecified, the process is assumed to be white (i.e., it 
% is assumed that maOrder = 0). If includeZeroLag is set (default), the plot 
% is from lag 0 to maxOrd, otherwise from lag 1.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function rho = acf( y, maxOrd, signLvl, plotIt, maOrder, includeZeroLag )


if nargin<3
    signLvl = 0.05;
end
if nargin<4
    plotIt = 0;
end
if nargin<5
    maOrder = 0;
end
if nargin<6,
    includeZeroLag = 1;
end
if maOrder>maxOrd
    error('ACF: call with maOrder larger than maxOrd.');
end     
if (signLvl < 0) | (signLvl>1 )
    error('ACF: not a valid level of significance.');
end
signScale = norminv( 1-signLvl/2, 0, 1 );

rho = xcorr( y-mean(y), maxOrd, 'biased' );
rho = rho/max(rho);
rho = rho(maxOrd+1:end);

if plotIt
    if includeZeroLag
        rangeLags = 0:maxOrd;
        maxRange  = 1.1; 
        startLag  = 0;
    else
        rangeLags = 1:maxOrd;
        rho = rho(2:end);
        maxRange = max( abs(rho) )*1.2;
        startLag  = 1;
    end
    stem(rangeLags,rho )
    xlabel('Lag')
    ylabel('Amplitude')
    condInt = signScale * ones(length(rangeLags),1)/sqrt(length(y));
    condInt = condInt * sqrt( 1 + 2*sum( rho(1:maOrder).^2 ) );
    hold on
    plot( rangeLags, condInt,'r--' )
    plot( rangeLags, -condInt,'r--' )
    hold off
    axis([startLag length(rangeLags)-1 -maxRange maxRange ])
end
