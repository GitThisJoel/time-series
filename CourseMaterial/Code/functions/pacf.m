% FUNCTION phi = pacf( y, maxOrd, [signLvl], [plotIt], [includeZeroLag] )
%
% The function estimates the PACF of y for lags up to maxOrd. If the
% optional parameter plotIt is given the PACF is plotted together with a
% (Gaussian) confidence interval with significance signLvl (default 0.05).
% If includeZeroLag is set (default), the plot is from lag 0 to maxOrd, 
% otherwise from lag 1.

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function phi = pacf( y, maxOrd, signLvl, plotIt, includeZeroLag )

if nargin<3
    signLvl = 0.05;
end
if nargin<4
    plotIt = 0;
end
if (signLvl < 0) || (signLvl>1 )
    error('PACF: not a valid level of significance.');
end
signScale = norminv( 1-signLvl/2, 0, 1 );
if nargin<5
    includeZeroLag = 1;
end

y = y - mean(y);
phi = zeros( maxOrd+1, 1 );
for k=0:maxOrd
    a = aryule( y, k );
    phi(k+1) = a(end);
end

if plotIt
    if includeZeroLag
        rangeLags = 0:maxOrd;
        maxRange  = 1.1; 
        startLag  = 0;
    else
        rangeLags = 1:maxOrd;
        phi = phi(2:end);
        maxRange = max( abs(phi) )*1.2;
        startLag  = 1;
    end  
    stem(rangeLags,phi )
    xlabel('Lag')
    ylabel('Amplitude')
    condInt = signScale*ones(length(rangeLags),1)/sqrt(length(y));
    hold on
    plot( rangeLags, condInt,'r--' )
    plot( rangeLags, -condInt,'r--' )
    hold off
    axis([startLag length(rangeLags)-1 -maxRange maxRange ])
end
