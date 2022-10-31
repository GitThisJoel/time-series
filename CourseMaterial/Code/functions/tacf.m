% FUNCTION rho = tacf( y, maxOrd, [alpha], [signLvl], [plotIt], [includeZeroLag] )
%
% The function estimates the alpha-trimmed ACF of y for lags up to maxOrd. 
% 
% As a guideline, alpha should be set to 1-2% for general time series, 3-5% 
% for medium contaminated series, and 6-10% for heavily contaminated series. 
% Default is 2%.
%
% If the optional parameter plotIt is given the ACF is plotted together with
% a (Gaussian) confidence interval with significance signLvl (default 0.05), 
% assuming that y can be modeled as a white noise. If includeZeroLag is set 
% (default), the plot is from lag 0 to maxOrd, otherwise from lag 1.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function rho = tacf( y, maxOrd, alpha, signLvl, plotIt, includeZeroLag )

if nargin<3
    alpha = 0.02;
end
if nargin<4
    signLvl = 0.05;
end
if nargin<5
    plotIt = 0;
end
if nargin<6
    includeZeroLag = 1;
end
if (signLvl < 0) || (signLvl>1 )
    error('TACF: not a valid level of significance.');
end
signScale = norminv( 1-signLvl/2, 0, 1 );

N  = length(y);
ys = sort( y );

% Form trimmed data.
g = N*alpha;
indY = y < ys( floor(g)+1 ) ;
y(indY) = 0;
indY = y > ys( floor(N-g+1)) ;
y(indY) = 0;

% Find kept indices. 
La =  (y > ys( floor(g)+1 )) & (y < ys( floor(N-g+1))) ;
Lg = zeros(N,1);
Lg(La) = 1;

% Form trimmed mean estimate.
my = sum(y) / sum(Lg);

% Estimate TACF.
rho = xcorr( Lg.*(y-my), maxOrd, 'none' );
g2 = xcorr( Lg, maxOrd, 'none' );
rho = rho(maxOrd+1:end) ./ g2(maxOrd+1:end);
rho = rho / rho(1);

% Display results.
if plotIt
    if includeZeroLag
        rangeLags = 0:maxOrd;
        maxRange  = 1.1; 
        startLag  = 0;
        g2 = g2(maxOrd+1:end);
    else
        rangeLags = 1:maxOrd;
        rho = rho(2:end);
        g2 = g2(maxOrd+2:end);
        maxRange = max( abs(rho) )*1.2;
        startLag  = 1;
    end
    stem(rangeLags,rho )
    xlabel('Lag')
    ylabel('Amplitude')
    condInt = signScale * ones(length(rangeLags),1) ./ sqrt( g2 );
    hold on
    plot( rangeLags, condInt,'--' )
    plot( rangeLags, -condInt,'--' )
    hold off
    axis([startLag length(rangeLags)-1 -maxRange maxRange ])
end
