% FUNCTION [ maxLam, bct, offsetValue ] = bcNormPlot( y, plotIt, lamRange )
%
% The function computes the Box-Cox normality plot, finding the lambda
% value withing lamRange (default [-2,2]) that maximize the plot. If plotIt
% is set (default), the Box-Cox normality curve is plotted.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ maxLam, bct, offsetValue ] = bcNormPlot( y, plotIt, lamRange )

if nargin<2
    plotIt = 1;
end
if nargin<3
    lamRange = linspace( -2, 2, 100 );
end

offsetValue = 0;
if min(y)<=0
    offsetValue = 1 - min(y);
    y = y + offsetValue;
end

N = length(y);
for k=1:length(lamRange)
    
    lambda = lamRange(k);
    if lambda == 0
        z = log( y );
    else
        z = ( y.^lambda - 1 )/lambda;
    end
    
    bct(k) = -(N/2).*log(std(z).^2) + (lambda-1)*(sum(log(y)));
end

[ ~, indC ] = max( bct );
maxLam = lamRange( indC );

if plotIt
    plot( lamRange, bct )
    xlabel('Lambda')
    ylabel('Log-likelihood')
end
