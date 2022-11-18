% FUNCTION [ foundModel, ey, acfEst, pacfEst ] = estimateARMA( data, Am, Cm, titleStr, noLags )
%
% The function estimates the specified ARMA model (Am, Cm) for data and
% returns the found model and the residual.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ foundModel, ey, acfEst, pacfEst ] = estimateARMA( data, Am, Cm, titleStr, noLags )

% It is not uncommon that PEM suffers instabilities for complicated model
% structures, getting stuck in a local minima. An often useful trick is to
% initiate the polynomials with smaller coefficients. If you gets very odd
% estimates, try reducing the multiple a bit further - if you are lucky, it
% might help!
Am(2:end) = Am(2:end) * 0.3;
Cm(2:end) = Cm(2:end) * 0.3;

dataContainer = iddata( data );
polyContainer = idpoly( Am,[],Cm );
if length(Cm)>1
    polyContainer.Structure.c.Free = Cm;        % Only estimate the specified C(z) parameters.
end
if length(Am)>1
    polyContainer.Structure.a.Free = Am;        % Only estimate the specified A(z) parameters.
end
foundModel = pem( dataContainer, polyContainer );
present( foundModel );

% Compute the residual and its variance.
ey = filter( foundModel.A, foundModel.C, data );  ey = ey(length(foundModel.A):end );
fprintf('\nThe variance of the residual is %5.3f.\n', var(ey) )

% Plot the ACF and PACF of the residual.
[acfEst, pacfEst] = plotACFnPACF( ey, noLags, titleStr );

% Is the residual white? Lets examine the Monti test.
checkIfNormal( ey, "" ); close; % closes the window that checkIfNormal produces... 
checkIfWhite( ey );
