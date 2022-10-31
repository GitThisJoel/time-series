% FUNCTION [ foundModel, ey, acfEst, pacfEst ] = estimateBJ( y, x, B, A2, C1, A1, titleStr, noLags )
%
% The function estimates the Box-Jenkins model
% 
%   y(t) = [ B(z) / A2(z)] x(t) + [ C1(z)/A1(z) ] e(t)
%
% and returns the found model and the residual.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%

% Note that the function idpoly use the notation: 
%
%       A(z) y(t) = [B(z)/F(z)] u(t) + [C(z)/D(z)] e(t)
%
% This means that:
%
%   A(z) = 1,       B(z) = B(z),    F(z) = A2(z)
%   C(z) = C1(z),   D(z) = A1(z)
%
function [ foundModel, ey, acfEst, pacfEst ] = estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )

% It is not uncommon that PEM suffers instabilities for complicated model
% structures, getting stuck in a local minima. An often useful trick is to
% initiate the polynomials with smaller coefficients. If you gets very odd
% estimates, try reducing the multiple a bit further - if you are lucky, it
% might help!
A1(2:end) = A1(2:end) * 0.3;
A2(2:end) = A2(2:end) * 0.3;
C1(2:end) = C1(2:end) * 0.3;
B(2:end)  = B(2:end)  * 0.3;

polyContainer = idpoly( 1, B, C1, A1, A2 );
dataContainer = iddata( y, x );
if length(B)>1                              % Only estimate the specified parameters.
    polyContainer.Structure.B.Free = B;
end
if length(A2)>1
    polyContainer.Structure.F.Free = A2;
end
if length(C1)>1
    polyContainer.Structure.C.Free = C1;
end
if length(A1)>1
    polyContainer.Structure.D.Free = A1;
end
foundModel = pem( dataContainer, polyContainer );
present( foundModel );

% Form the residual using all the polynomials - and remove the initial samples. 
ey = resid( foundModel, dataContainer );
ey = ey.y( length(foundModel.B):end );
fprintf('\nThe variance of the residual is %5.3f.\n', var(ey) )
 
% Plot the ACF and PACF of the residual.
[acfEst, pacfEst] = plotACFnPACF( ey, noLags, titleStr );
 
% Is the residual white? Lets examine the Monti test.
checkIfWhite( ey );
