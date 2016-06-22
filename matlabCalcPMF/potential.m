% ---------
% potential
% ---------
%
% Returns functions for the potential:
% BV(x,lambda) = 5*(x.^2 - 1).^2 + 6*(lambda - 1/2).*x
% which appears in 
% "Free Energy Calculations: Theory and Applications in Chemistry 
% and Biology" by Chipot and Pohorille (eds) on page 188.
%
% Outputs
% U:       the potential as a function of position 
%          and reaction coordinate lambda
% dUdX:    the derivative of the potential with respect to position
% P:       the normalized probability with respect to position
function [U,dUdX,P] = potential

U =    @(x,lambda) 5*(x.^2 - 1).^2 + 6*(lambda - 1/2).*x;
dUdX = @(x,lambda) 20*x.*(x.^2 - 1) + 6*(lambda - 1/2);

P = @(x,lambda) exp(-U(x,lambda));
sP = quadl(@(x)P(x,0),-5,5);
P = @(x,lambda) exp(-U(x,lambda))/sP;