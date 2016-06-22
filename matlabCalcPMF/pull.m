% ----
% pull
% ----
%
% Runs pulling simulations along 
% a one-dimensional surface specified by the 
% potential function.
%
% Inputs
% dUdX:    the derivative of the potential as a function
%          of position
% ks:      the bias spring constant
% iters:   the number of iterations
% lambda:  a 1xS matrix containing the schedule of the 
%          bias centers where S is the step
%
% Outputs
% Xs:      a AxS matrix containing X coordinate values,
% Ws:      a AxS matrix containing the accumulated work values
%
function [Xs,Ws] = pull(dUdX,ks,iters,lambda)

%%% Simulation parameters
% KbT = 1;
% D = 1;

steps = length(lambda);

dt = 1E-3;               %%% timestep size, seconds
equil_steps = 100;

%%% Collected data for all the simulations
Xs = zeros(iters,steps);    %%% positions
Ws = zeros(iters,steps);    %%% accumulated work

for n = 1:iters
  X = zeros(steps,1);    %%% positions
  W = zeros(steps,1);    %%% accumulated work
  
  %%% Equilibration
  if n>1
     X(1) = Xs(n-1,1);
  else
     X(1) = lambda(1);
  end
  R1 = randn(equil_steps,1)*sqrt(2*dt);   %%% random variables, variance 2D(dt)
  for t = 2:equil_steps
    dx = R1(t) + dt*(-ks*(X(1)-lambda(1)) - dUdX(X(1)));
    X(1) = X(1) + dx;
  end
  
  R1 = randn(steps,1)*sqrt(2*dt);   %%% random variables, variance 2D(dt)
  %%% Brownian Dynamics Engine
  for t = 2:steps
    dx = R1(t) + dt*(-ks*(X(t-1)-lambda(t-1)) - dUdX(X(t-1)));
    X(t) = X(t-1) + dx;
    W(t) = W(t-1) + ks/2*(X(t)-lambda(t)).^2 - ks/2*(X(t)-lambda(t-1)).^2;
  end

  Xs(n,:) = X;
  Ws(n,:) = W;
end