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
function [Xs,Ws] = loadtrajmlmr(ks,lambda)

%%% Simulation parameters
% KbT = 1;
% D = 1;
trajfs = dir('R-*');
iters = length(dir('R-*'));

steps = length(lambda);

%dt = 0.001;               %%% timestep size, seconds


%%% Collected data for all the simulations
Xs = zeros(iters,steps);    %%% positions
Ws = zeros(iters,steps);    %%% accumulated work

for n = 1:iters
  traj = importdata(trajfs(n).name);   %%%load trajectory(tt,xt) to traj
  X = traj(:,2);    %%% positions
  W = zeros(steps,1);    %%% accumulated work
  

  %%% Brownian Dynamics Engine
  for t = 2:steps
    W(t) = W(t-1) + ks/2*(X(t)-lambda(t)).^2 - ks/2*(X(t)-lambda(t-1)).^2;
  end

  Xs(n,:) = X;
  Ws(n,:) = W;
end