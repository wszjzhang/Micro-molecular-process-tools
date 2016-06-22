% ----
% loadtraj
% ----
%
% load trajectories and calculate accumulate work
% 
%
% Inputs
% trajecories: 
% xt:      the x position of particle
% tt:      the time
% parameters:
% ks:      the bias spring constant
% iters:   the number of iterations
% lambda:  a 1xS matrix containing the schedule of the 
%          bias centers where S is the step
% v:       pulling speed
%
% Outputs
% Xs:      a AxS matrix containing X coordinate values,
% Ws:      a AxS matrix containing the accumulated work values
%

    function [Xs,Ws] = loadtrajf(ks,v,lambda,dt)

%%% Simulation parameters
% KbT = 1;
% D = 1;

%%% Load trajectory files 
trajfs = dir('F-*');
iters = length(dir('F-*'));

steps = length(lambda);

%dt = 0.001;               %%% timestep size, seconds
%equil_steps = 100;

%%% Collected data for all the simulations
Xs = zeros(iters,steps);    %%% positions
Ws = zeros(iters,steps);    %%% accumulated work

for n = 1:iters

	  traj = importdata(trajfs(n).name);   %%%load trajectory(tt,xt) to traj
%	  X = zeros(steps,1);    %%% positions vt
%	  W = zeros(steps,1);    %%% accumulated work
    wi = 0;
  for t = 1:steps
	    X = traj(t,2);
            wi = wi + ks*v*(traj(t,2)-lambda(1))*dt; 
            W = 0.5*ks*v*v*traj(t,1)*traj(t,1) - wi;
  
  Xs(n,t) = X;
  Ws(n,t) = W;

  end

end
