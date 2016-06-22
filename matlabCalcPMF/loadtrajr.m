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

    function [Xs,Ws] = loadtrajr(ks,v,lambda,dt)

%%% Simulation parameters
% KbT = 1;
% D = 1;

%%% Load trajectory files 
trajfs = dir('R-*');
iters = length(dir('R-*'));

steps = length(lambda);

%dt = 0.001;               %%% timestep size, seconds
equil_steps = 100;

%%% Collected data for all the simulations
Xs = zeros(iters,steps);    %%% positions
Ws = zeros(iters,steps);    %%% accumulated work

for n = 1:iters
	  traj = importdata(trajfs(n).name);   %%%load trajectory(tt,xt) to traj
%	  X = zeros(steps,1);    %%% positions vt
%	  W = zeros(steps,1);    %%% accumulated work
    wi = 0;
 %calculate the total work in the time reverse process
  for t = 1:steps
            wi = wi + ks*v*(lambda(1)-traj(t,2))*dt;
  end
  wrt = 0.5*ks*v*v*traj(end,1)*traj(end,1) - wi;

 %calculate the work in process
  wi = 0;
  for t = 1:steps
	    X = lambda(steps-t+1);
            wi = wi + ks*v*(lambda(1)-traj(t,2))*dt; 
            W = wrt - (0.5*ks*v*v*traj(t,1)*traj(t,1) - wi);

   Ws(n,steps-t+1) = W;
   Xs(n,steps-t+1) = X;
  end 
 
end
