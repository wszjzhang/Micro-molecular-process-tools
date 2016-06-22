% -----------------
% FR method for PMF
% -----------------
%
% A implemetation of the FR method for calculating PMF as described in
% 
% Calculating potentials of mean force and diffusion coefficients
% from nonequilibrium processes without Jarzynski’s equality
%v
% Ioan Kosztin, Bogdan Barz, and Lorant Janosi   
% The Journal of Chemical Physics 124, 064106 (2006) 

%%% Simulation Parameters
ks = 1000;                                % Spring constant
bA = -1.5;                              % Bias position in state A
bB = 1.5;                               % Bias position in state B

iters = length(dir('F-*'));                            % Number of trajectories
traj_pars = dir('F-*');
% Read parameters from trajectory file

steps = length(importdata(traj_pars(1).name));  % Number of simulation steps
traj_in = importdata(traj_pars(1).name);
dt = traj_in(2,1) - traj_in(1,1);

lambdaF = linspace(bA,bB,steps);        % Forward pulling protocol
lambdaR = linspace(bB,bA,steps);        % Reverse pulling protocol
v = (bB - bA)/(steps*dt);               % pulling speed
xspace = bA:v*dt:bB;                 % Histogram edges

%%% Simulation
[XF,WF] = loadtrajf(ks,v,lambdaF,dt);
[XR,WR] = loadtrajr(ks,v,lambdaR,dt);

%%% Bidirectional Analysis

pmfFRctf = mean(WF)-(mean(WF.^2) - mean(WF).^2)/2;   %Core to calculate PMF using JE 
                                                     %cumulant  approximation
pmfFRctf = pmfFRctf - pmfFRctf(1);

pmfFRctr = mean(WR)-(mean(WR.^2) - mean(WR).^2)/2;   %Core to calculate PMF using JE 
                                                     %cumulant  approximation
pmfFRctr = pmfFRctr - pmfFRctr(1);


%%% Figure 1, JE PMF
figure(1);
clf

textpos = 17.5;

pts = 1:2:length(xspace);
pts2 = 2:2:length(xspace);

%save pmf to file
mkdir data
pmf_finalf = [lambdaF' pmfFRctf'];
save data/jef_pmf.dat pmf_finalf -ASCII;
pmf_finalr = [lambdaF' -pmfFRctr'];
save data/jer_pmf.dat pmf_finalr -ASCII;

hold on
plot(lambdaF, pmfFRctf)
plot(lambdaF, -pmfFRctr, 'r');

ylabel('U(x)');
xlabel('x');
text(-1.5,textpos,'JE cumulant approximation','FontSize',8)

