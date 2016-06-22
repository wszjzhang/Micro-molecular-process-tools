% -----------------
% FR method for PMF
% -----------------
%
% A implemetation of the FR method for calculating PMF as described in
% 
% Calculating potentials of mean force and diffusion coefficients
% from nonequilibrium processes without Jarzynski’s equality
%
% Ioan Kosztin, Bogdan Barz, and Lorant Janosi   
% The Journal of Chemical Physics 124, 064106 (2006) 

%%% Simulation Parameters
ks = 1000;                              % Spring constant
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

%save work distribution to file
workf = [lambdaF' WF'];
workr = [lambdaR' WR'];
mkdir data
save data/work_f.dat workf -ASCII
save data/work_r.dat workr -ASCII

%%% Bidirectional Analysis

pmfFRct = (mean(WF) - mean(WR))/2;   %Core to calculate PMF using MLM


%%% Figure 2, FR PMF
figure(2);
clf

textpos = 17.5;

pts = 1:1:length(xspace)-1;
pts2 = 2:2:length(xspace);

%save pmf to file

pmf_final = [lambdaF' pmfFRct'-pmfFRct(1)]
mkdir data
save data/fr_pmf.dat pmf_final -ASCII

%xspacev = xspace';
%PMF_exactv = PMF_exact';
%save hummer0pmfp.dat xspacev -ASCII
%save hummer0pmfminh.dat PMF_exactv -ASCII

hold on
plot(xspace(pts),pmfFRct(:,pts),'^b','MarkerSize',4);

ylabel('U(x)');
xlabel('x');
text(-1.5,textpos,'FR method PMF','FontSize',8)

