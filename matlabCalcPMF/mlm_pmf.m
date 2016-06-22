% -----------------
% mlm calculate pmf
% -----------------
%
% A demonstration of the "conjugate twin" method for combining
% bidirectionavl pulling data, as described in
% 
% D. Minh and A. Adib. Optimized free energies from bidirectional 
% single-molecule force spectroscopy.  s
% Physical Review Letters 100, 180602 (2008).

%%% Simulation Parameters

ks = 1000;                                % Spring constant
bA = -1.5;                              % Bias position in state A
bB = 1.5;                               % Bias position in state B

iters = length(dir('F-*'));                            % Number of trajectories
traj_pars = dir('F-*')
% Read parameters from trajectory file

steps = length(importdata(traj_pars(1).name));  % Number of simulation steps
traj_in = importdata(traj_pars(1).name);
dt = traj_in(2,1) - traj_in(1,1);



lambdaF = linspace(bA,bB,steps);        % Forward pulling protocol
lambdaR = linspace(bB,bA,steps);        % Reverse pulling protocol
v = (bB - bA)/(steps*dt);               % pulling speed
xspace = bA:v*dt:bB;                    % Histogram edges

%%% Simulation
[XF,WF] = loadtrajmlmf(ks,lambdaF);
[XR,WR] = loadtrajmlmr(ks,lambdaR);

%%% Bidirectional Analysis
deltaF = BAR(1,WF(:,end),WR(:,end));    % Estimates $\Delta F$ using BAR
[pmfFRct,Ft_ct] = pmf_pullFR(XF,XR,WF,WR,ks,lambdaF,xspace,deltaF);   %Core to calculate PMF using MLM
Ft_Chelli = Chelli(WF,WR);              % Estimates $\Detal F_t$ using Chellis method

pmfFRct = pmfFRct - pmfFRct(1);
%%% Figure 3, MLE PMF
%offset = 5;
%pmfF = pmfF - pmfF(1+offset) + PMF_exact(1+offset);
%pmfR = pmfR - pmfR(end-offset) + PMF_exact(end-offset);
%pmfFRct = pmfFRct - pmfFRct(1+offset) + PMF_exact(1+offset);

figure(3);
clf

textpos = 17.5;

pts = 2:1:length(xspace);
pts2 = 2:2:length(xspace);

%save pmf to file

pmf_final = [lambdaF' pmfFRct(:,pts)']
mkdir ../data
save ../data/mlm_pmf.dat pmf_final -ASCII

%xspacev = xspace';
%PMF_exactv = PMF_exact';
%save hummer0pmfp.dat xspacev -ASCII
%save hummer0pmfminh.dat PMF_exactv -ASCII
hold on
plot(lambdaF, pmfFRct(:,pts))
ylabel('U(x)');
xlabel('x');
text(-1.5,textpos,'MLE PMF','FontSize',8)

