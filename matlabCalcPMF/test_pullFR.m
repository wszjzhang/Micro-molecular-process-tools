% -----------
% test_pullFR
% -----------
%
% A demonstration of the "conjugate twin" method for combining
% bidirectional pulling data, as described in
% 
% D. Minh and A. Adib. Optimized free energies from bidirectional 
% single-molecule force spectroscopy.  
% Physical Review Letters 100, 180602 (2008).

%%% Simulation Parameters
steps = 750;                            % Number of simulation steps
iters = 500;                            % Number of trajectories

ks = 15;                                % Spring constant

bA = -1.5;                              % Bias position in state A
bB = 1.5;                               % Bias position in state B
lambdaF = linspace(bA,bB,steps);        % Forward pulling protocol
lambdaR = linspace(bB,bA,steps);        % Reverse pulling protocol

xspace = -1.5:0.05:1.5;                 % Histogram edges

%%% Analytical Results
[U,dUdX] = potential;                   % Original potential energy surface
V = @(x,b) ks/2*(x-b).^2;               % Harmonic biasing potential

PMF_exact = U(xspace,1)-min(U(xspace,1));

%%% Exact $\Delta F_t$
qoff = 5;
Ft_exact = zeros(1,steps);
for t = 1:steps
    Ft_exact(t) = -log(quadgk(@(x) exp(-(U(x,1)+V(x,lambdaF(t)))),...
        lambdaF(t)-qoff,lambdaF(t)+qoff));
end
Ft_exact = Ft_exact - Ft_exact(1);
deltaF_exact = Ft_exact(end)-Ft_exact(1);

%%% Simulation
[XF,WF] = pull(@(x) dUdX(x,1),ks,iters,lambdaF);
[XR,WR] = pull(@(x) dUdX(x,1),ks,iters,lambdaR);

%%% Unidirectional Analysis
[pmfF,FtF] = pmf_pullF(XF,WF,ks,lambdaF,xspace);
[pmfR,FtR] = pmf_pullF(XR,WR,ks,lambdaR,xspace);

%%% Removes half the data for a fair comparison
XF = XF(1:2:end,:);
XR = XR(1:2:end,:);
WF = WF(1:2:end,:);
WR = WR(1:2:end,:);

%%% Bidirectional Analysis
deltaF = BAR(1,WF(:,end),WR(:,end));    % Estimates $\Delta F$ using BAR
[pmfFRct,Ft_ct] = pmf_pullFR(XF,XR,WF,WR,ks,lambdaF,xspace,deltaF);
Ft_Chelli = Chelli(WF,WR);              % Estimates $\Detal F_t$ using Chelli's method

%%% Output

%%% Figure 1, $\Delta F_t$
FtF = FtF - FtF(1);
FtR = FtR - FtR(end) + deltaF;
Ft_ct = Ft_ct - Ft_ct(1);

time = (1:steps)*1E-3;
skip = 24;
tskip = 1:skip:steps;
tskip2 = skip/3:skip:steps;
tskip3 = skip*2/3:skip:steps;

figure(1)
clf
ax = axes;
hold on
plot(time,Ft_exact,'Color',[0.8 0.8 0.8],'LineWidth',2);
plot(time(tskip),FtF(:,tskip),'>g','MarkerSize',4);
plot(time(tskip2),FtR(:,tskip2),'<r','MarkerSize',4);
plot(time(tskip3),Ft_ct(:,tskip3),'^b','MarkerSize',4);
plot(time,Ft_Chelli,':m' ,'LineWidth',2);
hold off

xlabel('Time (s)')
ylabel('\Delta F_t (k_BT)')
set(ax,'XLim',[0 steps*1E-3]);

%%% Figure 2, PMF
offset = 5;
pmfF = pmfF - pmfF(1+offset) + PMF_exact(1+offset);
pmfR = pmfR - pmfR(end-offset) + PMF_exact(end-offset);
pmfFRct = pmfFRct - pmfFRct(1+offset) + PMF_exact(1+offset);

figure(2);
clf

xmin = -1.6;
xmax = 1.6;
ymin = -1;
ymax = 20;
textpos = 17.5;

pts = 1:2:length(xspace);
pts2 = 2:2:length(xspace);


posi = xspace(pts);
pmfminh = pmfFRct(:,pts);
posiv = posi';
pmfminhv = pmfminh';
save hummerminhpmfp.dat posiv -ASCII
save hummerminhpmfminh.dat pmfminhv -ASCII

xspacev = xspace';
PMF_exactv = PMF_exact';
save hummer0pmfp.dat xspacev -ASCII
save hummer0pmfminh.dat PMF_exactv -ASCII

pa = subplot('Position',[0.14 0.25 0.4 0.70]);
set(pa,'FontSize',8);
hold on
plot(xspace,PMF_exact,'Color',[0.8 0.8 0.8],'Linewidth',2)

plot(xspace(pts),pmfF(:,pts),'>g','MarkerSize',4);
plot(xspace(pts2),pmfR(:,pts2),'<r','MarkerSize',4);
axis([xmin xmax ymin ymax]);
ylabel('PMF (k_BT)')
text(-1.5,textpos,'(a) Unidirectional','FontSize',8)

pb = subplot('Position',[0.55 0.25 0.4 0.70]);
set(pb,'FontSize',8);
hold on

plot(xspace,PMF_exact,'Color',[0.8 0.8 0.8],'Linewidth',2)
plot(xspace(pts),pmfFRct(:,pts),'^b','MarkerSize',4);
axis([xmin xmax ymin ymax]);
set(pb,'YTick',[]);
xlabel('position');
text(-1.5,textpos,'(b) Bidirectional','FontSize',8)

