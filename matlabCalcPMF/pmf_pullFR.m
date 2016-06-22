% ------------
% pmf_pullFRct
% ------------
%
% Calculates the PMF from a set of bidirectional pulling experiments 
% using the method of Hummer and Szabo with conjugate twins, as described
% in:
%
% D. Minh and A. Adib. Optimized free energies from bidirectional 
% single-molecule force spectroscopy.  
% Physical Review Letters 100, 180602 (2008).
%
%
% [pmfFR,deltaF] =  pmf_pullF(XF,XR,WF,WR,ks,iters,lambda,lambdaR,xspace);
%
% Parameters
% XF:      a AxS matrix containing X coordinate values,
%          where A is the iteration and S is the step (nm)
% XR:      a AxS matrix containing X coordinate values,
%          where A is the iteration and S is the step (nm)
% WF:      a AxS matrix containing the accumulated work values
% WR:      a AxS matrix containing the accumulated work values
% ks:      the bias spring constant
% lambda: a 1xS matrix containing the schedule of the 
%          bias centers where S is the step
% lambdaR: a 1xS matrix containing the schedule of the 
%          bias centers where S is the step
% xspace:  a 1xB vector marking the centers of the bins along X,
%          where B is the number of bins (nm)
%
function [pmfFRct,Ftct] =  pmf_pullFRct(XF,XR,WF,WR,...
  ks,lambda,xspace,deltaF)

xsp = xspace(2)-xspace(1);

iters = size(XF,1);
steps = size(XF,2);

NF = size(XF,1);
NR = size(XR,1);

reverse = linspace(steps,1,steps);

%%% Forward time slices with reverse congugate twins
%%% Calculate conjugate twin position and work
XRct = XR(:,reverse);
WRi = WR(:,2:end)-WR(:,1:end-1);
WRct = [zeros(iters,1) -cumsum(WRi(:,reverse(2:end)),2)];

%%% Optimal weights for trajectories
wFft = NF./(NF+NR*exp(-(WF(:,end)-deltaF)))*ones(1,steps);
wRct = NR./(NF+NR*exp(WR(:,end)+deltaF))*ones(1,steps);

%%% Time slices
eWminusVFft = exp(-(WF-(ks/2*(XF-ones(iters,1)*lambda).^2)));
eWminusVRct = exp(-(WRct-(ks/2*(XRct-ones(iters,1)*lambda).^2)));
timesliceF = histw(XF,wFft.*eWminusVFft,xspace-xsp/2);
timesliceR = histw(XRct,wRct.*eWminusVRct,xspace-xsp/2);
timeslice = timesliceF/NF + timesliceR/NR;

Vt = ks/2*((xspace'*ones(1,steps))...
     -(ones(length(xspace),1)*lambda)).^2;

%%% Free energies for states
Ftct = -log(sum(wFft.*exp(-WF))/NF + sum(wRct.*exp(-WRct))/NR);

%%% WHAM weights
wF = exp(-(Vt-ones(length(xspace),1)*Ftct));

%%% PMF using conjugate twin estimator for Ft
rho = sum(timeslice.*wF,2)./sum(wF,2);
rho = rho/sum(rho);
pmfFRct = -log(rho');