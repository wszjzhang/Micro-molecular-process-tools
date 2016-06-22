% ---------
% pmf_pullF
% ---------
%
% Calculates the PMF from a set of unidirectional pulling experiments 
% using the method of Hummer and Szabo.
%
% G. Hummer and A. Szabo. Free energy reconstruction from nonequilibrium
% single-molecule pulling experiments. Proceedings of the National Academy
% of Sciences U.S.A. 98, 3658-3661 (2001). 
%
% [pmfF,Ft] =  pmf_pullF(Xs,Ws,ks,iters,lambda,xspace)
%
% Parameters
% Xs:      a AxS matrix containing X coordinate values,
%          where A is the iteration and S is the step (nm)
% Ws:      a AxS matrix containing the accumulated work values
% ks:      the bias spring constant
% lambda:  a 1xS matrix containing the schedule of the 
%          bias centers where S is the step
% xspace:  a 1xB vector marking the centers of the bins along X,
%          where B is the number of bins (nm)
%
% Outputs
% pmfF:    Potential of mean force estimate
% Ft:      Free energy estimate (Jarzynski's equality)
%
function [pmfF,Ft] =  pmf_pullF(Xs,Ws,ks,lambda,xspace)

xsp = xspace(2)-xspace(1);

iters = size(Xs,1);
steps = size(Xs,2);

%%% Weights for each time slice
eWminusVt = exp(-(Ws-(ks/2*(Xs-ones(iters,1)*lambda).^2)));
%%% Time slices
timeslice = histw(Xs,eWminusVt,xspace-xsp/2);

Ft = -log(mean(exp(-Ws)));

V = ks/2*((xspace'*ones(1,steps))...
     -(ones(length(xspace),1)*lambda)).^2;

w = exp(-(V-ones(length(xspace),1)*Ft));
rw = sum(timeslice.*w,2)./sum(w,2);
pmfF = -log(rw');