% ------
% Chelli
% ------
% 
% Calculates the free energy difference as a function of time for a
% nonequilibrium switching process between two states using the method of
% Chelli et. al.  This is solved by finding the zero of the implicit
% formula.
%
% R. Chelli, S. Marsili, and P. Procacci. Calculation of the potential of
% mean force from nonequilibrium measurements via maximum likelihood 
% estimators. Physical Review E 77, 031104 (2008).
%
% Here, we assume that $\beta = 1% and NF = NR
%
% Parameters
% WF:       A NxT matrix of work values in the forward direction, 
%           where N is the iteration and T is the time step.
% WR:       A NxT matrix of work values in the reverse direction.
%
% Output
% Ft:       The free energy difference estimate
function Ft = Chelli(WF,WR)

steps = size(WF,2);
Ft = zeros(1,steps);
fermi = @(x) 1./(1+exp(x));

deltaF = BAR(1,WF(:,end),WR(:,end));

for t = 1:steps
 WFAq = WF(:,t);
 WFqb = WF(:,steps)-WFAq;
 WRBq = WR(:,1+steps-t);
 WRqa = WR(:,steps)-WRBq;

 dFAq = -log(mean(exp(-WFAq)));
 dFBQ = -log(mean(exp(-WRBq)));
 
 % Equation 8, for reference
%dlogL = @(dF) sum(fermi(WFAq*ones(size(dF))-ones(size(WFAq))*dF)) - ...
%    exp(dFBQ)*sum((exp(-WRBq)*ones(size(dF))).*fermi(WRqa*ones(size(dF))+ones(size(WRqa))*dF));

 % Equation 16
 dlogL = @(dF) sum(fermi(WFAq*ones(size(dF))-ones(size(WFAq))*dF)) - ...
     exp(dFBQ)*sum((exp(-WRBq)*ones(size(dF)))...
     .*fermi(WRqa*ones(size(dF))+ones(size(WRqa))*dF)) - ...
     exp(dFAq)*sum((exp(-WFAq)*ones(size(dF)))...
     .*fermi(WFqb*ones(size(dF))+ones(size(WFqb))*(dF-deltaF))) + ...
     sum(fermi(WRBq*ones(size(dF))-ones(size(WRBq))*(dF-deltaF)));

 Ft(t) = fzero(dlogL,dFAq);
end