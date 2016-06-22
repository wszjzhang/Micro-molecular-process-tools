% ---v
% BAR
% ---
%
% Calculates the free energy difference between two states using the
% classic Bennett Acceptance Ratio method, as extended by Crooks.  The
% equation is solved by finding the zero of the implicit function.
%
% C. Bennett. Efficient Estimation of Free Energy Differences from Monte
% Carlo Data. Journal of Computational Physics 22, 245-268 (1976).
% G. Crooks. Path-ensemble averages in systems driven far from equilibrium.
% Physical Review E 61, 2361-2366 (2000). 
% 
% deltaF = BAR_sc(beta,WF,WR)
%
% Parameters
% beta:     The inverse temperature
% WF:       The total work done in forward processes
% WR:       The total work done in reverse processes
%
% Output
% deltaF:   The free energy difference estimate
function deltaF =  BAR(beta,WF,WR)

%%% Use Jarzynski's equality in the forward direction as
%%% an initial estimate of deltaF
deltaF = -log(mean(exp(-beta*WF)))/beta;

fermi = @(x) 1./(1+exp(x));
dlogL = @(dF) sum(fermi(beta*WF*ones(size(dF))-ones(size(WF))*dF)) - ...
    sum(fermi(beta*WR*ones(size(dF))+ones(size(WR))*dF));

deltaF = fzero(dlogL,deltaF);