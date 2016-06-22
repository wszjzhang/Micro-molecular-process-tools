% -----
% histw
% -----
%
% Calculates the weighted histogram over columns.
%
% [z, xspace] = histw(X,w,xspace)
%
% Parameters
% X:        a matrix containing the coordinate values
% w:        a corresponding matrix containing the weights
% xspace:   the histogram edges
%
% Outputs
% z:        a matrix of weighted histogram counts
% xspace:   the histogram edges
function [z, xspace] = histw(X,w,xspace)

%%% Default histogram edges
if nargin < 3
  allX = X(1:end);
  minX = min(allX);
  maxX = max(allX);
  xspace = linspace(minX,maxX,50);
  xsp = xspace(2) - xspace(1);
else
  minX = min(xspace);
  xsp = xspace(2) - xspace(1);
end

xind = floor((X-minX)/xsp)+1;

z = zeros(length(xspace),size(X,2));
maxInd = length(xspace) + 1;

for a = 1:size(X,1)
  for b = 1:size(X,2)
    if ((xind(a,b)>0)&&(xind(a,b)<maxInd))
      z(xind(a,b),b) = z(xind(a,b),b) + w(a,b);
    end
  end
end