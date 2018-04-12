%[px xbin atxbin dx] = get_margin_prob(x,nbin)
%return marginal probability based on histogram calculation
%given x, and number of bin 'nbin'
%px is marginal prob. 
%xbin is bin range of x
% atxbin is number of x assign to each x
% since this based on histogram, the total area under histogram
% must be equal to 1. the total area is sum(px*dx)
function [px xbin atxbin dx] = get_margin_prob(x,nbin,xbin)

if nargin == 1
    nbin=round(sqrt(length(x)));
end
if nargin < 3
    x0 = min(x); xn = max(x); dx = (xn-x0)/(nbin-1);
%     xbin = (x0-dx/2:dx:xn+dx/2)';
    xbin=linspace(x0-dx/2,xn+dx/2,nbin+1);
end
[px atxbin]= histc(x,xbin);
px = px(1:end-1)/length(x); %normalize so that the 'area (height*width)' of the histogram equal to 1
end
