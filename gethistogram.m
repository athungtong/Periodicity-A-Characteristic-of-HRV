function [f edges]=gethistogram(data,dx,minx,maxx)
% [f edges]=gethistogram(data,dx) return histogram and bin start from min
% to max of the data.

if nargin<2
    minx=min(data);
    maxx=max(data);
    dx= (maxx-minx)/20;
end

if nargin < 3
    minx=min(data);
    maxx=max(data);
end

edges = minx:dx:maxx;
f=histc(data,edges);
%f=f/length(data);
