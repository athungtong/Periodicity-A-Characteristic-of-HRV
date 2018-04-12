function [IXY HX HY HXY xbin px ybin py pxy] = mutual_info(x,y,normal)
% [IXY HX HY HXY] = mutual_info(x,y) computes mutual information IXY of x
% and y with #bin=round(sqrt(length(x))) and normalize the result based on the 'normal' option (default=1). 
% normal = 0: No normalized
%
% normal = 1: Normalized mutual information to be [0 1] interval
% IXY_norm = (HX+HY-HXY)/(0.5*(HX+HY));
% ref: Zhi-Yuan Su; Tzuyin Wu, Yeng-Tseng Wang, Hsin-Yi Huang. An
% investigation into the linear and non linear correlation of two music
% walk sequences
% Physica D 237 (2008) 1815-1824
%
% normal = 2: Normalized mutual information to be [0 1] interval method 2
% basic idea is to use log base equal to number of bin
% Q. Wang et al. Physica D 200 (2005) 287-295
%important note: not gaurantee mi=1 for the identical signal

if nargin == 2
    normal = 1;
end

% x = (x-mean(x))/std(x);
% y = (y-mean(y))/std(y);

% nxbin=GSSEntropy(@bar_height_entropy,x,0,length(x)/8,.1);
% nybin=GSSEntropy(@bar_height_entropy,y,0,length(y)/8,.1);
% nxbin=round(sqrt(length(x)));
nxbin=16;
nybin=nxbin;
[px xbin atxbin] = get_margin_eq_prob(x,nxbin);
[py ybin atybin]= get_margin_eq_prob(y,nybin);
% pxy = get_joint_prob(x,y,xbin,ybin,atybin,nxbin,nybin);
[Hxy pxy] = jointprob([],[],[],length(x),atxbin,nxbin,atybin,nybin);
% surf(xbin,ybin,pxy

[IXY HX HY HXY] = get_shannon(px,py,pxy,normal);


xbin(end)=[]; ybin(end)=[];

end
