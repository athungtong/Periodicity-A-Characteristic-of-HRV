function [mstar norm norm2 X r]=fnn(s,whatnorm,t,thres,isplot,r,M)
% fnn=fnn(s,r,M) try to find best embedding dimension of the time series s
% based on the FNN method suggested by Kantz. 
% 
% Input Arguments
% s is the time series, 
% whatnorm is type of norm (max or Euc)
% t is the time delay. Default value is 1
% thres is the threshold used (percent) used to make decision if m is good
% enough. Default value is 5 percent.
% r is the tolorent vector the unit is "times" for example, r=1 mean equal,
% r=2 means double and so on. Default value is r=2:0.05:4
% M is the maximum dimension one wants to test (start from m=2 to M-1).
% Default is 10
% 
% Output Arguments
% mstar is the optimal dimension m
% norm is the norm of the reconstructed phase plan with mstar
% norm2 is the norm of mstar+1
% X is M-2 x length(r) matrix. 
% r is the again the tolorence.
% 
% The idea of this method is the following. for a given dimension m, find
% distance to the nearest neighbor of all points i. Then change m to m+1
% and find the distance to the nearest neighbor of all point i again (note:
% number of points decrease by t when increase m to m+1). Then for each
% point, compare the distance to NN at m and m+1. If the ratio
% RNN(m+1)/RNN(m) > r, we consider the neighbor at point i at m is FNN.
% Then we count all points that the ratio is greater than r and divided by
% all point times 100 to find the percent FNN. We do this for
% r=r(1)....r(end). In other word, we compute
% X(i,j)=sum( (r(i)<=rnn2(m(j+1)./rnn))/length(rnn(m(j)))*100;
% if X(i,j)< thres, then we stop and m(j) is the optimum dimension
%
%Another way to determine the embedded dimension is to find corre dimen or
%entropy of different m starting from m=2....M. then see which m that
%doesn't change the corre dimen or entropy. 
% 
% Reference: (not sure this is matched) Improved false nearest neighbor method to detect determinism in time series data 
% Hegger, Rainer and Kantz, Holger
% Physical Review E Statistical Physics Plasmas Fluids And Related Interdisciplinary Topics (1999)
% Volume: 60, Issue: 4 Pt B, Pages: 4970-4973

switch nargin 
    case{1}
        whatnorm = 'max'; t = 1; thres=5; r=2:0.05:4; M = 10; isplot=1;
    case{2}
        t= 1; thres=5; M=10; r= 2:0.05:4;    isplot=1;
    case{3}
        thres=5; M=10; r= 2:0.05:4;  isplot=1;
    case{4}
        M=10; r= 2:0.05:4;  isplot=1;
    case{5}
        M=10; r= 2:0.05:4; 
end

s=center(s);
m = 2:M+1;
X = zeros(length(r),length(m)-1);

if strcmp(whatnorm,'max')
    norm = getmaxnorm(s,m(1),t); %start from m=2;
else
    norm = getnorm(s,m(1),t);
end

rnn = getNN2(norm,t);
mstar = m(end); %default value
for j=1:length(m)-1   
    norm2 = get_norm_mplusone(s,m(j+1),t,norm,whatnorm);
    rnn2 = getNN2(norm2,0);   
    for i=1:length(r)
%         X(i,j)=sum( (r(i)<=rnn2./rnn) .* (rnn<=1/r(i)) ) /
%         sum(rnn<=1/r(i)) * 100; %This method is suffer from denom close
%         to 0 so X will swing when r increase
          X(i,j)=sum( (r(i)<=rnn2./rnn))/length(rnn)*100;
        if isnan(X(i,j)),X(i,j)=0;end %set to 0 if denom is 0
    end
    if mean(X(:,j))<thres, mstar = j+1; break;end  % stop if reach optimal m             
    norm=norm2;
    rnn=rnn2(1:end-t);
end
X(:,j+1:end)=[];

if isplot
    subplot(211);
    for i=1:size(X,2)
        plot(r,X(:,i),'.-'); hold on;
        text(r(2+2*i),X(2+2*i,i),['m=' num2str(i+1)]);
    end
    grid on; hold off;
    meaning = ['When m increase from 2 to 3, the distance from point i to its nearest neighbor increases twice (double) for ' num2str(X(1,1)) ' % of all points'];
    text(r(1),X(1,1),meaning);
    xlabel('tolorence r');  ylabel('% False Nearest Neighbors');

    subplot(212)
    plot(m(1:size(X,2)),mean(X),'*-');hold on;
    plot(m(mstar-1),mean(X(:,mstar-1)),'*r'); hold off;
    xlabel('m');ylabel('average % FNN');
    grid on;hold off;
    text(m(mstar-1),mean(X(:,mstar-1)),'optimal m');
end
 
% %  Junk
% DX = zeros(length(r),length(m)-2);
% for j=1:size(X,2)-1
%     DX(:,j)=((X(:,j+1)-X(:,j)));
% end
% 
% subplot(212);
% for i=1:size(DX,2)
%     plot(r,DX(:,i),'.-'); hold on;
%     text(r(2+2*i),DX(2+2*i,i),['m=' num2str(i+1)]);
% end
% hold off;
