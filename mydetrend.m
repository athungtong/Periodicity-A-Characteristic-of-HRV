function [yout yhat]=mydetrend(y,d,x,pf,m)
%%
if nargin <2,d=1;x=(1:length(y))'; end
if nargin <3,x=(1:length(y))'; end

if nargin<4
    pf=0.1;m=3;
end

if size(y,2)>size(y,1),y=y';end
if size(x,2)>size(x,1),x=x';end

if length(y)==1,yout=0;yhat=y;return;end

%
[out] = Palarm(y,pf,m);
Kout=(find(diff(out)))';

PlatStarts=[1 Kout+1];

PlatEnds=[Kout length(y)];

yout=zeros(size(y));
yhat = yout;
for i = 1:length(PlatStarts) 
    yi = y(PlatStarts(i):PlatEnds(i));
    if length(yi)==1, yhat(PlatStarts(i):PlatEnds(i))=yi;continue;end
    
    if d==1
    %This is for straigth line fit
        xi = x(PlatStarts(i):PlatEnds(i));   
        meanx=mean(xi); xdat=xi-meanx;   
        m = sum( xdat.*yi )/sum( (xdat).*xi);
        yhat(PlatStarts(i):PlatEnds(i))=m*xdat+mean(yi);
    else 
        temp = mypolyfit(yi,d);
        yhat(PlatStarts(i):PlatEnds(i)) = temp(:,d);
    end
%     plot(PlatStarts(i):PlatEnds(i),yhat(:,d),'g');
end

yout=y-yhat;

