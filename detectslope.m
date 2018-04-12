function [bestSlope SlopeEst]=detectslope(Y,X,p,isplot,label)
if nargin <3
    p=0.1;isplot=0; label='Correlation Dimension Slope';
end
if nargin <4
    isplot=0;label='Correlation Dimension Slope';
end
if nargin < 5
    label='Correlation Dimension Slope';
end
Len=length(Y);

Slopes =zeros(1,Len);

% Create a vector of pointwise slopes
for i=2:Len-1
    Slopes(i)=(Y(i+1)-Y(i-1))/(X(i+1)-X(i-1));
end
Slopes(1) = (Y(2)-Y(1))/(X(2)-X(1));
Slopes(Len) = (Y(Len)-Y(Len-1))/(X(Len)-X(Len-1));

Slopes(isnan(Slopes))=[];
if isplot
    subplot(121)
    plot(X,Y,'.'); set(gca,'fontsize',5);
    hold on;
    title(label,'FontSize',8)
    ylabel('C(r)','FontSize',8);   
    
    subplot(122)
    plot(X,Slopes,'o','markersize',2);set(gca,'fontsize',5);
    hold on
    plot(X,Slopes,'-c','markersize',2);set(gca,'fontsize',5);
end



[Out]=Palarm(Slopes,p,round(length(Slopes)/100));

Kout=(find(diff(Out)))';

PlatStarts=[1 Kout+1];

PlatEnds=[Kout Len];
if isplot
    subplot(122)
    for i=1:length(PlatStarts)
        plot(X(PlatStarts(i):PlatEnds(i)),Out(PlatStarts(i):PlatEnds(i)),'color','r','LineWidth',1)
        hold on;
    end
    hold off;
    xlabel('Log(r)','FontSize',8)
    ylabel('Slopes','FontSize',8)
    h=legend('Pointwise Slopes','Detected Slopes');
    set(h,'fontsize',6);
end

if size(X,1)> size(X,2)
    X=X'; Y=Y';
end

% Linear regression
SlopeEst = zeros(size(PlatStarts));
Intercept=SlopeEst;
for i = 1:length(PlatStarts)    
    x = X(PlatStarts(i):PlatEnds(i)); mx=mean(x);
    y = Y(PlatStarts(i):PlatEnds(i));
    SlopeEst(i) = sum( (x-mx).*y )/sum( (x-mx).*x);   
    if isplot
        Intercept(i) = mean(y)- SlopeEst(i)*mean(x);        
        subplot(121)
        plot(x,SlopeEst(i)*x+Intercept(i),'r','linewidth',1);hold on;
    end   
end
if isplot
    subplot(121)
    h=legend('D2','Estimated Slopes');
    set(h,'fontsize',6);
    hold off;
end

lenFit = PlatEnds-PlatStarts+1;
[~, i]=max(lenFit);
bestSlope=SlopeEst(i);

